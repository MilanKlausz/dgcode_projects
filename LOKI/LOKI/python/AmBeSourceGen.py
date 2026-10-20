"""AmBe neutron source in a cylindrical HDPE moderator (event_gen=ambe of the LOKI simulation).

The neutrons start on the surface of the moderator cylinder (vertical axis, centred at x/y/z_meters): the curved side
(uniform in height and around the axis), the top and the bottom (uniform on the disks), in proportion to the neutrons
leaving through each surface in the generated energy range. The energy is sampled from the spectrum of the neutrons
leaving through that surface, and the direction from their angular distribution relative to the outward surface
normal (in energy groups; or isotropic in the outward hemisphere: angular_distribution=isotropic). The spectra and
angular distributions are tables in LOKI/data from a Geant4 simulation of the AmBe source in the moderator (see the
header of a table). The moderator itself is not in the geometry. The tables of LOKI/data are for moderator diameters
of 8 to 34 cm (ambe_hdpe_d<diameter>cm_spectrum.txt, height 27 cm); the moderator radius and height are taken from
the header of the table unless given.

Normalisation: one generated neutron corresponds to 1 / (neutrons leaving the generated surfaces in the generated
energy range per source neutron) source neutrons; this yield is printed at the start (an AmBe source of 1 Ci emits
~2.2e6 neutrons per second).
"""
import bisect
import math
import re

import G4CustomPyGen
from Core import FindData
from Units import units

DEFAULT_SPECTRUM = "ambe_hdpe_d23cm_spectrum.txt"
SURFACES = ("side", "top", "bottom")


def read_table(path):
    """The table of a moderator: (lower edges, upper edges, {surface: neutrons per source neutron in the bins}) [eV],
    and the angular distributions {surface: [(lower energy, upper energy, fractions in equal cos bins)]}. Tables with
    a single spectrum column are the spectrum through the side (without angular distributions)."""
    low, high, spectra, angular = [], [], {}, {}
    with open(path) as table:
        for line in table:
            if not line.strip() or line.startswith("#"):
                continue
            words = line.split()
            if words[0] == "angular":
                angular.setdefault(words[1], []).append((float(words[2]), float(words[3]), [float(x) for x in words[4:]]))
                continue
            values = [float(x) for x in words]
            low.append(values[0]), high.append(values[1])
            for surface, value in zip(SURFACES, values[2:]):
                spectra.setdefault(surface, []).append(value)
    return low, high, spectra, angular


def read_moderator_size(path):
    """(radius, height) [cm] of the moderator of a table (its header line "# Moderator: diameter D cm, height H cm");
    None if the table has no such line."""
    with open(path) as table:
        for line in table:
            match = re.match(r"#\s*Moderator: diameter ([0-9.eE+-]+) cm, height ([0-9.eE+-]+) cm", line)
            if match:
                return 0.5 * float(match.group(1)), float(match.group(2))
    return None


class EnergySampler:
    """Samples a binned spectrum (log-uniform within a bin), optionally restricted to [emin, emax] eV."""

    def __init__(self, low, high, value, emin=0.0, emax=float("inf")):
        self.bins = []
        weights = []
        for a, b, v in zip(low, high, value):
            lo, hi = max(a, emin), min(b, emax)
            if v > 0 and hi > lo:
                # the part of the bin in the range (log-uniform density within the bin)
                weights.append(v * math.log(hi / lo) / math.log(b / a))
                self.bins.append((lo, hi))
        self.total = sum(value)  # neutrons per source neutron (all energies)
        self.selected = sum(weights)  # in the selected energy range
        self.cdf = []
        acc = 0.0
        for w in weights:
            acc += w
            self.cdf.append(acc / self.selected if self.selected > 0 else 1.0)

    def sample(self, u1, u2):
        """An energy [eV] for two uniform random numbers."""
        i = min(bisect.bisect_left(self.cdf, u1), len(self.cdf) - 1)
        lo, hi = self.bins[i]
        return lo * (hi / lo) ** u2


class AngularSampler:
    """cos(angle to the outward surface normal) from tabulated distributions (equal cos bins from 0 to 1, uniform
    within a bin) in energy groups."""

    def __init__(self, groups):
        self.groups = []
        for lo, hi, fractions in groups:
            total = sum(fractions)
            cdf, acc = [], 0.0
            for f in fractions:
                acc += f
                cdf.append(acc / total)
            self.groups.append((lo, hi, cdf))

    def sample(self, energy, u1, u2):
        cdf = next((c for lo, hi, c in self.groups if lo <= energy < hi), self.groups[-1][2])
        i = min(bisect.bisect_left(cdf, u1), len(cdf) - 1)
        return (i + u2) / len(cdf)


class AmBeSourceGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        # the centre of the moderator cylinder [m] (vertical axis)
        self.addParameterDouble("x_meters", 0.0)
        self.addParameterDouble("y_meters", 0.0)
        self.addParameterDouble("z_meters", 0.0)
        # the moderator [cm] (-1: from the header of the table)
        self.addParameterDouble("moderator_radius_cm", -1.0)
        self.addParameterDouble("moderator_height_cm", -1.0)
        # the table (in the LOKI data, or a path) and the energy range to generate [eV]
        self.addParameterString("spectrum_file", DEFAULT_SPECTRUM)
        self.addParameterDouble("min_energy_eV", 0.0)
        self.addParameterDouble("max_energy_eV", 1e9)
        # the surfaces to start from (comma separated: side, top, bottom)
        self.addParameterString("surfaces", "side,top,bottom")
        # the direction relative to the outward surface normal: "moderator" (the angular distributions of the table)
        # or "isotropic" (in the outward hemisphere)
        self.addParameterString("angular_distribution", "moderator")

    def init_generator(self, gun):
        gun.set_type("neutron")
        name = self.spectrum_file
        path = name if "/" in name else FindData("LOKI", name)
        low, high, spectra, angular = read_table(path)
        surfaces = [s.strip() for s in self.surfaces.split(",") if s.strip()]
        for surface in surfaces:
            if surface not in SURFACES:
                raise ValueError(f"AmBeSourceGen: unknown surface {surface!r} (side, top, bottom)")
            if surface not in spectra:
                raise ValueError(f"AmBeSourceGen: {name} has no spectrum for the {surface} (a side-only table)")
        if self.angular_distribution not in ("moderator", "isotropic"):
            raise ValueError("AmBeSourceGen: angular_distribution must be moderator or isotropic")
        if self.angular_distribution == "moderator" and any(s not in angular for s in surfaces):
            raise ValueError(f"AmBeSourceGen: {name} has no angular distributions (use angular_distribution=isotropic)")
        self._samplers = {s: EnergySampler(low, high, spectra[s], self.min_energy_eV, self.max_energy_eV) for s in surfaces}
        self._angular = ({s: AngularSampler(angular[s]) for s in surfaces}
                         if self.angular_distribution == "moderator" else None)
        self._surfaces = [s for s in surfaces if self._samplers[s].selected > 0]
        if not self._surfaces:
            raise ValueError(f"AmBeSourceGen: no neutrons between {self.min_energy_eV:g} and {self.max_energy_eV:g} eV")
        self.yield_per_source_neutron = sum(self._samplers[s].selected for s in self._surfaces)
        self._surface_cdf, acc = [], 0.0
        for s in self._surfaces:
            acc += self._samplers[s].selected
            self._surface_cdf.append(acc / self.yield_per_source_neutron)
        size = read_moderator_size(path)
        self._radius_cm = self.moderator_radius_cm if self.moderator_radius_cm > 0 else (size or (None, None))[0]
        self._height_cm = self.moderator_height_cm if self.moderator_height_cm > 0 else (size or (None, None))[1]
        if self._radius_cm is None or self._height_cm is None:
            raise ValueError(f"AmBeSourceGen: give moderator_radius_cm and moderator_height_cm (not in {name})")
        parts = ", ".join(f"{s} {self._samplers[s].selected:.4f}" for s in self._surfaces)
        print(f"AmBeSourceGen: {name}: {self.yield_per_source_neutron:.4f} neutrons per source neutron in "
              f"{self.min_energy_eV:g} - {self.max_energy_eV:g} eV leaving through the generated surfaces ({parts}); "
              f"one generated neutron = {1.0 / self.yield_per_source_neutron:.4f} source neutrons; directions: "
              f"{self.angular_distribution}; moderator radius {self._radius_cm:g} cm, height {self._height_cm:g} cm")

    def generate_event(self, gun):
        surface = self._surfaces[min(bisect.bisect_left(self._surface_cdf, self.rand()), len(self._surfaces) - 1)]
        energy = self._samplers[surface].sample(self.rand(), self.rand())
        r, half = self._radius_cm, 0.5 * self._height_cm
        if surface == "side":
            phi = 2.0 * math.pi * self.rand()
            normal = (math.cos(phi), 0.0, math.sin(phi))
            point = (r * normal[0], (self.rand() - 0.5) * self._height_cm, r * normal[2])
            tangents = ((-math.sin(phi), 0.0, math.cos(phi)), (0.0, 1.0, 0.0))
        else:
            sign = 1.0 if surface == "top" else -1.0
            rho, phi = r * math.sqrt(self.rand()), 2.0 * math.pi * self.rand()
            normal = (0.0, sign, 0.0)
            point = (rho * math.cos(phi), sign * half, rho * math.sin(phi))
            tangents = ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0))
        gun.set_position(self.x_meters * units.m + point[0] * units.cm,
                         self.y_meters * units.m + point[1] * units.cm,
                         self.z_meters * units.m + point[2] * units.cm)
        if self._angular is None:
            cos_t = self.rand()  # isotropic in the outward hemisphere
        else:
            cos_t = self._angular[surface].sample(energy, self.rand(), self.rand())
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        psi = 2.0 * math.pi * self.rand()
        a, b = sin_t * math.cos(psi), sin_t * math.sin(psi)
        gun.set_direction(*(cos_t * n + a * t1 + b * t2 for n, t1, t2 in zip(normal, *tangents)))
        gun.set_energy(energy * units.eV)
