"""AmBe neutron source in a cylindrical HDPE moderator (event_gen=ambe of sb_lokiambe_sim).

The neutrons start on the surface of the moderator cylinder (vertical axis, centred at x/y/z_meters): the curved side
(uniform in height and around the axis), the top and the bottom (uniform on the disks), in proportion to the neutrons
leaving through each surface in the generated energy range. The energy is sampled from the spectrum of the neutrons
leaving through that surface, and the direction from their angular distribution relative to the outward surface
normal (in energy groups). The spectra and angular distributions are a table in LokiAmBe/data, made with the moderator
simulation of this package (see the header of the table and LokiAmBe/scripts/makespectrum). The moderator itself is not
in the geometry. The table of the source of the LoKI detector test (October 2026) is
ambe_hdpe_d26cm_h26cm_spectrum.txt (diameter 26 cm, height 26 cm); the moderator radius and height are taken from the
header of the table.

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

DEFAULT_SPECTRUM = "ambe_hdpe_d26cm_h26cm_spectrum.txt"
SURFACES = ("side", "top", "bottom")


def _cdf(values):
    """The normalised cumulative sum of the values."""
    total, acc, cdf = sum(values), 0.0, []
    for value in values:
        acc += value
        cdf.append(acc / total)
    return cdf


def _pick(cdf, u):
    """The index of the bin of a uniform random number."""
    return min(bisect.bisect_left(cdf, u), len(cdf) - 1)


def read_table(path):
    """The table of a moderator: (lower edges, upper edges, {surface: neutrons per source neutron in the bins}) [eV],
    and the angular distributions {surface: [(lower energy, upper energy, fractions in equal cos bins)]}."""
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
            low.append(values[0])
            high.append(values[1])
            for surface, value in zip(SURFACES, values[2:]):
                spectra.setdefault(surface, []).append(value)
    return low, high, spectra, angular


def read_moderator_size(path):
    """(radius, height) [cm] of the moderator of a table (its header line "# Moderator: diameter D cm, height H cm")."""
    with open(path) as table:
        for line in table:
            match = re.match(r"#\s*Moderator: diameter ([0-9.eE+-]+) cm, height ([0-9.eE+-]+) cm", line)
            if match:
                return 0.5 * float(match.group(1)), float(match.group(2))
    raise ValueError(f"AmBeSourceGen: {path} has no line with the moderator size")


class EnergySampler:
    """Samples a binned spectrum (log-uniform within a bin), optionally restricted to energies up to emax [eV]."""

    def __init__(self, low, high, value, emax=float("inf")):
        self.bins = []
        weights = []
        for a, b, v in zip(low, high, value):
            hi = min(b, emax)
            if v > 0 and hi > a:
                # the part of the bin in the range (log-uniform density within the bin)
                weights.append(v * math.log(hi / a) / math.log(b / a))
                self.bins.append((a, hi))
        self.total = sum(value)  # neutrons per source neutron (all energies)
        self.selected = sum(weights)  # in the selected energy range
        self.cdf = _cdf(weights)

    def sample(self, u1, u2):
        """An energy [eV] for two uniform random numbers."""
        lo, hi = self.bins[_pick(self.cdf, u1)]
        return lo * (hi / lo) ** u2


class AngularSampler:
    """cos(angle to the outward surface normal) from tabulated distributions (equal cos bins from 0 to 1, uniform
    within a bin) in energy groups."""

    def __init__(self, groups):
        self.groups = [(lo, hi, _cdf(fractions)) for lo, hi, fractions in groups]

    def _cdf_at(self, energy):
        return next((c for lo, hi, c in self.groups if lo <= energy < hi), self.groups[-1][2])

    def sample(self, energy, u1, u2):
        cdf = self._cdf_at(energy)
        return (_pick(cdf, u1) + u2) / len(cdf)


class AmBeSourceGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        # the centre of the moderator cylinder [m] (vertical axis)
        self.addParameterDouble("x_meters", 0.0)
        self.addParameterDouble("y_meters", 0.0)
        self.addParameterDouble("z_meters", 0.0)
        # the table (in the LokiAmBe data, or a path) and the highest energy to generate [eV]
        self.addParameterString("spectrum_file", DEFAULT_SPECTRUM)
        self.addParameterDouble("max_energy_eV", 1e9)
        # the particle: "neutron", or "geantino" for line-of-sight (geometry) studies with the same sampling
        self.addParameterString("particle", "neutron")

    def init_generator(self, gun):
        gun.set_type(self.particle)
        name = self.spectrum_file
        path = name if "/" in name else FindData("LokiAmBe", name)
        low, high, spectra, angular = read_table(path)
        if any(s not in spectra or s not in angular for s in SURFACES):
            raise ValueError(f"AmBeSourceGen: {name} needs the spectra and angular distributions of the side, top and bottom")
        self._samplers = {s: EnergySampler(low, high, spectra[s], self.max_energy_eV) for s in SURFACES}
        self._angular = {s: AngularSampler(angular[s]) for s in SURFACES}
        self._surfaces = [s for s in SURFACES if self._samplers[s].selected > 0]
        if not self._surfaces:
            raise ValueError(f"AmBeSourceGen: no neutrons below {self.max_energy_eV:g} eV")
        selected = [self._samplers[s].selected for s in self._surfaces]
        self.yield_per_source_neutron = sum(selected)
        self._surface_cdf = _cdf(selected)
        self._radius_cm, self._height_cm = read_moderator_size(path)
        # (the yield in this line is read by LokiOutput/python/AmBeRate.py)
        parts = ", ".join(f"{s} {self._samplers[s].selected:.4f}" for s in self._surfaces)
        print(f"AmBeSourceGen: {name}: {self.yield_per_source_neutron:.4f} neutrons per source neutron below "
              f"{self.max_energy_eV:g} eV leaving the moderator ({parts}); one generated neutron = "
              f"{1.0 / self.yield_per_source_neutron:.4f} source neutrons; moderator radius {self._radius_cm:g} cm, "
              f"height {self._height_cm:g} cm")

    def _start_on(self, surface):
        """A random start point on the surface: (point [Geant4 units], outward normal, two tangents)."""
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
        start = (self.x_meters * units.m + point[0] * units.cm,
                 self.y_meters * units.m + point[1] * units.cm,
                 self.z_meters * units.m + point[2] * units.cm)
        return start, normal, tangents

    def _emitted_direction(self, normal, tangents, energy, surface):
        """A direction from the angular distribution of the surface."""
        cos_t = self._angular[surface].sample(energy, self.rand(), self.rand())
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        psi = 2.0 * math.pi * self.rand()
        a, b = sin_t * math.cos(psi), sin_t * math.sin(psi)
        return tuple(cos_t * n + a * t1 + b * t2 for n, t1, t2 in zip(normal, *tangents))

    def generate_event(self, gun):
        surface = self._surfaces[_pick(self._surface_cdf, self.rand())]
        energy = self._samplers[surface].sample(self.rand(), self.rand())
        start, normal, tangents = self._start_on(surface)
        gun.set_position(*start)
        gun.set_direction(*self._emitted_direction(normal, tangents, energy, surface))
        gun.set_energy(energy * units.eV)
