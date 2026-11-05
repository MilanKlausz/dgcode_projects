"""AmBe neutron source in a cylindrical HDPE moderator (event_gen=ambe of sb_lokiambe_sim).

The neutrons start on the surface of the moderator cylinder (vertical axis, centred at x/y/z_meters): the curved side
(uniform in height and around the axis), the top and the bottom (uniform on the disks), in proportion to the neutrons
leaving through each surface in the generated energy range. The energy is sampled from the spectrum of the neutrons
leaving through that surface, and the direction from their angular distribution relative to the outward surface
normal (in energy groups; or isotropic in the outward hemisphere: angular_distribution=isotropic). The spectra and
angular distributions are a table in LokiAmBe/data, made with the moderator simulation of this package (see the header
of the table and LokiAmBe/scripts/makespectrum). The moderator itself is not in the geometry. The table of the source
of the LoKI detector test (October 2026) is ambe_hdpe_d26cm_h26cm_spectrum.txt (diameter 26 cm, height 26 cm); the
moderator radius and height are taken from the header of the table unless given.

Direction biasing (aim_half_angle_deg > 0): the directions are sampled uniformly within a cone of that half angle
around the direction from the start point to the aim point (aim_x/y/z_meters), with the weight (true angular density
x solid angle of the cone), e.g. to simulate a far bank with enough statistics; the neutrons emitted outside the cone
(which could reach the aim region only by scattering) are not simulated.

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


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(v):
    norm = math.sqrt(sum(c * c for c in v))
    return tuple(c / norm for c in v)


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
            low.append(values[0])
            high.append(values[1])
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

    def density(self, energy, mu):
        """probability density per unit cos(angle to the normal) at mu (0 <= mu <= 1)"""
        cdf = self._cdf_at(energy)
        n = len(cdf)
        i = min(int(mu * n), n - 1)
        return (cdf[i] - (cdf[i - 1] if i > 0 else 0.0)) * n


class AmBeSourceGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        # the centre of the moderator cylinder [m] (vertical axis)
        self.addParameterDouble("x_meters", 0.0)
        self.addParameterDouble("y_meters", 0.0)
        self.addParameterDouble("z_meters", 0.0)
        # the moderator [cm] (-1: from the header of the table)
        self.addParameterDouble("moderator_radius_cm", -1.0)
        self.addParameterDouble("moderator_height_cm", -1.0)
        # the table (in the LokiAmBe data, or a path) and the energy range to generate [eV]
        self.addParameterString("spectrum_file", DEFAULT_SPECTRUM)
        self.addParameterDouble("min_energy_eV", 0.0)
        self.addParameterDouble("max_energy_eV", 1e9)
        # the surfaces to start from (comma separated: side, top, bottom)
        self.addParameterString("surfaces", "side,top,bottom")
        # the direction relative to the outward surface normal: "moderator" (the angular distributions of the table)
        # or "isotropic" (in the outward hemisphere)
        self.addParameterString("angular_distribution", "moderator")
        # direction biasing: a cone of this half angle around the direction to the aim point (0: no biasing)
        self.addParameterDouble("aim_half_angle_deg", 0.0)
        self.addParameterDouble("aim_x_meters", 0.0)
        self.addParameterDouble("aim_y_meters", 0.0)
        self.addParameterDouble("aim_z_meters", 10.0)
        # the particle: "neutron", or "geantino" for line-of-sight (geometry) studies with the same sampling
        self.addParameterString("particle", "neutron")

    def init_generator(self, gun):
        gun.set_type(self.particle)
        name = self.spectrum_file
        path = name if "/" in name else FindData("LokiAmBe", name)
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
        selected = [self._samplers[s].selected for s in self._surfaces]
        self.yield_per_source_neutron = sum(selected)
        self._surface_cdf = _cdf(selected)
        radius, height = read_moderator_size(path) or (None, None)
        self._radius_cm = self.moderator_radius_cm if self.moderator_radius_cm > 0 else radius
        self._height_cm = self.moderator_height_cm if self.moderator_height_cm > 0 else height
        if self._radius_cm is None or self._height_cm is None:
            raise ValueError(f"AmBeSourceGen: give moderator_radius_cm and moderator_height_cm (not in {name})")
        if not 0.0 <= self.aim_half_angle_deg < 90.0:
            raise ValueError("AmBeSourceGen: aim_half_angle_deg must be in [0, 90)")
        self._cos_aim = math.cos(math.radians(self.aim_half_angle_deg))
        self._cone_fraction = (1.0 - self._cos_aim)  # solid angle of the cone / (2 pi)
        self._aim = (self.aim_x_meters * units.m, self.aim_y_meters * units.m, self.aim_z_meters * units.m)
        if self.aim_half_angle_deg > 0:
            print(f"AmBeSourceGen: directions biased into a cone of {self.aim_half_angle_deg:g} deg half angle towards "
                  f"({self.aim_x_meters:g}, {self.aim_y_meters:g}, {self.aim_z_meters:g}) m (weighted)")
        # (the yield in this line is read by LokiOutput/python/AmBeRate.py)
        parts = ", ".join(f"{s} {self._samplers[s].selected:.4f}" for s in self._surfaces)
        print(f"AmBeSourceGen: {name}: {self.yield_per_source_neutron:.4f} neutrons per source neutron in "
              f"{self.min_energy_eV:g} - {self.max_energy_eV:g} eV leaving through the generated surfaces ({parts}); "
              f"one generated neutron = {1.0 / self.yield_per_source_neutron:.4f} source neutrons; directions: "
              f"{self.angular_distribution}; moderator radius {self._radius_cm:g} cm, height {self._height_cm:g} cm")

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
        """A direction from the angular distribution of the surface (or isotropic in the outward hemisphere)."""
        if self._angular is None:
            cos_t = self.rand()
        else:
            cos_t = self._angular[surface].sample(energy, self.rand(), self.rand())
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        psi = 2.0 * math.pi * self.rand()
        a, b = sin_t * math.cos(psi), sin_t * math.sin(psi)
        return tuple(cos_t * n + a * t1 + b * t2 for n, t1, t2 in zip(normal, *tangents))

    def _aimed_direction(self, start, normal, energy, surface):
        """(direction uniform in the cone towards the aim point, weight = true density / cone density)"""
        axis = _unit([a - b for a, b in zip(self._aim, start)])
        helper = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.9 else (0.0, 1.0, 0.0)
        e1 = _unit(_cross(axis, helper))  # an orthonormal basis around the axis
        e2 = _cross(axis, e1)
        cos_c = 1.0 - self.rand() * self._cone_fraction
        sin_c = math.sqrt(max(0.0, 1.0 - cos_c * cos_c))
        psi = 2.0 * math.pi * self.rand()
        direction = tuple(cos_c * a + sin_c * (math.cos(psi) * b + math.sin(psi) * c) for a, b, c in zip(axis, e1, e2))
        mu = sum(a * b for a, b in zip(direction, normal))  # cos of the angle to the outward normal
        if mu <= 0.0:
            return direction, 0.0  # into the moderator: not emitted
        density = 1.0 if self._angular is None else self._angular[surface].density(energy, min(mu, 1.0))
        return direction, density * self._cone_fraction

    def generate_event(self, gun):
        surface = self._surfaces[_pick(self._surface_cdf, self.rand())]
        energy = self._samplers[surface].sample(self.rand(), self.rand())
        start, normal, tangents = self._start_on(surface)
        gun.set_position(*start)
        if self.aim_half_angle_deg > 0:
            direction, weight = self._aimed_direction(start, normal, energy, surface)
            gun.set_weight(weight)
        else:
            direction = self._emitted_direction(normal, tangents, energy, surface)
        gun.set_direction(*direction)
        gun.set_energy(energy * units.eV)
