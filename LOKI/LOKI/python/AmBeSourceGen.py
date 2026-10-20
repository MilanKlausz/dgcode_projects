"""AmBe neutron source in a cylindrical HDPE moderator (event_gen=ambe of the LOKI simulation).

The neutrons start on the curved side surface of the moderator cylinder (vertical axis, centred at x/y/z_meters;
uniform in height and around the axis), with an isotropic direction in the outward hemisphere, and an energy sampled
from the spectrum of the neutrons leaving the moderator through its side (a table in LOKI/data, from a Geant4
simulation of the AmBe source in the moderator: see the header of the table). The moderator itself is not in the
geometry. Neutrons leaving through the top and bottom are not generated. The tables of LOKI/data are for moderator
diameters of 8 to 34 cm (ambe_hdpe_d<diameter>cm_side_spectrum.txt, height 27 cm); the moderator radius and height
are taken from the header of the table unless given.

Normalisation: one generated neutron corresponds to 1 / (side yield x sampled fraction) source neutrons, with the side
yield (neutrons through the side per source neutron) and the fraction of it in the sampled energy range printed at
the start (an AmBe source of 1 Ci emits ~2.2e6 neutrons per second).
"""
import bisect
import math
import re

import G4CustomPyGen
from Core import FindData
from Units import units

DEFAULT_SPECTRUM = "ambe_hdpe_d23cm_side_spectrum.txt"


def read_spectrum(path):
    """(lower edges, upper edges, neutrons per source neutron) [eV] of a spectrum table."""
    low, high, value = [], [], []
    with open(path) as table:
        for line in table:
            if line.strip() and not line.startswith("#"):
                a, b, v = (float(x) for x in line.split())
                low.append(a), high.append(b), value.append(v)
    return low, high, value


def read_moderator_size(path):
    """(radius, height) [cm] of the moderator of a spectrum table (its header line "# Moderator: diameter D cm,
    height H cm"); None if the table has no such line."""
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
        if not weights:
            raise ValueError(f"the spectrum has no neutrons between {emin} and {emax} eV")
        self.total = sum(value)  # neutrons per source neutron (all energies)
        self.selected = sum(weights)  # in the selected energy range
        self.cdf = []
        acc = 0.0
        for w in weights:
            acc += w
            self.cdf.append(acc / self.selected)

    def sample(self, u1, u2):
        """An energy [eV] for two uniform random numbers."""
        i = min(bisect.bisect_left(self.cdf, u1), len(self.cdf) - 1)
        lo, hi = self.bins[i]
        return lo * (hi / lo) ** u2


class AmBeSourceGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        # the centre of the moderator cylinder [m] (vertical axis)
        self.addParameterDouble("x_meters", 0.0)
        self.addParameterDouble("y_meters", 0.0)
        self.addParameterDouble("z_meters", 0.0)
        # the moderator [cm] (-1: from the header of the spectrum table)
        self.addParameterDouble("moderator_radius_cm", -1.0)
        self.addParameterDouble("moderator_height_cm", -1.0)
        # the spectrum table (in the LOKI data, or a path) and the energy range to generate [eV]
        self.addParameterString("spectrum_file", DEFAULT_SPECTRUM)
        self.addParameterDouble("min_energy_eV", 0.0)
        self.addParameterDouble("max_energy_eV", 1e9)

    def init_generator(self, gun):
        gun.set_type("neutron")
        name = self.spectrum_file
        path = name if "/" in name else FindData("LOKI", name)
        self._sampler = EnergySampler(*read_spectrum(path), self.min_energy_eV, self.max_energy_eV)
        size = read_moderator_size(path)
        self._radius_cm = self.moderator_radius_cm if self.moderator_radius_cm > 0 else (size or (None, None))[0]
        self._height_cm = self.moderator_height_cm if self.moderator_height_cm > 0 else (size or (None, None))[1]
        if self._radius_cm is None or self._height_cm is None:
            raise ValueError(f"AmBeSourceGen: give moderator_radius_cm and moderator_height_cm (not in {name})")
        print(f"AmBeSourceGen: spectrum {name}: {self._sampler.total:.4f} neutrons through the moderator side per "
              f"source neutron, {self._sampler.selected / self._sampler.total * 100:.2f}% of them in "
              f"{self.min_energy_eV:g} - {self.max_energy_eV:g} eV (generated); moderator radius {self._radius_cm:g} cm, "
              f"height {self._height_cm:g} cm")

    def generate_event(self, gun):
        phi = 2.0 * math.pi * self.rand()
        normal = (math.cos(phi), 0.0, math.sin(phi))  # outward normal of the side surface
        y = (self.rand() - 0.5) * self._height_cm
        r = self._radius_cm
        gun.set_position(self.x_meters * units.m + r * normal[0] * units.cm,
                         self.y_meters * units.m + y * units.cm,
                         self.z_meters * units.m + r * normal[2] * units.cm)
        # isotropic in the outward hemisphere: cos(angle to the normal) uniform in [0, 1]
        cos_t = self.rand()
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        psi = 2.0 * math.pi * self.rand()
        # an orthonormal basis (normal, vertical, tangent) of the surface point
        tangent = (-math.sin(phi), 0.0, math.cos(phi))
        gun.set_direction(cos_t * normal[0] + sin_t * math.cos(psi) * tangent[0],
                          sin_t * math.sin(psi),
                          cos_t * normal[2] + sin_t * math.cos(psi) * tangent[2])
        gun.set_energy(self._sampler.sample(self.rand(), self.rand()) * units.eV)
