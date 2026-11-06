"""The bare AmBe source of the moderator simulation (sb_lokiambe_moderatorsim): isotropic neutrons at the origin, with
an approximate ISO 8529-1 AmBe energy spectrum."""
import numpy as np
import G4CustomPyGen
from Units import units as Units

# Approximate ISO 8529-1 241Am-Be spectrum: dN/dE (relative, per MeV) at E [MeV], linearly interpolated (mean
# energy 4.19 MeV, structure near 3, 5 and 8 MeV, endpoint ~11 MeV). Replace with the official table if
# available: after >~10 cm of HDPE the leakage spectrum depends little on these details.
AMBE_E_MEV = np.array([0.0, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5,
                       8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.2])
AMBE_DNDE = np.array([0.0, 0.042, 0.046, 0.045, 0.039, 0.028, 0.032, 0.037, 0.036, 0.031, 0.032, 0.034, 0.031,
                      0.025, 0.023, 0.024, 0.025, 0.022, 0.016, 0.012, 0.011, 0.008, 0.004, 0.001, 0.0])


class AmBeModeratorGen(G4CustomPyGen.GenBase):
    def declare_parameters(self):
        self.addParameterInt("table_points", 20000)

    def init_generator(self, gun):
        gun.set_type("neutron")
        gun.set_position(0, 0, 0)
        e = np.linspace(AMBE_E_MEV[0], AMBE_E_MEV[-1], self.table_points)
        pdf = np.interp(e, AMBE_E_MEV, AMBE_DNDE)
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(e))])
        self._e, self._cdf = e, cdf / cdf[-1]

    def generate_event(self, gun):
        energy = np.interp(self.rand(), self._cdf, self._e)
        cos_t = 2 * self.rand() - 1
        phi = 2 * np.pi * self.rand()
        sin_t = np.sqrt(1 - cos_t * cos_t)
        gun.set_direction(sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t)
        gun.set_energy(max(energy, 1e-6) * Units.MeV)
