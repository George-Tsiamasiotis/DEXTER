"""Defines a Machine's `UnitRegistry` with defined normalized units."""

import numpy as np
from pint import UnitRegistry
from dexter.types import ParticleSpecies

PROTON_MASS = 1.672621923e-27
PROTON_CHARGE = 1.602176634e-19


class _Registry(UnitRegistry[float | np.ndarray]):

    def define_normalizations(
        self, raxis: float, baxis: float, species: ParticleSpecies
    ) -> None:
        self.define(f"Proton_mass   = {PROTON_MASS} kilogram")
        self.define(f"Proton_charge = {PROTON_CHARGE} coulomb")

        M = get_mass_number(species)
        Z = get_charge(species)

        CHARGE = Z * PROTON_CHARGE
        MASS = M * PROTON_MASS

        onaxis_gyrofrequency = CHARGE / MASS * baxis  # s^-1
        energy_unit = MASS * onaxis_gyrofrequency**2 * raxis**2  # Joule

        self.define(f"NormTime       = {1/onaxis_gyrofrequency} second")
        self.define(f"NormFrequency  = {onaxis_gyrofrequency} Hz")
        self.define(f"NormBfield     = {baxis} tesla")
        self.define(f"NormEnergy     = {energy_unit/PROTON_CHARGE} electron_volt")
        self.define(f"NormDistance   = {raxis} meter")
        self.define(f"NormPotential  = NormEnergy / ({PROTON_CHARGE} coulomb)")


def get_mass_number(species: ParticleSpecies) -> float:
    """Returns the mass number of a particle species."""
    match species:
        case "Electron":
            return 0.0005446623
        case "Proton":
            return 1
        case "Deuterium":
            return 2
        case "Tritium":
            return 3
        case "Alpha":
            return 4
        case "He3":
            return 3


def get_charge(species: ParticleSpecies) -> float:
    """Returns the charge of a particle species."""
    match species:
        case "Electron":
            return -1
        case "Proton":
            return 1
        case "Deuterium":
            return 1
        case "Tritium":
            return 1
        case "Alpha":
            return 2
        case "He3":
            return 2
        case _:
            raise TypeError("Invalid particle species")
