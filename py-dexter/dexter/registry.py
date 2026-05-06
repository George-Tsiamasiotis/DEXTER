"""Defines an equilibrium-specific UnitRegistry."""

from pint import UnitRegistry
from pint.facets.plain import PlainQuantity

from dexter.types import ArrayLike, ParticleSpecies

PROTON_MASS = 1.672621923e-27
PROTON_CHARGE = 1.602176634e-19


class _Registry(UnitRegistry):
    r"""Defines a custom UnitRegistry with all the additional quantities."""

    raxis: float
    baxis: float
    species: ParticleSpecies
    frequency_unit: float
    energy_unit: float

    def __after_init__(
        self,
        raxis: float,
        baxis: float,
        species: ParticleSpecies = "Proton",
    ):
        r"""Must be called after super().__init__() *and* super()._after_init()."""

        self.define(f"Proton_mass   = {PROTON_MASS} kilogram = NUkilogram")
        self.define(f"Proton_charge = {PROTON_CHARGE} coulomb = NUCoulomb")

        M = _get_mass_number(species)
        Z = _get_charge(species)

        frequency_unit = (Z / M) * PROTON_CHARGE / PROTON_MASS * baxis  # s^-1
        energy_unit = PROTON_MASS * frequency_unit**2 * raxis**2  # Joule

        self.define(f"time_unit         = {frequency_unit} Hz")
        self.define(f"frequency_unit    = {frequency_unit} Hz")
        self.define(f"energy_unit       = {energy_unit/PROTON_CHARGE} electron_volt")
        self.define(f"bfield_unit       = {baxis} Tesla")

        self.raxis = raxis
        self.baxis = baxis
        self.species = species
        self.frequency_unit = frequency_unit
        self.energy_unit = energy_unit

    def quantity(self, value: float | ArrayLike, units: str) -> PlainQuantity:
        """Creates a PlainQuantity."""
        return self.Quantity(value, units)


def _get_mass_number(species: ParticleSpecies) -> float:
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
        case "He3":
            return 3
        case "He":
            return 4
        case _:
            raise TypeError("Invalid particle species")


def _get_charge(species: ParticleSpecies) -> float:
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
        case "He3":
            return 2
        case "He":
            return 2
        case _:
            raise TypeError("Invalid particle species")
