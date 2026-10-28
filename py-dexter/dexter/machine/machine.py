"""Defines `Machine`, a container for all information a device."""

import numpy as np
from typing import override
from pint.facets.plain import PlainQuantity

from dexter.machine.flux import MagneticFlux
from dexter.types import (
    Interpolation1dType,
    Interpolation2dType,
    ParticleSpecies,
    Unit,
)

from dexter.machine.geometries import GeometryObject, NcGeometry
from dexter.machine.qfactors import NcQfactor, QfactorObject
from dexter.machine.currents import CurrentObject, NcCurrent
from dexter.machine.bfields import BfieldObject, NcBfield
from dexter.machine.efields import EfieldObject, NoEfield
from dexter.machine.perturbation import Perturbation

from dexter.machine._registry import _Registry


class Machine:
    r"""A Machine.

    Container for all information about the magnetic field, geometry and perturbations of a device.

    Parameters
    ----------
    qfactor
        The machine's q-factor.
    current
        The machine's current.
    bfield
        The machine's bfield.
    efield
        The machine's efield.
    perturbation
        The machine's perturbation.
    geometry
        The machine's geometry. This field is not required for routines, but is required for
        conversion to laboratory coordinates and SI units.
    species
        The particle species under study. This only affects unit conversions to SI.

    Attributes
    ----------
    qfactor
        The machine's q-factor.
    current
        The machine's current.
    bfield
        The machine's bfield.
    efield
        The machine's efield.
    perturbation
        The machine's perturbation.
    geometry
        The machine's geometry. This field is not required for routines, but is required for
        conversion to laboratory coordinates and SI units.
    species
        The particle species under study. This only affects unit conversions to SI.
    psit_last
        The value of the last closed toroidal flux $\psi_{t,LCFS}$.
    psip_last
        The value of the last closed poloidal flux $\psi_{p,LCFS}$.
    baxis
        The magnetic field strength on the axis, $B_{axis}$.
    raxis
        The device's major radius, $R_{axis}$.
    zaxis
        The vertical position of the magnetic axis in $[m]$.
    rgeo
        The horizontal position of the geometric axis (device major radius) in $[m]$.
    rlast
        The radial coordinate's value at the last closed flux surface, $r_{LCFS}$.

    Example
    -------

    ```python title="Machine creation"
    >>> lcfs = dex.MagneticFlux.Toroidal(0.05)
    >>> machine = dex.Machine(
    ...     qfactor=dex.UnityQfactor(lcfs),
    ...     current=dex.LarCurrent(),
    ...     bfield=dex.LarBfield(),
    ...     perturbation=dex.Perturbation(
    ...          [
    ...             dex.FluteMode(1e-5, lcfs, 2, 1, 0),
    ...             dex.FluteMode(1e-5, lcfs, 3, 1, 0),
    ...          ]
    ...     ),
    ... )

    ```

    """

    geometry: GeometryObject | None
    qfactor: QfactorObject
    current: CurrentObject
    bfield: BfieldObject
    efield: EfieldObject
    perturbation: Perturbation
    _reg = _Registry | None

    def __init__(
        self,
        qfactor: QfactorObject,
        current: CurrentObject,
        bfield: BfieldObject,
        *,
        geometry: GeometryObject | None = None,
        efield: EfieldObject | None = None,
        perturbation: Perturbation | None = None,
        species: ParticleSpecies = "Proton",
    ) -> None:
        self.qfactor = qfactor
        self.current = current
        self.bfield = bfield
        self.geometry = geometry

        if efield is not None:
            self.efield = efield
        else:
            self.efield = NoEfield()

        if perturbation is not None:
            self.perturbation = perturbation
        else:
            self.perturbation = Perturbation()

        if geometry is not None:
            self._reg = _Registry(case_sensitive=False, cache_folder=":auto:")
            self._reg.define_normalizations(geometry.raxis, geometry.baxis, species)
        else:
            self._reg = None

    @classmethod
    def FromNetcdf(
        cls,
        path: str,
        interp1d_type: Interpolation1dType,
        interp2d_type: Interpolation2dType,
    ) -> Machine:
        r"""Constructs a `Machine` from a netCDF file with no perturbations and default settings.

        Parameters
        ----------
        path
            The path to the netCDF file.
        interp1d_type
            The 1D interpolation type for the 1D quantities.
        interp2d_type
            The 2D interpolation type for the 2D quantities.

        Example
        -------

        ```python title="Machine creation from netCDF file"
        >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")

        ```

        """
        return Machine(
            geometry=NcGeometry(path, interp1d_type, interp2d_type),
            qfactor=NcQfactor(path, interp1d_type),
            current=NcCurrent(path, interp1d_type),
            bfield=NcBfield(path, interp2d_type),
        )

    def quantity[T: float | np.ndarray](
        self,
        value: T,
        units: Unit,
    ) -> PlainQuantity[T]:
        r"""Defines a [`Quantity`](https://pint.readthedocs.io/en/stable/user/defining-quantities.html).

        This constructor is aware of the specific [`machine's normalizations`][dexter.types.Unit] on top of
        all units defined by `pint`.

        !!! important

            The constructor is defined only if the `machine`'s `geometry` is not `None`.

        Parameters
        ----------
        value
            The numeric value(s) of the quantity.
        units
            A string of the quantity's units. Can be any unit that `pint` supports or one of the
            [`normalized units`][dexter.types.Unit].

        Example
        -------

        ``` py
        >>> geometry=dex.LarGeometry(baxis=2, raxis=2.75, rlast=0.5)
        >>> lcfs = geometry.psit_last
        >>> machine = dex.Machine(
        ...     geometry=geometry,
        ...     qfactor=dex.UnityQfactor(lcfs),
        ...     current=dex.LarCurrent(),
        ...     bfield=dex.LarBfield(),
        ... )
        >>>
        >>> machine.quantity(geometry.raxis, "meter")
        Quantity(2.75, "meter")
        >>> machine.quantity(geometry.raxis, "meter").to("NormMeter").magnitude
        1.0
        >>>
        >>> # Electric field strength in normalized units
        >>> efield_strength = machine.quantity(50, "kiloVolt / meter")
        >>> efield_strength.to("NormVolt / NormMeter")
        Quantity(4.7453113229756864e-05, "NormVolt / NormMeter")

        ```
        """
        if not isinstance(self._reg, _Registry):
            raise RuntimeError("UnitRegistry has not been defined")

        return self._reg.Quantity(value, units)

    @property
    def psit_last(self) -> MagneticFlux:
        return self.qfactor.psit_last

    @property
    def psip_last(self) -> MagneticFlux:
        return self.qfactor.psip_last

    @property
    def baxis(self) -> float:
        if self.geometry is not None:
            return self.geometry.baxis
        else:
            raise AttributeError("`baxis` is not defined")

    @property
    def raxis(self) -> float:
        if self.geometry is not None:
            return self.geometry.raxis
        else:
            raise AttributeError("`raxis` is not defined")

    @property
    def zaxis(self) -> float:
        if self.geometry is not None:
            return self.geometry.zaxis
        else:
            raise AttributeError("`zaxis` is not defined")

    @property
    def rgeo(self) -> float:
        if self.geometry is not None:
            return self.geometry.rgeo
        else:
            raise AttributeError("`rgeo` is not defined")

    @property
    def rlast(self) -> float:
        if self.geometry is not None:
            return self.geometry.rlast
        else:
            raise AttributeError("`rlast` is not defined")

    @override
    def __str__(self) -> str:
        string = "Machine:\n"
        if self.geometry is not None:
            string += str(getattr(self, "geometry", "")) + "\n"
        string += str(getattr(self, "qfactor", "")) + "\n"
        string += str(getattr(self, "current", "")) + "\n"
        string += str(getattr(self, "bfield", "")) + "\n"
        string += str(getattr(self, "efield", "")) + "\n"
        string += str(getattr(self, "perturbation", "")) + "\n"
        return string

    @override
    def __repr__(self) -> str:
        return self.__str__()
