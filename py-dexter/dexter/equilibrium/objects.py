"""Final wrappers of the exported rust types.

Note
----

Some objects have a `_dyn` class attribute. This is an attempt at reflection for methods that are generic on
the corresponding wrapped types. The string is used to compose the final (manually) monomorphized method and
call it.

Also see `.traits` module.
"""

import numpy as np
from typing import TypeAlias
from numpy.typing import ArrayLike, NDArray

from dexter._core import _PyLarGeometry, _PyNcGeometry
from dexter._core import _PyUnityQfactor, _PyParabolicQfactor, _PyNcQfactor
from dexter._core import _PyLarCurrent, _PyNcCurrent
from dexter._core import _PyLarBfield, _PyNcBfield
from dexter._core import _PyCosHarmonic, _PyNcHarmonic
from dexter._core import _PyCosPerturbation, _PyNcPerturbation
from .traits import (
    _FluxCommuteTrait,
    _GeometryTrait,
    _QfactorTrait,
    _CurrentTrait,
    _BfieldTrait,
    _HarmonicTrait,
)
from ._plotters import (
    _FluxPlotter,
    _GeometryPlotter,
    _NumericalGeometryPlotter,
    _QfactorPlotter,
    _CurrentPlotter,
    _HarmonicPlotter,
)
from dexter.types import (
    ArrayShape,
    Array1,
    Array2,
    FluxWall,
    NetCDFVersion,
    EquilibriumType,
    Interp1DType,
    Interp2DType,
    FluxState,
    PhaseMethod,
)


class LarGeometry(_GeometryTrait, _GeometryPlotter):
    r"""Analytical Large Aspect Ratio Geometry of a circular cross section device.

    Parameters
    ----------
    baxis
        The magnetic field strength on the magnetic axis $B_0$ in $[T]$.
    raxis
        The horizontal position of the magnetic axis $R_0$ in $[m]$.
    rwall
        The value of the $r_{wall}$ coordinate at the wall in $[m]$.

    Example
    -------
    ```python title="LarGeometry creation"
    >>> geometry = dex.LarGeometry(baxis=2, raxis=1.75, rwall=0.5)

    ```
    """

    _dyn: str = "larG"
    _rust: _PyLarGeometry  # type: ignore[assignment]

    def __init__(self, baxis: float, raxis: float, rwall: float) -> None:
        setattr(self, "_rust", _PyLarGeometry(baxis=baxis, raxis=raxis, rwall=rwall))
        _GeometryTrait.__init__(self)

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def baxis(self) -> float:
        """The magnetic field strength on the magnetic axis $B_0$ in $[T]$."""
        return self._rust.baxis

    @property
    def raxis(self) -> float:
        """The horizontal position of the magnetic axis $R_0$ in $[m]$."""
        return self._rust.raxis

    @property
    def rwall(self) -> float:
        """The value of the $r_{wall}$ coordinate at the wall in $[m]$."""
        return self._rust.rwall

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units."""
        return self._rust.psi_wall

    @property
    def rlab_wall(self) -> Array1:
        """Last $R_{lab}$ values that correspond to the device's wall (last closed flux surface)."""
        return self._rust.rlab_wall

    @property
    def zlab_wall(self) -> Array1:
        """Last $Z_{lab}$ values that correspond to the device's wall (last closed flux surface)."""
        return self._rust.zlab_wall

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class NcGeometry(
    _GeometryTrait,
    _FluxCommuteTrait,
    _FluxPlotter,
    _GeometryPlotter,
    _NumericalGeometryPlotter,
):
    r"""Object describing the general geometry of a numerical equilibrium.

    Stores relates scalars and arrays, and provides interpolation methods for converting
    between different variables and coordinate systems.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp1d_type
        The type of 1D Interpolation.
    interp2d_type
        The type of 2D Interpolation.

    Example
    -------
    ```python title="NcGeometry creation"
    >>> geometry = dex.NcGeometry(path, "Cubic", "Bicubic")

    ```
    """

    _dyn: str = "ncdG"
    _rust: _PyNcGeometry  # type: ignore[assignment]

    def __init__(
        self,
        path: str,
        interp1d_type: Interp1DType,
        interp2d_type: Interp2DType,
    ) -> None:
        setattr(
            self,
            "_rust",
            _PyNcGeometry(
                path=path,
                interp1d_type=interp1d_type,
                interp2d_type=interp2d_type,
            ),
        )
        _FluxCommuteTrait.__init__(self)
        _GeometryTrait.__init__(self)

    @property
    def path(self) -> str:
        """The path of the netCDF file."""
        return self._rust.path

    @property
    def netcdf_version(self) -> NetCDFVersion:
        """The netCDF convention version (SemVer)."""
        return self._rust.netcdf_version

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def interp1d_type(self) -> str:
        """The 1D Interpolation type."""
        return self._rust.interp1d_type

    @property
    def interp2d_type(self) -> str:
        """The 2D Interpolation type."""
        return self._rust.interp2d_type

    @property
    def baxis(self) -> float:
        """The magnetic field strength on the magnetic axis $B_0$ in $[T]$."""
        return self._rust.baxis

    @property
    def raxis(self) -> float:
        """The horizontal position of the magnetic axis $R_0$ in $[m]$."""
        return self._rust.raxis

    @property
    def zaxis(self) -> float:
        """The vertical position of the magnetic axis in $[m]$."""
        return self._rust.zaxis

    @property
    def rgeo(self) -> float:
        """The geometrical axis (device major radius) in $[m]$."""
        return self._rust.rgeo

    @property
    def rwall(self) -> float:
        """The value of the $r_{wall}$ coordinate at the wall in $[m]$."""
        return self._rust.rwall

    @property
    def shape(self) -> ArrayShape:
        r"""The shape of the 2 dimensional data arrays, as in $(len(\psi/\psi_p), len(\theta))$.

        If both coordinates are “good”, they are guaranteed to be of the same length.
        """
        return self._rust.shape

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units."""
        return self._rust.psi_wall

    @property
    def psip_wall(self) -> float:
        r"""The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units."""
        return self._rust.psip_wall

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    @property
    def psi_array(self) -> Array1:
        r"""The NetCDF $\psi$ data."""
        return self._rust.psi_array

    @property
    def psip_array(self) -> Array1:
        r"""The NetCDF $\psi_p$ data."""
        return self._rust.psip_array

    @property
    def theta_array(self) -> Array1:
        r"""The NetCDF $\theta$ data."""
        return self._rust.theta_array

    @property
    def r_array(self) -> Array1:
        r"""The NetCDF $r$ data."""
        return self._rust.r_array

    @property
    def rlab_array(self) -> Array2:
        """The NetCDF $R_{lab}$ data array."""
        return self._rust.rlab_array

    @property
    def zlab_array(self) -> Array2:
        """The NetCDF $Z_{lab}$ data array."""
        return self._rust.zlab_array

    @property
    def jacobian_array(self) -> Array2:
        """The Jacobian $J$ data array."""
        return self._rust.jacobian_array

    @property
    def rlab_wall(self) -> Array1:
        """Last $R_{lab}$ values that correspond to the device's wall (last closed flux surface)."""
        return self._rust.rlab_wall

    @property
    def zlab_wall(self) -> Array1:
        """Last $Z_{lab}$ values that correspond to the device's wall (last closed flux surface)."""
        return self._rust.zlab_wall

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


# ================================================================================================


class UnityQfactor(_QfactorTrait, _FluxCommuteTrait, _FluxPlotter, _QfactorPlotter):
    r"""Analytical q-factor profile of $q=1$ and $\psi=\psi_p$.

    Example
    -------
    ```python title="UnityQfactor creation"
    >>> qfactor = dex.UnityQfactor()

    ```
    """

    _dyn: str = "uniQ"
    _rust: _PyUnityQfactor  # type: ignore[assignment]

    def __init__(self) -> None:
        setattr(self, "_rust", _PyUnityQfactor())
        _FluxCommuteTrait.__init__(self)
        _QfactorTrait.__init__(self)

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class ParabolicQfactor(_QfactorTrait, _FluxCommuteTrait, _FluxPlotter, _QfactorPlotter):
    r"""Analytical q-factor of parabolic q(ψ) profile.

    Described by the following formulas:

    Note
    ----
    A ParabolicQfactor is defined with the help of the [`FluxWall`][dexter.FluxWall]
    helper type, which changes the position where qwall is met.

    Parameters
    ----------
    qaxis
        The value of $q$ on the magnetic axis.
    qwall
        The value of $q$ on the wall.
    flux_wall
        Helper type to define a ParabolicQfactor with respect to one of the two fluxes’ values at the wall.

    Example
    -------
    ```python title="UnityQfactor creation"
    >>> # Define q(ψ=ψwall=0.45) = qwall = 3.8
    >>> flux_wall: dex.types.FluxWall = ("Toroidal", 0.45)
    >>> qfactor = dex.ParabolicQfactor(qaxis=1.1, qwall=3.8, flux_wall=flux_wall)

    ```
    """

    _dyn: str = "parQ"
    _rust: _PyParabolicQfactor  # type: ignore[assignment]

    def __init__(
        self,
        qaxis: float,
        qwall: float,
        flux_wall: FluxWall,
    ) -> None:
        setattr(
            self,
            "_rust",
            _PyParabolicQfactor(
                qaxis=qaxis,
                qwall=qwall,
                flux_wall=flux_wall,
            ),
        )
        _FluxCommuteTrait.__init__(self)
        _QfactorTrait.__init__(self)

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def qaxis(self) -> float:
        """The value of $q$ on the magnetic axis."""
        return self._rust.qaxis

    @property
    def qwall(self) -> float:
        """The value of $q$ on the wall."""
        return self._rust.qwall

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units."""
        return self._rust.psi_wall

    @property
    def psip_wall(self) -> float:
        r"""The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units."""
        return self._rust.psip_wall

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class NcQfactor(_FluxCommuteTrait, _QfactorTrait, _FluxPlotter, _QfactorPlotter):
    r"""Numerical q-factor reconstructed from a NetCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D Interpolation type.

    Example
    -------
    ```python title="NcQfactor creation"
    >>> qfactor = dex.NcQfactor(path, "Steffen")

    ```
    """

    _dyn: str = "ncdQ"
    _rust: _PyNcQfactor  # type: ignore[assignment]

    def __init__(
        self,
        path: str,
        interp_type: Interp1DType,
    ) -> None:
        setattr(
            self,
            "_rust",
            _PyNcQfactor(
                path=path,
                interp_type=interp_type,
            ),
        )
        _FluxCommuteTrait.__init__(self)
        _QfactorTrait.__init__(self)

    @property
    def path(self) -> str:
        """The path of the netCDF file."""
        return self._rust.path

    @property
    def netcdf_version(self) -> NetCDFVersion:
        """The netCDF convention version (SemVer)."""
        return self._rust.netcdf_version

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def interp_type(self) -> str:
        """The Interpolation type."""
        return self._rust.interp_type

    @property
    def qaxis(self) -> float:
        """The value of $q$ on the magnetic axis."""
        return self._rust.qaxis

    @property
    def qwall(self) -> float:
        """The value of $q$ on the wall."""
        return self._rust.qwall

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units."""
        return self._rust.psi_wall

    @property
    def psip_wall(self) -> float:
        r"""The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units."""
        return self._rust.psip_wall

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    @property
    def psi_array(self) -> Array1:
        r"""The NetCDF $\psi$ data."""
        return self._rust.psi_array

    @property
    def psip_array(self) -> Array1:
        r"""The NetCDF $\psi_p$ data."""
        return self._rust.psip_array

    @property
    def q_array(self) -> Array1:
        r"""The NetCDF $q$ data."""
        return self._rust.q_array

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


# ================================================================================================


class LarCurrent(_CurrentTrait, _CurrentPlotter):
    """Analytical Large Aspect Ratio Current with $g=1$ and $I=0$.

    Example
    -------
    ```python title="LarCurrent creation"
    >>> current = dex.LarCurrent()

    ```
    """

    _dyn: str = "larC"
    _rust: _PyLarCurrent  # type: ignore[assignment]

    def __init__(self) -> None:
        setattr(self, "_rust", _PyLarCurrent())
        _CurrentTrait.__init__(self)

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class NcCurrent(_CurrentTrait, _CurrentPlotter):
    r"""Numerical plasma current reconstructed from a NetCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The 1D Interpolation type.

    Example
    -------
    ```python title="NcCurrent creation"
    >>> current = dex.NcCurrent(path, "Steffen")

    ```
    """

    _dyn: str = "ncdC"
    _rust: _PyNcCurrent  # type: ignore[assignment]

    def __init__(
        self,
        path: str,
        interp_type: Interp1DType,
    ) -> None:
        setattr(
            self,
            "_rust",
            _PyNcCurrent(
                path=path,
                interp_type=interp_type,
            ),
        )
        _CurrentTrait.__init__(self)

    @property
    def path(self) -> str:
        """The path of the netCDF file."""
        return self._rust.path

    @property
    def netcdf_version(self) -> NetCDFVersion:
        """The netCDF convention version (SemVer)."""
        return self._rust.netcdf_version

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def interp_type(self) -> str:
        """The Interpolation type."""
        return self._rust.interp_type

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units."""
        return self._rust.psi_wall

    @property
    def psip_wall(self) -> float:
        r"""The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units."""
        return self._rust.psip_wall

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    @property
    def psi_array(self) -> Array1:
        r"""The NetCDF $\psi$ data."""
        return self._rust.psi_array

    @property
    def psip_array(self) -> Array1:
        r"""The NetCDF $\psi_p$ data."""
        return self._rust.psip_array

    @property
    def g_array(self) -> Array1:
        r"""The NetCDF $g$ data."""
        return self._rust.g_array

    @property
    def i_array(self) -> Array1:
        r"""The NetCDF $I$ data."""
        return self._rust.i_array

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


# ================================================================================================


class LarBfield(_BfieldTrait):
    r"""Analytical Large Aspect Ratio magnetic field with
    $B(\psi, \theta) = 1 - \sqrt{2\psi}\cos(\theta)$.

    Example
    -------
    ```python title="LarBfield creation"
    >>> bfield = dex.LarBfield()

    ```
    """

    _dyn: str = "larB"
    _rust: _PyLarBfield  # type: ignore[assignment]

    def __init__(self) -> None:
        setattr(self, "_rust", _PyLarBfield())
        _BfieldTrait.__init__(self)

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class NcBfield(_BfieldTrait):
    r"""Numerical magnetic field profile reconstructed from a netCDF file.

    Related quantities are computed by interpolating over the data arrays.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The type of 2D Interpolation.

    Example
    -------
    ```python title="NcBfield creation"
    >>> bfield = dex.NcBfield(path, "Bicubic")

    ```
    """

    _dyn: str = "ncdB"
    _rust: _PyNcBfield  # type: ignore[assignment]

    def __init__(
        self,
        path: str,
        interp_type: Interp2DType,
    ) -> None:
        setattr(
            self,
            "_rust",
            _PyNcBfield(
                path=path,
                interp_type=interp_type,
            ),
        )
        _BfieldTrait.__init__(self)

    @property
    def path(self) -> str:
        """The path of the netCDF file."""
        return self._rust.path

    @property
    def netcdf_version(self) -> NetCDFVersion:
        """The netCDF convention version (SemVer)."""
        return self._rust.netcdf_version

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def interp_type(self) -> str:
        """The Interpolation type."""
        return self._rust.interp_type

    @property
    def baxis(self) -> float:
        """The magnetic field strength on the magnetic axis $B_0$ in $[T]$."""
        return self._rust.baxis

    @property
    def shape(self) -> ArrayShape:
        r"""The shape of the 2 dimensional data arrays, as in $(len(\psi/\psi_p), len(\theta))$.

        If both coordinates are “good”, they are guaranteed to be of the same length.
        """
        return self._rust.shape

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units."""
        return self._rust.psi_wall

    @property
    def psip_wall(self) -> float:
        r"""The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units."""
        return self._rust.psip_wall

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    @property
    def psi_array(self) -> Array1:
        r"""The NetCDF $\psi$ data."""
        return self._rust.psi_array

    @property
    def psip_array(self) -> Array1:
        r"""The NetCDF $\psi_p$ data."""
        return self._rust.psip_array

    @property
    def theta_array(self) -> Array1:
        r"""The NetCDF $\theta$ data."""
        return self._rust.theta_array

    @property
    def b_array(self) -> Array2:
        r"""The NetCDF $g$ data."""
        return self._rust.b_array

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


# ================================================================================================


class CosHarmonic(_HarmonicTrait, _HarmonicPlotter):
    r"""A simple analytical Harmonic of the form `αcos(mθ-nζ+φ)`, where `α`, `φ` constants.

    Parameters
    ----------
    alpha
        The harmonic's constant amplitude $\alpha$ in Normalized Units.
    m
        The poloidal mode number $m$.
    n
        The poloidal mode number $n$.
    phase
        The harmonic's constant phase $\phi$ in $[rads]$.

    Example
    -------
    ```python title="CosHarmonic creation"
    >>> harmonic = dex.CosHarmonic(alpha=1e-3, m=3, n=2, phase=0)

    ```
    """

    _dyn: str = "cosH"
    _rust: _PyCosHarmonic  # type: ignore[assignment]

    def __init__(
        self,
        alpha: float,
        m: int,
        n: int,
        phase: float,
    ) -> None:
        setattr(
            self,
            "_rust",
            _PyCosHarmonic(
                alpha=alpha,
                m=m,
                n=n,
                phase=phase,
            ),
        )
        _HarmonicTrait.__init__(self)

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def alpha(self) -> float:
        r"""The harmonic's constant amplitude $\alpha$ in Normalized Units."""
        return self._rust.alpha

    @property
    def m(self) -> int:
        """The harmonic's poloidal mode number $m$."""
        return self._rust.m

    @property
    def n(self) -> int:
        """The harmonic's toroidal mode number $n$."""
        return self._rust.n

    @property
    def phase(self) -> float:
        r"""The harmonic's constant phase $\phi$ in $[rads]$."""
        return self._rust.phase

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


class NcHarmonic(_HarmonicTrait, _HarmonicPlotter):
    r"""Single perturbation harmonic from a netCDF file.

    !!! info "Phase calculation configuration"

        $\phi$ calculation can be further configured with the `phase_method` optional parameter, which
        defaults to `Resonance`, meaning that a constant value equal to the value of $\phi$ at the
        resonance is used. If no valid value can be found, it falls back to `Zero`. See
        [`PhaseMethod`][dexter.PhaseMethod]) for available configurations.

        !!! example

            ```python title="Phase calculation configuration"
            >>> harmonic = dex.NcHarmonic(path, "cubic", 3, 2, phase_method = "Average")
            >>> harmonic = dex.NcHarmonic(path, "cubic", 3, 2, phase_method = "Interpolation")
            >>> harmonic = dex.NcHarmonic(path, "cubic", 3, 2, phase_method = ("Custom", 3.1415))

            ```

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp_type
        The type of 1D Interpolation.
    m
        The poloidal mode number $m$.
    n
        The poloidal mode number $n$.
    phase_method
        The method of calculation of the phase $\phi(\psi/\psi_p)$. Defaults to "Resonance".
    """

    _dyn: str = "ncdH"
    _rust: _PyNcHarmonic  # type: ignore[assignment]

    def __init__(
        self,
        path: str,
        interp_type: Interp1DType,
        m: int,
        n: int,
        phase_method: PhaseMethod = "Zero",
    ) -> None:
        setattr(
            self,
            "_rust",
            _PyNcHarmonic(
                path=path,
                interp_type=interp_type,
                m=m,
                n=n,
                phase_method=phase_method,
            ),
        )
        _HarmonicTrait.__init__(self)

    @property
    def path(self) -> str:
        """The path of the netCDF file."""
        return self._rust.path

    @property
    def netcdf_version(self) -> NetCDFVersion:
        """The netCDF convention version (SemVer)."""
        return self._rust.netcdf_version

    @property
    def equilibrium_type(self) -> EquilibriumType:
        """The object's equilibrium's type."""
        return self._rust.equilibrium_type

    @property
    def interp_type(self) -> str:
        """The Interpolation type."""
        return self._rust.interp_type

    @property
    def m(self) -> int:
        """The harmonic's poloidal mode number $m$."""
        return self._rust.m

    @property
    def n(self) -> int:
        """The harmonic's toroidal mode number $n$."""
        return self._rust.n

    @property
    def phase_method(self) -> PhaseMethod:
        r"""The method of calculation of the phase $\phi(\psi/\psi_p)$."""
        return self._rust.phase_method

    @property
    def phase_average(self) -> float:
        """The average value of the phase arrays, if `phase_method` is `Average`."""
        return self._rust.phase_average

    @property
    def psi_phase_resonance(self) -> float:
        """The toroidal flux’s value where the resonance is met, if `phase_method` is `Resonance` and
        the resonance is in bounds."""
        return self._rust.psi_phase_resonance

    @property
    def psip_phase_resonance(self) -> float:
        """The poloidal flux’s value where the resonance is met, if `phase_method` is `Resonance` and
        the resonance is in bounds."""
        return self._rust.psip_phase_resonance

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux value at the wall $\psi_{wall}$ in Normalized Units."""
        return self._rust.psi_wall

    @property
    def psip_wall(self) -> float:
        r"""The poloidal flux value at the wall $\psi_{p,wall}$ in Normalized Units."""
        return self._rust.psip_wall

    @property
    def psi_state(self) -> FluxState:
        """The state of the toroidal flux coordinate."""
        return self._rust.psi_state

    @property
    def psip_state(self) -> FluxState:
        """The state of the poloidal flux coordinate."""
        return self._rust.psip_state

    @property
    def psi_array(self) -> Array1:
        r"""The NetCDF $\psi$ data."""
        return self._rust.psi_array

    @property
    def psip_array(self) -> Array1:
        r"""The NetCDF $\psi_p$ data."""
        return self._rust.psip_array

    @property
    def alpha_array(self) -> Array1:
        r"""The NetCDF $\alpha$ data."""
        return self._rust.alpha_array

    @property
    def phase_array(self) -> Array1:
        r"""The NetCDF $\phi$ data."""
        return self._rust.phase_array

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


# ================================================================================================


class Perturbation:
    r"""A sum of an arbitrary number of Harmonics.

    Note
    ----
    All the harmonics must be of the same type.

    Parameters
    ----------
    harmonics
        List of the contained harmonics.

    Example
    -------
    ```python title="Create a Perturbation of specific cosine harmonics"
    >>> perturbation = dex.Perturbation(
    ...     [
    ...         dex.CosHarmonic(1e-3, 1, 2, 0.0),
    ...         dex.CosHarmonic(1e-3, 1, 3, 0.0),
    ...         dex.CosHarmonic(1e-3, 1, 4, 0.0),
    ...         dex.CosHarmonic(1e-3, 1, 5, 0.0),
    ...     ]
    ... )

    ```

    ```python title="CosPerturbation creation with list iteration"
    >>> perturbation = dex.Perturbation(
    ...     [dex.CosHarmonic(1e-3, 1, n, 0.0) for n in range(1, 8)] # modes with m=1 and n=1-7
    ... )

    ```

    ```python title="Create a Perturbation from numerical harmonics"
    >>> perturbation = dex.Perturbation(
    ...     [
    ...         dex.NcHarmonic(path, "cubic", 2, 1, phase_method="Interpolation"),
    ...         dex.NcHarmonic(path, "cubic", 2, 2, phase_method="Interpolation"),
    ...         dex.NcHarmonic(path, "cubic", 3, 2, phase_method="Interpolation"),
    ...     ]
    ... )

    ```
    """

    _dyn: str = "cosP"  # treat it as zero CosHarmonic
    _rust: _PyCosPerturbation | _PyNcPerturbation

    def __init__(self, harmonics: list[CosHarmonic] | list[NcHarmonic]) -> None:
        self._harmonics = harmonics
        match harmonics:
            case []:
                setattr(self, "_rust", _PyCosPerturbation(harmonics=[]))
            case [*cos] if all([isinstance(harmonic, CosHarmonic) for harmonic in cos]):
                _harmonics = [harmonic._rust for harmonic in harmonics]
                self._dyn = "cosP"
                setattr(self, "_rust", _PyCosPerturbation(harmonics=_harmonics))  # type: ignore
            case [*nc] if all([isinstance(harmonic, NcHarmonic) for harmonic in nc]):
                _harmonics = [harmonic._rust for harmonic in harmonics]
                self._dyn = "ncdP"
                setattr(self, "_rust", _PyNcPerturbation(harmonics=_harmonics))  # type: ignore
            case _:
                raise TypeError("All harmonics must be of the same type")

        self._p_of_psi = np.vectorize(self._rust.p_of_psi)
        self._p_of_psip = np.vectorize(self._rust.p_of_psip)
        self._dp_dpsi = np.vectorize(self._rust.dp_dpsi)
        self._dp_dpsip = np.vectorize(self._rust.dp_dpsip)
        self._dp_of_psi_dtheta = np.vectorize(self._rust.dp_of_psi_dtheta)
        self._dp_of_psip_dtheta = np.vectorize(self._rust.dp_of_psip_dtheta)
        self._dp_of_psi_dzeta = np.vectorize(self._rust.dp_of_psi_dzeta)
        self._dp_of_psip_dzeta = np.vectorize(self._rust.dp_of_psip_dzeta)
        self._dp_of_psi_dt = np.vectorize(self._rust.dp_of_psi_dt)
        self._dp_of_psip_dt = np.vectorize(self._rust.dp_of_psip_dt)

    @property
    def harmonics(self) -> list[CosHarmonic] | list[NcHarmonic]:
        return self._harmonics

    def p_of_psi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's value $p(\psi, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._p_of_psi(psi, theta, zeta, t)[()]

    def p_of_psip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's value $p(\psi_p, \theta, \zeta, t)$ value in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._p_of_psip(psip, theta, zeta, t)[()]

    def dp_dpsi(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to $\psi$, $\partial p(\psi, \theta, \zeta, t)/\partial\psi$
        in Normalized Units.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_dpsi(psi, theta, zeta, t)[()]

    def dp_dpsip(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to $\psi_p$, $\partial p(\psi_p, \theta, \zeta, t)/\partial \psi_p$
        in Normalized Units.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_dpsip(psip, theta, zeta, t)[()]

    def dp_of_psi_dtheta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to $\theta$, $\partial p(\psi, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_of_psi_dtheta(psi, theta, zeta, t)[()]

    def dp_of_psip_dtheta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to $\theta$, $\partial p(\psi_p, \theta, \zeta, t)/\partial \theta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_of_psip_dtheta(psip, theta, zeta, t)[()]

    def dp_of_psi_dzeta(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to $\zeta$, $\partial p(\psi, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_of_psi_dzeta(psi, theta, zeta, t)[()]

    def dp_of_psip_dzeta(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to $\zeta$, $\partial p(\psi_p, \theta, \zeta, t)/\partial \zeta$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_of_psip_dzeta(psip, theta, zeta, t)[()]

    def dp_of_psi_dt(
        self,
        psi: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to the time $t$, $\partial p(\psi, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi$.

        Parameters
        ----------
        psi
            The toroidal flux $\psi$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_of_psi_dt(psi, theta, zeta, t)[()]

    def dp_of_psip_dt(
        self,
        psip: ArrayLike,
        theta: ArrayLike,
        zeta: ArrayLike,
        t: ArrayLike,
    ) -> NDArray:
        r"""The perturbation's derivative with respect to the time $t$, $\partial p(\psi_p, \theta, \zeta, t)/\partial t$
        in Normalized Units, as a function of $\psi_p$.

        Parameters
        ----------
        psip
            The poloidal flux $\psi_p$ in Normalized Units.
        theta
            The $\theta$ angle in $[rads]$.
        zeta
            The $\zeta$ angle in $[rads]$.
        t
            The time in Normalized Units
        """
        return self._dp_of_psip_dt(psip, theta, zeta, t)[()]

    def __getitem__(self, index: int):
        """Makes the object indexable."""
        return self._rust[index]

    def __len__(self) -> int:
        """Returns the number of the contained harmonics."""
        return self._rust.__len__()

    def __str__(self) -> str:
        return self._rust.__str__()

    def __repr__(self) -> str:
        return self.__str__()


# ================================================================================================

Geometry: TypeAlias = LarGeometry | NcGeometry
"""Available 'Geometry' Objects"""

Qfactor: TypeAlias = UnityQfactor | ParabolicQfactor | NcQfactor
"""Available 'Qfactor' Objects"""

Current: TypeAlias = LarCurrent | NcCurrent
"""Available 'Current' Objects"""

Bfield: TypeAlias = LarBfield | NcBfield
"""Available 'Bfield' Objects"""

Harmonic: TypeAlias = CosHarmonic | NcHarmonic
"""Available 'Harmonic' Objects"""
