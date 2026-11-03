"""Defines the different q-factor objects as wrappers over `_PyQfactor`."""

import numpy as np
from semver import Version
from typing import TypeAlias
from functools import cached_property

from dexter._core import _PyQfactor
from dexter._utils import _into_pyarray1f64

from dexter.machine.flux import MagneticFlux, MagneticFluxArray
from dexter.machine.base import MachineObject, Qfactor
from dexter.types import (
    Array1,
    ArrayLike,
    Interpolation1dType,
    MagneticFluxKind,
    NetCDFVersion,
)


class UnityQfactor(MachineObject, Qfactor):
    r"""Analytical q-factor profile of $q=1$ and $\psi_t=\psi_p$.

    Parameters
    ----------
    lcfs
        The last closed flux surface. Only used for bounds checking.

    Example
    -------
    ``` py
    >>> lcfs = dex.MagneticFlux.Toroidal(0.05)
    >>> qfactor = dex.UnityQfactor(lcfs)

    ```
    """

    _r: _PyQfactor

    def __init__(self, lcfs: MagneticFlux) -> None:
        self._r = _PyQfactor.build_unity(lcfs._r)
        MachineObject.__init__(self)
        Qfactor.__init__(self)


class ParabolicQfactor(MachineObject, Qfactor):
    r"""Analytical parabolic q-factor profile.

    Parameters
    ----------
    qaxis
        The q-factor's value at the magnetic axis, $q_{axis}$.
    qlast
        The q-factor's value at the last closed flux surface, $q_{LCFS}$.

    Example
    -------
    ``` py
    >>> lcfs = dex.MagneticFlux.Toroidal(0.05)
    >>> qfactor = dex.ParabolicQfactor(1.1, 3.9, lcfs)

    ```
    """

    _r: _PyQfactor

    def __init__(self, qaxis: float, qlast: float, lcfs: MagneticFlux) -> None:
        self._r = _PyQfactor.build_parabolic(qaxis, qlast, lcfs._r)
        MachineObject.__init__(self)
        Qfactor.__init__(self)


class ArrayQfactor(MachineObject, Qfactor):
    r"""Numerical q-factor profile reconstructed from $\psi-q$ pairs.

    !!! note "Consistency of evaluations"

        Regardless of the kind of the passed `flux_array`, two interpolators are built: a
        $\iota(\psi_t)$ and a $q(\psi_p)$ interpolator, where the values of the other flux are
        obtained by integrating the corresponding interpolator.

        When calling an evaluation method, the kind of the passed `flux` argument determines the
        interpolator that will be used; an evaluation with respect to a **toroidal** flux always uses
        the $\iota(\psi_t)$ interpolator while an evaluation with resepct to a **poloidal** flux
        always uses the $q(\psi_p)$ interpolator. For example, `#!python qfactor.eval_q(psi)` returns
        `#!python 1/qfactor.eval_iota(psi)` and `#!python qfactor.eval_iota_deriv(psip)` returns
        `#!python -qfactor.eval_q_deriv(psip)/qfactor.eval_q(psip)^2`.

        The above rule guarantees that no evaluation will pass through both interpolators (with the
        exception of `eval_deriv_wrt_other`, which is a shorthand of `eval_other -> eval_deriv_of_other`).
        This is important as passing through both interpolators introduces numerical inconsistencies.

        However, when `#!python flux_array.kind() == Toroidal`, the *very construction* of the $q(\psi_p)$
        interpolator is inconsistent , as the extra step of calculating $\psi_p$ from $\psi_t$ is needed
        before constructing the $q(\psi_p)$ interpolator (similarly when `#!python flux_array.kind() ==
        Poloidal` for $\iota(\psi_t)$). Therefore, while the evaluation methods guarantee that the
        correct interpolator will be used, it is up to you to ensure you evaluate the `qfactor` with
        respect to the flux you constructed it with.

        !!! tip

            The `defining_flux` attribute stores the kind of flux w.r.t. which the `ArrayQfactor` was defined.

    Parameters
    ----------
    flux_array
        The magnetic flux array.
    q_array
        The $q$ values array.
    interp_type
        The 1D interpolation type.

    Attributes
    ----------
    interp_type
        The 1D interpolation type.
    defining_flux
        The flux w.r.t. which the `qfactor` was defined.
    psit_array
        The toroidal flux' values.
    psip_array
        The poloidal flux' values.
    q_array
        The q-factor's values.

    Raises
    ------
    Exception
        An exception is raised if:

        + Data arrays are empty,
        + Data arrays have different sizes,
        + The first element of `flux_array` is not `0.0`.
        + `flux_array` values are non-rising.
        + The interpolators fail to build.


    Example
    -------
    ``` py
    >>> psip_values = [0, 0.1, 0.2, 0.3, 0.4]
    >>> q_array = [0.9, 0.7, 0.5, 0.1, -0.2]
    >>> qfactor = dex.ArrayQfactor(
    ...     flux_array=dex.MagneticFluxArray.Poloidal(psip_values),
    ...     q_array=q_array,
    ...     interp_type="Cubic",
    ...     )

    ```
    """

    _r: _PyQfactor

    def __init__(
        self,
        flux_array: MagneticFluxArray,
        q_array: ArrayLike,
        interp_type: Interpolation1dType,
    ) -> None:
        self._r = _PyQfactor.build_array(
            flux_array=flux_array._r,
            q_array=_into_pyarray1f64(np.asarray(q_array)),
            interp_type=interp_type,
        )
        MachineObject.__init__(self)
        Qfactor.__init__(self)

    @cached_property  # immutable type
    def psit_array(self) -> Array1:
        return self._r.get_array("psit_array")

    @cached_property  # immutable type
    def psip_array(self) -> Array1:
        return self._r.get_array("psip_array")

    @cached_property  # immutable type
    def q_array(self) -> Array1:
        return self._r.get_array("q_array")

    @cached_property  # immutable type
    def interp_type(self) -> Interpolation1dType:
        return self._r.interp_type

    @cached_property  # immutable type
    def defining_flux(self) -> MagneticFluxKind:
        return self._r.defining_flux


class NcQfactor(MachineObject, Qfactor):
    r"""Numerical q-factor profile reconstructed from a netCDF file.

    This type is a wrapper type of [`ArrayQfactor`][dexter.ArrayQfactor] where the $\psi$ and $q$
    data are extracted from the netCDF file.

    Parameters
    ----------
    path
        The path to the netCDF file.
    interp_type
        The 1D interpolation type.

    Attributes
    ----------
    path
        The path to the netCDF file.
    netcdf_version
        The netCDF file's version (SemVer).
    interp_type
        The 1D interpolation type.
    psit_array
        The toroidal flux' values.
    psip_array
        The poloidal flux' values.
    q_array
        The q-factor's values.

    Example
    -------
    ``` py
    >>> qfactor = dex.NcQfactor(path, "Cubic")

    ```
    """

    _r: _PyQfactor

    def __init__(self, path: str, interp_type: Interpolation1dType) -> None:
        self._r = _PyQfactor.build_nc(path, interp_type)
        MachineObject.__init__(self)
        Qfactor.__init__(self)

    @cached_property  # immutable type
    def path(self) -> str:
        return self._r.path

    @cached_property  # immutable type
    def netcdf_version(self) -> NetCDFVersion:
        return Version.parse(self._r.netcdf_version)

    @cached_property  # immutable type
    def interp_type(self) -> Interpolation1dType:
        return self._r.interp_type

    @cached_property  # immutable type
    def defining_flux(self) -> MagneticFluxKind:
        return self._r.defining_flux

    @cached_property  # immutable type
    def psit_array(self) -> Array1:
        return self._r.get_array("psit_array")

    @cached_property  # immutable type
    def psip_array(self) -> Array1:
        return self._r.get_array("psip_array")

    @cached_property  # immutable type
    def q_array(self) -> Array1:
        return self._r.get_array("q_array")


QfactorObject: TypeAlias = UnityQfactor | ParabolicQfactor | ArrayQfactor | NcQfactor
r"""Available [`Qfactor`][dexter.machine.base.Qfactor] objects.

+ [`UnityQfactor`][dexter.UnityQfactor]: Analytical q-factor profile of $q=1$ and $\psi_t=\psi_p$.
+ [`ParabolicQfactor`][dexter.ParabolicQfactor]: Analytical parabolic q-factor profile.
+ [`ArrayQfactor`][dexter.ArrayQfactor]: Numerical q-factor profile reconstructed from $\psi-q$ pairs.
+ [`NcQfactor`][dexter.NcQfactor]: Numerical q-factor profile reconstructed from a netCDF file.
"""
