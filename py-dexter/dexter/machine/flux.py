import numpy as np

from collections.abc import Callable
from typing import override
from functools import cached_property

from dexter.types import Array1, MagneticFluxKind, ArrayLike

from dexter._core import _PyMagneticFlux, _PyMagneticFluxArray
from dexter._utils import _ReprStrImpl, _RustTypeWrapper, _into_pyarray1f64


class MagneticFlux(_ReprStrImpl, _RustTypeWrapper):
    r"""Helper type to define the magnetic flux $\psi_t$ or $\psi_p$.

    This type is instantiated through the [`Toroidal`][dexter.MagneticFlux.Toroidal] and
    [`Poloidal`][dexter.MagneticFlux.Poloidal] class methods.

    Attributes
    ----------
    value
        The MagneticFlux' value, regardless of its kind.
    kind
        The MagneticFlux' kind.
    """

    _r: _PyMagneticFlux
    value: float
    kind: MagneticFluxKind

    def __init__(self) -> None:
        raise RuntimeError("Cannot instantiate class")

    @classmethod
    def Toroidal(cls, value: float) -> MagneticFlux:
        r"""Defines a toroidal magnetic flux $\psi_t$.

        Parameters
        ----------
        value
            The value of the toroidal magnetic flux.

        Example
        -------
        ```python title="Toroidal MagneticFlux creation"
        >>> psit = dex.MagneticFlux.Toroidal(0.05)

        ```
        """
        return cls._init(value, _PyMagneticFlux.toroidal)

    @classmethod
    def Poloidal(cls, value: float) -> MagneticFlux:
        r"""Defines a poloidal magnetic flux $\psi_p$.

        Parameters
        ----------
        value
            The value of the poloidal magnetic flux.

        Example
        -------
        ```python title="Poloidal MagneticFlux creation"
        >>> psip = dex.MagneticFlux.Poloidal(0.04)

        ```
        """
        return cls._init(value, _PyMagneticFlux.poloidal)

    @classmethod
    def _init(
        cls,
        value: float,
        _method: Callable[[float], _PyMagneticFlux],
    ) -> MagneticFlux:
        """Creates a `MagneticFlux` by calling the appropriate constructor."""
        obj = MagneticFlux.__new__(MagneticFlux)
        obj._r = _method(value)
        obj.value = obj._r.value
        obj.kind = obj._r.kind

        return obj

    def __add__(self, scalar: float) -> MagneticFlux:
        return MagneticFlux._wrap(self._r.__add__(scalar))

    def __sub__(self, scalar: float) -> MagneticFlux:
        return MagneticFlux._wrap(self._r.__sub__(scalar))

    def __neg__(self) -> MagneticFlux:
        return MagneticFlux._wrap(self._r.__neg__())

    def __mul__(self, scalar: float) -> MagneticFlux:
        return MagneticFlux._wrap(self._r.__mul__(scalar))

    def __truediv__(self, scalar: float) -> MagneticFlux:
        return MagneticFlux._wrap(self._r.__truediv__(scalar))

    def __radd__(self, scalar: float) -> MagneticFlux:
        return self.__add__(scalar)

    def __rsub__(self, scalar: float) -> MagneticFlux:
        return self.__sub__(scalar).__neg__()

    def __rmul__(self, scalar: float) -> MagneticFlux:
        return self.__mul__(scalar)

    @override
    def __eq__(self, other: object) -> bool:
        """Returns `true` if both `kind` and `value` of `other` are equal to `self`."""
        if isinstance(other, MagneticFlux):
            if self.value == other.value and self.kind == other.kind:
                return True
        return False

    @override
    @classmethod
    def _wrap(cls, _r: _PyMagneticFlux) -> MagneticFlux:
        """Wraps the `_r` type to a `MagneticFlux`."""
        match _r.kind:
            case "Toroidal":
                return MagneticFlux.Toroidal(_r.value)
            case "Poloidal":
                return MagneticFlux.Poloidal(_r.value)
            case _:
                raise RuntimeError("unreachable")


class MagneticFluxArray(_ReprStrImpl, _RustTypeWrapper):
    r"""A container of multiple magnetic flux values and their kind.

    This type is instantiated through the [`Toroidal`][dexter.MagneticFluxArray.Toroidal] and
    [`Poloidal`][dexter.MagneticFluxArray.Poloidal] class methods.

    Attributes
    ----------
    values
        The contained values as a 1D Array.
    kind
        The magnetic flux' kind.

    """

    _r: _PyMagneticFluxArray

    def __init__(self) -> None:
        raise RuntimeError("Cannot instantiate class")

    @classmethod
    def Toroidal(cls, array: ArrayLike) -> MagneticFluxArray:
        r"""Defines an array of toroidal fluxes $\psi_t$.

        Parameters
        ----------
        array
            The values of the toroidal magnetic fluxes.

        Example
        -------
        ```python title="Toroidal MagneticFluxArray creation"
        >>> psit_array = dex.MagneticFluxArray.Toroidal(np.linspace(0, 0.05, 30))

        ```
        """
        return cls._init(array, _PyMagneticFluxArray.toroidal)

    @classmethod
    def Poloidal(cls, array: ArrayLike) -> MagneticFluxArray:
        r"""Defines an array of poloidal fluxes $\psi_p$.

        Parameters
        ----------
        array
            The values of the poloidal magnetic fluxes.

        Example
        -------
        ```python title="Poloidal MagneticFluxArray creation"
        >>> psip_array = dex.MagneticFluxArray.Poloidal([0.1, 0.2, 0.3])

        ```
        """
        return cls._init(array, _PyMagneticFluxArray.poloidal)

    @classmethod
    def _init(
        cls, array: ArrayLike, _method: Callable[[np.ndarray], _PyMagneticFluxArray]
    ) -> MagneticFluxArray:
        """Creates a `MagneticFluxArray` by calling the appropriate constructor."""
        array = np.atleast_1d(np.asarray(array, dtype=np.float64))
        array = _into_pyarray1f64(array)
        if not np.all(np.isfinite(array)):
            raise ValueError("NaN values encounter in MagneticFluxArray")

        obj = MagneticFluxArray.__new__(MagneticFluxArray)
        obj._r = _method(array)
        return obj

    @property
    def kind(self) -> MagneticFluxKind:
        return self._r.kind

    @cached_property  # immutable type
    def values(self) -> Array1:
        return self._r.values

    def __len__(self) -> int:
        """Returns the number of contained values."""
        return self._r.__len__()

    def __getitem__(self, index: int) -> MagneticFlux:
        """Returns a [`MagneticFlux`][dexter.MagneticFlux] of the same kind and value of
        `#!python self.values[index]`."""
        _r = self._r[index]
        if _r is None:
            raise IndexError
        else:
            return MagneticFlux._wrap(_r)
