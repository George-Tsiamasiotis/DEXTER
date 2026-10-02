from collections.abc import Callable
from typing import override

from dexter.types import MagneticFluxKind

from dexter._core import _PyMagneticFlux
from dexter._utils import _ReprStrImpl, _RustTypeWrapper


class MagneticFlux(_ReprStrImpl, _RustTypeWrapper):
    r"""Helper type to define the magnetic flux $\psi$ or $\psi_p$.

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
        r"""Defines a toroidal `MagneticFlux` $\psi$.

        Parameters
        ----------
        value
            The value of the toroidal magnetic flux.

        Example
        -------
        ```python title="Toroidal MagneticFlux creation"
        >>> flux = dex.MagneticFlux.Toroidal(0.05)

        ```
        """
        return cls._init(value, _PyMagneticFlux.toroidal)

    @classmethod
    def Poloidal(cls, value: float) -> MagneticFlux:
        r"""Defines a poloidal `MagneticFlux` $\psi_p$.

        Parameters
        ----------
        value
            The value of the poloidal magnetic flux.

        Example
        -------
        ```python title="Poloidal MagneticFlux creation"
        >>> flux = dex.MagneticFlux.Poloidal(0.04)

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

    def __mul__(self, scalar: float) -> MagneticFlux:
        """Multiplies the inner value with `scalar` without changing the `kind`."""
        if self.kind == "Toroidal":
            return MagneticFlux.Toroidal(self.value * scalar)
        elif self.kind == "Poloidal":
            return MagneticFlux.Poloidal(self.value * scalar)
        else:
            raise RuntimeError("unreachable")

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
        if _r.kind == "Toroidal":
            return MagneticFlux.Toroidal(_r.value)
        elif _r.kind == "Poloidal":
            return MagneticFlux.Poloidal(_r.value)
        else:
            raise RuntimeError("unreachable")
