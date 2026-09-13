import numpy as np

from dexter._core import _PyMagneticFlux

from dexter.types import MagneticFluxKind, Array1
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
        obj = MagneticFlux.__new__(MagneticFlux)
        obj._r = _PyMagneticFlux.toroidal(value)
        obj.value = obj._r.value
        obj.kind = obj._r.kind

        return obj

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
        obj = MagneticFlux.__new__(MagneticFlux)
        obj._r = _PyMagneticFlux.poloidal(value)
        obj.value = obj._r.value
        obj.kind = obj._r.kind

        return obj

    def __eq__(self, other: object) -> bool:
        """Returns `true` if both `kind` and `value` of `other` are equal to `self`."""
        if isinstance(other, MagneticFlux):
            if self.value == other.value and self.kind == other.kind:
                return True
        return False

    @classmethod
    def _wrap(cls, _r: _PyMagneticFlux) -> MagneticFlux:
        """Wraps the `_r` type to a `MagneticFlux`."""
        if _r.kind == "Toroidal":
            return MagneticFlux.Toroidal(_r.value)
        elif _r.kind == "Poloidal":
            return MagneticFlux.Poloidal(_r.value)
        else:
            raise RuntimeError("unreachable")
