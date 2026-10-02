"""Common helper types to be used across all submodules."""

import inspect
import numpy as np
from collections.abc import Callable
from typing import Any, TypeAlias, Protocol, override
from pint.facets.plain import PlainQuantity

from dexter.types import Array1

_PyAny: TypeAlias = Any
"""The exported type."""


class _ReprStrImpl(Protocol):
    """Exports the wrapped type's `__str__` and `__repr__`."""

    _r: _PyAny

    @override
    def __repr__(self) -> str:
        return self._r.__repr__()

    @override
    def __str__(self) -> str:
        return self._r.__repr__()


class _RustTypeWrapper(Protocol):
    """Wrappers around the exported Rust types."""

    _r: _PyAny

    @classmethod
    def _wrap(cls, _r: _PyAny) -> Any:
        """Creates a wrapped type from the corresponding exported type."""
        obj = cls.__new__(cls)
        obj._r = _r
        return obj


def get_default_args(func: Callable[[Any], Any]):
    """Returns a `Signature` with the optional parameters' names and default values.

    Useful in providing `argparse` the default values.
    """
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }


def tex_unit(pint_unit: PlainQuantity[Any]) -> str:
    r"""Converts a Quantity's units to a Latex-printable format."""
    units = pint_unit.units
    if pint_unit.is_compatible_with("second"):
        match str(units):
            case "femtosecond":
                return r"fm"
            case "picosecond":
                return r"ps"
            case "nanosecond":
                return r"ns"
            case "microsecond":
                return r"\mu s"
            case "millisecond":
                return r"ms"
            case "second":
                return r"s"
            case _:
                return str(units)
    else:
        assert False, "unimplemented"


def _into_pyarray1f64(array: Array1) -> Array1:
    """Makes sure that a numpy array is C-contiguous and of the correct type.

    Use this function whenever Rust expects a `PyArray1<f64>`.
    """
    if array.ndim != 1:
        raise TypeError("Array must be 1-dimensional")
    if not array.flags.c_contiguous:
        array = np.ascontiguousarray(array)

    return array.astype(np.float64)
