"""Common helper types to be used across all submodules."""

import inspect
from typing import Any, TypeAlias

from pint.util import UnitsContainer
from pint.facets.plain import PlainQuantity

_PyAny: TypeAlias = Any
"""The exported type."""


class _ReprStrImpl:
    """Exports the wrapped type's `__str__` and `__repr__`."""

    _r: _PyAny

    def __repr__(self) -> str:
        return self._r.__repr__()

    def __str__(self) -> str:
        return self._r.__repr__()


class _RustTypeWrapper:
    """Wrappers around the exported Rust types."""

    _r: _PyAny

    @classmethod
    def _wrap(cls, _r: _PyAny) -> Any:
        """Creates a wrapped type from the corresponding exported type."""
        obj = cls.__new__(cls)
        obj._r = _r
        return obj


def _get_default_args(func):
    """Returns a `Signature` with the optional parameters' names and default values.

    Useful in providing `argparse` the default values.
    """
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }


def _tex_unit(pint_unit: PlainQuantity) -> str:
    r"""Converts a Quantity's units to a Latex-printable format."""
    units = pint_unit.units
    if pint_unit.is_compatible_with("second"):
        match units:
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
