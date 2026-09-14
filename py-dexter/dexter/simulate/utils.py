from pint.facets.plain import PlainQuantity
from pint.util import UnitsContainer


def tex_unit(pint_unit: PlainQuantity) -> str:
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
            case "milisecond":
                return r"ms"
            case "second":
                return r"s"
            case _:
                return str(units)
    else:
        assert False, "unimplemented"
