"""Utilities for colors used in plots."""

from matplotlib.patches import Patch

from dexter.types import OrbitType


def orbit_color(orbit_type: OrbitType) -> str:
    r"""Returns each orbit type's color string.

    Confined = 'bright', Lost = 'deep'
    """
    match orbit_type:
        case "Undefined":
            return "xkcd:coral"
        case "TrappedConfined":
            return "xkcd:bright red"
        case "TrappedLost":
            return "xkcd:deep red"
        case "CoPassingConfined":
            return "xkcd:bright blue"
        case "CoPassingLost":
            return "xkcd:deep blue"
        case "CuPassingConfined":
            return "xkcd:bright green"
        case "CuPassingLost":
            return "xkcd:deep green"
        case "Potato":
            return "xkcd:tan"
        case "Stagnated":
            return "xkcd:sky"
        case "Unclassified":
            return "xkcd:bright purple"
        case _:  # Failed(..), also update legend handlers
            return "xkcd:indigo"


def _orbit_color_legend_handles(counter: Counter) -> list[Patch]:
    """Constructs a legend handle for each found orbit type.

    To be used in `ax.legend(handles=<>)`/
    """
    res = []
    for orbit_type, counts in counter.items():
        if orbit_type.startswith("Failed"):
            continue
        res.append(
            Patch(color=orbit_color(orbit_type), label=f"{orbit_type} ({counts})")
        )
    return res
