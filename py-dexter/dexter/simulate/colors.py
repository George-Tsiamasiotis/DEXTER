"""Utilities for colors used in plots."""

from dexter.types import OrbitType


def orbit_color(orbit_type: OrbitType) -> str:
    r"""Confined = 'bright', Lost = 'deep'"""
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
