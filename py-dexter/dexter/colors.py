from dexter.types import BoundaryExtremumKind, CriticalPointKind, OrbitType


def particle_color(index: int) -> str:
    COLORS = "cbrmk"
    return COLORS[index % len(COLORS)]


def orbit_color(orbit_type: OrbitType) -> str:
    r"""Returns an orbit type's color string.

    Helps with being consistent with the coloring on different plots.
    """
    # In general, Confined = 'bright', Lost = 'deep'
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
            return "xkcd:azure"
        case _:
            return "xkcd:indigo"


def critical_point_color(kind: CriticalPointKind) -> str:
    r"""Returns a critical point type's color string.

    Helps with being consistent with the coloring on different plots.
    """
    match kind:
        case "OpointMin":
            return "xkcd:royal blue"
        case "OpointMax":
            return "xkcd:bright blue"
        case "Xpoint":
            return "r"
        case _:
            return "k"


def boundary_extremum_color(kind: BoundaryExtremumKind) -> str:
    r"""Returns a boundary extremum type's color string.

    Helps with being consistent with the coloring on different plots.
    """
    match kind:
        case "Minimum":
            return "xkcd:royal purple"
        case "Maximum":
            return "xkcd:bright purple"
        case _:
            return "k"
