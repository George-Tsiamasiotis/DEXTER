from dexter.types import (
    BoundaryExtremumKind,
    CriticalPointKind,
    OrbitType,
)


def particle_color(index: int) -> str:
    """Returns a color string by indexing the `"cbrmk"` color cycle.

    This method is necessary to ensure that each particle is assigned the same color across different
    plots, as its color only depends on its index inside the [`Queue`][dexter.Queue] (e.g. its order
    of instantiation).

    !!! note

        Unlike [`orbit_color`][dexter.colors.orbit_color], the coloring of the particle is arbitrary as it is
        only decided by its position inside the `queue`, however it is guaranteed that the same particle will
        have the same color in every plot.

    !!! tip

        When iterating over the particles to be plotted, use:
        ``` py
        for i, p in enumerate(queue.particles()):
            color = particle_color(i)
            ...
        ```

    Parameters
    ----------
    index
        The particle's index
    """
    COLORS = "cbrmk"
    return COLORS[index % len(COLORS)]


def orbit_color(orbit_type: OrbitType) -> str:
    r"""Returns an orbit type's color string, so that all particles of the same orbit type have the same color.

    Helps with being consistent with the coloring on different plots.

    !!! tip

        When iterating over the particles to be plotted, use:
        ``` py
        for i, p in enumerate(queue.particles()):
            color = orbit_color(p.orbit_type)
            ...
        ```

    Parameters
    ----------
    orbit_type
        The orbit type of the particle.
    """
    # In general, Confined = 'bright', Lost = 'deep'
    match orbit_type:
        case "Undefined":
            return "xkcd:coral"
        case "TrappedConfined":
            return "xkcd:saffron"
        case "TrappedLost":
            return "xkcd:burnt orange"
        case "CoPassingConfined":
            return "xkcd:ocean blue"
        case "CoPassingLost":
            return "xkcd:marine blue"
        case "CuPassingConfined":
            return "xkcd:apple green"
        case "CuPassingLost":
            return "xkcd:deep green"
        case "Potato":
            return "xkcd:tan"
        case "Stagnated":
            return "xkcd:sky"
        case "Unclassified":
            return "xkcd:charcoal"
        case _:
            return "xkcd:indigo"


def critical_point_color(kind: CriticalPointKind) -> str:
    r"""Returns a critical point type's color string.

    Helps with being consistent with the coloring on different plots.

    Parameters
    ----------
    kind
        The kind of the critical point
    """
    match kind:
        case "OpointMin":
            return "xkcd:royal blue"
        case "OpointMax":
            return "xkcd:vivid blue"
        case "Xpoint":
            return "r"
        case _:
            return "k"


def boundary_extremum_color(kind: BoundaryExtremumKind) -> str:
    r"""Returns a boundary extremum type's color string.

    Helps with being consistent with the coloring on different plots.

    Parameters
    ----------
    kind
        The kind of the extremum point
    """
    match kind:
        case "Minimum":
            return "xkcd:purply"
        case "Maximum":
            return "xkcd:bright purple"
        case _:
            return "k"
