import matplotlib
import matplotlib.pyplot

matplotlib.use("gtk3agg")
matplotlib.pyplot.rcParams["text.usetex"] = True

from dexter.equilibrium import LastClosedFluxSurface
from dexter.equilibrium import Geometry, LarGeometry, NcGeometry
from dexter.equilibrium import Qfactor, UnityQfactor, ParabolicQfactor, NcQfactor
from dexter.equilibrium import Current, LarCurrent, NcCurrent
from dexter.equilibrium import Bfield, LarBfield, NcBfield
from dexter.equilibrium import Harmonic, CosHarmonic, NcHarmonic
from dexter.equilibrium import Perturbation
from dexter.equilibrium import Equilibrium, numerical_equilibrium

from dexter._core import (
    DEFAULT_THETA_PADDING_WIDTH,
    DEFAULT_STEPPING_METHOD,
    DEFAULT_FIRST_STEP,
    DEFAULT_MAX_STEPS,
    DEFAULT_SAFETY_FACTOR,
    DEFAULT_ENERGY_REL_TOL,
    DEFAULT_ENERGY_ABS_TOL,
    DEFAULT_ERROR_REL_TOL,
    DEFAULT_ERROR_ABS_TOL,
)

from dexter.types import (
    ArrayLike,
    Array,
    Array1,
    Array2,
    ArrayShape,
    Canvas,
    MultiCanvas,
    Locator,
    NetCDFVersion,
    EquilibriumType,
    Interp1DType,
    Interp2DType,
    FluxCoordinate,
    FluxState,
    PhaseMethod,
    CoordinateSet,
    Intersection,
    IntegrationStatus,
    OrbitType,
    SteppingMethod,
    Routine,
)

from dexter.simulate import (
    COMs,
    InitialFlux,
    InitialConditions,
    IntersectParams,
    Particle,
    InitialFluxArray,
    QueueInitialConditions,
    Queue,
    plot_energy_contour,
    plot_particle_poloidal_drift,
)

from dexter.common import get_max_threads, set_num_threads

__all__ = [
    # Configuration Constants
    "DEFAULT_THETA_PADDING_WIDTH",
    "DEFAULT_STEPPING_METHOD",
    "DEFAULT_MAX_STEPS",
    "DEFAULT_FIRST_STEP",
    "DEFAULT_SAFETY_FACTOR",
    "DEFAULT_ENERGY_REL_TOL",
    "DEFAULT_ENERGY_ABS_TOL",
    "DEFAULT_ERROR_REL_TOL",
    "DEFAULT_ERROR_ABS_TOL",
    # Free functions
    "get_max_threads",
    "set_num_threads",
    # Types
    "ArrayLike",
    "Array",
    "Array1",
    "Array2",
    "ArrayShape",
    "Canvas",
    "MultiCanvas",
    "Locator",
    "NetCDFVersion",
    "EquilibriumType",
    "Interp1DType",
    "Interp2DType",
    "LastClosedFluxSurface",
    "FluxCoordinate",
    "FluxState",
    "PhaseMethod",
    "CoordinateSet",
    "Intersection",
    "IntegrationStatus",
    "OrbitType",
    "SteppingMethod",
    "Routine",
    # Equilibrium
    "Geometry",
    "Qfactor",
    "Current",
    "Bfield",
    "Harmonic",
    "Perturbation",
    "LarGeometry",
    "NcGeometry",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
    "LarBfield",
    "NcBfield",
    "CosHarmonic",
    "NcHarmonic",
    "Perturbation",
    "Equilibrium",
    "numerical_equilibrium",
    # Simulate
    "COMs",
    "InitialFlux",
    "InitialConditions",
    "IntersectParams",
    "Particle",
    "InitialFluxArray",
    "QueueInitialConditions",
    "Queue",
    "plot_energy_contour",
    "plot_particle_poloidal_drift",
]
