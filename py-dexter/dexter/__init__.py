import matplotlib
import matplotlib.pyplot

matplotlib.use("gtk3agg")
matplotlib.pyplot.rcParams["text.usetex"] = True


from dexter.equilibrium.objects import (
    LastClosedFluxSurface,
    Geometry,
    LarGeometry,
    NcGeometry,
    Qfactor,
    UnityQfactor,
    ParabolicQfactor,
    NcQfactor,
    Current,
    LarCurrent,
    NcCurrent,
    Bfield,
    LarBfield,
    NcBfield,
    Harmonic,
    CosHarmonic,
    NcHarmonic,
    Perturbation,
)

from dexter.equilibrium.equilibrium import (
    Equilibrium,
    numerical_equilibrium,
)

from dexter.types import (
    ArrayLike,
    Array,
    Array1,
    Array2,
    ArrayShape,
    Canvas,
    Canvas3d,
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
    EnergyPzetaPosition,
    OrbitType,
    SteppingMethod,
    Routine,
)

from dexter.simulate.objects import (
    COMs,
    InitialFlux,
    InitialConditions,
    IntersectParams,
    Particle,
    InitialFluxArray,
    QueueInitialConditions,
    Queue,
    Parabola,
    EnergyPzetaPlane,
)

from dexter.simulate.energy_contour import (
    plot_energy_contour,
    plot_particle_poloidal_drift,
)

from dexter.common import get_max_threads, set_num_threads

from dexter.compound_plots.plot_parabolas import (
    plot_parabolas,
    plot_qkinetic_tricontour,
)

__all__ = [
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
    "Canvas3d",
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
    "EnergyPzetaPosition",
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
    "Parabola",
    "EnergyPzetaPlane",
    # Plotting methods
    "plot_parabolas",
    "plot_qkinetic_tricontour",
]
