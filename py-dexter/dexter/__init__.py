import matplotlib
import matplotlib.pyplot

matplotlib.use("gtk3agg")
matplotlib.pyplot.rcParams["text.usetex"] = True

from dexter.equilibrium import Geometry, LarGeometry, NcGeometry
from dexter.equilibrium import Qfactor, UnityQfactor, ParabolicQfactor, NcQfactor
from dexter.equilibrium import Current, LarCurrent, NcCurrent
from dexter.equilibrium import Bfield, LarBfield, NcBfield
from dexter.equilibrium import Harmonic, CosHarmonic, NcHarmonic
from dexter.equilibrium import Perturbation
from dexter.equilibrium import Equilibrium, numerical_equilibrium

from dexter.types import (
    ArrayLike,
    Array,
    Array1,
    Array2,
    ArrayShape,
    NetCDFVersion,
    EquilibriumType,
    Interp1DType,
    Interp2DType,
    FluxWall,
    FluxCoordinate,
    FluxState,
    PhaseMethod,
    CoordinateSet,
    Intersection,
    IntegrationStatus,
    SteppingMethod,
    Routine,
)

from dexter.simulate import (
    InitialFlux,
    InitialConditions,
    IntersectParams,
    Particle,
    InitialFluxArray,
    QueueInitialConditions,
    Queue,
)

__all__ = [
    # Types
    "ArrayLike",
    "Array",
    "Array1",
    "Array2",
    "ArrayShape",
    "NetCDFVersion",
    "EquilibriumType",
    "Interp1DType",
    "Interp2DType",
    "FluxWall",
    "FluxCoordinate",
    "FluxState",
    "PhaseMethod",
    "CoordinateSet",
    "Intersection",
    "IntegrationStatus",
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
    "InitialFlux",
    "InitialConditions",
    "IntersectParams",
    "Particle",
    "InitialFluxArray",
    "QueueInitialConditions",
    "Queue",
]
