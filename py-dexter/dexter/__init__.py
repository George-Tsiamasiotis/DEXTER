import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("gtk3agg")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.dpi"] = 180
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.autolayout"] = False
plt.rcParams["figure.constrained_layout.use"] = True

from dexter.types import (
    Array,
    Array1,
    Array2,
    ArrayLike,
    ArrayShape,
    MachineType,
    NetCDFVersion,
    MagneticFluxKind,
    FluxCoordinateState,
    Interpolation1dType,
    Interpolation2dType,
    PhaseMethod,
    CoordinateSet,
    ParticleSpecies,
    IntegrationStatus,
    EnergyPzetaPosition,
    OrbitType,
)

from dexter.utils import get_max_threads, set_num_threads

from dexter.machine.flux import MagneticFlux
from dexter.machine.geometries import GeometryObject, LarGeometry, NcGeometry
from dexter.machine.qfactors import (
    QfactorObject,
    UnityQfactor,
    ParabolicQfactor,
    NcQfactor,
)
from dexter.machine.currents import CurrentObject, LarCurrent, NcCurrent
from dexter.machine.bfields import BfieldObject, LarBfield, NcBfield
from dexter.machine.modes import ModeObject, FluteMode, NcFluteMode
from dexter.machine.perturbation import Perturbation
from dexter.machine.machine import Machine

from dexter.machine.plot import plot_current, plot_qfactor, plot_bfield, plot_mode
from dexter.simulate.plot import (
    plot_evolution,
    plot_poloidal_drift,
    plot_pzeta_poincare,
    plot_rz_poincare,
)

from dexter.simulate.initial import (
    InitialConditions,
    MagneticFluxArray,
    QueueInitialConditions,
)
from dexter.simulate.params import SteppingMethod, SolverParams, IntersectParams
from dexter.simulate.particle import Particle
from dexter.simulate.queue import Queue

from dexter.simulate.energy import (
    create_poloidal_grid,
    energy_of_psi_grid,
    energy_of_psip_grid,
)

from dexter.comspace.energy_pzeta_plane import EnergyPzetaPlane

__all__ = [
    # Type Aliases
    "Array",
    "Array1",
    "Array2",
    "ArrayLike",
    "ArrayShape",
    "MachineType",
    "NetCDFVersion",
    "MagneticFluxKind",
    "FluxCoordinateState",
    "Interpolation1dType",
    "Interpolation2dType",
    "PhaseMethod",
    "CoordinateSet",
    "ParticleSpecies",
    "IntegrationStatus",
    "EnergyPzetaPosition",
    "OrbitType",
    # Utilities
    "get_max_threads",
    "set_num_threads",
    # Machine
    "GeometryObject",
    "QfactorObject",
    "CurrentObject",
    "BfieldObject",
    "ModeObject",
    "MagneticFlux",
    "LarGeometry",
    "NcGeometry",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
    "LarBfield",
    "NcBfield",
    "FluteMode",
    "NcFluteMode",
    "Perturbation",
    "Machine",
    # Machine (plot)
    "plot_current",
    "plot_qfactor",
    "plot_bfield",
    "plot_mode",
    # Simulate
    "InitialConditions",
    "MagneticFluxArray",
    "QueueInitialConditions",
    "SteppingMethod",
    "SolverParams",
    "IntersectParams",
    "Particle",
    "Queue",
    "create_poloidal_grid",
    "energy_of_psi_grid",
    "energy_of_psip_grid",
    # Simulate (plot)
    "plot_evolution",
    "plot_poloidal_drift",
    "plot_pzeta_poincare",
    "plot_rz_poincare",
    # Comspace
    "EnergyPzetaPlane",
]
