import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("gtk3agg")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.dpi"] = 180
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.autolayout"] = False
plt.rcParams["figure.constrained_layout.use"] = True

from typing import TypeAlias

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
    SteppingMethod,
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
from dexter.simulate.plot import plot_evolution

from dexter.simulate.initial import (
    InitialConditions,
    MagneticFluxArray,
    QueueInitialConditions,
)
from dexter.simulate.particle import Particle
from dexter.simulate.queue import Queue

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
    "SteppingMethod",
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
    "Particle",
    "Queue",
    # Simulate (plot)
    "plot_evolution",
]
