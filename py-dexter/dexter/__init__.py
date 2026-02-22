import matplotlib.pyplot

matplotlib.pyplot.rcParams["text.usetex"] = True

from dexter.equilibrium import Geometry, LarGeometry, NcGeometry
from dexter.equilibrium import Qfactor, UnityQfactor, ParabolicQfactor, NcQfactor
from dexter.equilibrium import Current, LarCurrent, NcCurrent
from dexter.equilibrium import Bfield, LarBfield, NcBfield
from dexter.equilibrium import Harmonic, CosHarmonic, NcHarmonic
from dexter.equilibrium import Perturbation, CosPerturbation, NcPerturbation
from dexter.equilibrium import NcEquilibrium, perturbation

from dexter.simulate import InitialConditions, Particle

__all__ = [
    # equilibrium
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
    "perturbation",
    "CosPerturbation",
    "NcPerturbation",
    "NcEquilibrium",
    # Simulate
    "InitialConditions",
    "Particle",
]
