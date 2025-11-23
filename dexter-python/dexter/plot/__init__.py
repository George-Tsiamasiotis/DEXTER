from .qfactor_plots import q_plot, psi_plot
from .currents_plots import g_plot, i_plot
from .bfield_plots import b_plot, db_plot
from .harmonic_plots import alpha_plot, phase_plot
from .particle_plots import evolution_plot

from importlib.util import find_spec

# gtk3agg backend needs PyGObject(gi), which needs a C compiler to be installed.
# gkt4agg spams warnings for no reason
if find_spec("gi") is not None:
    import matplotlib
    import matplotlib.pyplot
    import matplotlib

    matplotlib.use("gtk3agg")
    matplotlib.pyplot.rcParams["text.usetex"] = True

__all__ = [
    "q_plot",
    "psi_plot",
    "g_plot",
    "i_plot",
    "b_plot",
    "db_plot",
    "alpha_plot",
    "phase_plot",
    "evolution_plot",
]
