from dexter._core import NcGeometry
from dexter._core import UnityQfactor, ParabolicQfactor, NcQfactor
from dexter._core import LarCurrent, NcCurrent

from dexter.plotters.flux_plotter import _FluxPlotter
from dexter.plotters.qfactor_plotter import _QfactorPlotter
from dexter.plotters.current_plotter import _CurrentPlotter
import matplotlib.pyplot

_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/test_netcdf.nc"
_TOROIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/toroidal_test_netcdf.nc"
_POLOIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/poloidal_test_netcdf.nc"

matplotlib.pyplot.rcParams["text.usetex"] = True

# =============== Attach plotting methods

for geometry in [NcGeometry]:
    setattr(geometry, "plot_psip_of_psi", _FluxPlotter.plot_psip_of_psi)
    setattr(geometry, "plot_psi_of_psip", _FluxPlotter.plot_psi_of_psip)

for qfactor in [UnityQfactor, ParabolicQfactor, NcQfactor]:
    setattr(qfactor, "plot_q_of_psi", _QfactorPlotter.plot_q_of_psi)
    setattr(qfactor, "plot_q_of_psip", _QfactorPlotter.plot_q_of_psip)
    setattr(qfactor, "plot_psip_of_psi", _FluxPlotter.plot_psip_of_psi)
    setattr(qfactor, "plot_psi_of_psip", _FluxPlotter.plot_psi_of_psip)
    setattr(qfactor, "plot_dpsip_dpsi", _QfactorPlotter.plot_dpsip_dpsi)
    setattr(qfactor, "plot_dpsi_dpsip", _QfactorPlotter.plot_dpsi_dpsip)
    setattr(qfactor, "plot_iota_of_psi", _QfactorPlotter.plot_iota_of_psi)
    setattr(qfactor, "plot_iota_of_psip", _QfactorPlotter.plot_iota_of_psip)

for current in [LarCurrent, NcCurrent]:
    setattr(current, "plot_g_of_psi", _CurrentPlotter.plot_g_of_psi)
    setattr(current, "plot_g_of_psip", _CurrentPlotter.plot_g_of_psip)
    setattr(current, "plot_i_of_psi", _CurrentPlotter.plot_i_of_psi)
    setattr(current, "plot_i_of_psip", _CurrentPlotter.plot_i_of_psip)
    setattr(current, "plot_dg_dpsi", _CurrentPlotter.plot_dg_dpsi)
    setattr(current, "plot_dg_dpsip", _CurrentPlotter.plot_dg_dpsip)
    setattr(current, "plot_di_dpsi", _CurrentPlotter.plot_di_dpsi)
    setattr(current, "plot_di_dpsip", _CurrentPlotter.plot_di_dpsip)


__all__ = [
    "_TOROIDAL_TEST_NETCDF_PATH",
    "_POLOIDAL_TEST_NETCDF_PATH",
    "_TEST_NETCDF_PATH",
    "NcGeometry",
    "UnityQfactor",
    "ParabolicQfactor",
    "NcQfactor",
    "LarCurrent",
    "NcCurrent",
]
