from dexter.equilibrium import _TEST_NETCDF_PATH as netcdf_path
from dexter import (
    UnityQfactor,
    ParabolicQfactor,
    NcQfactor,
    LarCurrent,
    NcCurrent,
    # LarGeometry,
    NcGeometry,
)
from dexter import Current, Qfactor, Geometry
from dexter.types import FluxWall

# Unsure
# def test_lar_geometry():
#     geometry = LarGeometry(baxis=2, raxis=1.75, rwall=0.5)
#     _test_all_geometry_plots(geometry)


def test_nc_geometry():
    geometry = NcGeometry(netcdf_path, "Steffen", "Bicubic")
    _test_all_flux_plots(geometry)
    _test_all_geometry_plots(geometry)


def test_unity_qfactor():
    qfactor = UnityQfactor()
    _test_all_qfactor_plots(qfactor)
    _test_all_flux_plots(qfactor)


def test_parabolic_qfactor():
    flux_wall: FluxWall = ("Toroidal", 0.45)
    qfactor = ParabolicQfactor(1.1, 3.8, flux_wall)
    _test_all_qfactor_plots(qfactor)
    _test_all_flux_plots(qfactor)


def test_nc_qfactor():
    qfactor = NcQfactor(netcdf_path, "Steffen")
    _test_all_qfactor_plots(qfactor)
    _test_all_flux_plots(qfactor)


def test_lar_current():
    current = LarCurrent()
    _test_all_current_plots(current)


def test_nc_current():
    current = NcCurrent(netcdf_path, "Steffen")
    _test_all_current_plots(current)


def _test_all_flux_plots(obj: Qfactor | NcGeometry):
    obj.plot_psip_of_psi(points=50, data=True)
    obj.plot_psip_of_psi(points=50, data=False)
    obj.plot_psi_of_psip(points=50, data=True)
    obj.plot_psi_of_psip(points=50, data=False)


def _test_all_qfactor_plots(qfactor: Qfactor):
    qfactor.plot_q_of_psi(points=50, data=True)
    qfactor.plot_q_of_psi(points=50, data=False)
    qfactor.plot_q_of_psip(points=50, data=True)
    qfactor.plot_q_of_psip(points=50, data=False)
    qfactor.plot_dpsip_dpsi(points=50)
    qfactor.plot_dpsi_dpsip(points=50)
    qfactor.plot_iota_of_psi(points=50)
    qfactor.plot_iota_of_psip(points=50)


def _test_all_current_plots(current: Current):
    current.plot_g_of_psi(points=50, data=True)
    current.plot_g_of_psi(points=50, data=False)
    current.plot_g_of_psip(points=50, data=True)
    current.plot_g_of_psip(points=50, data=False)
    current.plot_i_of_psi(points=50, data=True)
    current.plot_i_of_psi(points=50, data=False)
    current.plot_i_of_psip(points=50, data=True)
    current.plot_i_of_psip(points=50, data=False)

    current.plot_dg_dpsi(points=50)
    current.plot_dg_dpsip(points=50)
    current.plot_di_dpsi(points=50)
    current.plot_di_dpsip(points=50)


def _test_all_geometry_plots(geometry: Geometry):
    geometry.plot_r_of_psi(points=50, data=True)
    geometry.plot_r_of_psi(points=50, data=False)
    geometry.plot_r_of_psip(points=50, data=True)
    geometry.plot_r_of_psip(points=50, data=False)
    geometry.plot_psi_of_r(points=50, data=True)
    geometry.plot_psi_of_r(points=50, data=False)
    geometry.plot_psip_of_r(points=50, data=True)
    geometry.plot_psip_of_r(points=50, data=False)
    geometry.plot_flux_surfaces()
    geometry.plot_jacobian()
