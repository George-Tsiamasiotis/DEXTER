from dexter import _TEST_NETCDF_PATH as netcdf_path
from dexter import UnityQfactor, ParabolicQfactor, NcQfactor, LarCurrent, NcCurrent
from dexter.types import Current, Qfactor


def test_unity_qfactor():
    qfactor = UnityQfactor()
    _test_all_qfactor_plots(qfactor)


def test_parabolic_qfactor():
    qfactor = ParabolicQfactor(1.1, 3.9, 0.45)
    _test_all_qfactor_plots(qfactor)


def test_nc_qfactor():
    qfactor = NcQfactor(netcdf_path, "Steffen")
    _test_all_qfactor_plots(qfactor)


def test_lar_current():
    current = LarCurrent()
    _test_all_current_plots(current)


def test_nc_current():
    current = NcCurrent(netcdf_path, "Steffen")
    _test_all_current_plots(current)


def _test_all_qfactor_plots(qfactor: Qfactor):
    qfactor.plot_q_of_psi(points=50, data=True)
    qfactor.plot_q_of_psi(points=50, data=False)
    qfactor.plot_q_of_psip(points=50, data=True)
    qfactor.plot_q_of_psip(points=50, data=False)
    qfactor.plot_psip_of_psi(points=50, data=True)
    qfactor.plot_psip_of_psi(points=50, data=False)
    qfactor.plot_psi_of_psip(points=50, data=True)
    qfactor.plot_psi_of_psip(points=50, data=False)
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
