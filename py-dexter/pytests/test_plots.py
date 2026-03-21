import numpy as np

from dexter.equilibrium import _TEST_NETCDF_PATH as netcdf_path
from dexter import (
    Equilibrium,
    Geometry,
    Qfactor,
    Current,
    Harmonic,
    UnityQfactor,
    ParabolicQfactor,
    NcQfactor,
    LarCurrent,
    NcCurrent,
    NcGeometry,
    CosHarmonic,
    NcHarmonic,
    InitialFlux,
    BoozerInitialConditions,
    Particle,
    InitialFluxArray1,
    IntersectParams,
    Queue,
    QueueInitialConditions,
)
from dexter.types import FluxWall

# Unsure
# def test_lar_geometry():
#     geometry = LarGeometry(baxis=2, raxis=1.75, rwall=0.5)
#     _test_all_geometry_plots(geometry)


def test_nc_geometry():
    geometry = NcGeometry(netcdf_path, "Steffen", "Bicubic")
    _test_all_flux_plots(geometry)
    _test_all_geometry_plots(geometry)
    _test_all_numerical_geometry_plots(geometry)


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


def test_cos_harmonic():
    harmonic = CosHarmonic(1e-3, 3, 2, 0)
    _test_all_harmonic_plots(harmonic)


def test_nc_harmonic(nc_harmonic: NcHarmonic):
    _test_all_harmonic_plots(nc_harmonic)


def test_particle_plots(nc_equilibrium: Equilibrium):
    initial = BoozerInitialConditions(
        t0=0,
        flux0=InitialFlux("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        rho0=1e-4,
        mu0=7e-6,
    )
    particle = Particle(initial)
    particle.integrate(
        nc_equilibrium,
        (0, 100),
        stepping_method=("FixedStep", 0.5),
    )
    particle.plot_evolution()
    particle.plot_evolution(percentage=2)
    particle.plot_evolution(downsample=False)
    particle.plot_poloidal_drift(percentage=2)
    particle.plot_db_drift(nc_equilibrium)


def test_queue_plots(nc_equilibrium: Equilibrium):
    initials = QueueInitialConditions(
        t0=np.zeros(3),
        flux0=InitialFluxArray1("Toroidal", np.linspace(0, 0.001, 3)),
        theta0=np.zeros(3),
        zeta0=np.zeros(3),
        rho0=np.full(3, 1e-6),
        mu0=np.zeros(3),
    )
    queue = Queue(initials)
    intersect_params = IntersectParams("ConstTheta", 0.0, 5)
    queue.intersect(nc_equilibrium, intersect_params)
    queue.plot_const_theta_cartesian_poincare()

    intersect_params = IntersectParams("ConstZeta", 0.0, 5)
    queue.intersect(nc_equilibrium, intersect_params)
    queue.plot_const_zeta_cartesian_poincare(initial=True)


# ================================================================================================


def _test_all_flux_plots(obj: Qfactor | NcGeometry):
    obj.plot_psip_of_psi(points=50, data=True)
    obj.plot_psip_of_psi(points=50, data=False)
    obj.plot_psi_of_psip(points=50, data=True)
    obj.plot_psi_of_psip(points=50, data=False)


def _test_all_geometry_plots(geometry: Geometry):
    geometry.plot_r_of_psi(points=50, data=True)
    geometry.plot_r_of_psi(points=50, data=False)
    geometry.plot_r_of_psip(points=50, data=True)
    geometry.plot_r_of_psip(points=50, data=False)
    geometry.plot_psi_of_r(points=50, data=True)
    geometry.plot_psi_of_r(points=50, data=False)
    geometry.plot_psip_of_r(points=50, data=True)
    geometry.plot_psip_of_r(points=50, data=False)


def _test_all_numerical_geometry_plots(geometry: NcGeometry):
    geometry.plot_flux_surfaces()
    geometry.plot_jacobian()


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


def _test_all_harmonic_plots(harmonic: Harmonic):
    harmonic.plot_alpha_of_psi(points=50, data=True)
    harmonic.plot_alpha_of_psi(points=50, data=False)
    harmonic.plot_alpha_of_psip(points=50, data=True)
    harmonic.plot_alpha_of_psip(points=50, data=False)
    harmonic.plot_dalpha_of_psi(points=50)
    harmonic.plot_dalpha_of_psi(points=50)
    harmonic.plot_dalpha_of_psip(points=50)
    harmonic.plot_dalpha_of_psip(points=50)
    harmonic.plot_phase_of_psi(points=50, data=True)
    harmonic.plot_phase_of_psi(points=50, data=False)
    harmonic.plot_phase_of_psip(points=50, data=True)
    harmonic.plot_phase_of_psip(points=50, data=False)
