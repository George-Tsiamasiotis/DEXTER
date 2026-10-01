import dexter as dex
import numpy as np


def test_plot_machine_analytical(lar_machine_perturbed: dex.Machine):
    machine = lar_machine_perturbed
    dex.plot_qfactor(machine, points=50, data=True, show=False)
    dex.plot_current(machine, points=50, data=True, show=False)
    dex.plot_bfield(machine, levels=30, show=False)


def test_plot_machine_numerical(nc_machine_perturbed: dex.Machine):
    machine = nc_machine_perturbed
    dex.plot_qfactor(machine, points=50, data=True, show=False)
    dex.plot_current(machine, points=50, data=True, show=False)
    dex.plot_bfield(machine, levels=30, show=False)


def test_plot_mode_analytical():
    LCFS = dex.MagneticFlux.Toroidal(0.03)
    flute_mode = dex.FluteMode(1e-5, LCFS, 5, 2, 0)
    dex.plot_mode(flute_mode, points=50, show=False)

    LCFS = dex.MagneticFlux.Poloidal(0.03)
    flute_mode = dex.FluteMode(1e-5, LCFS, 5, 2, 0)
    dex.plot_mode(flute_mode, points=50, data=True, show=False)


def test_plot_mode_numerical(nc_flute_mode: dex.NcFluteMode):
    dex.plot_mode(nc_flute_mode, points=50, show=False)


def test_plot_particle_analytical(lar_machine: dex.Machine):
    psi_last = lar_machine.psi_last.value
    flux0 = dex.MagneticFlux.Toroidal(psi_last * 0.5)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-6, 1e-7)
    particle = dex.Particle(initial)
    particle.close(lar_machine)
    assert particle.integration_status == "ClosedPeriods(1)"
    assert 100 < particle.steps_taken < 10_000
    dex.plot_evolution(lar_machine, particle, show=False)
    particle.plot_evolution(lar_machine, downsample=True, show=False)
    dex.plot_poloidal_drift(lar_machine, particle, array_shape=(10, 10), show=False)
    particle.plot_poloidal_drift(lar_machine, array_shape=(10, 10), show=False)


def test_plot_particle_numerical(nc_machine: dex.Machine):
    psi_last = nc_machine.psi_last.value
    flux0 = dex.MagneticFlux.Toroidal(psi_last * 0.5)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-6, 1e-7)
    particle = dex.Particle(initial)
    particle.close(nc_machine)
    assert particle.integration_status == "ClosedPeriods(1)"
    assert 100 < particle.steps_taken < 10_000
    dex.plot_evolution(nc_machine, particle, show=False)
    particle.plot_evolution(nc_machine, downsample=True, show=False)
    dex.plot_poloidal_drift(nc_machine, particle, array_shape=(10, 10), show=False)
    particle.plot_poloidal_drift(nc_machine, array_shape=(10, 10), show=False)


def test_plot_queue_analytical(lar_machine: dex.Machine):
    num = 4
    psi0s = dex.MagneticFluxArray.Toroidal(
        np.linspace(1e-6, lar_machine.psi_last.value, num)
    )
    initial = dex.QueueInitialConditions.Boozer(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        rho0=np.full(num, 1e-5),
        mu0=np.full(num, 1e-6),
    )
    queue = dex.Queue(initial)
    intersect_params = dex.IntersectParams("ConstZeta", 0, 10, "Both")
    queue.intersect(lar_machine, intersect_params)

    queue.plot_pzeta_poincare(lar_machine, intersect_params, initial=True)
    queue.plot_rz_poincare(lar_machine, intersect_params, initial=True)


def test_plot_queue_numerical(nc_machine: dex.Machine):
    num = 4
    psi0s = dex.MagneticFluxArray.Toroidal(
        np.linspace(1e-6, nc_machine.psi_last.value, num)
    )
    initial = dex.QueueInitialConditions.Boozer(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        rho0=np.full(num, 1e-5),
        mu0=np.full(num, 1e-6),
    )
    queue = dex.Queue(initial)
    intersect_params = dex.IntersectParams("ConstZeta", 0, 10, "Both")
    queue.intersect(nc_machine, intersect_params)

    queue.plot_pzeta_poincare(nc_machine, intersect_params, initial=True)
    queue.plot_rz_poincare(nc_machine, intersect_params, initial=True)


def test_plot_energy_pzeta_plane_analytical(lar_machine: dex.Machine):
    plane = dex.EnergyPzetaPlane(lar_machine, 1e-5)
    plane.plot()


def test_plot_energy_pzeta_plane_numerical(nc_machine: dex.Machine):
    plane = dex.EnergyPzetaPlane(nc_machine, 1e-5)
    plane.plot()
