import numpy as np
import dexter as dex


def test_energy_pzeta_plane_analytical(lar_machine: dex.Machine):
    plane = dex.EnergyPzetaPlane(lar_machine, 3e-6)
    plane.show(show=False)

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
    queue.classify(lar_machine)

    plane.add_particles(queue[2])
    plane.show(show=False)
    plane.clear_particles()
    plane.add_particles(queue)
    plane.show(show=False)


def test_energy_pzeta_plane_numerical(nc_machine: dex.Machine):
    plane = dex.EnergyPzetaPlane(nc_machine, 3e-6)
    plane.show(show=False)

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
    queue.classify(nc_machine)

    plane.add_particles(queue[2])
    plane.show(show=False)
    plane.clear_particles()
    plane.add_particles(queue)
    plane.show(show=False)
