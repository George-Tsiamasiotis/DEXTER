import numpy as np
import dexter as dex


def test_queue_routines_analytical(lar_machine: dex.Machine):
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

    queue.integrate(lar_machine, (0, 100))
    queue.intersect(lar_machine, intersection="ConstZeta", angle=0, turns=4)
    queue.close(lar_machine)
    queue.classify(lar_machine)


import numpy as np
import dexter as dex


def test_queue_routines_numerical(nc_machine: dex.Machine):
    num = 4
    psi0s = dex.MagneticFluxArray.Toroidal(
        np.linspace(1e-6, nc_machine.psi_last.value, num)
    )
    initial = dex.QueueInitialConditions.Boozer(
        t0=np.zeros(num),
        flux0=psi0s,
        theta0=np.zeros(num),
        zeta0=np.zeros(num),
        rho0=np.full(num, 1e-3),
        mu0=np.full(num, 1e-5),
    )
    queue = dex.Queue(initial)

    queue.integrate(nc_machine, (0, 100))
    queue.intersect(nc_machine, intersection="ConstZeta", angle=0, turns=4)
    queue.close(nc_machine)
    queue.classify(nc_machine)
