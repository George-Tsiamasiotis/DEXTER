import numpy as np
import dexter as dex


def test_queue_from_particles(lar_machine: dex.Machine):
    psi0 = dex.MagneticFlux.Toroidal(0.01)
    i = dex.InitialConditions.Boozer(0, psi0, 0, 0, 1e-4, 1e-6)
    p1 = dex.Particle(i)
    p2 = dex.Particle(i)
    p1.close(lar_machine)
    queue = dex.Queue.FromParticles([p1, p2])
    particles = queue.particles
    assert len(particles) == 2
    assert isinstance(particles[0], dex.Particle)
    assert isinstance(particles[1], dex.Particle)
    assert "Closed" in particles[0].integration_status
    assert particles[1].integration_status == "Initialized"


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
