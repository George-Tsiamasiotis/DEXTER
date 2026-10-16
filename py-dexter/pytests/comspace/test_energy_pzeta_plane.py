import numpy as np
import dexter as dex


def test_energy_pzeta_plane(lar_machine: dex.Machine):
    plane = dex.EnergyPzetaPlane(lar_machine, 1e-5)
    plane.show(show=False)

    flux0 = lar_machine.psi_last / 2
    initial = dex.InitialConditions.Mixed(0, flux0, 3, 2, 1e-6, 1e-7)
    particle = dex.Particle(initial)
    particle.classify(lar_machine)  # calculate energy
    plane.add_particles(particle)
    plane.show(show=False)

    plane.clear_particles()
    plane.show(show=False)


def test_energy_pzeta_plane_critical_points(lar_machine: dex.Machine):
    plane = dex.EnergyPzetaPlane(lar_machine, 1e-5)
    pzeta_array = np.linspace(-1, 0, 10) * lar_machine.psip_last.value
    plane.calculate_critical_points(pzeta_array)

    points = plane.critical_points
    assert len(points) > 0
    plane.show()

    particle = dex.Particle.FromCriticalPoint(points[-1])
    particle.close(lar_machine)
    assert "ClosedPeriods" in particle.integration_status
