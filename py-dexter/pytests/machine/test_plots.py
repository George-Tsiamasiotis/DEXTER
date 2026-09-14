import dexter as dex


def test_plots_analytical(lar_machine: dex.Machine):
    dex.plot_qfactor(lar_machine, points=50, data=True, show=False)
    dex.plot_current(lar_machine, points=50, data=True, show=False)
    dex.plot_bfield(lar_machine, levels=30, show=False)


def test_plots_nc(nc_machine: dex.Machine):
    dex.plot_qfactor(nc_machine, points=50, data=True, show=False)
    dex.plot_current(nc_machine, points=50, data=True, show=False)
    dex.plot_bfield(nc_machine, levels=30, show=False)


def test_plot_mode(nc_flute_mode: dex.NcFluteMode):
    LCFS = dex.MagneticFlux.Toroidal(0.03)
    flute_mode = dex.FluteMode(1e-5, LCFS, 5, 2, 0)
    dex.plot_mode(flute_mode, points=50, show=False)
    LCFS = dex.MagneticFlux.Poloidal(0.03)
    flute_mode = dex.FluteMode(1e-5, LCFS, 5, 2, 0)
    dex.plot_mode(flute_mode, points=50, show=False)
    dex.plot_mode(flute_mode, points=50, data=True, show=False)


def test_plot_particle(lar_machine: dex.Machine):
    psi_last = lar_machine.psi_last.value
    flux0 = dex.MagneticFlux.Toroidal(psi_last * 0.5)
    initial = dex.InitialConditions.Boozer(0, flux0, 1, 2, 1e-6, 1e-7)
    particle = dex.Particle(initial)
    particle.close(lar_machine)
    assert particle.integration_status == "ClosedPeriods(1)"
    assert 100 < particle.steps_taken < 10_000
    dex.plot_evolution(lar_machine, particle, show=False)
    particle.plot_evolution(lar_machine, downsample=True, show=False)
