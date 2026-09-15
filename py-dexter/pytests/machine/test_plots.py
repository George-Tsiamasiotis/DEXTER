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
