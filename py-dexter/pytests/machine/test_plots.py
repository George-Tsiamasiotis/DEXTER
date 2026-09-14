import dexter as dex


def test_plots_analytical(lar_machine: dex.Machine):
    dex.plot_qfactor(lar_machine, points=50, data=True, show=False)
    dex.plot_current(lar_machine, points=50, data=True, show=False)


def test_plots_nc(nc_machine: dex.Machine):
    dex.plot_qfactor(nc_machine, points=50, data=True, show=False)
    dex.plot_current(nc_machine, points=50, data=True, show=False)
