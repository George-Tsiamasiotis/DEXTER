"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space in a LAR equilibrium."""

import dexter as dex

LCFS = dex.LastClosedFluxSurface(kind="Poloidal", value=0.01)
equilibrium = dex.Equilibrium(
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

coms = dex.COMs(mu=7e-6)

coms.plot_parabolas(equilibrium, ymax=3)
