"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space in a LAR equilibrium."""

import dexter as dex

LCFS = dex.LastClosedFluxSurface(kind="Poloidal", value=0.04)
equilibrium = dex.Equilibrium(
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

dex.plot_parabolas(equilibrium, mu=2e-4, xlim=(-2.1, 1), ymax=3)
