"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space in a LAR equilibrium."""

import dexter as dex

geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = dex.LastClosedFluxSurface("Toroidal", geometry.psi_last)
equilibrium = dex.Equilibrium(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

dex.plot_parabolas(equilibrium, mu=6e-5, xlim=(-2.3, 1), ymax=3)
