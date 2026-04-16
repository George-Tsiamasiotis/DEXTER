"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space in a numerical equilibrium."""

import dexter as dex

equilibrium = dex.numerical_equilibrium("data.nc", "Steffen", "Bicubic")

coms = dex.COMs(mu=7e-3)

coms.plot_parabolas(equilibrium, xlim=(-1.2, 0.3), ymax=3, density=3000)
