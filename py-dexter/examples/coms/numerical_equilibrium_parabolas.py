"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space in a numerical equilibrium."""

import dexter as dex

equilibrium = dex.numerical_equilibrium("data.nc", "Steffen", "Bicubic")

dex.plot_parabolas(equilibrium, mu=7e-3, xlim=(-1.2, 0.3), ymax=3, density=3000)
