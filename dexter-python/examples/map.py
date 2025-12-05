"""Calculates a Particle's intersections with θ=const plane."""

import matplotlib.pyplot as plt
import dexter as dx
from dexter.plot import evolution_plot


geometry = dx.Geometry("./data.nc", "steffen", "bicubic")
qfactor = dx.Qfactor("./data.nc", "steffen")
currents = dx.Currents("./data.nc", "steffen")
bfield = dx.Bfield("./data.nc", "bicubic")
perturbation = dx.Perturbation(
    [
        dx.Harmonic("./data.nc", "steffen", m=1, n=8),
        dx.Harmonic("./data.nc", "steffen", m=1, n=9),
    ]
)

initial = dx.InitialConditions(
    time0=0,
    theta0=0,
    psip0=0.4 * geometry.psip_wall,
    rho0=0.01,
    zeta0=0.0,
    mu=0,
)

particle = dx.Particle(initial)

params = dx.MappingParameters("ConstTheta", 3.14, 100)

particle.map(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
    params=params,
)
print(particle)

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
evolution_plot(fig, particle)
