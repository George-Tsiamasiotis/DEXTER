"""Integrates a Particle and plots its time evolution."""

import matplotlib.pyplot as plt
import dexter as dx
from dexter.plot import evolution_plot


qfactor = dx.Qfactor("./data.nc", "steffen")
currents = dx.Currents("./data.nc", "steffen")
bfield = dx.Bfield("./data.nc", "bicubic")
harmonics = [dx.Harmonic("./data.nc", "steffen", m=0, n=n) for n in range(1, 24)]
harmonics.pop(6)
perturbation = dx.Perturbation(harmonics)


initial = dx.InitialConditions(
    time0=0,
    theta0=3.14,
    psip0=0.9 * qfactor.psip_wall,
    rho0=0.01,
    zeta0=0.0,
    mu=0,
)

particle = dx.Particle(initial)

particle.integrate(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
    t_eval=(0, 20000),
)
print(particle)

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
evolution_plot(fig, particle)
