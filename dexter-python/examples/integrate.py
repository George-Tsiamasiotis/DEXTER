# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
# ]
# ///
import matplotlib.pyplot as plt
import dexter as dx
from dexter.plot import evolution_plot


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
    theta0=3.14,
    psip0=0.32 * qfactor.psip_wall,
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
    t_eval=(0, 15000),
)
print(particle)

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
evolution_plot(fig, particle)
