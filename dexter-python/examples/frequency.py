"""Integrates a Particle for a single θ-ψp period to calculate ωθ, ωζ and qkinetic."""

import matplotlib.pyplot as plt
import dexter as dx
from dexter.plot import evolution_plot


qfactor = dx.Qfactor("./data.nc", "steffen")
currents = dx.Currents("./data.nc", "steffen")
bfield = dx.Bfield("./data.nc", "bicubic")
perturbation = dx.Perturbation([])

initial = dx.InitialConditions(
    time0=0,
    theta0=2,
    psip0=0.5 * qfactor.psip_wall,
    rho0=0.000008,
    zeta0=0.0,
    mu=0,
)

particle = dx.Particle(initial)

particle.calculate_frequencies(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
)
print(particle)

print(f"Particle's qkinetic = {particle.frequencies.qkinetic}")
print(f"Particle's q(ψ0)    = {qfactor.q(particle.initial_conditions.psip0)}")

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
evolution_plot(fig, particle)
