#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "matplotlib",
#   "pyncare"
# ]
# ///
import matplotlib
import pyncare as pc
import matplotlib.pyplot as plt
from math import pi

matplotlib.use("gtk3agg")

qfactor = pc.Qfactor("./data.nc", "akima")
currents = pc.Currents("./data.nc", "akima")
bfield = pc.Bfield("./data.nc", "bicubic")
perturbation = pc.Perturbation([])

initial = pc.InitialConditions(
    time0=0,
    theta0=1.0,
    psip0=0.5 * qfactor.psip_wall,
    rho0=0.01,
    zeta0=0.1,
    mu=0,
)

particle = pc.Particle(initial)

particle.calculate_omega_theta(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
)

print(particle)
print(f"t-span = [{particle.evolution.time[0]}, {particle.evolution.time[-1]}]")
print(
    f"θ-span = [{particle.evolution.theta[0]}, {particle.evolution.theta[-1]%(2*pi)}]"
)
print(f"ψp-span = [{particle.evolution.psip[0]}, {particle.evolution.psip[-1]}]")
print(f"Pθ-span = [{particle.evolution.ptheta[0]}, {particle.evolution.ptheta[-1]}]")

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
pc.orbit_plot(fig, particle)
