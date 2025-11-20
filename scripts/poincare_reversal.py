#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "matplotlib",
#   "numpy",
#   "pyncare"
# ]
# ///

"""30840_103_axi_rev_m48_n0_nub257_nvb1_allmodes"""

import matplotlib.pyplot as plt
import pyncare as pc
import numpy as np

qfactor = pc.Qfactor("./data.nc", "akima")
currents = pc.Currents("./data.nc", "akima")
bfield = pc.Bfield("./data.nc", "bicubic")
perturbation = pc.Perturbation(
    [
        pc.Harmonic("./data.nc", "steffen", m=0, n=1),
    ]
)

psip_lo = 0.084
psip_hi = qfactor.psip_wall

psips = np.concat(
    (
        np.linspace(psip_lo, psip_hi, 50),
        psip_hi * np.ones(50),
    )
)
zetas = np.concat(
    (
        -np.pi * np.ones(50),
        np.linspace(-np.pi, np.pi, 50),
    )
)

init = pc.PoincareInit(
    psips=psips,
    zetas=zetas,
    thetas=np.zeros(len(psips)),
    rhos=0.001 * np.ones(len(psips)),
    mus=np.zeros(len(psips)),
)
params = pc.MappingParameters(section="theta", alpha=np.pi, intersections=1000)

poincare = pc.Poincare(init=init, params=params)
poincare.run(
    qfactor=qfactor,
    currents=currents,
    bfield=bfield,
    perturbation=perturbation,
)
print(poincare)

fig = plt.figure(figsize=(10, 5), layout="constrained", dpi=120)
ax = fig.add_subplot()
pc.poincare_plot(ax, poincare, qfactor, radial=True, color=False)
