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
        pc.Harmonic("./data.nc", "akima", m=0, n=1),
    ]
)

num = 200
psip_lo = 0.093
psip_mi = 0.094
psip_hi = qfactor.psip_wall

"Phase Average"
init = pc.PoincareInit(
    thetas=np.zeros(num),
    psips=np.concat(
        (
            psip_lo * np.ones(50),
            psip_mi * np.ones(50),
            psip_hi * np.ones(100),
        )
    ),
    rhos=0.001 * np.ones(num),
    zetas=np.concat(
        (
            np.linspace(-np.pi, np.pi, 50),
            np.linspace(-np.pi, np.pi, 50),
            np.linspace(-np.pi, np.pi, 100),
        )
    ),
    mus=np.zeros(num),
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
pc.poincare_plot(ax, poincare, wall=qfactor.psip_wall)
