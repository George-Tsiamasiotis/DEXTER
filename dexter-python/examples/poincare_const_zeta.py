"""Calculates the Poincare map of a Heap of Particles."""

import matplotlib.pyplot as plt
import numpy as np
import dexter as dx
from dexter.plot import poincare_plot


geometry = dx.Geometry("./data.nc", "steffen", "bicubic")
qfactor = dx.Qfactor("./data.nc", "steffen")
currents = dx.Currents("./data.nc", "steffen")
bfield = dx.Bfield("./data.nc", "bicubic")
perturbation = dx.Perturbation([dx.Harmonic("./data.nc", "steffen", m=1, n=7)])

rs = np.linspace(0.13, geometry.r_wall, 40)
rs = np.concat((rs, [0.1843]))  # Add the separatrix manually
rs = np.concat((rs, [0.17976, 0.186]))  # Add some extra trapped orbits
psips = np.asarray([geometry.psip(r) for r in rs])
zetas = np.pi * np.ones(len(psips))

initials = dx.HeapInitialConditions(
    psips=psips,
    zetas=zetas,
    thetas=np.zeros(len(psips)),
    rhos=0.0001 * np.ones(len(psips)),
    mus=np.zeros(len(psips)),
)

heap = dx.Heap(initials)

params = dx.MappingParameters("ConstZeta", np.pi, 2000)

heap.poincare(qfactor, currents, bfield, perturbation, params)
print(heap)

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
ax = fig.add_subplot()
poincare_plot(
    ax,
    heap,
    params,
    geometry,
    radial=True,
    lab=True,
    color=True,
    wall=True,
)
