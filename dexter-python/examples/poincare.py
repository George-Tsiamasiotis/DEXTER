"""Calculates the Poincare map of a Heap of Particles."""

import matplotlib.pyplot as plt
import numpy as np
import dexter as dx
from dexter.plot import poincare_plot


qfactor = dx.Qfactor("./data.nc", "steffen")
currents = dx.Currents("./data.nc", "steffen")
bfield = dx.Bfield("./data.nc", "bicubic")
harmonics = [dx.Harmonic("./data.nc", "steffen", m=0, n=n) for n in range(1, 24)]
harmonics.pop(6)
perturbation = dx.Perturbation(harmonics)

psip = 0.09378

num1 = 31
num2 = 21
psips = np.concat(
    (
        psip * np.ones(num1),
        psip * np.ones(num2),
    )
)
zetas = np.concat(
    (
        np.linspace(-np.pi, np.pi, num1),
        np.linspace(-np.pi / 2, np.pi / 2, num2),
    ),
)

initials = dx.HeapInitialConditions(
    psips=psips,
    zetas=zetas,
    thetas=np.zeros(len(psips)),
    rhos=0.0001 * np.ones(len(psips)),
    mus=np.zeros(len(psips)),
)

heap = dx.Heap(initials)

params = dx.MappingParameters("ConstTheta", 3.14, 1000)

heap.poincare(qfactor, currents, bfield, perturbation, params)
print(heap)

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
ax = fig.add_subplot()
poincare_plot(ax, heap, params, qfactor, radial=False, color=False, wall=True)
