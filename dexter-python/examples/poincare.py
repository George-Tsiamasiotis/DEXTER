# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "numpy",
# ]
# ///
import matplotlib.pyplot as plt
import numpy as np
import dexter as dx


qfactor = dx.Qfactor("./data.nc", "steffen")
currents = dx.Currents("./data.nc", "steffen")
bfield = dx.Bfield("./data.nc", "bicubic")
perturbation = dx.Perturbation(
    [
        dx.Harmonic("./data.nc", "steffen", m=1, n=8),
        dx.Harmonic("./data.nc", "steffen", m=1, n=9),
    ]
)

num = 20
initials = dx.HeapInitialConditions(
    thetas=np.zeros(num),
    psips=0.0944 * np.ones(num),
    rhos=0.001 * np.ones(num),
    zetas=np.linspace(-np.pi, np.pi, num),
    mus=np.zeros(num),
)

heap = dx.Heap(initials)

params = dx.MappingParameters("ConstTheta", 3.14, 1000)

heap.poincare(qfactor, currents, bfield, perturbation, params)
print(heap)

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
