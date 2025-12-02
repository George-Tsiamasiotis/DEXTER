"""Calculates the ωθ, ωζ and qkinetic of a Heap of Particles."""

import matplotlib.pyplot as plt
import numpy as np
import dexter as dx


geometry = dx.Geometry("./data.nc", "steffen", "bicubic")
qfactor = dx.Qfactor("./data.nc", "steffen")
currents = dx.Currents("./data.nc", "steffen")
bfield = dx.Bfield("./data.nc", "bicubic")
perturbation = dx.Perturbation([])

num = 5000
rho0 = np.log10(1e-5)
rho1 = np.log10(1e-1)
initials = dx.HeapInitialConditions(
    psips=0.7 * geometry.psip_wall * np.ones(num),
    zetas=0 * np.ones(num),
    thetas=np.ones(num),
    rhos=np.logspace(rho0, rho1, num),
    mus=7e-6 * np.zeros(num),
)

heap = dx.Heap(initials)

heap.calculate_frequencies(qfactor, currents, bfield, perturbation)
print(heap)

fig = plt.figure(figsize=(9, 6), layout="constrained", dpi=120)
ax = fig.add_subplot()

qs = np.asarray([heap[i].frequencies.qkinetic for i in range(len(heap))])
energies = np.asarray([heap[i].initial_energy for i in range(len(heap))])

ax.scatter(energies, qs, c="b", s=1)

ax.set_title(r"$q_{kinetic}-E$")
ax.set_xlabel("$E$ $[Normalized]$")
ax.set_ylabel("$q_{kinetic}$")
ax.grid(True)

plt.show()
