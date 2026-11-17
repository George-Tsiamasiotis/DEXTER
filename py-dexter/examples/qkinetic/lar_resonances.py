import numpy as np
import dexter as dex
import matplotlib.pyplot as plt
from alpha_shapes import Alpha_Shaper  # type: ignore
from matplotlib.axes import Axes

from dexter.utils import FractionFormatter, prepare_contour_levels

RNG = np.random.default_rng(42)

# ========================= Machine setup

geometry = dex.LarGeometry(3, 1.75, 0.5)
LCFS = geometry.psit_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

# ========================= Initial conditions and Queue setup

mu = 5e-5
pzeta_lim = (-2.03, 0.71)

num = 200000

# Avoid going too close to the axis
psit0_values = RNG.uniform(low=1e-3, high=1, size=num) * LCFS.value
psit0s = dex.MagneticFluxArray.Toroidal(psit0_values)
theta0s = RNG.random(num) * 2 * np.pi

pzeta0s = np.linspace(*pzeta_lim, num) * machine.psip_last.value

initial_conditions = dex.QueueInitialConditions.Mixed(
    t0=np.zeros(num),
    flux0=psit0s,
    theta0=theta0s,
    zeta0=np.zeros(num),
    pzeta0=pzeta0s,
    mu0=np.full(num, mu),
)

queue = dex.Queue(initial_conditions)
queue.finalize_initial_conditions(machine)
queue.retain_energy(0, 3.05 * mu)

# ========================= Bifurcation analysis

plane = dex.EnergyPzetaPlane(machine, mu)
bifurcation_analysis_pzetas = np.linspace(*pzeta_lim, 2000) * machine.psip_last.value
plane.bifurcation_analysis(bifurcation_analysis_pzetas)
plane.print_durations()

# ========================= Particle classification and filtering

queue.classify(machine)
plane.add_particles(queue)
plane.show(colors="Confinement")
plane.clear_particles()

queue.discard_lost()
plane.add_particles(queue)
plane.show(colors="Confinement")
plane.clear_particles()


co_queue = dex.Queue.FromParticles(
    [
        particle
        for particle in queue.particles()
        if particle.initial_conditions.theta_dot0 > 0
    ]
)
plane.add_particles(co_queue)
plane.show(colors="Confinement")
plane.clear_particles()

cu_queue = dex.Queue.FromParticles(
    [
        particle
        for particle in queue.particles()
        if particle.initial_conditions.theta_dot0 < 0
    ]
)
cu_queue.discard_orbit_types(["TrappedConfined", "Potato", "Stagnated"])
plane.add_particles(cu_queue)
plane.show(colors="Confinement")
plane.clear_particles()

# ========================= Computation of qkinetic

solver_params = dex.SolverParams(
    method=dex.SteppingMethod.EnergyAdaptiveStep(rel_tol=1e-7, abs_tol=1e-8),
    max_steps=50_000,
)

co_queue.close(machine, solver_params=solver_params)
co_queue.retain_closed()

cu_queue.close(machine, solver_params=solver_params)
cu_queue.retain_closed()

# ========================= Contour setup

levels = [
    2,
    1 / 6,
    1 / 3,
    4 / 3,
    8 / 7,
    5 / 3,
]
levels = np.asarray(levels)
levels = np.concat((levels, -levels))
levels = prepare_contour_levels(levels)


def qkinetic_contour(queue: dex.Queue, ax: Axes):
    """Plots a Tricontour plot on `ax` of the particles contained in `queue`."""

    # Extract data
    points = np.asarray(
        (
            queue.pzeta_array / machine.psip_last.value,
            queue.energy_array / mu,
        )
    ).T
    qkinetic_array = queue.qkinetic_array

    # Avoids contour artifacts on the convex Trapped-Lost region
    triang = Alpha_Shaper(points, normalize=False)
    triang.set_mask_at_alpha(10)

    # Plot and add labels for each found qkinetic resonance
    tri = ax.tricontour(
        triang,
        qkinetic_array,
        levels=levels,
        cmap="tab20b",
        linewidths=1.5,
    )
    ax.clabel(
        tri,
        fmt=FractionFormatter(),
        fontsize=10,
        zorder=100,
    )


_, ax = plane.show(show=False)
qkinetic_contour(co_queue, ax)
plt.show()

_, ax = plane.show(show=False)
qkinetic_contour(cu_queue, ax)
plt.show()

_, ax = plane.show(show=False)
qkinetic_contour(co_queue, ax)
qkinetic_contour(cu_queue, ax)
plt.show()
