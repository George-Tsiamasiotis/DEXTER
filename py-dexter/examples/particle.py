import dexter as dex

geometry = dex.LarGeometry(baxis=3, raxis=1.75, rlast=0.5)
LCFS = geometry.psi_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

initial = dex.InitialConditions.Mixed(
    t0=0,
    flux0=0.6 * machine.psi_last,
    theta0=0,
    zeta0=0,
    pzeta0=-0.4 * machine.psip_last.value,
    mu0=8e-5,
)

particle = dex.Particle(initial)
particle.classify(machine)
particle.close(machine, 4)
print(particle)

dex.plot_evolution(machine, particle)
particle.plot_poloidal_drift(machine)
assert particle.orbit_type == "TrappedConfined"
assert particle.integration_status == "ClosedPeriods(4)"

# =================================

machine.perturbation = dex.Perturbation(
    [
        dex.FluteMode(1e-4, LCFS, 1, 3, 0),
        dex.FluteMode(1e-4, LCFS, 3, 2, 0),
        dex.FluteMode(1e-4, LCFS, 5, 2, 0),
        dex.FluteMode(1e-4, LCFS, 7, 3, 0),
    ]
)

# =================================

initial = dex.InitialConditions.Mixed(
    t0=0,
    flux0=0.15 * machine.psi_last,
    theta0=0,
    zeta0=0,
    pzeta0=-0.5 * machine.psip_last.value,
    mu0=8e-5,
)

particle = dex.Particle(initial)
particle.integrate(machine, (0, 2e4))
print(particle)

dex.plot_evolution(machine, particle)
particle.plot_poloidal_drift(machine)
assert particle.integration_status == "Integrated"


# =================================

initial = dex.InitialConditions.Mixed(
    t0=0,
    flux0=0.5 * machine.psi_last,
    theta0=0,
    zeta0=0,
    pzeta0=0.1 * machine.psip_last.value,
    mu0=4e-5,
)

machine.perturbation = dex.Perturbation(
    [
        dex.FluteMode(6e-5, LCFS, 4, 3, 0),
    ]
)

particle = dex.Particle(initial)
intersect_params = dex.IntersectParams("ConstZeta", 0, 1000)
solver_params = dex.SolverParams(
    method=dex.SteppingMethod.EnergyAdaptiveStep(1e-7, 1e-8),
    max_steps=10_000_000,
)
particle.intersect(machine, intersect_params, solver_params)
print(particle)

dex.plot_evolution(machine, particle)
particle.plot_poloidal_drift(machine)
assert particle.integration_status == "Intersected"
