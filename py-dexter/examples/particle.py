import dexter as dex

rlast = 0.5
raxis = 2.75
LCFS = dex.MagneticFlux.Toroidal((rlast / raxis) ** 2 / 2)
machine = dex.Machine(
    geometry=dex.LarGeometry(baxis=3, raxis=raxis, rlast=rlast),
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
    mu0=2e-5,
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
    mu0=2e-5,
)

particle = dex.Particle(initial)
particle.integrate(machine, (0, 5e4))
print(particle)

dex.plot_evolution(machine, particle)
particle.plot_poloidal_drift(machine)
assert particle.integration_status == "Integrated"


# =================================

initial = dex.InitialConditions.Mixed(
    t0=0,
    flux0=0.1 * machine.psi_last,
    theta0=0,
    zeta0=0,
    pzeta0=-0.7 * machine.psip_last.value,
    mu0=4e-5,
)

machine.perturbation = dex.Perturbation(
    [
        dex.FluteMode(1e-4, LCFS, 4, 3, 0),
    ]
)

particle = dex.Particle(initial)
particle.intersect(
    machine, intersection="ConstZeta", angle=0, turns=400, max_steps=10_000_000
)
print(particle)

dex.plot_evolution(machine, particle)
particle.plot_poloidal_drift(machine)
assert particle.integration_status == "Intersected"
