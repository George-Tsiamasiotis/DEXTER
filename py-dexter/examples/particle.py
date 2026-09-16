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

flux0 = dex.MagneticFlux.Toroidal(0.5 * machine.psi_last.value)
initial = dex.InitialConditions.Mixed(
    t0=0,
    flux0=flux0,
    theta0=0,
    zeta0=0,
    pzeta0=-0.4 * machine.psip_last.value,
    mu0=2e-5,
)

particle = dex.Particle(initial)
particle.classify(machine)

particle.close(machine, 4)
assert particle.integration_status == "ClosedPeriods(4)"
print(particle)
dex.plot_evolution(machine, particle)


machine.perturbation = dex.Perturbation(
    [
        dex.FluteMode(1e-5, LCFS, 3, 2, 0),
        dex.FluteMode(1e-5, LCFS, 5, 3, 0),
    ]
)

particle.integrate(machine, (0, 5e4))
assert particle.integration_status == "Integrated"
print(particle)
dex.plot_evolution(machine, particle)

particle.intersect(machine, intersection="ConstZeta", angle=0, turns=100)
assert particle.integration_status == "Intersected"
print(particle)
particle.plot_evolution(machine)
