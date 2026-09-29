//! Comparison of the different integration `SteppingMethods`.

use dexter_machine::*;
use dexter_simulate::*;

fn main() {
    use MagneticFlux::*;
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = UnityQfactor::new(lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-2, lcfs, 1, 1, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    let energy_params = SolverParams {
        method: SteppingMethod::EnergyAdaptiveStep {
            rel_tol: 1e-8,
            abs_tol: 1e-10,
        },
        ..Default::default()
    };
    let error_params = SolverParams {
        method: SteppingMethod::ErrorAdaptiveStep {
            rel_tol: 1e-16,
            abs_tol: 1e-19,
        },
        ..Default::default()
    };

    let fixed_params = SolverParams {
        method: SteppingMethod::FixedStep(60.0),
        ..Default::default()
    };

    let initial = InitialConditions::boozer(0.0, Toroidal(0.3), 0.0, 0.0, 1e-4, 1e-6);

    let mut energy_particle = Particle::new(&initial);
    let mut error_particle = Particle::new(&initial);
    let mut fixed_particle = Particle::new(&initial);

    let teval = (0.0, 1e7);
    energy_particle.integrate(machine, teval, &energy_params);
    error_particle.integrate(machine, teval, &error_params);
    fixed_particle.integrate(machine, teval, &fixed_params);

    use IntegrationStatus::*;
    assert_eq!(energy_particle.integration_status(), Integrated);
    assert_eq!(error_particle.integration_status(), Integrated);
    assert_eq!(fixed_particle.integration_status(), Integrated);

    println!("Energy adaptive step:");
    print_results(&energy_particle);
    println!("Error adaptive step:");
    print_results(&error_particle);
    println!("Fixed step:");
    print_results(&fixed_particle);
}

fn print_results(particle: &Particle) {
    println!("\tSteps taken: {}", particle.steps_taken());
    println!("\tEnergy variance: {:?}", particle.energy_var());
    println!("\tDuration: {:?}", particle.duration());
}
