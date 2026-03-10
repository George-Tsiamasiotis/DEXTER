//! Comparison of the different integration `SteppingMethods`.

use dexter_equilibrium::*;
use dexter_simulate::*;

fn main() {
    use InitialFlux::*;
    let qfactor = UnityQfactor::new();
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-2, 1, 1, 0.0),
        CosHarmonic::new(1e-3, 1, 2, 0.0),
        CosHarmonic::new(1e-3, 1, 4, 0.0),
    ]);

    let energy_params = SolverParams {
        method: SteppingMethod::EnergyAdaptiveStep,
        energy_rel_tol: 1e-14,
        energy_abs_tol: 1e-15,
        ..Default::default()
    };
    let error_params = SolverParams {
        method: SteppingMethod::ErrorAdaptiveStep,
        error_rel_tol: 1e-22,
        error_abs_tol: 1e-23,
        ..Default::default()
    };

    let fixed_params = SolverParams {
        method: SteppingMethod::FixedStep(2.0),
        ..Default::default()
    };

    let initial = InitialConditions {
        t0: 0.0,
        flux0: Toroidal(0.3),
        theta0: 0.0,
        zeta0: 0.0,
        rho0: 1e-4,
        mu0: 1e-6,
    };

    let mut energy_particle = Particle::new(&initial);
    let mut error_particle = Particle::new(&initial);
    let mut fixed_particle = Particle::new(&initial);

    let teval = (0.0, 1e5);
    energy_particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        teval,
        &energy_params,
    );
    error_particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        teval,
        &error_params,
    );
    fixed_particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        teval,
        &fixed_params,
    );

    println!("Energy adaptive step:");
    print_results(&energy_particle);
    println!("Error adaptive step:");
    print_results(&error_particle);
    println!("Fixed step:");
    print_results(&fixed_particle);
}

fn print_results(particle: &Particle) {
    println!("\tSteps taken: {}", particle.steps_taken());
    println!("\tEnergy variance: {:?}", particle.energy_var().unwrap());
}
