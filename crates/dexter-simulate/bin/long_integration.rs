use dexter_equilibrium::*;
use dexter_simulate::*;

fn main() {
    // Equilibrium setup
    use InitialFlux::*;
    let qfactor = ParabolicQfactor::new(1.1, 3.9, FluxWall::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-3, 1, 2, 0.0),
        CosHarmonic::new(1e-3, 1, 4, 0.0),
    ]);

    // Particle setup
    let initial = InitialConditions {
        t0: 0.0,
        flux0: Toroidal(0.2),
        theta0: 0.0,
        zeta0: 0.0,
        rho0: 1e-4,
        mu0: 1e-6,
    };
    let mut particle = Particle::new(&initial);

    // Integrate
    let teval = (0.0, 10.0);
    let solver_params = SolverParams {
        method: SteppingMethod::EnergyAdaptiveStep,
        energy_rel_tol: 1e-17,
        energy_abs_tol: 1e-18,
        max_steps: 10_000_000,
        ..Default::default()
    };
    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        teval,
        &solver_params,
    );
    dbg!(&particle);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::TimedOut(..)
    ));
}
