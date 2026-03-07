use dexter_equilibrium::*;
use dexter_simulate::*;

fn main() {
    // TODO: numerical_equilibrium_intersect();
    analytical_equilibrium_intersect();
}

fn analytical_equilibrium_intersect() {
    // Equilibrium setup
    use InitialFlux::*;
    let qfactor = ParabolicQfactor::new(1.1, 3.8, FluxWall::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(3e-3, 3, 1, 0.0),
        CosHarmonic::new(2e-3, 7, 2, 0.0),
        CosHarmonic::new(1e-3, 15, 4, 0.0),
    ]);

    // Particle setup
    let initial = InitialConditions {
        t0: 0.0,
        flux0: Toroidal(0.02),
        theta0: 0.0,
        zeta0: 0.0,
        rho0: 8e-3,
        mu0: 1e-6,
    };
    let intersect_params = IntersectParams::new(Intersection::ConstZeta, 0.0, 20);
    let mut particle = Particle::new(&initial);

    // Calculate intersections
    let solver_params = SolverParams {
        method: SteppingMethod::FixedStep(0.4),
        ..Default::default()
    };
    particle.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &solver_params,
    );
    dbg!(&particle);
}
