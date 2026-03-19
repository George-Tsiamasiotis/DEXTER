//! Tests corner cases for unwanted particle integration results

mod common;
use crate::common::*;
use dexter_simulate::*;

#[test]
fn time_out() {
    let (qfactor, current, bfield, perturbation) = lar_equilibrium();

    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.1), 0.0, 0.0, 1e-4, 1e-6);
    let solver_params = SolverParams {
        max_steps: 10,
        ..Default::default()
    };

    // Particle integration
    let mut particle = Particle::new(&initial);
    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        (0.0, 1e5),
        &solver_params,
    );
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::TimedOut(..)
    ));

    // Particle intersection
    let mut particle = Particle::new(&initial);
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 1.0, 10);
    particle.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &solver_params,
    );
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::TimedOut(..)
    ));
}

#[test]
fn out_of_bounds_initialization() {
    use IntegrationStatus::OutOfBoundsInitialization;
    let (qfactor, current, bfield, perturbation) = lar_equilibrium();

    let initial = InitialConditions::boozer(0.0, InitialFlux::Poloidal(1e10), 0.0, 0.0, 1e-4, 1e-6);

    // Particle integration
    let mut particle = Particle::new(&initial);
    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        (0.0, 1e5),
        &SolverParams::default(),
    );
    assert!(matches!(
        particle.integration_status(),
        OutOfBoundsInitialization
    ));

    // Particle intersection
    let mut particle = Particle::new(&initial);
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 0.0, 10);
    particle.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &SolverParams::default(),
    );
    assert!(matches!(
        particle.integration_status(),
        OutOfBoundsInitialization
    ));
}

#[test]
fn intersected_time_out() {
    let (qfactor, current, bfield, perturbation) = lar_equilibrium();

    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.1), 0.0, 0.0, 1e-4, 0.0);
    let solver_params = SolverParams {
        max_steps: 10000,
        ..Default::default()
    };

    // Particle intersection
    let mut particle = Particle::new(&initial);
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 1.0, 10000);
    particle.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &solver_params,
    );
    dbg!(&particle);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::IntersectedTimedOut
    ));
}
