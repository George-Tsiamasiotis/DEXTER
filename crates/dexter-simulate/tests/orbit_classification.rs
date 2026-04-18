//! Test particle orbit classification.
//!
//! Theses tests are accompanied by Python scripts to visualize the `(E, Pζ)` space and the orbits.

use dexter_equilibrium::*;
use dexter_simulate::*;

fn create_equilibrium() -> (ParabolicQfactor, LarCurrent, LarBfield) {
    let lcfs = LastClosedFluxSurface::Toroidal(0.03);
    (
        ParabolicQfactor::new(1.1, 3.9, lcfs),
        LarCurrent::new(),
        LarBfield::new(),
    )
}

/// CuPassing-Confined inside the left wall parabola.
#[test]
#[rustfmt::skip]
fn orbit1() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.001);
    let pzeta0 = - 0.8 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CoPassingConfined inside the magnetic axis parabola.
#[test]
#[rustfmt::skip]
fn orbit2() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.015);
    let pzeta0 = - 0.1 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::CoPassingConfined);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CuPassing-Lost inside the right wall parabola, outside the left wall parabola and left from the
/// `Pζ/ψp_last = -1` line.
#[test]
#[rustfmt::skip]
fn orbit3() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.02);
    let pzeta0 = - 1.5 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::CuPassingLost);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// Trapped-Lost inside the right wall parabola.
#[test]
#[rustfmt::skip]
fn orbit4() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.01);
    let pzeta0 = - 0.8 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::TrappedLost);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::Escaped);
}

/// Trapped-Confined outside of all parabolas.
#[test]
#[rustfmt::skip]
fn orbit5() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.018);
    let pzeta0 = - 0.5 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::TrappedConfined);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Potato
#[test]
#[rustfmt::skip]
fn orbit6() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.0045);
    let pzeta0 = - 0.0448 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::Potato);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// Stagnated
#[test]
#[rustfmt::skip]
fn orbit7() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.0014);
    let pzeta0 = - 0.0 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::Stagnated);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}

/// CuPassingConfined outside of all parabolas and above the tp boundary.
#[test]
#[rustfmt::skip]
fn orbit8() {
    let (qfactor, current, bfield) = create_equilibrium();

    let psi0 = InitialFlux::Toroidal(0.0014);
    let pzeta0 = - 0.4 * qfactor.psip_last();
    let initial = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, pzeta0, 6e-5);

    let mut particle = Particle::new(&initial);
    assert_eq!(particle.orbit_type(), OrbitType::Undefined);

    particle.classify(&qfactor, &current, &bfield);

    assert_eq!(particle.orbit_type(), OrbitType::CuPassingConfined);

    particle.close(&qfactor, &current, &bfield, &Perturbation::zero(), 1,&SolverParams::default());
    assert_eq!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1));
}
