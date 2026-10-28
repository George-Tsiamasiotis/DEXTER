//! `θ-ψ` orbit closing and qkinetic calculation checks.

#![allow(non_snake_case)]

mod common;

use approx::assert_relative_eq as ar;
use dexter_machine::*;
use dexter_simulate::*;
use ndarray::Array1;
use std::f64::consts::TAU;

/// Simplest case
/// Use a very low `ρ` and `μ=0` to track a magnetic field line, which should have
/// `qkinetic = qmagnetic = 1`
///
/// Note that this orbit has essentially Δψ=0. The high order hamiltonian terms are essentially
/// zero, so the `first` and `final` values should be almost equal, up to the solver's tolerance.
///
/// Results calculated manually with the `integrate()` routine.
#[test]
#[rustfmt::skip]
fn field_line_single_period_uniQ() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.5));
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.3), 1.0, 0.0, 1e-12, 0.0);
    let mut particle = Particle::new(&initial);
    particle.close(machine, 1, &SolverParams::default());
    assert!(matches!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1)));

    let expected_dt = 17084889680236.904;
    let expected_omega_theta = TAU / expected_dt;
    let expected_omega_zeta = TAU / expected_dt;

    let eps=1e-10;
    ar!(particle.t_values().last().copied().unwrap(), expected_dt, epsilon=eps);
    ar!(particle.psit_values().first().copied().unwrap(), particle.psit_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.psip_values().first().copied().unwrap(), particle.psip_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.theta_values().first().copied().unwrap(), particle.theta_values().last().copied().unwrap() - TAU, epsilon=eps);
    ar!(particle.zeta_values().first().copied().unwrap(), particle.zeta_values().last().copied().unwrap() - TAU, epsilon=eps);
    ar!(particle.rho_values().first().copied().unwrap(), particle.rho_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.pzeta_values().first().copied().unwrap(), particle.pzeta_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.ptheta_values().first().copied().unwrap(), particle.ptheta_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.frequencies().omega_theta.unwrap(), expected_omega_theta, epsilon=eps);
    ar!(particle.frequencies().omega_zeta.unwrap(), expected_omega_zeta, epsilon=eps);
    ar!(particle.frequencies().qkinetic.unwrap(), 1.0, epsilon=eps);
}

#[test]
#[rustfmt::skip]
fn trapped_particle_single_period_uniQ() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.1));
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.02), 1.0, 0.0, 1e-6, 1e-6);
    let mut particle = Particle::new(&initial);
    particle.close(machine, 1, &SolverParams::default());
    particle.classify(machine);
    assert!(matches!(particle.integration_status(), IntegrationStatus::ClosedPeriods(1)));
    assert_eq!(particle.orbit_type(), OrbitType::TrappedConfined);

    let expected_dt = 17703.530171220227;
    let expected_dzeta = 0.07681124462891854; // from manual integration
    let expected_omega_theta = TAU / expected_dt;
    let expected_omega_zeta = expected_dzeta / expected_dt;

    let eps = 1e-6;
    ar!(particle.t_values().last().copied().unwrap(), expected_dt, epsilon=eps);
    ar!(particle.psit_values().first().copied().unwrap(), particle.psit_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.theta_values().first().copied().unwrap(), particle.theta_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.ptheta_values().first().copied().unwrap(), particle.ptheta_values().last().copied().unwrap(), epsilon=eps);
    ar!(particle.frequencies().omega_theta.unwrap(), expected_omega_theta, epsilon=eps);
    ar!(particle.frequencies().omega_zeta.unwrap(), expected_omega_zeta, epsilon=eps);
    ar!(particle.frequencies().qkinetic.unwrap(), 0.012222584887591207, epsilon=eps); // derived
}

#[test]
#[rustfmt::skip]
fn multiple_periods_frequencies_calculation() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.1));
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    // Taken from `trapped_particle_single_period_uniQ`
    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.02), 1.0, 0.0, 1e-6, 1e-6);

    let mut omega_thetas = Vec::new();
    let mut omega_zetas = Vec::new();
    let mut qkinetics = Vec::new();

    let PERIODS = 5;

    for periods in 1..=PERIODS{
        let mut particle = Particle::new(&initial);
        particle.close(machine, periods, &SolverParams::default());
        particle.classify(machine);
        match particle.integration_status(){
            IntegrationStatus::ClosedPeriods(closed_periods) => assert_eq!(periods, *closed_periods),
            _ => panic!("wrong integration status")
        }
        assert_eq!(particle.orbit_type(), OrbitType::TrappedConfined);
        omega_thetas.push(particle.omega_theta().unwrap());
        omega_zetas.push(particle.omega_zeta().unwrap());
        qkinetics.push(particle.qkinetic().unwrap());
    }


    // Taken from `trapped_particle_single_period_uniQ`
    let expected_dt = 17703.530171220227;
    let expected_dzeta = 0.07681124462891854; // from manual integration
    let expected_omega_theta = TAU / expected_dt;
    let expected_omega_zeta = expected_dzeta / expected_dt;
    let expected_qkinetic = expected_omega_zeta / expected_omega_theta;

    let omega_thetas = Array1::from_vec(omega_thetas);
    let omega_zetas = Array1::from_vec(omega_zetas);
    let qkinetics = Array1::from_vec(qkinetics);
    let expected_omega_thetas = Array1::from_elem(PERIODS, expected_omega_theta);
    let expected_omega_zetas = Array1::from_elem(PERIODS, expected_omega_zeta);
    let expected_qkinetics = Array1::from_elem(PERIODS, expected_qkinetic);

    let eps = 1e-5;
    ar!(omega_thetas, expected_omega_thetas, epsilon=eps);
    ar!(omega_zetas, expected_omega_zetas, epsilon=eps);
    ar!(qkinetics, expected_qkinetics, epsilon=eps);
}
