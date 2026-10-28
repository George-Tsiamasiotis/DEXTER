//! Guiding Center fixed interval integration.

#![allow(non_snake_case)]

mod common;

use crate::common::check_integrated_particle_time_series;
use approx::{RelativeEq, assert_abs_diff_eq, assert_relative_eq};
use dexter_machine::{
    Interpolation1dType::Akima,
    Interpolation2dType::Bicubic,
    extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH, TOROIDAL_TEST_NETCDF_PATH},
    *,
};
use dexter_simulate::*;
use ndarray::{Array1, ArrayView1};
use std::{f64::consts::TAU, path::PathBuf};

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_uniQ_larC_larB_cosP() {
    use MagneticFlux::*;
    let lcfs = Toroidal(0.45);
    let qfactor = UnityQfactor::new(lcfs);
    let current = LarCurrent;
    let bfield = LarBfield;
    let perturbation = Perturbation::new(vec![
            Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
        ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    let solver_params = SolverParams::default();

    let initial = InitialConditions::boozer(0.0, Toroidal(0.3), 0.0, 0.0, 1e-4, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(machine, teval, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_values().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_rsd().unwrap() < 1e-7);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_time_series(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_ncdQ_ncdC_ncdB_ncdP() {
    use MagneticFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(TOROIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, Akima).build().unwrap();
    let current = NcCurrentBuilder::new(&path, Akima).build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
    let perturbation = Perturbation::new(vec![
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 2).with_phase_method(Interpolation).build().unwrap()),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    let solver_params = SolverParams::default();

    let initial = InitialConditions::boozer(0.0, Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(machine, teval, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_values().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_rsd().unwrap() < 1e-7);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_time_series(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_poloidal_integration_ncdQ_ncdC_ncdB_ncdP() {
    use MagneticFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, Akima).build().unwrap();
    let current = NcCurrentBuilder::new(&path, Akima).build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
    let perturbation = Perturbation::new(vec![
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 2).with_phase_method(Interpolation).build().unwrap()),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    let solver_params = SolverParams::default();

    let initial = InitialConditions::boozer(0.0, Poloidal(0.2), 0.0, 0.0, 1e-4, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(machine, teval, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_values().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_rsd().unwrap() < 1e-7);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_time_series(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_gcmotion_check_uniQ_larC_larB_cosP() {
    use MagneticFlux::*;
    let lcfs = Toroidal(0.45);
    let qfactor = UnityQfactor::new(lcfs);
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let solver_params = SolverParams::default();

    let psit0 = Toroidal(0.018365472910927463);
    let initial = InitialConditions::boozer(0.0, psit0, 0.0, 0.0, -0.006634527089072539, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Initialized
    ));

    let teval = (0.0, 3e3);
    particle.integrate(machine, teval, &solver_params);

    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Integrated
    ));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_values().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 20000);
    assert!(particle.energy_rsd().unwrap() < 1e-6);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-8);

    assert_relative_eq!(particle.initial_energy().unwrap(), 1.5189224863170239e-05, epsilon = 1e-20);
    assert_relative_eq!(initial.mu0() * machine.bfield().eval_b(psit0, 0.0, &mut Accelerator2d::new()).unwrap(),
        8.083468084746437e-07,
        epsilon = 1e-15
    );
    assert_relative_eq!(*particle.psit_values().last().unwrap(), 0.02059090128879661, epsilon = 1e-5);
    assert_abs_diff_eq!(*particle.theta_values().last().unwrap(), -15.887142479457443, epsilon = 0.1);
    assert_abs_diff_eq!(*particle.zeta_values().last().unwrap(), -15.80106331080741, epsilon = 0.1);
    assert_relative_eq!(*particle.rho_values().last().unwrap(), -0.004409098711254563, epsilon = 1e-5);
    assert_relative_eq!(*particle.ptheta_values().last().unwrap(), 0.020590901288745446, epsilon = 1e-5);
    assert_relative_eq!(ArrayView1::from(particle.pzeta_values()).var(1.0), 0.0, epsilon = 1e-10);
    assert_relative_eq!(*particle.pzeta_values().first().unwrap(), -0.025, epsilon = 1e-4);
    assert_relative_eq!(*particle.pzeta_values().last().unwrap(), -0.025, epsilon = 1e-4);

    check_integrated_particle_time_series(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_poloidal_equivalence() {
    use MagneticFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, Akima).build().unwrap();
    let current = NcCurrentBuilder::new(&path, Akima).build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
    let perturbation = Perturbation::new(vec![
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 2).with_phase_method(Interpolation).build().unwrap()),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    let solver_params = SolverParams::default();

    let psit0 = Toroidal(0.2);
    let psip0 = machine.qfactor().eval_other(psit0, &mut Accelerator::new()).unwrap();
    let tor_initial = InitialConditions::boozer(0.0, psit0, 0.0, 0.0, 1e-4, 1e-6);
    let pol_initial = InitialConditions::boozer(0.0, psip0, 0.0, 0.0, 1e-4, 1e-6);

    let mut tor_particle = Particle::new(&tor_initial);
    let mut pol_particle = Particle::new(&pol_initial);
    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Initialized));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 2.56e5);
    tor_particle.integrate(machine, teval, &solver_params);
    pol_particle.integrate(machine, teval, &solver_params);
    dbg!(&tor_particle);
    dbg!(&pol_particle);

    assert!(matches!(tor_particle.integration_status(), IntegrationStatus::Integrated));
    assert!(matches!(pol_particle.integration_status(), IntegrationStatus::Integrated));

    // Make sure we dont integrate the same particle twice. The energy should have a small
    // arithmetic error due to the ψ interpolation.
    assert_ne!(tor_particle.initial_energy(), pol_particle.initial_energy());

    assert!(tor_particle.steps_stored() > 10_000);
    assert!(pol_particle.steps_stored() > 10_000);
    assert!(tor_particle.steps_stored() < 10_100);
    assert!(pol_particle.steps_stored() < 10_100);

    check_integrated_particle_time_series(&tor_particle);
    check_integrated_particle_time_series(&pol_particle);

    assert_relative_eq!(tor_particle.initial_energy().unwrap(), pol_particle.initial_energy().unwrap(), epsilon=1e-12);
    assert_relative_eq!(tor_particle.final_energy().unwrap(), pol_particle.final_energy().unwrap(), epsilon=1e-11);

    assert_abs_diff_eq!(
        tor_particle.t_values().last().copied().unwrap(),
        pol_particle.t_values().last().copied().unwrap(),
        epsilon = teval.1 * 1e-3
    );
    assert_abs_diff_eq!(
        tor_particle.psit_values().last().copied().unwrap(),
        pol_particle.psit_values().last().copied().unwrap(),
        epsilon = 1e-3*machine.qfactor().psit_last().value()
    );
    assert_abs_diff_eq!(
        tor_particle.psip_values().last().copied().unwrap(),
        pol_particle.psip_values().last().copied().unwrap(),
        epsilon = 1e-3*machine.qfactor().psip_last().value()
    );
    assert_abs_diff_eq!(
        tor_particle.rho_values().last().copied().unwrap(),
        pol_particle.rho_values().last().copied().unwrap(),
        epsilon = 1e-2*tor_initial.rho0().unwrap() // This deviation might be because rho varies quite fast
    );
    assert_abs_diff_eq!(
        tor_particle.theta_values().last().copied().unwrap(),
        pol_particle.theta_values().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_abs_diff_eq!(
        tor_particle.zeta_values().last().copied().unwrap(),
        pol_particle.zeta_values().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_relative_eq!(
        tor_particle.ptheta_values().last().copied().unwrap(),
        pol_particle.ptheta_values().last().copied().unwrap(),
        epsilon = 1e-5
    );
    assert_relative_eq!(
        tor_particle.pzeta_values().last().copied().unwrap(),
        pol_particle.pzeta_values().last().copied().unwrap(),
        epsilon = 1e-5
    );
}

#[test]
#[rustfmt::skip]
fn gc_mixed_boozer_equivalence() {
    use MagneticFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, Akima).build().unwrap();
    let current = NcCurrentBuilder::new(&path, Akima).build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
    let perturbation = Perturbation::new(vec![
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 1).with_phase_method(Interpolation).build().unwrap()),
        Box::new(NcFluteModeBuilder::new(&path, Akima, 2, 2).with_phase_method(Interpolation).build().unwrap()),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    let solver_params = SolverParams::default();

    let rho0 = 1e-4;
    let psip0 = Poloidal(0.2);
    let teval = (0.0, 3.2e5);

    let boozer_initial = InitialConditions::boozer(0.0, psip0, 0.0, 0.0, rho0, 1e-6);
    let mut boozer_particle = Particle::new(&boozer_initial);
    assert!(matches!(boozer_particle.integration_status(), IntegrationStatus::Initialized));
    boozer_particle.integrate(machine, teval, &solver_params);
    dbg!(&boozer_particle);

    // use the `pzeta` calculated by the boozer particle
    let pzeta = boozer_particle.initial_conditions().pzeta0().unwrap();

    let mixed_initial = InitialConditions::mixed(0.0, psip0, 0.0, 0.0, pzeta, 1e-6);
    let mut mixed_particle = Particle::new(&mixed_initial);
    assert!(matches!(mixed_particle.integration_status(), IntegrationStatus::PartiallyInitialized));
    mixed_particle.integrate(machine, teval, &solver_params);
    dbg!(&mixed_particle);

    assert!(matches!(boozer_particle.integration_status(), IntegrationStatus::Integrated));
    assert!(matches!(mixed_particle.integration_status(), IntegrationStatus::Integrated));

    // Make sure we dont integrate the same particle twice. The energy should have a small
    // arithmetic error due to the ψ interpolation.
    assert_ne!(boozer_particle.initial_energy(), mixed_particle.initial_energy());

    assert!(boozer_particle.steps_stored() > 10_000);
    assert!(mixed_particle.steps_stored() > 10_000);
    assert!(boozer_particle.steps_stored() < 10_100);
    assert!(mixed_particle.steps_stored() < 10_100);

    check_integrated_particle_time_series(&boozer_particle);
    check_integrated_particle_time_series(&mixed_particle);

    assert_relative_eq!(boozer_particle.initial_energy().unwrap(), mixed_particle.initial_energy().unwrap(), epsilon=1e-12);
    assert_relative_eq!(boozer_particle.final_energy().unwrap(), mixed_particle.final_energy().unwrap(), epsilon=1e-11);

    assert_abs_diff_eq!(
        boozer_particle.t_values().last().copied().unwrap(),
        mixed_particle.t_values().last().copied().unwrap(),
        epsilon = teval.1 * 1e-3
    );
    assert_abs_diff_eq!(
        boozer_particle.psit_values().last().copied().unwrap(),
        mixed_particle.psit_values().last().copied().unwrap(),
        epsilon = 1e-3*machine.qfactor().psit_last().value()
    );
    assert_abs_diff_eq!(
        boozer_particle.psip_values().last().copied().unwrap(),
        mixed_particle.psip_values().last().copied().unwrap(),
        epsilon = 1e-3*machine.qfactor().psip_last().value()
    );
    assert_abs_diff_eq!(
        boozer_particle.rho_values().last().copied().unwrap(),
        mixed_particle.rho_values().last().copied().unwrap(),
        epsilon = 1e-3*boozer_initial.rho0().unwrap()
    );
    assert_abs_diff_eq!(
        boozer_particle.theta_values().last().copied().unwrap(),
        mixed_particle.theta_values().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_abs_diff_eq!(
        boozer_particle.zeta_values().last().copied().unwrap(),
        mixed_particle.zeta_values().last().copied().unwrap(),
        epsilon = 1e-3*TAU
    );
    assert_relative_eq!(
        boozer_particle.ptheta_values().last().copied().unwrap(),
        mixed_particle.ptheta_values().last().copied().unwrap(),
        epsilon = 1e-5
    );
    assert_relative_eq!(
        boozer_particle.pzeta_values().last().copied().unwrap(),
        mixed_particle.pzeta_values().last().copied().unwrap(),
        epsilon = 1e-5
    );
}

#[test]
#[rustfmt::skip]
fn gc_const_pzeta_toroidal_integration_uniQ_larC_larB_cosP() {
    use MagneticFlux::*;
    let lcfs = Toroidal(0.45);
    let qfactor = UnityQfactor::new(lcfs);
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let solver_params = SolverParams::default();
    let initial = InitialConditions::mixed(0.0, Toroidal(0.03), 0.0, 0.0, -0.03, 1e-6);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::PartiallyInitialized));

    let teval = (0.0, 6e5);
    particle.integrate(machine, teval, &solver_params);
    let pzeta_array = ArrayView1::from(particle.pzeta_values());
    dbg!(&particle);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_values().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_rsd().unwrap() < 1e-7);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);
    assert!(
        pzeta_array.relative_eq(
            &Array1::from_elem(particle.steps_stored(), initial.pzeta0().unwrap()),
            1e-13,
            1e-15
        )
    );
    assert!(pzeta_array.std(1.0) <= 1e-15);

    check_integrated_particle_time_series(&particle);
}
