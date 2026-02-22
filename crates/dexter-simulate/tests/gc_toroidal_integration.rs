#![allow(non_snake_case)]

mod common;

use crate::common::check_integrated_particle_arrays;
use approx::{assert_abs_diff_eq, assert_relative_eq};
use dexter_equilibrium::{extract::TEST_NETCDF_PATH, *};
use dexter_simulate::*;
use rsl_interpolation::{Accelerator, Cache};
use std::path::PathBuf;

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_uniQ_larC_larB_cosP() {
    use InitialFlux::*;
    let qfactor = UnityQfactor::new();
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-3, 1, 2, 0.0),
        CosHarmonic::new(1e-3, 1, 4, 0.0),
    ]);

    let params = IntegrationParams::default();

    let initial = InitialConditions {t0: 0.0, flux0: Toroidal(0.3), theta0: 0.0, zeta0: 0.0, rho0: 1e-4, mu0: 1e-6};

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(&qfactor, &current, &bfield, &perturbation, teval, &params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_ncdQ_ncdC_ncdB_ncdP() {
    use InitialFlux::*;
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();
    let perturbation = Perturbation::new(&[
        NcHarmonicBuilder::new(&path, "steffen", 2, 1).with_phase_method(Interpolation).build().unwrap(),
        NcHarmonicBuilder::new(&path, "steffen", 3, 2).with_phase_method(Interpolation).build().unwrap(),
    ]);

    let params = IntegrationParams::default();

    let initial = InitialConditions {t0: 0.0, flux0: Toroidal(0.2), theta0: 0.0, zeta0: 0.0, rho0: 1e-4, mu0: 1e-6};

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e-1);
    particle.integrate(&qfactor, &current, &bfield, &perturbation, teval, &params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_gcmotion_check_uniQ_larC_larB_cosP() {
    use InitialFlux::*;
    let qfactor = UnityQfactor::new();
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::zero();

    let params = IntegrationParams::default();

    let psi0 = 0.018365472910927463;
    let initial = InitialConditions {
        t0: 0.0,
        flux0: Toroidal(psi0),
        theta0: 0.0,
        zeta0: 0.0,
        rho0: -0.006634527089072539,
        mu0: 1e-6,
    };

    let mut particle = Particle::new(&initial);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Initialized
    ));

    let teval = (0.0, 3e3);
    particle.integrate(&qfactor, &current, &bfield, &perturbation, teval, &params);

    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::Integrated
    ));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-18);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-8);

    assert_relative_eq!(particle.initial_energy().unwrap(), 1.5189224863170239e-05, epsilon = 1e-20);
    assert_relative_eq!(initial.mu0 * bfield.b_of_psi(psi0, 0.0, &mut Accelerator::new(), &mut Accelerator::new(), &mut Cache::new()).unwrap(),
        8.083468084746437e-07,
        epsilon = 1e-15
    );
    assert_relative_eq!(*particle.psi_array().last().unwrap(), 0.02059090128879661, epsilon = 1e-5);
    assert_abs_diff_eq!(*particle.theta_array().last().unwrap(), -15.887142479457443, epsilon = 0.1);
    assert_abs_diff_eq!(*particle.zeta_array().last().unwrap(), -15.80106331080741, epsilon = 0.1);
    assert_relative_eq!(*particle.rho_array().last().unwrap(), -0.004409098711254563, epsilon = 1e-5);
    assert_relative_eq!(*particle.ptheta_array().last().unwrap(), 0.020590901288745446, epsilon = 1e-5);
    assert_relative_eq!(particle.pzeta_array().var(1.0), 0.0, epsilon = 1e-10);
    assert_relative_eq!(*particle.pzeta_array().first().unwrap(), -0.025, epsilon = 1e-4);
    assert_relative_eq!(*particle.pzeta_array().last().unwrap(), -0.025, epsilon = 1e-4);

    check_integrated_particle_arrays(&particle);
}
