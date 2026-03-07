#![allow(non_snake_case)]

mod common;
use approx::assert_relative_eq;
use dexter_equilibrium::{extract::TEST_NETCDF_PATH, *};
use dexter_simulate::*;
use std::path::PathBuf;

use crate::common::check_integrated_particle_arrays;

#[test]
#[rustfmt::skip]
fn gc_toroidal_integration_uniQ_larC_larB_cosP() {
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

    let solver_params = SolverParams::default();

    let initial = InitialConditions {t0: 0.0, flux0: Poloidal(0.2), theta0: 0.0, zeta0: 0.0, rho0: 1e-4, mu0: 1e-6};

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    let teval = (0.0, 1e5);
    particle.integrate(&qfactor, &current, &bfield, &perturbation, teval, &solver_params);
    dbg!(&particle);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Integrated));
    assert_eq!(particle.steps_stored(), particle.steps_taken());
    assert!(particle.t_array().last().copied().unwrap() >= teval.1);
    assert!(particle.steps_stored() > 1000);
    assert!(particle.steps_stored() < 10000);
    assert!(particle.energy_var().unwrap() < 1e-20);
    assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10);

    check_integrated_particle_arrays(&particle);
}
