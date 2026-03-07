#![allow(non_snake_case)]

mod common;

use crate::common::check_integrated_particle_arrays;
// use approx::assert_relative_eq;
use dexter_equilibrium::extract::POLOIDAL_TEST_NETCDF_PATH;
use dexter_equilibrium::*;
use dexter_simulate::*;
use std::path::PathBuf;

#[test]
#[rustfmt::skip]
fn gc_toroidal_intersect_uniQ_larC_larB_noP() {
    use InitialFlux::*;
    let qfactor = UnityQfactor::new();
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::zero();

    let solver_params = SolverParams::default();

    let initial = InitialConditions {t0: 0.0, flux0: Toroidal(0.02), theta0: 3.14, zeta0: 0.0, rho0: 1e-4, mu0: 1e-6};
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 3.14, 10);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    particle.intersect(&qfactor, &current, &bfield, &perturbation, &intersect_params, &solver_params);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Intersected));
    assert_eq!(particle.steps_stored(), 11);
    assert!(particle.steps_taken() > 10);
    assert!(particle.energy_var().unwrap() < 1e-20);
    // assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10); // FIXME:

    check_integrated_particle_arrays(&particle);
}

#[test]
#[rustfmt::skip]
fn gc_poloidal_intersect_ncdQ_ncdC_ncdB_noP() {
    use InitialFlux::*;
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();
    let perturbation = Perturbation::zero();

    let solver_params = SolverParams::default();

    let initial = InitialConditions {t0: 0.0, flux0: Poloidal(0.1), theta0: 3.14, zeta0: 0.0, rho0: 1e-4, mu0: 1e-6};
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 3.14, 3);

    let mut particle = Particle::new(&initial);
    assert!(matches!(particle.integration_status(), IntegrationStatus::Initialized));

    particle.intersect(&qfactor, &current, &bfield, &perturbation, &intersect_params, &solver_params);
    dbg!(&particle);

    assert!(matches!(particle.integration_status(), IntegrationStatus::Intersected));
    assert!(particle.steps_stored() == 4);
    assert!(particle.steps_taken() > 3);
    assert!(particle.energy_var().unwrap() < 1e-20);
    // assert_relative_eq!(particle.initial_energy().unwrap(), particle.final_energy().unwrap(), epsilon = 1e-10); // FIXME:

    check_integrated_particle_arrays(&particle);
}
