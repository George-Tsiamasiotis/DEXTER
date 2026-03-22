//! Test `Queue::intersect` routine.

#![allow(non_snake_case)]

use dexter_equilibrium::extract::POLOIDAL_TEST_NETCDF_PATH;
use dexter_equilibrium::*;
use dexter_simulate::*;
use ndarray::Array1;
use std::path::PathBuf;

#[test]
fn queue_intersect_const_theta_parQ_larC_larB_cosP() -> Result<()> {
    // Equilibrium setup
    let qfactor = ParabolicQfactor::new(1.1, 1.9, FluxWall::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-3, 1, 2, 0.0),
        CosHarmonic::new(1e-3, 1, 4, 0.0),
    ]);

    // Initial Conditions setup
    let particle_count = 10;
    let psis = qfactor.psi_wall() * Array1::linspace(0.0, 1.0, particle_count);
    let psis = toroidal_fluxes(&psis.to_vec());
    let initial_conditions = QueueInitialConditions::boozer(
        &vec![0.0; particle_count],
        &psis.to_vec(),
        &vec![0.0; particle_count],
        &vec![0.0; particle_count],
        &vec![1e-6; particle_count],
        &vec![2e-6; particle_count],
    )?;

    let mut queue = Queue::new(&initial_conditions);

    assert_eq!(queue.particle_count(), particle_count);
    assert_eq!(queue.routine(), Routine::None);
    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::Initialized)
            .count(),
        particle_count
    );

    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 0.0, 10);
    queue.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &SolverParams::default(),
    );

    // All but the first at `ψ=0` should be intersected.
    assert!(queue[0].integration_status() == IntegrationStatus::OutOfBoundsInitialization);
    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::Intersected)
            .count(),
        particle_count - 1
    );

    println!("{queue:#?}");
    Ok(())
}

#[test]
fn queue_poloidal_intersect_const_zeta_ncdQ_ncdC_ncdB_ncdP() -> Result<()> {
    // Equilibrium setup
    use PhaseMethod::Interpolation;
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();
    let perturbation = Perturbation::new(&[
        NcHarmonicBuilder::new(&path, "steffen", 2, 1)
            .with_phase_method(Interpolation)
            .build()
            .unwrap(),
        NcHarmonicBuilder::new(&path, "steffen", 3, 2)
            .with_phase_method(Interpolation)
            .build()
            .unwrap(),
    ]);

    // Initial Conditions setup
    let particle_count = 10;
    let psips = qfactor.psip_wall().unwrap() * Array1::linspace(0.0, 0.97, particle_count);
    let psips = poloidal_fluxes(&psips.to_vec());
    let initial_conditions = QueueInitialConditions::boozer(
        &vec![0.0; particle_count],
        &psips.to_vec(),
        &vec![0.0; particle_count],
        &vec![0.0; particle_count],
        &vec![1e-6; particle_count],
        &vec![2e-6; particle_count],
    )?;

    let mut queue = Queue::new(&initial_conditions);

    assert_eq!(queue.particle_count(), particle_count);
    assert_eq!(queue.routine(), Routine::None);
    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::Initialized)
            .count(),
        particle_count
    );

    let intersect_params = IntersectParams::new(Intersection::ConstZeta, 0.0, 10);
    queue.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &SolverParams::default(),
    );

    // All but the first at `ψ=0` should be intersected.
    assert!(queue[0].integration_status() == IntegrationStatus::OutOfBoundsInitialization);
    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::Intersected)
            .count(),
        particle_count - 1
    );

    println!("{queue:#?}");
    Ok(())
}
