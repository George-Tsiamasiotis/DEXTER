//! Test `Queue::integrate` routine.

#![allow(non_snake_case)]

use dexter_equilibrium::*;
use dexter_simulate::*;
use ndarray::Array1;

#[test]
fn queue_close_parQ_larC_larB_cosP() -> Result<()> {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = ParabolicQfactor::new(1.1, 1.9, LastClosedFluxSurface::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-3, lcfs, 1, 2, 0.0),
        CosHarmonic::new(1e-3, lcfs, 1, 4, 0.0),
    ]);

    // Initial Conditions setup
    let particle_count = 10;
    let psis = qfactor.psi_last() * Array1::linspace(0.0, 1.0, particle_count);
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

    queue.close(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        1,
        &SolverParams::default(),
    );

    // All but the first at `ψ=0` should be closed.
    assert!(queue[0].integration_status() == IntegrationStatus::OutOfBoundsInitialization);
    assert_eq!(
        queue
            .iter()
            .filter(|particle| particle.integration_status() == IntegrationStatus::ClosedPeriods(1))
            .count(),
        particle_count - 1
    );

    println!("{queue:#?}");
    Ok(())
}
