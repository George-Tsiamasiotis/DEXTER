//! `Queue::intersect` routine.

use dexter_equilibrium::*;
use dexter_simulate::*;
use ndarray::Array1;

fn main() -> Result<()> {
    // Equilibrium setup
    let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-3, 1, 2, 0.0),
        CosHarmonic::new(1e-3, 1, 4, 0.0),
    ]);

    // Initial Conditions setup
    let particle_count = 200;
    let psis = qfactor.psi_last() * Array1::linspace(0.1, 0.9, particle_count);
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
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 0.0, 100);
    queue.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &SolverParams::default(),
    );
    println!("{queue:#?}");
    Ok(())
}
