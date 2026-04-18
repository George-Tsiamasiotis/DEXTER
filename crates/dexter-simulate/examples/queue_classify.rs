//! Compare `Queue::classify` and `Queue::classify_common_mu` routines

#![allow(non_snake_case)]
#![allow(unused_variables)]

use std::time::Instant;

use dexter_equilibrium::*;
use dexter_simulate::*;
use ndarray::Array1;

fn main() -> Result<(), SimulationError> {
    // Equilibrium setup
    let lcfs = LastClosedFluxSurface::Toroidal(0.05);
    let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();

    // Initial Conditions setup
    let particle_count = 1000000;
    let pzetas = qfactor.psip_last() * Array1::linspace(-1.4, 0.2, particle_count);
    let psis = qfactor.psi_last() * Array1::linspace(0.001, 0.5, particle_count);
    let psis = toroidal_fluxes(&psis.to_vec());
    let initial_conditions = QueueInitialConditions::mixed(
        Array1::zeros(particle_count).as_slice().unwrap(),
        &psis.to_vec(),
        Array1::ones(particle_count).as_slice().unwrap(),
        Array1::zeros(particle_count).as_slice().unwrap(),
        pzetas.as_slice().unwrap(),
        Array1::from_elem(particle_count, 7e-6).as_slice().unwrap(),
    )?;

    // ===========================================================================================

    let mut no_common_queue = Queue::new(&initial_conditions);
    let no_common_start = Instant::now();
    no_common_queue.classify(&qfactor, &current, &bfield);
    let no_common_elapsed = no_common_start.elapsed();

    let mut common_queue = Queue::new(&initial_conditions);
    let common_start = Instant::now();
    common_queue.classify_common_mu(&qfactor, &current, &bfield);
    let common_elapsed = common_start.elapsed();

    // Sanity check
    for i in 0..no_common_queue.particle_count() {
        assert_eq!(
            no_common_queue[i].orbit_type(),
            common_queue[i].orbit_type()
        );
    }

    println!("===============================================================");
    println!("Number of particles: {particle_count}");
    println!("`Queue::classify` duration: {no_common_elapsed:?}");
    println!("`Queue::classify_common_mu` duration: {common_elapsed:?}");
    println!("===============================================================");

    Ok(())
}
