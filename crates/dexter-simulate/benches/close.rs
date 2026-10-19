//! [`Particle::close`] benchmark

#![allow(unused_results)]

use std::{hint::black_box, time::Duration};

use dexter_machine::*;
use dexter_simulate::*;

use criterion::{Criterion, criterion_group, criterion_main};

fn close_benchmark(c: &mut Criterion) {
    let lcfs = MagneticFlux::Toroidal(0.45);
    let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 4, 0.0)),
    ]);
    let machine = MachineBuilder::new(&qfactor, &current, &bfield)
        .with_perturbation(&perturbation)
        .build();

    // Particle setup
    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);
    let mut particle = Particle::new(&initial);

    // Close
    let solver_params = SolverParams {
        method: SteppingMethod::EnergyAdaptiveStep {
            rel_tol: 1e-9,
            abs_tol: 1e-10,
        },
        max_steps: 10_000,
        ..Default::default()
    };

    // ===========================================================================================

    let mut group = c.benchmark_group("Particle::close");
    group.measurement_time(Duration::from_secs(10));
    let _ = group.bench_function("Particle::close", |b| {
        b.iter(|| {
            black_box(particle.close(machine, 1, &solver_params));
            assert!(particle.steps_taken() > 1_000);
            assert!(matches!(
                particle.integration_status(),
                IntegrationStatus::ClosedPeriods(1)
            ),);
        })
    });

    // ===========================================================================================
}

criterion_group!(benches, close_benchmark);
criterion_main!(benches);
