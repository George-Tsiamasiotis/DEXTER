//! Integration of a [`Particle`] for a specific time interval.

use std::time::Instant;

use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, HarmonicCache, Qfactor};

use crate::particle::{Caches, EqObjects, Particle, ParticleCacheStats};
use crate::state::GCState;
use crate::system::{SolverParams, Stepper};

use super::IntegrationStatus;

// ===============================================================================================

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`IntegrationStatus`] variant for each possible error.
pub(super) fn integrate<Q, C, B, H>(
    particle: &mut Particle,
    objects: &EqObjects<Q, C, B, H>,
    teval: (f64, f64),
    solver_params: &SolverParams,
) where
    Q: Qfactor + FluxCommute,
    C: Current,
    B: Bfield,
    H: Harmonic,
{
    // =============== Setup

    let start = Instant::now();
    particle.evolution.reset();
    let mut caches = Caches::<H::Cache> {
        harmonic_caches: objects.perturbation.generate_caches(),
        ..Default::default()
    };

    let mut state1 = match GCState::new(&particle.initial, objects, &mut caches) {
        Ok(state) => state,
        Err(_) => {
            particle.status = IntegrationStatus::OutOfBoundsInitialization;
            return;
        }
    };
    particle.initial_energy = Some(state1.energy());
    let mut state2: GCState;
    let mut dt = solver_params.first_step;

    // =============== Main loop

    loop {
        if particle.evolution.steps_taken == solver_params.max_steps {
            particle.status = IntegrationStatus::TimedOut(start.elapsed());
            break;
        }
        if particle.evolution.tf().is_some_and(|t| t > teval.1) {
            particle.status = IntegrationStatus::Integrated;
            break;
        }

        // Perform a step
        let mut stepper = Stepper::new(&state1);
        state2 = match stepper
            .start(dt, objects, &mut caches)
            .inspect(|_| dt = stepper.calculate_optimal_step(dt, solver_params))
            .and_then(|_| stepper.next_state(dt, objects, &mut caches))
        {
            Ok(state) => state,
            Err(_) => {
                // `start()` and `next_state()` can only fail if an evaluation is out of bounds.
                particle.status = IntegrationStatus::Escaped;
                break;
            }
        };

        // Store and continue
        particle.evolution.push_state(&state1);
        particle.evolution.steps_taken += 1;
        state1 = state2;
    }

    // =============== Finalize

    particle.evolution.duration = start.elapsed();
    particle.final_energy = Some(state1.energy());
    particle.evolution.finish();
    particle.stats = ParticleCacheStats {
        psi_acc: caches.psi_acc,
        psip_acc: caches.psip_acc,
        theta_acc: caches.theta_acc,
        harmonic_cache_hits: caches.harmonic_caches.iter().map(|c| c.hits()).sum(),
        harmonic_cache_misses: caches.harmonic_caches.iter().map(|c| c.misses()).sum(),
    }
}
