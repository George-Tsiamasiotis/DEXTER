//! Integration of a [`Particle`] for a specific time interval.

use std::time::Instant;

use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, HarmonicCache, Qfactor};

use crate::particle::{Caches, EqObjects, Particle, ParticleStats};
use crate::state::{GCState, SteppingMethod};
use crate::stepper_params_trait_impl;
use crate::system::{Stepper, StepperParams};

use super::IntegrationStatus;

/// Defines the parameters of the [`Particle::integrate`] routine.
///
/// See [`IntegrationParams::default`] for the default values.
#[derive(Debug, Clone)]
pub struct IntegrationParams {
    /// The optimal step calculation method.
    pub method: SteppingMethod,
    /// The maximum amount of steps a particle can make before terminating its integration.
    pub max_steps: usize,
    /// The initial time step for the RKF45 adaptive step method. The value is empirical.
    pub first_step: f64,
    /// The safety factor of the solver. Should be less than 1.0
    pub safety_factor: f64,
    /// The relative tolerance of the energy difference in every step.
    pub energy_rel_tol: f64,
    /// The absolute tolerance of the energy difference in every step.
    pub energy_abs_tol: f64,
    /// The relative tolerance of the local truncation error in every step.
    pub error_rel_tol: f64,
    /// The absolute tolerance of the local truncation error in every step.
    pub error_abs_tol: f64,
}

stepper_params_trait_impl!(&IntegrationParams);

impl Default for IntegrationParams {
    fn default() -> Self {
        Self {
            method: SteppingMethod::EnergyAdaptiveStep,
            max_steps: 1_000_000,
            first_step: 1e-1,
            safety_factor: 0.9,
            energy_rel_tol: 1e-12,
            energy_abs_tol: 1e-14,
            error_rel_tol: 1e-12,
            error_abs_tol: 1e-14,
        }
    }
}

// ===============================================================================================

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`IntegrationStatus`] variant for each possible error.
pub(super) fn integrate<Q, C, B, H>(
    particle: &mut Particle,
    objects: &EqObjects<Q, C, B, H>,
    teval: (f64, f64),
    params: &IntegrationParams,
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
    let mut dt = params.first_step;

    // =============== Main loop

    loop {
        if particle.evolution.steps_taken == params.max_steps {
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
            .inspect(|_| dt = stepper.calculate_optimal_step(dt, &params))
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
    particle.stats = Some(ParticleStats {
        psi_acc: caches.psi_acc,
        psip_acc: caches.psip_acc,
        theta_acc: caches.theta_acc,
        harmonic_cache_hits: caches.harmonic_caches.iter().map(|c| c.hits()).sum(),
        harmonic_cache_misses: caches.harmonic_caches.iter().map(|c| c.misses()).sum(),
    })
}
