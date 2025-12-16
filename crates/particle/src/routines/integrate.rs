//! Integration of a [`Particle`] for a specific time interval.

use std::time::Instant;

use equilibrium::{Bfield, Current, Perturbation, Qfactor};

use crate::IntegrationConfig;
use crate::{Evolution, Particle, State, Stepper};
use crate::{ParticleError, Result};

/// Integrates the particle for a specific time interval.
pub(crate) fn integrate(
    particle: &mut Particle,
    qfactor: &impl Qfactor,
    currents: &impl Current,
    bfield: &impl Bfield,
    perturbation: &impl Perturbation,
    t_eval: (f64, f64),
    config: &IntegrationConfig,
) -> Result<()> {
    // ==================== Setup

    let res: Result<()>;
    let start = Instant::now();
    particle.evolution = Evolution::default();
    particle
        .initial_state
        .evaluate(qfactor, currents, bfield, perturbation)?;
    let mut state1 = particle.initial_state.clone();
    let mut state2: State;
    let mut dt = config.first_step;

    // ==================== Main loop

    loop {
        if particle.evolution.steps_taken() == config.max_steps {
            res = Err(ParticleError::TimedOut(start.elapsed()));
            break;
        }
        if particle.evolution.final_time().unwrap_or(t_eval.0) > t_eval.1 {
            res = Ok(());
            break;
        }

        // Perform a step
        let mut stepper = Stepper::new(&state1);
        stepper.start(dt, qfactor, currents, bfield, perturbation)?;
        dt = stepper.calculate_optimal_step(dt, config)?;
        state2 = stepper
            .next_state(dt)
            .into_evaluated(qfactor, currents, bfield, perturbation)?;
        particle.evolution.steps_taken += 1;

        // Store and continue
        particle.evolution.push_state(&state2);
        state1 = state2;
    }

    // ==================== Finalization

    particle.final_state = state1.into_evaluated(qfactor, currents, bfield, perturbation)?;
    particle.evolution.finish();
    particle.calculate_orbit_type();
    particle.evolution.duration = start.elapsed();
    res
}
