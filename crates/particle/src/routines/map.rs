//! Integration of a [`Particle`] and calculation of its exact intersections with a
//! constant θ/ζ surface.

use std::f64::consts::TAU;
use std::time::Instant;

use crate::routines::henon::{
    calculate_intersection_state, calculate_mod_state1, calculate_mod_state2, calculate_mod_step,
    intersected,
};
use crate::{Evolution, MappingConfig, Particle, State, Stepper};
use crate::{ParticleError, Result};

use equilibrium::{Bfield, Current, Perturbation, Qfactor};

/// Defines the surface of the Poincare section.
#[derive(Debug, Clone, Copy)]
pub enum PoincareSection {
    /// Defines a surface of xᵢ= θ.
    ConstTheta,
    /// Defines a surface of xᵢ= ζ.
    ConstZeta,
}

/// Defines all the necessary parameters of a Poincare Map.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub struct MappingParameters {
    /// The surface of section Σ, defined by an equation xᵢ= α, where xᵢ= θ or ζ.
    pub section: PoincareSection,
    /// The constant that defines the surface of section.
    pub alpha: f64,
    /// The number of interections to calculate.
    pub intersections: usize,
}

impl MappingParameters {
    /// Creates a new [`MappingParameters`].
    pub fn new(section: PoincareSection, alpha: f64, intersections: usize) -> Self {
        // mod `alpha` to avoid modding it in every step
        Self {
            section,
            alpha: alpha.rem_euclid(TAU),
            intersections,
        }
    }
}

/// Calculates the PoincareSection=const intersections.
pub(crate) fn map_integrate(
    particle: &mut Particle,
    qfactor: &impl Qfactor,
    current: &impl Current,
    bfield: &impl Bfield,
    perturbation: &impl Perturbation,
    params: &MappingParameters,
    config: &MappingConfig,
) -> Result<()> {
    // ==================== Setup

    let res: Result<_>;
    let start = Instant::now();
    particle.evolution = Evolution::default();
    particle
        .initial_state
        .evaluate(qfactor, current, bfield, perturbation)?;
    let mut state1 = particle.initial_state.clone();
    let mut state2: State;
    let mut dt = config.first_step;

    // ==================== Main loop

    loop {
        if particle.evolution.steps_taken() == config.max_steps {
            res = Err(ParticleError::TimedOut(start.elapsed()));
            break;
        }
        if particle.evolution.steps_stored() > params.intersections {
            res = Ok(());
            break;
        }

        // Perform a step on the normal system.
        let mut stepper = Stepper::new(&state1);
        stepper.start(dt, qfactor, current, bfield, perturbation)?;
        dt = stepper.calculate_optimal_step(dt, config)?;
        state2 = stepper
            .next_state(dt)
            .into_evaluated(qfactor, current, bfield, perturbation)?;
        particle.evolution.steps_taken += 1;

        // Hénon's trick.
        // Depending on the PoincareSection, the independent variable becomes either `zeta` or
        // `theta`. Checking its value in every function and every loop has negligible performance
        // impact and produces much more readable code, instead of rewritting the same function
        // twice.
        let (old_angle, new_angle) = match params.section {
            PoincareSection::ConstTheta => (state1.theta, state2.theta),
            PoincareSection::ConstZeta => (state1.zeta, state2.zeta),
        };
        if intersected(old_angle, new_angle, params.alpha) {
            let mod_state1 = calculate_mod_state1(&state1, &params.section);
            let dtau = calculate_mod_step(&state1, &state2, params);
            let mod_state2 =
                calculate_mod_state2(qfactor, current, bfield, perturbation, mod_state1, dtau)?;
            let intersection_state = calculate_intersection_state(
                qfactor,
                current,
                bfield,
                perturbation,
                params,
                mod_state2,
            )?;

            particle.evolution.push_state(&intersection_state);

            // NOTE: Even after landing on the intersection, we must continue the integration from
            // state2. If we continue from the intersection state, a wrong sign change detection
            // will most likely occur, which causes the particle to get stuck. However, starting
            // from state2 is not a problem, since state2 was calculated from the solver, so it is
            // a valid state with a valid step size within the solver's tolerance.
        }
        // In both cases, continue from the next state.
        state1 = state2;
    }

    // ==================== Finalization

    check_mapping_accuracy(&particle.evolution, &params.section, config)?;
    particle.final_state = state1.into_evaluated(qfactor, current, bfield, perturbation)?;
    particle.evolution.finish();
    particle.calculate_orbit_type();
    particle.evolution.duration = start.elapsed();
    res
}

/// Checks if all the value diffs in the array are within the threshold.
fn check_mapping_accuracy(
    evolution: &Evolution,
    section: &PoincareSection,
    config: &MappingConfig,
) -> Result<()> {
    let intersections_array = match section {
        PoincareSection::ConstZeta => &evolution.zeta,
        PoincareSection::ConstTheta => &evolution.theta,
    };
    _check_mapping_accuracy(intersections_array, config)
}

/// Extract for testing
fn _check_mapping_accuracy(intersections_array: &[f64], config: &MappingConfig) -> Result<()> {
    match intersections_array
        .windows(2)
        .skip(1) // Skip the starting point, since it is often not on the intersection
        .all(|v| (v[1] - v[0]).abs() - TAU < config.map_threshold)
    {
        true => Ok(()),
        false => Err(crate::ParticleError::IntersectionError),
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn test_accuracy_check() {
        let config = MappingConfig::default();
        let ok1 = [
            0.0 * TAU,
            1.0 * TAU,
            2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let ok2 = [
            100.0,
            0.0 * TAU,
            1.0 * TAU,
            2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let ok3 = [
            1.0 + 0.0 * TAU,
            1.0 + 1.0 * TAU,
            1.0 + 2.0 * TAU + 1e-12,
            1.0 + 3.0 * TAU - 1e-12,
            1.0 + 4.0 * TAU,
        ];

        assert!(_check_mapping_accuracy(&ok1, &config).is_ok());
        assert!(_check_mapping_accuracy(&ok2, &config).is_ok());
        assert!(_check_mapping_accuracy(&ok3, &config).is_ok());

        let not_ok1 = [
            0.0 * TAU,
            1.0 * TAU,
            // 2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let not_ok2 = [
            100.0,
            0.0 * TAU,
            1.0 * TAU,
            // 2.0 * TAU + 1e-12,
            3.0 * TAU - 1e-12,
            4.0 * TAU,
        ];
        let not_ok3 = [
            100.0,
            1.0 * TAU,
            2.0 * TAU,
            3.0 * TAU + 1.0,
            4.0 * TAU,
            5.0 * TAU,
        ];

        assert!(_check_mapping_accuracy(&not_ok1, &config).is_err());
        assert!(_check_mapping_accuracy(&not_ok2, &config).is_err());
        assert!(_check_mapping_accuracy(&not_ok3, &config).is_err());
    }
}
