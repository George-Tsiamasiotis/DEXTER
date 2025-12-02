use crate::Particle;
use crate::ParticleError;
use crate::Radians;
use crate::Result;
use crate::State;
use crate::Stepper;
use crate::henon::{
    calculate_intersection_state, calculate_mod_state1, calculate_mod_state2, calculate_mod_step,
    intersected,
};
use config::*;

use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use std::f64::consts::TAU;
use std::time::Duration;

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
    pub alpha: Radians,
    /// The number of interections to calculate.
    pub intersections: usize,
}

impl MappingParameters {
    /// Creates a new [`MappingParameters`].
    pub fn new(section: PoincareSection, alpha: Radians, intersections: usize) -> Self {
        // mod `alpha` to avoid modding it in every step
        Self {
            section,
            alpha: alpha % TAU,
            intersections,
        }
    }
}

/// Calculates the PoincareSection=const intersections.
pub(crate) fn map_integrate(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    currents: &Currents,
    perturbation: &Perturbation,
    params: &MappingParameters,
) -> Result<()> {
    // Last two states of the particle.
    let mut state1 = particle.initial_state.clone(); // Already evaluated
    let mut state2: State;

    let mut dt = RKF45_FIRST_STEP;

    while particle.evolution.steps_stored() <= params.intersections {
        // Perform a step on the normal system.
        let mut stepper = Stepper::default();
        stepper.init(&state1);
        stepper.start(dt, qfactor, bfield, currents, perturbation)?;
        dt = stepper.calculate_optimal_step(dt)?;
        state2 = stepper.next_state(dt);
        state2.evaluate(qfactor, currents, bfield, perturbation)?;

        if particle.evolution.steps_taken() >= MAX_STEPS {
            return Err(ParticleError::TimedOut(Duration::default()));
        }

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
                calculate_mod_state2(qfactor, bfield, currents, perturbation, mod_state1, dtau)?;
            let intersection_state = calculate_intersection_state(
                qfactor,
                bfield,
                currents,
                perturbation,
                params,
                mod_state2,
            )?;

            // Store the intersection state.
            particle.evolution.push_state(&intersection_state);

            // NOTE: Even after landing on the intersection, we must continue the integration from
            // state2. If we continue from the intersection state, a wrong sign change detection
            // will most likely occur, which causes the particle to get stuck. However, starting
            // from state2 is not a problem, since state2 was calculated from the solver, so it is
            // a valid state with a valid step size within the solver's tolerance.
        }
        // In both cases, continue from the next state.
        state1 = state2;
        particle.evolution.steps_taken += 1;
    }
    particle.final_state = state1.into_evaluated(qfactor, currents, bfield, perturbation)?;
    Ok(())
}

/// Checks if all the value diffs in the array are within the threshold.
pub(crate) fn check_accuracy(array: &[Radians], threshold: Radians) -> Result<()> {
    // array.iter().skip(1).for_each(|v| {
    //     dbg!(v % TAU);
    // });
    match array
        .windows(2)
        .skip(1) // Skip the starting point, since it is usually not on the intersection
        .all(|v| (v[1] - v[0]).abs() - TAU < threshold)
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

        assert!(check_accuracy(&ok1, MAP_THRESHOLD).is_ok());
        assert!(check_accuracy(&ok2, MAP_THRESHOLD).is_ok());
        assert!(check_accuracy(&ok3, MAP_THRESHOLD).is_ok());

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

        assert!(check_accuracy(&not_ok1, MAP_THRESHOLD).is_err());
        assert!(check_accuracy(&not_ok2, MAP_THRESHOLD).is_err());
        assert!(check_accuracy(&not_ok3, MAP_THRESHOLD).is_err());
    }
}
