//! Classification of a particle's orbit by projection on the `(E, Pζ, μ)` space.

use dexter_comspace::{EnergyPzetaPlane, EnergyPzetaPosition, OrbitType};
use dexter_machine::Machine;

use crate::particle::Particle;

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`OrbitType`] variant for each valid orbit calculation (including "erroneous" orbits),
/// and set [`OrbitType::Failed`] for hard errors.
///
/// If a `_plane` is passed, then that plane is used for the calculations, avoiding the
/// need for generating a new one. See [`Particle::_classify`].
pub(super) fn classify(
    particle: &mut Particle,
    machine: Machine,
    _plane: Option<&EnergyPzetaPlane>,
) {
    // =============== Particle Setup

    // Do not alter the `evolution` or `integration_status`
    particle.orbit_type = OrbitType::Undefined; // Fallback

    // =============== Energy-Pζ plane Setup

    let mu = particle.initial_conditions.mu0();
    let plane = match _plane {
        Some(plane) => {
            #[expect(clippy::float_cmp, reason = "we need bit-to-bit equivalence")]
            if mu != plane.mu() {
                unreachable!("New EnergyPzetaPlanes must be generated");
            }
            plane
        }
        None => &EnergyPzetaPlane::from_mu(machine, mu),
    };

    check_parabola_alphas(plane);

    // If any of the two failed to evaluate
    let (Some(point), Some(theta_dot)) = (
        particle.initial_conditions.energy_pzeta_point(),
        particle.initial_conditions.theta_dot0(),
    ) else {
        particle.energy_pzeta_position = EnergyPzetaPosition::Unclassified;
        particle.orbit_type = particle.energy_pzeta_position.orbit_type();
        return;
    };

    // =============== Routine

    particle.energy_pzeta_position = EnergyPzetaPosition::new(point, plane, theta_dot);
    particle.orbit_type = particle.energy_pzeta_position.orbit_type();
}

// ===============================================================================================

/// Checks that all parabolas have nonzero `α` coefficients.
///
/// This is a corner case that should be fatal.
fn check_parabola_alphas(plane: &EnergyPzetaPlane) {
    assert!(
        plane.axis_parabola().a != 0.0,
        "Encountered zero 'a' coefficient in magnetic axis parabola"
    );
    assert!(
        plane.left_wall_parabola().a != 0.0,
        "Encountered zero 'a' coefficient in left wall parabola"
    );
    assert!(
        plane.right_wall_parabola().a != 0.0,
        "Encountered zero 'a' coefficient in right wall parabola"
    );
}
