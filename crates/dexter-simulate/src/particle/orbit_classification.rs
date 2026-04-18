//! Classification of a particle's orbit by projection on the `(E, Pζ, μ)` space.

use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, Qfactor};
use parabola::Point;

use crate::particle::{EqObjects, IntegrationCaches, Particle};
use crate::state::GCState;
use crate::{EnergyPzetaPlane, OrbitType};

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`OrbitType`] variant for each valid orbit calculation (including "erroneous" orbits),
/// and set [`OrbitType::Failed`] for hard errors.
///
/// If a `_plane` is passed, then that plane is used for the calculations, avoiding the
/// need for generating a new one. See [`Particle::_classify`].
pub(super) fn classify<Q, C, B, H>(
    particle: &mut Particle,
    objects: &EqObjects<Q, C, B, H>,
    _plane: Option<&EnergyPzetaPlane>,
) where
    Q: Qfactor + FluxCommute,
    C: Current,
    B: Bfield,
    H: Harmonic,
{
    // =============== Particle Setup

    // Do not alter the `evolution` or `integration_status`
    let mut caches = IntegrationCaches::<H::Cache> {
        harmonic_caches: objects.perturbation.generate_caches(),
        ..Default::default()
    };

    // Return early if the initial flux happens to be exactly 0.0 or out of bounds.
    if particle.initial_conditions().flux0.value() == 0.0 {
        particle.orbit_type = OrbitType::Failed("Out of bounds initialization".into());
        return;
    }
    if particle.initial_conditions.finalize(objects).is_err() {
        particle.orbit_type = OrbitType::Failed("Invalid initial conditions".into());
        return;
    }
    let Ok(initial_state) = GCState::new(&particle.initial_conditions, objects, &mut caches) else {
        particle.orbit_type = OrbitType::Failed("Invalid initial conditions".into());
        return;
    };

    particle.initial_energy = Some(initial_state.energy());
    particle.orbit_type = OrbitType::Unclassified; // Fallback

    // =============== Energy-Pζ plane Setup

    let mu = particle.initial_conditions.mu0;
    let plane = match _plane {
        Some(plane) => {
            #[expect(clippy::float_cmp, reason = "we need bit-to-bit equivalence")]
            if mu != plane.mu() {
                unreachable!("New EnergyPzetaPlanes must be generated");
            }
            plane
        }
        None => &EnergyPzetaPlane::from_mu(objects.qfactor, objects.current, objects.bfield, mu),
    };

    check_parabola_alphas(plane);

    let pzeta = particle
        .initial_conditions
        .pzeta0
        .expect("Initial conditions have been finalized");
    let energy = particle
        .initial_energy
        .expect("initial energy has been calculated");
    let point = Point {
        x: pzeta,
        y: energy,
    };
    let psip_last = -plane
        .left_wall_parabola()
        .axis()
        .expect("parabola's 'a' checked");
    let tp = plane.tp_boundary();

    // =============== Routine

    let is_in_axis = plane.axis_parabola().contains(point);
    let is_in_left_wall = plane.left_wall_parabola().contains(point);
    let is_in_right_wall = plane.right_wall_parabola().contains(point);
    let is_in_psip = (-psip_last..=0.0).contains(&pzeta);
    let is_trapped = is_in_psip && tp.contains(energy, pzeta); // Short-circuit
    let is_above_tp = tp.is_above(energy, pzeta);
    let is_below_tp = tp.is_below(energy, pzeta);

    if is_in_left_wall {
        particle.orbit_type = OrbitType::CuPassingConfined;
    }
    if is_in_axis && !is_trapped {
        particle.orbit_type = OrbitType::CoPassingConfined;
    }

    if !is_in_axis && (pzeta.is_sign_positive()) {
        particle.orbit_type = OrbitType::Stagnated;
    }

    if is_in_right_wall && !is_in_left_wall {
        if pzeta < -psip_last {
            particle.orbit_type = OrbitType::CuPassingLost;
        } else {
            if is_trapped {
                particle.orbit_type = OrbitType::TrappedLost;
            } else {
                particle.orbit_type = OrbitType::Failed(
                        "Unimplemented(h, i, j orbits): is_in_right_wall, !is_in_left_wall, !is_trapped, !is_in_psip".into(),
                    );
            }
        }
    }

    if is_trapped && !is_in_right_wall && !is_in_axis {
        particle.orbit_type = OrbitType::TrappedConfined;
    }

    if is_trapped && is_in_axis {
        particle.orbit_type = OrbitType::Potato;
    }

    if !is_trapped && !is_in_axis && !is_in_right_wall && is_in_psip {
        if is_above_tp {
            particle.orbit_type = OrbitType::CuPassingConfined;
        } else if is_below_tp {
            // WARN: not sure
            particle.orbit_type = OrbitType::Stagnated;
        } else {
            unreachable!("One of the above checks should pass")
        }
    }
}

// ===============================================================================================

/// Checks that all parabolas have nonzero `α` coefficients.
///
/// This is a corner case that should fatal.
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
