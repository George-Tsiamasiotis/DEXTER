//! Definition of [`EnergyPzetaPosition`], and [`OrbitType`] structs.

use parabola::Point;
use rsl_interpolation::Accelerator;

use crate::EnergyPzetaPlane;

/// The position of an `(E-Pζ)` point on the `(E-Pζ)` plane, relative to the orbit
/// classification curves.
///
/// See the [`diagram`][EnergyPzetaPlane] for the different (E-Pζ) regions.
#[non_exhaustive]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[expect(missing_docs, reason = "see diagram")]
pub enum EnergyPzetaPosition {
    /// No classification has been attempted.
    #[default]
    Undefined,
    Alpha,
    Beta,
    Gamma,
    Delta,
    Epsilon,
    Zeta,
    Eta,
    Theta,
    Iota,
    Kappa,
    Lambda,
    Mu,
    /// Not falling under any of the above categories.
    Unclassified,
}

/// A particle's orbit type.
///
/// As described by [`R. B. White`], an orbit is classified depending on its location relative
/// to the well-defined `(E, Pζ, μ=const)`.
///
/// [`R. B. White`]: https://doi.org/10.1142/P440
#[derive(Default, Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum OrbitType {
    /// Particle has not been classified.
    #[default]
    Undefined,
    /// A Trapped-Lost particle.
    ///
    /// # Definition
    ///
    /// A particle is called *trapped* if there exists a mirror point where `ρ=0`.
    TrappedLost,
    /// A Trapped-Confined particle.
    ///
    /// # Definition
    ///
    /// A particle is called *trapped* if there exists a mirror point where `ρ=0`.
    TrappedConfined,
    /// A CoPassing-Lost particle.
    ///
    /// # Definition
    ///
    /// A particle is called *co-passing* if it is not trapped and it holds that `dot(θ)>0`.
    CoPassingLost,
    /// A CoPassing-Confined particle.
    ///
    /// # Definition
    ///
    /// A particle is called *co-passing* if it is not trapped and it holds that `dot(θ)>0`.
    CoPassingConfined,
    /// A CounterPassing-Lost particle.
    ///
    /// # Definition
    ///
    /// A particle is called *counter-passing* if it is not trapped and it holds that `dot(θ)<0`.
    CuPassingLost,
    /// A CounterPassing-Confined particle.
    ///
    /// # Definition
    ///
    /// A particle is called *counter-passing* if it is not trapped and it holds that `dot(θ)<0`.
    CuPassingConfined,
    /// A Potato particle.
    ///
    /// # Definition
    ///
    /// A particle's orbit is called a *potato* orbit if it is trapped but still circles the
    /// magnetic axis due to its drift. In the `(E, Pζ)` plane, those lie inside the intersection
    /// of the trapped-passing boundary and the magnetic axis parabola.
    Potato,
    /// A Stagnated particle.
    ///
    /// # Definition
    ///
    /// A particle is called *stagnated* if it always has positive parallel velocity but does
    /// not circle the magnetic axis. In the `(E, Pζ)` plane, those lie to the right of the
    /// trapped-passing boundary and above the magnetic axis parabola.
    Stagnated,
    /// Not falling under any of the other categories.
    Unclassified,
}

impl EnergyPzetaPosition {
    /// Creates a new `EnergyPzetaPosition` from an `(E, Pζ)` plane and a point.
    #[expect(clippy::redundant_else, reason = "hopeless")]
    #[must_use]
    pub fn new(point: Point, plane: &EnergyPzetaPlane, theta0_dot: f64) -> Self {
        #[expect(clippy::missing_panics_doc, reason = "checked")]
        let psip_last = -plane
            .left_wall_parabola()
            .axis()
            .expect("parabola's 'a' checked");

        let pzeta = point.x;
        let energy = point.y;

        let mut acc = Accelerator::new();
        let is_in_axis = plane.axis_parabola().contains(point);
        let is_in_left_wall = plane.left_wall_parabola().contains(point);
        let is_in_right_wall = plane.right_wall_parabola().contains(point);
        let is_in_psip = (-psip_last..=0.0).contains(&pzeta);

        // Short-circuit condition to avoid the more expensive interpolations
        let might_be_trapped = !is_in_left_wall && is_in_psip;

        let tpb = plane.tp_boundary();
        let is_above_tp = might_be_trapped && tpb.is_above(energy, pzeta, &mut acc);
        let is_below_tp = might_be_trapped && tpb.is_below(energy, pzeta, &mut acc);
        let is_in_tpb = is_in_psip && !is_above_tp && !is_below_tp;

        // Inside left wall

        if is_in_left_wall {
            if is_in_axis {
                return Self::Unclassified; // No allowed orbits here.
            } else {
                return Self::Alpha; // CounterPassing-Confined
            }
        }

        // Inside right wall, outside left wall. No need to check for left wall again

        if is_in_right_wall {
            if pzeta <= -psip_last {
                return Self::Beta; // CounterPassing-Lost
            } else if is_in_tpb {
                return Self::Gamma; // Trapped-Lost
            } else if is_in_axis {
                return Self::Zeta; // CoPassing-Lost
            } else {
                if theta0_dot.is_sign_positive() {
                    return Self::Delta; // CoPassing-Lost
                } else {
                    return Self::Epsilon; // CounterPassing-Confined
                }
            }
        }

        // Inside magnetic axis only. No need to check for walls again

        if is_in_axis {
            if is_in_tpb {
                return Self::Theta; // Potato
            } else {
                return Self::Eta; // CoPassing-Confined
            }
        }

        // Outside all parabolas. No need to check for any of them.

        if is_in_tpb {
            return Self::Iota; // Trapped-Confined
        }

        if is_above_tp {
            if theta0_dot.is_sign_positive() {
                return Self::Kappa; // CoPassing-Confined
            } else {
                return Self::Lambda; // CounterPassing-Lost
            }
        } else {
            if pzeta > -psip_last {
                return Self::Mu; // Stagnated
            }
        }

        Self::Unclassified
    }

    /// Calculates the [`OrbitType`] from the resolved [`EnergyPzetaPosition`].
    #[expect(clippy::match_same_arms, reason = "clearer this way")]
    #[must_use]
    pub fn orbit_type(&self) -> OrbitType {
        match *self {
            Self::Undefined => OrbitType::Undefined,
            Self::Alpha => OrbitType::CuPassingConfined,
            Self::Beta => OrbitType::CuPassingLost,
            Self::Gamma => OrbitType::TrappedLost,
            Self::Delta => OrbitType::CoPassingLost,
            Self::Epsilon => OrbitType::CuPassingConfined,
            Self::Zeta => OrbitType::CoPassingLost,
            Self::Eta => OrbitType::CoPassingConfined,
            Self::Theta => OrbitType::Potato,
            Self::Iota => OrbitType::TrappedConfined,
            Self::Kappa => OrbitType::CoPassingConfined,
            Self::Lambda => OrbitType::CuPassingConfined,
            Self::Mu => OrbitType::Stagnated,
            Self::Unclassified => OrbitType::Unclassified,
        }
    }
}
