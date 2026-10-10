//! Calculations on the Constants of Motion (COMs) space.
//!
//! ### Energy grid calculation
//!
//! + [`create_poloidal_grid`]: Creates a `(θ, ψ/ψp)` meshgrid from the two 1D arrays.
//! + [`energy_of_psi_grid`]: Calculates the energy on a 2D meshgrid of the `θ` and `ψ` arrays.
//! + [`energy_of_psip_grid`]: Calculates the energy on a 2D meshgrid of the `θ` and `ψp` arrays.
//!
//! ### Fixed Point Analysis
//!
//! + [`StationaryCurve`]: Representation of a configuration's stationary curve.
//!     - [`StationaryCurveSegment`]: Segments comprising the stationary curve.
//!
//! ### Calculations on the E - Pζ plane
//!
//! + [`EnergyPzetaPlane`]: Representation of the COM space `(E, Pζ, μ=const)`. Contains the three
//!   classification parabolas and the Trapped-Passing Boundary:
//!     - [`EnergyPzetaPlane::axis_parabola`]: Magnetic Axis Parabola (MA).
//!     - [`EnergyPzetaPlane::left_wall_parabola`]: Left Wall Parabola (LW).
//!     - [`EnergyPzetaPlane::right_wall_parabola`]: Right Wall Parabola (RW).
//!     - [`EnergyPzetaPlane::tp_boundary`][TrappedPassingBoundary]: Representation of the
//!       Trapped-Passing Boundary.
//! + [`EnergyPzetaPosition`]: A region of the E - Pζ plane defined by the classification curves.
//! + [`OrbitType`]: A particle's orbit type. Can be calculated with [`EnergyPzetaPosition::orbit_type`].

mod bifurcation;
mod energy_grid;
mod energy_pzeta_plane;

// ============== Re-exports

pub use parabola::Point;

// ============== Public API

pub use energy_pzeta_plane::{
    EnergyPzetaPlane, EnergyPzetaPosition, OrbitType, TrappedPassingBoundary,
};

pub use energy_grid::{create_poloidal_grid, energy_of_psi_grid, energy_of_psip_grid};

pub use bifurcation::{StationaryCurve, StationaryCurveSegment};

/// Crate configuration constants.
pub mod constants {
    use std::f64::consts::TAU;
    use std::ops::Range;

    /// The density of the trapped-passing boundary curves' points.
    ///
    /// A higher number is needed to better classify Potato and Stagnated orbits.
    pub const TRAPPED_PASSING_BOUNDARY_DENSITY: usize = 500;

    // ============== Stationary Curve

    /// The density of the `θ` array when constructing the grid for the calculation of the
    /// [`StationaryCurve`][crate::StationaryCurve].
    pub const SC_CONTOUR_THETA_POINTS: usize = 300;

    /// The density of the flux array when constructing the grid for the calculation of the
    /// [`StationaryCurve`][crate::StationaryCurve].
    pub const SC_CONTOUR_FLUX_POINTS: usize = 300;

    /// The bounds of the `theta` array in which a stationary curve segment is considered valid.
    pub const SC_THETA_BOUND: Range<f64> = (-1e-4 * TAU)..((1.0 - 1e-4) * TAU);

    /// The bounds of the `flux` array in which a stationary curve segment is considered valid.
    /// Units are normalized w.r.t. the flux value at the last closed surface.
    pub const SC_FLUX_BOUND: Range<f64> = (1e-8)..(1.0 - 1e-10);

    /// The `theta` range in which slightly negative values are clamped to 0.
    pub const SC_THETA_CLAMP_RANGE: Range<f64> = -1e-4 * TAU..0.0;

    /// The maximum number of gradient decent steps each segment is allowed to take.
    pub const SCGD_MAX_ITERATIONS: usize = 20;

    /// The maximum allowed deviation of a point from the true stationary curve. If a point has a value
    /// of `𝜕B/𝜕𝜃` less than `MAX_DEVIATION`, then it is allowed to terminate its gradient decent.
    pub const SCGD_MAX_DEVIATION: f64 = 1e-7;

    /// The maximum step allowed during gradient descent, with respect to `(2π, LCFS)`.
    pub const SCGD_MAX_STEP: (f64, f64) = (1e-4, 1e-4);

    /// The minimum step allowed during gradient descent, with respect to `(2π, LCFS)`.
    pub const SCGD_MIN_STEP: (f64, f64) = (1e-9, 1e-9);

    /// The first step of the gradient descent, with respect to `(2π, LCFS)`.
    pub const SCGD_FIRST_STEP: (f64, f64) = (1e-4, 1e-4);
}
