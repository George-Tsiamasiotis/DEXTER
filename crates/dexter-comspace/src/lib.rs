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
    /// The density of the trapped-passing boundary curves' points.
    ///
    /// A higher number is needed to better classify Potato and Stagnated orbits.
    pub const TRAPPED_PASSING_BOUNDARY_DENSITY: usize = 500;

    /// The density of the `θ` array when constructing the grid for the calculation of the
    /// [`StationaryCurve`][crate::StationaryCurve].
    pub const SC_CONTOUR_THETA_POINTS: usize = 200;

    /// The density of the flux array when constructing the grid for the calculation of the
    /// [`StationaryCurve`][crate::StationaryCurve].
    pub const SC_CONTOUR_FLUX_POINTS: usize = 100;
}
