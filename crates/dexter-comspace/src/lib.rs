//! Calculations on the Constants of Motion (COMs) space.
//!
//! ### [`Guiding Center Hamiltonian`][GCHamiltonian]
//!
//! The guiding center Hamiltonian.
//!
//! #### Evaluations
//!
//! + [`GCHamiltonian::eval`]: Evaluates the Hamiltonian in a `(θ, ψ/ψp; Pζ, μ)` point.
//! + [`GCHamiltonian::eval_array`]: Evaluates the Hamiltonian in a `([θ], [ψ/ψp]; Pζ, μ)` 1D array.
//! + [`GCHamiltonian::eval_grid`]: Evaluates the Hamiltonian in a `([θ], [ψ/ψp]; Pζ, μ)` 2D grid.
//! + [`GCHamiltonian::eval_theta_dot_array`]: Evaluates `dot(θ)` in a `([θ], [ψ/ψp]; Pζ, μ)` 1D array.
//! + [`GCHamiltonian::build_stationary_curve`]: Builds the [`StationaryCurve`] `𝜕H/𝜕𝜃 = 0` of `machine`.
//!
//! #### Stationary Curve
//!
//! + [`StationaryCurve`]: Representation of a configuration's stationary curve `𝜕H/𝜕𝜃 = 0`.
//! + [`ParametrizedSegment`]: A 2D parametrized curve.
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
//!
//! + [`EnergyPzetaPosition`]: A region of the E - Pζ plane defined by the classification curves.
//! + [`OrbitType`]: A particle's orbit type. Can be calculated with [`EnergyPzetaPosition::orbit_type`].

mod energy_pzeta_plane;
mod error;
mod gc_hamiltonian;

// ============== Re-exports

pub use parabola::Point;

// ============== Public API

pub use error::COMError;

pub use gc_hamiltonian::GCHamiltonian;
pub use gc_hamiltonian::{ParametrizedSegment, StationaryCurve};

pub use energy_pzeta_plane::{
    EnergyPzetaPlane, EnergyPzetaPosition, OrbitType, TrappedPassingBoundary,
};

/// Crate configuration constants.
pub mod constants {
    /// The density of the trapped-passing boundary curves' points.
    ///
    /// A higher number is needed to better classify Potato and Stagnated orbits.
    pub const TRAPPED_PASSING_BOUNDARY_DENSITY: usize = 500;
}
