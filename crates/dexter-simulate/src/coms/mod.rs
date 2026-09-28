//! Calculations on the Constants of Motion (COMs) space.

mod energy;
mod energy_pzeta_plane;
mod tp_boundary;

pub(crate) use energy_pzeta_plane::EnergyPzetaPlane;
pub(crate) use tp_boundary::TrappedPassingBoundary;

pub use energy::{create_poloidal_grid, energy_of_psi_grid, energy_of_psip_grid};
