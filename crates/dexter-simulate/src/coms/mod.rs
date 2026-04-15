//! Calculations on the Constants of Motion (COMs) space.

mod energy;

/// The constants of motion in an unperturbed equilibrium.
#[derive(Debug)]
pub struct COMs {
    /// The Energy in Normalized Units.
    pub energy: Option<f64>,
    /// The canonical momentum `Pζ` in Normalized Units.
    pub pzeta: Option<f64>,
    /// The magnetic moment `μ` in Normalized Units.
    pub mu: Option<f64>,
}
