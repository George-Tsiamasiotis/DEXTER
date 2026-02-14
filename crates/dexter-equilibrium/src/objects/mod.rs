pub(crate) mod getters;

pub(crate) mod bfield;
pub(crate) mod currents;
pub(crate) mod geometries;
pub(crate) mod harmonics;
pub(crate) mod nc_harmonic;
pub(crate) mod perturbation;
pub(crate) mod qfactors;

/// Describes the type of equilibrium the object represents.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumType {
    /// Equilibrium reconstructed from numerical data.
    ///
    /// Evaluations are calculated by interpolation.
    Numerical,
    /// Equilibrium described by analytical formulas.
    ///
    /// Evaluations are calculated by simply evaluating the formulas.
    Analytical,
}
