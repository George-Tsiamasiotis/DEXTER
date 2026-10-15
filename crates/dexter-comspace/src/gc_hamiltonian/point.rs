//! Definition of a Critical Point of the Hamiltonian.

use dexter_machine::{Machine, MagneticFlux};

use crate::{COMError, GCHamiltonian};

/// The kind of critical point.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum CriticalPointKind {
    /// Unknown kind.
    Unknown,
}

/// Representation of a Critical Point of the Hamiltonian.
#[derive(Debug, Clone)]
pub struct CriticalPoint {
    /// The point's `θ` position.
    theta: f64,
    /// The point's `ψ/ψp` position.
    flux: MagneticFlux,
    /// The point's energy.
    energy: f64,
    /// The point's `Pζ`.
    pzeta: f64,
    /// The point's `μ`.
    mu: f64,
    /// The point's kind.
    kind: CriticalPointKind,
}

impl CriticalPoint {
    /// Creates a new `CriticalPoint` from its `(θ, ψ/ψp)` position and its constants `Pζ` and `μ`.
    ///
    /// The point's energy is calculated from the passed arguments.
    ///
    /// The point should be constructed with parameters obtained by bifurcation analysis.
    ///
    /// # Errors
    ///
    /// This method returns an [`EvalError`][dexter_machine::EvalError] if the calculation of its
    /// energy fails.
    pub fn new(
        machine: Machine,
        theta: f64,
        flux: MagneticFlux,
        pzeta: f64,
        mu: f64,
    ) -> Result<Self, COMError> {
        let energy = GCHamiltonian::eval(machine, theta, flux, pzeta, mu)?;
        Ok(Self {
            theta,
            flux,
            energy,
            pzeta,
            mu,
            kind: CriticalPointKind::Unknown,
        })
    }

    /// Classifies the critical point.
    ///
    /// # Panics
    ///
    /// This method panics as it is not yet implemented.
    pub fn classify(&mut self) {
        unimplemented!()
    }
}

/// Getters.
impl CriticalPoint {
    /// Returns the point's `θ` position.
    #[must_use]
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// Returns the point's `ψ/ψp` position.
    #[must_use]
    pub fn flux(&self) -> MagneticFlux {
        self.flux
    }

    /// Returns the point's `Pζ`.
    #[must_use]
    pub fn pzeta(&self) -> f64 {
        self.pzeta
    }

    /// Returns the point's `μ`.
    #[must_use]
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// Returns the point's energy.
    #[must_use]
    pub fn energy(&self) -> f64 {
        self.energy
    }

    /// Returns the point's kind.
    #[must_use]
    pub fn kind(&self) -> CriticalPointKind {
        self.kind
    }
}
