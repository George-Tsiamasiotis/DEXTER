//! Definition of a Critical Point of the Hamiltonian.

use dexter_machine::{Machine, MagneticFlux};

use crate::gc_hamiltonian::Hessian;
use crate::{COMError, GCHamiltonian};

/// The kind of a [`CriticalPoint`].
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum CriticalPointKind {
    /// Unknown kind.
    Unknown,
    /// O-point and local minimum of the Hamiltonian.
    OpointMin,
    /// O-point and local maximum of the Hamiltonian.
    OpointMax,
    /// X-point of the Hamiltonian.
    Xpoint,
}

/// Representation of a Critical Point of the Hamiltonian.
#[derive(Debug, Clone)]
pub struct CriticalPoint {
    /// The point's `θ` position.
    theta: f64,
    /// The point's `ψ` position.
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
    /// Creates a new `CriticalPoint` from its `(θ, ψ)` position and its constants `Pζ` and `μ`.
    ///
    /// The point's energy is calculated from the passed arguments.
    ///
    /// The point should be constructed with parameters obtained by bifurcation analysis.
    ///
    /// # Errors
    ///
    /// This method returns an [`EvalError`][dexter_machine::EvalError] if the calculation of its
    /// energy fails.
    pub(crate) fn new(
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

    /// Calculates the point's [`CriticalPointKind`] by evaluating the Hessian's determinant and
    /// sign of `𝜕²H/𝜕θ²`.
    pub fn classify(&mut self, machine: Machine) {
        let hessian = match Hessian::eval(machine, self.theta, self.flux, self.pzeta, self.mu) {
            Ok(hessian) => hessian,
            Err(err) => {
                unreachable!("Hessian::eval cannot fail as all critical points are valid ({err})")
            }
        };

        self.kind = if hessian.det().is_sign_positive() {
            if hessian.dtheta2().is_sign_positive() {
                CriticalPointKind::OpointMin
            } else {
                CriticalPointKind::OpointMax
            }
        } else {
            CriticalPointKind::Xpoint
        }
    }
}

/// Getters.
impl CriticalPoint {
    /// Returns the point's `θ` position.
    #[must_use]
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// Returns the point's `ψ` position.
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
