//! Definition of a Boundary Extremum of the Hamiltonian.

use std::f64::consts::TAU;

use dexter_machine::{Accelerator2d, EvalError, Machine, MagneticFlux};

use crate::{COMError, GCHamiltonian};

/// The kind of boundary extremum.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum BoundaryExtremumKind {
    /// Unknown kind.
    Unknown,
    /// Extremum point of minimum energy.
    Minimum,
    /// Extremum point of maximum energy.
    Maximum,
}

/// Representation of a Boundary extremum of the Hamiltonian.
#[derive(Debug, Clone)]
pub struct BoundaryExtremum {
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
    kind: BoundaryExtremumKind,
}

impl BoundaryExtremum {
    /// Creates a new `BoundaryExtremum` from its `(θ, ψ)` position and its constants `Pζ` and `μ`.
    ///
    /// The point's energy is calculated from the passed arguments.
    ///
    /// The point should be constructed with parameters obtained by boundary extremum analysis.
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
            kind: BoundaryExtremumKind::Unknown,
        })
    }

    /// Classifies the extremum point by evaluating the Hamiltonian's second derivative with respect
    /// to `θ`.
    pub fn classify(&mut self, machine: Machine) {
        match self._classify(machine) {
            Ok(_) => (),
            Err(err) => {
                unreachable!("Classification cannot fail as all boundary extrema are valid ({err})")
            }
        };
    }

    /// Does the classification.
    ///
    /// The [`Self::classify`] method is called internally after a valid computation of a boundary
    /// extremum, which by construction is infallible, therefore it can never err. However we return
    /// a Result to keep the code clearer.
    fn _classify(&mut self, machine: Machine) -> Result<(), EvalError> {
        use dexter_machine::MagneticFluxKind::*;

        // Follow the same methodology as in the Hessian's evaluation. This case is much simpler as
        // no derivatives of `ρ` appear.

        let theta = self.theta;
        let flux = self.flux;
        let acc = &mut Accelerator2d::new();
        let theta = theta.rem_euclid(TAU);
        let psip = match flux.kind() {
            Toroidal => machine.qfactor().eval_other(flux, acc.xacc())?,
            Poloidal => flux,
        };

        let g = machine.current().eval_g(flux, acc.xacc())?;
        let rho = (self.pzeta + psip.value()) / g;

        let b = machine.bfield().eval_b(flux, theta, acc)?;
        let db_dtheta = machine.bfield().eval_deriv_theta(flux, theta, acc)?;
        let db_dtheta2 = machine.bfield().eval_deriv_theta2(flux, theta, acc)?;

        let rb = rho * b;
        let mu_par = self.mu + rho * rb;

        let dtheta2 = mu_par * db_dtheta2 + (rho * db_dtheta).powi(2);

        self.kind = if dtheta2.is_sign_positive() {
            BoundaryExtremumKind::Minimum
        } else {
            BoundaryExtremumKind::Maximum
        };
        Ok(())
    }
}

/// Getters.
impl BoundaryExtremum {
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
    pub fn kind(&self) -> BoundaryExtremumKind {
        self.kind
    }
}
