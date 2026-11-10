//! Defines the `Hessian` struct.

use std::f64::consts::TAU;

use dexter_machine::{Accelerator2d, EvalError, Machine, MagneticFlux};

/// The Hamiltonian's Hessian evaluated at a specific `(θ, ψ, Pζ, μ)` point.
///
/// The kind of flux depends on the `flux` parameter of the [`Self::eval`] method.
#[non_exhaustive]
pub(crate) struct Hessian {
    /// The `𝜕²H/𝜕𝜃²` value.
    dtheta2: f64,
    /// The `𝜕²H/𝜕ψ²` value.
    dflux2: f64,
    /// The `𝜕²H/𝜕θ𝜕ψ` value.
    mixed: f64,
}

impl Hessian {
    /// Evaluates the Hessian at a specific `(θ, ψ, Pζ, μ)` point.
    ///
    /// # Errors
    ///
    /// This method returns an error if any of the machine's evaluation methods fails. This should
    /// however be unreachable as this method is called by properly calculated critical points.
    pub(crate) fn eval(
        machine: Machine,
        theta: f64,
        flux: MagneticFlux,
        pzeta: f64,
        mu: f64,
    ) -> Result<Self, EvalError> {
        use dexter_machine::MagneticFluxKind::*;

        // NOTE: In this scope, all quantities and their derivatives are expressed with respect to
        // the passed `flux` by properly taking its `kind` into account. Note that in all evaluation
        // methods, only `flux` is passed as an argument while `psip` is only used once in the
        // calculation of `ρ`. This ensures that the kind of flux is always 'good' and no imprecise
        // evaluations with respect to the other flux take place.

        let acc = &mut Accelerator2d::new();
        let theta = theta.rem_euclid(TAU);
        let psip = match flux.kind() {
            Toroidal => machine.qfactor().eval_other(flux, acc.xacc())?,
            Poloidal => flux,
        };

        // Define the 'intermediate' quantity `ρ=(Pζ+ψp)/g`.
        //
        // This term 'absorbs' the Hamiltonian's explicit dependence on `ψp`. This means that any
        // necessary 'manual' chain rules that depend on the kind of `flux` can be absorbed by `ρ`
        // and its derivatives, while all other quantities that appear in the Hamiltonian's second
        // derivatives are automatically calculated with respect to the correct kind of `flux`.
        //
        // For simplicity, `ρ` and its derivatives are expressed as functions of `ψp`, meaning that
        // in the case where `flux.kind == Toroidal`, the proper chain rules must be performed.
        let g = machine.current().eval_g(flux, acc.xacc())?;
        let rho = (pzeta + psip.value()) / g;

        // Calculate `dρ/dψ` and `d²ρ/dψ²`, and perform the chain rules if needed.
        let mut dg = machine.current().eval_g_deriv(flux, acc.xacc())?;
        let mut dg2 = machine.current().eval_g_deriv2(flux, acc.xacc())?;
        let (drho, drho2) = match flux.kind() {
            Toroidal => {
                // Use the first and second order chain rules to calculate the correct derivative.
                //
                // Important: At this point, `dg` and `dg2` are calculated with respect to `psit`.
                // The second order chain rule *must be performed first* before `dg` is converted,
                // otherwise an erroneous double chain rule takes place.
                let iota = machine.qfactor().eval_iota(flux, acc.xacc())?;
                let diota = machine.qfactor().eval_iota_deriv(flux, acc.xacc())?;
                dg2 = iota.powi(2) * dg2 + diota * dg;
                dg *= iota;

                let (drho, drho2) = Self::rho_derivs(rho, g, dg, dg2);
                let drho2 = iota.powi(2) * drho2 + diota * drho;
                let drho = iota * drho;
                (drho, drho2)
            }
            Poloidal => Self::rho_derivs(rho, g, dg, dg2),
        };

        // =============================

        // Here, `flux.kind()` automatically ensures that the quantities are calculated with respect
        // to the correct flux, and can therefore be inserted in the Hamiltonian's derivatives as
        // they are without further chain rules.

        let b = machine.bfield().eval_b(flux, theta, acc)?;
        let db_dflux = machine.bfield().eval_deriv_flux(flux, theta, acc)?;
        let db_dflux2 = machine.bfield().eval_deriv_flux2(flux, theta, acc)?;
        let db_dtheta = machine.bfield().eval_deriv_theta(flux, theta, acc)?;
        let db_dtheta2 = machine.bfield().eval_deriv_theta2(flux, theta, acc)?;
        let db_mixed = machine.bfield().eval_deriv_mixed(flux, theta, acc)?;
        let dphi2 = machine.efield().eval_phi_deriv2(flux)?;

        let rb = rho * b;
        let mu_par = mu + rho * rb;
        let rb_par = drho * b + rho * db_dflux;

        let dtheta2 = mu_par * db_dtheta2 + (rho * db_dtheta).powi(2);

        let dflux2 = rb_par.powi(2)
            + rb * (drho2 * b + rho * db_dflux2 + 2.0 * drho * db_dflux)
            + mu * db_dflux2
            + dphi2;

        let mixed = 2.0 * rho * drho * b * db_dtheta
            + rho.powi(2) * db_dtheta * db_dflux
            + mu_par * db_mixed;

        Ok(Self {
            dtheta2,
            dflux2,
            mixed,
        })
    }

    /// Calculates the first and second derivatives of `ρ=(Pζ+ψp)/g` with respect to the poloidal
    /// flux `ψp`, with `dg` and `dg2` being the derivatives of `g` with respect to `ψp`.
    fn rho_derivs(rho: f64, g: f64, dg: f64, dg2: f64) -> (f64, f64) {
        let drho = (1.0 - rho * dg) / g;
        let drho2 = -2.0 * dg / g.powi(2) - rho * (dg2 / g - 2.0 * (dg / g).powi(2));
        (drho, drho2)
    }

    /// Returns the Hessian's determinant.
    pub(crate) fn det(&self) -> f64 {
        self.dtheta2 * self.dflux2 - self.mixed.powi(2)
    }

    /// Returns the Hessian's `𝜕²H/𝜕θ²` term.
    pub(crate) fn dtheta2(&self) -> f64 {
        self.dtheta2
    }
}
