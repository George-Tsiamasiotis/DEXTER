//! Definitions of evaluation methods of equilibrium objects.
//!
//! For analytical equilibria, this is achieved by evaluation of analytical formulas, while for
//! numerical equilibria by interpolation over the reconstructed data arrays.

use rsl_interpolation::{Accelerator, Cache};

use crate::HarmonicCache;
use crate::Result;
use crate::{Flux, Length, Radians};

// TODO: (maybe) add doctests

/// Equilibrium geometry related quantities computation
pub trait Geometry {
    /// Calculates the radial coordinate `r(ψp)` **in \[m\]**.
    fn r(&self, psip: Flux) -> Result<Length>;

    /// Calculates the poloidal flux `ψp(r)`, where r is **in \[m\]**.
    fn psip(&self, r: Length) -> Result<Flux>;

    /// Calculates the toroidal flux `ψ(ψp)`.
    fn psi(&self, r: Length) -> Result<Flux>;

    /// Calculates `R(ψp, θ)`,
    fn rlab(&self, psip: Flux, theta: Radians) -> Result<f64>;

    /// Calculates `Z(ψp, θ)`,
    fn zlab(&self, psip: Flux, theta: Radians) -> Result<f64>;

    /// Calculates the Jacobian `J(ψp, θ)`,
    fn jacobian(&self, psip: Flux, theta: Radians) -> Result<f64>;
}

/// q-factor related quantities computation.
pub trait Qfactor {
    /// Calculates the q-factor `q(ψp)`.
    fn q(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the toroidal flux `ψ(ψp)`.
    fn psi(&self, psip: Flux, acc: &mut Accelerator) -> Result<Flux>;

    /// Calculates the derivative `dψ/dψp`.
    fn dpsi_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;
}

/// Plasma current related quantities computation.
pub trait Current {
    /// Calculates `g(ψp)`
    fn g(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `I(ψp)`
    fn i(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕g(ψp)/𝜕ψp`
    fn dg_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕I(ψp)/𝜕ψp`
    fn di_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;
}

/// Magnetic field related quantities computation.
pub trait Bfield {
    /// Calculates `B(ψp, θ)`,
    fn b(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕B(ψp, θ) /𝜕ψp`.
    fn db_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕B(ψp, θ) /𝜕𝜃`.
    fn db_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;
}

/// Single Harmonic related quantities computation
pub trait Harmonic {
    /// Calculates the harmonic `α(ψp) * cos(mθ-nζ+φ(ψp))`.
    fn h(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the harmonic derivative `𝜕h/𝜕ψp`.
    fn dh_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the harmonic derivative `𝜕h/𝜕θ`.
    fn dh_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the perturbation derivative `𝜕h/𝜕ζ`.
    fn dh_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the perturbation derivative `𝜕h/𝜕t`.
    #[allow(unused_variables)]
    fn dh_dt(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        // Time-independent perturbations at the moment.
        Ok(0.0)
    }

    /// Calculates the harmonic's *amplitude* `α(ψp)`.
    fn a(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the harmonic's *amplitude* derivative `dα(ψp)/dψp`.
    fn da_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the harmonic's *phase* `φ(ψp)`.
    fn phase(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;
}

/// Perturbation related quantities computation
pub trait Perturbation {
    /// Calculates the Perturbation `Σ{ α(n,m)(ψp) * cos(mθ-nζ+φ0) }`.
    fn p(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `ψp`,
    fn dp_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `θ`.
    fn dp_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `ζ`.
    fn dp_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `t`.
    fn dp_dt(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;
}
