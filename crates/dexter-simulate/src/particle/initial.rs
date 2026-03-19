//! Definition of a Particle's initial conditions in various coordinate sets.

use crate::InitialFlux;

/// A set of initial conditions.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct InitialConditions {
    /// The initial time.
    pub t0: f64,
    /// The initial toroidal/poloidal flux, depending on the value of [`InitialFlux`].
    pub flux0: InitialFlux,
    /// The initial `θ` angle.
    pub theta0: f64,
    /// The initial `ζ` angle.
    pub zeta0: f64,
    /// The initial parallel gyro radius `ρ`.
    pub rho0: f64,
    /// The magnetic moment `μ`.
    pub mu0: f64,
}

impl InitialConditions {
    /// Constructs an [`InitialConditions`] from a set of Boozer variables.
    ///
    /// The set can be in either the `(t, ψ, θ, ζ, ρ, μ)` or `(t, ψp, θ, ζ, ρ, μ)` configuration space.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);
    /// ```
    #[must_use]
    pub fn boozer(
        t0: f64,
        flux0: InitialFlux,
        theta0: f64,
        zeta0: f64,
        rho0: f64,
        mu0: f64,
    ) -> Self {
        Self {
            t0,
            flux0,
            theta0,
            zeta0,
            rho0,
            mu0,
        }
    }
}
