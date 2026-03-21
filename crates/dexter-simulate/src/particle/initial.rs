//! Definition of a Particle's initial conditions in various coordinate sets.

use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, Qfactor};
use rsl_interpolation::Accelerator;

use crate::particle::EqObjects;
use crate::{InitialFlux, Result};

/// The kind of [`InitialConditions`] set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoordinateSet {
    /// Initial conditions set in the `(t, ψ, θ, ζ, ρ, μ)` space.
    BoozerToroidal,
    /// Initial conditions set in the `(t, ψp, θ, ζ, ρ, μ)` space.
    BoozerPoloidal,
    /// Initial conditions set in the `(t, Pζ, ψ, θ, ζ, μ)` space.
    MixedToroidal,
    /// Initial conditions set in the `(t, Pζ, ψp, θ, ζ, μ)` space.
    MixedPoloidal,
}

/// A set of initial conditions.
///
/// The set can be created from one the following combinations:
///
/// + Boozer-Toroidal: Initial conditions set in the `(t, ψ, θ, ζ, ρ, μ)` space.
/// + Boozer-Poloidal: Initial conditions set in the `(t, ψp, θ, ζ, ρ, μ)` space.
/// + Mixed-Toroidal: Initial conditions set in the `(t, Pζ, ψ, θ, ζ, μ)` space.
/// + Mixed-Poloidal: Initial conditions set in the `(t, Pζ, ψp, θ, ζ, μ)` space.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct InitialConditions {
    /// The initial time.
    pub(crate) t0: f64,
    /// The initial toroidal/poloidal flux, depending on the value of [`InitialFlux`].
    pub(crate) flux0: InitialFlux,
    /// The initial `θ` angle.
    pub(crate) theta0: f64,
    /// The initial `ζ` angle.
    pub(crate) zeta0: f64,
    /// The initial parallel gyro radius `ρ`.
    pub(crate) rho0: f64,
    /// The magnetic moment `μ`.
    pub(crate) mu0: f64,
    /// The canonical momentum `Pζ`.
    pub(crate) pzeta0: f64,
    /// Which coordinate set were the initial conditions defined.
    coordinate_set: CoordinateSet,
}

impl InitialConditions {
    /// Constructs an [`InitialConditions`] from a set of Boozer variables.
    ///
    /// The set can be in either the `(t, ψ, θ, ζ, ρ, μ)` or `(t, ψp, θ, ζ, ρ, μ)` configuration space.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// let psi0 = InitialFlux::Toroidal(0.2);
    /// let initial = InitialConditions::boozer(0.0, psi0, 0.0, 0.0, 1e-4, 1e-6);
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
        let coordinate_set = match flux0 {
            InitialFlux::Toroidal(_) => CoordinateSet::BoozerToroidal,
            InitialFlux::Poloidal(_) => CoordinateSet::BoozerPoloidal,
        };
        Self {
            t0,
            flux0,
            theta0,
            zeta0,
            rho0,
            mu0,
            coordinate_set,
            pzeta0: f64::NAN,
        }
    }

    /// Constructs an [`InitialConditions`] from a set of Boozer variables.
    ///
    /// The set can be in either the `(t, ψ, θ, ζ, ρ, μ)` or `(t, ψp, θ, ζ, ρ, μ)` configuration space.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// # use std::f64::consts::PI;
    /// let psip0 = InitialFlux::Poloidal(0.01);
    /// let initial = InitialConditions::mixed(0.0, -0.027, psip0, PI, PI/2.0, 1e-6);
    /// ```
    #[must_use]
    pub fn mixed(
        t0: f64,
        pzeta0: f64,
        flux0: InitialFlux,
        theta0: f64,
        zeta0: f64,
        mu0: f64,
    ) -> Self {
        let coordinate_set = match flux0 {
            InitialFlux::Toroidal(_) => CoordinateSet::MixedToroidal,
            InitialFlux::Poloidal(_) => CoordinateSet::MixedPoloidal,
        };
        Self {
            t0,
            pzeta0,
            flux0,
            theta0,
            zeta0,
            mu0,
            coordinate_set,
            rho0: f64::NAN,
        }
    }

    /// Calculates the missing coordinates needed for the integration routines.
    ///
    /// # Errors
    ///
    /// Returns a [`SimulationError`] if the missing coordinates cannot be calculated. This can
    /// occur if the [`Current`] object specifically does not define g(ψ) (`MixedToroidal` case) or
    /// g(ψp) (`MixedPoloidal` case), which are necessary for calculating the fluxes.
    pub(crate) fn finalize<Q, C, B, H>(&mut self, objects: &EqObjects<Q, C, B, H>) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let acc = &mut Accelerator::new();
        // Calculate `rho0`
        match self.coordinate_set {
            CoordinateSet::MixedToroidal => {
                let psi0 = self.flux0.value();
                let g_of_psi0 = objects.current.g_of_psi(psi0, acc)?;
                let psip0 = objects.qfactor.psip_of_psi(psi0, acc)?;
                self.rho0 = (self.pzeta0 + psip0) / g_of_psi0;
            }
            CoordinateSet::MixedPoloidal => {
                let psip0 = self.flux0.value();
                let g_of_psip0 = objects.current.g_of_psip(psip0, acc)?;
                self.rho0 = (self.pzeta0 + psip0) / g_of_psip0;
            }
            _ => (), // Boozer is ready for integration
        };
        self.assert_integration_ready();
        Ok(())
    }

    /// Asserts the necessary coordinates are initialized.
    fn assert_integration_ready(&self) {
        assert!(
            self.t0.is_finite()
                && self.flux0.value().is_finite()
                && self.theta0.is_finite()
                && self.zeta0.is_finite()
                && self.rho0.is_finite()
                && self.mu0.is_finite(),
            "Invalid InitialConditions"
        )
    }
}

/// Getters.
impl InitialConditions {
    /// Returns the initial time `t0`.
    #[must_use]
    pub fn t0(&self) -> f64 {
        self.t0
    }

    /// Returns the initial flux `flux0`.
    #[must_use]
    pub fn flux0(&self) -> InitialFlux {
        self.flux0
    }

    /// Returns the initial `θ`.
    #[must_use]
    pub fn theta0(&self) -> f64 {
        self.theta0
    }

    /// Returns the initial `ζ`.
    #[must_use]
    pub fn zeta0(&self) -> f64 {
        self.zeta0
    }

    /// Returns the initial `ρ`.
    #[must_use]
    pub fn rho0(&self) -> f64 {
        self.rho0
    }

    /// Returns the initial `μ`.
    #[must_use]
    pub fn mu0(&self) -> f64 {
        self.mu0
    }

    /// Returns the initial `Pζ`.
    #[must_use]
    pub fn pzeta0(&self) -> f64 {
        self.pzeta0
    }

    /// Returns the [`CoordinateSet`].
    #[must_use]
    pub fn coordinate_set(&self) -> CoordinateSet {
        self.coordinate_set.clone()
    }
}

#[cfg(test)]
mod test {
    use dexter_equilibrium::extract::POLOIDAL_TEST_NETCDF_PATH;
    use dexter_equilibrium::*;
    use std::{f64::consts::PI, path::PathBuf};

    use super::*;
    use crate::*;
    use InitialFlux::*;

    #[test]
    fn boozer_initial_conditions() {
        let objects = EqObjects {
            qfactor: &UnityQfactor::new(),
            current: &LarCurrent::new(),
            bfield: &LarBfield::new(),
            perturbation: &Perturbation::zero(),
        };
        let mut i1 = InitialConditions::boozer(0.0, Toroidal(0.02), PI, PI, 1e-4, 1e-6);
        i1.finalize(&objects).unwrap();
        assert_eq!(i1.coordinate_set, CoordinateSet::BoozerToroidal);

        let mut i2 = InitialConditions::boozer(0.0, Poloidal(0.02), PI, PI, 1e-4, 1e-6);
        i2.finalize(&objects).unwrap();
        assert_eq!(i2.coordinate_set, CoordinateSet::BoozerPoloidal);

        let mut initial = InitialConditions::boozer(0.0, Toroidal(0.02), PI, PI, 1e-4, 1e-6);
        initial.finalize(&objects).unwrap();
        let _ = initial.t0();
        let _ = initial.flux0();
        let _ = initial.theta0();
        let _ = initial.zeta0();
        let _ = initial.rho0();
        let _ = initial.mu0();
        let _ = initial.pzeta0();
    }

    #[test]
    fn mixed_toroidal_initial_conditions() {
        let objects = EqObjects {
            qfactor: &UnityQfactor::new(),
            current: &LarCurrent::new(),
            bfield: &LarBfield::new(),
            perturbation: &Perturbation::zero(),
        };
        let mut initial = InitialConditions::mixed(0.0, -0.027, Toroidal(0.01), PI, PI, 1e-6);
        initial.finalize(&objects).unwrap();
        assert_eq!(initial.coordinate_set, CoordinateSet::MixedToroidal);

        let mut particle = Particle::new(&initial);
        particle.integrate(
            objects.qfactor,
            objects.current,
            objects.bfield,
            objects.perturbation,
            (0.0, 1e2),
            &SolverParams::default(),
        );
        assert!(particle.steps_taken() > 10);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        // LarCurrent defines ψp evaluations but LarBfield cannot integrate with respect to ψp.
        InitialConditions::mixed(0.0, -0.027, Poloidal(0.01), PI, PI, 1e-6)
            .finalize(&objects)
            .unwrap();
    }

    #[test]
    fn mixed_poloidal_initial_conditions() {
        let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
        let qfactor = &NcQfactorBuilder::new(&path, "steffen").build().unwrap();
        let current = &NcCurrentBuilder::new(&path, "steffen").build().unwrap();
        let bfield = &NcBfieldBuilder::new(&path, "bicubic").build().unwrap();
        let perturbation = &Perturbation::zero();

        let objects = EqObjects {
            qfactor,
            current,
            bfield,
            perturbation,
        };

        let mut initial = InitialConditions::mixed(0.0, -0.027, Poloidal(0.01), PI, PI, 1e-6);
        initial.finalize(&objects).unwrap();

        let mut particle = Particle::new(&initial);
        particle.integrate(
            objects.qfactor,
            objects.current,
            objects.bfield,
            objects.perturbation,
            (0.0, 1e2),
            &SolverParams::default(),
        );
        assert!(particle.steps_taken() > 10);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        let _ = InitialConditions::mixed(0.0, -0.027, Toroidal(0.01), PI, PI, 1e-6)
            .finalize(&objects)
            .unwrap_err();
    }
}
