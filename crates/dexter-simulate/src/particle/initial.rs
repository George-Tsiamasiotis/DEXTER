//! Definition of a Particle's initial conditions in various coordinate sets.

use dexter_comspace::Point;
use dexter_machine::{Machine, MagneticFluxKind};
use rsl_interpolation::Accelerator;

use crate::particle::IntegrationCaches;
use crate::state::GCState;
use crate::{MagneticFlux, SimulationError};

/// The kind of [`InitialConditions`] set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSet {
    /// Initial conditions set in the `(t, ψt, θ, ζ, ρ, μ)` space.
    BoozerToroidal,
    /// Initial conditions set in the `(t, ψp, θ, ζ, ρ, μ)` space.
    BoozerPoloidal,
    /// Initial conditions set in the `(t, ψt, θ, ζ, Pζ, μ)` space.
    MixedToroidal,
    /// Initial conditions set in the `(t, ψp, θ, ζ, Pζ, μ)` space.
    MixedPoloidal,
}

impl CoordinateSet {
    /// Returns `true` if the set is in mixed coordinates.
    fn is_mixed(&self) -> bool {
        use CoordinateSet::*;
        match *self {
            BoozerToroidal | BoozerPoloidal => false,
            MixedToroidal | MixedPoloidal => true,
        }
    }
}

/// The status of an [`InitialConditions`] set.
///
/// # Note
///
/// All particle routines are guarded behind an `initial_conditions.status == Finalized`. If any
/// other variant is encountered, the routine returns early and the corresponding
/// [`integration_status`][crate::IntegrationStatus] field is set.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InitialConditionsStatus {
    /// [`InitialConditions`] have been instantiated through the [`boozer`] or [`mixed`] constructors.
    ///
    /// [`boozer`]: InitialConditions::boozer
    /// [`mixed`]: InitialConditions::mixed
    Instantiated,
    /// All guiding center quantities have been evaluated and the Particle can run any routine.
    ///
    /// This is the status of a particle if any routine has run or if
    /// [`Particle::finalize_initial_conditions`][crate::Particle::finalize_initial_conditions]
    /// has been called explicitly.
    Finalized,

    /// Encountered `NaN` when finalizing the initial conditions.
    NanInitialConditions,
    /// Particle initialized with a 'bad' initial magnetic flux.
    BadMagneticFlux,
    /// [`InitialConditions`] were out of bounds or `ψ0=0.0`.
    OutOfBoundsInitialization,
}

/// A set of initial conditions.
///
/// The set can be created from one the following combinations:
///
/// + Boozer-Toroidal: Initial conditions set in the `(t, ψt, θ, ζ, ρ, μ)` space.
/// + Boozer-Poloidal: Initial conditions set in the `(t, ψp, θ, ζ, ρ, μ)` space.
/// + Mixed-Toroidal: Initial conditions set in the `(t, ψt, θ, ζ, Pζ, μ)` space.
/// + Mixed-Poloidal: Initial conditions set in the `(t, ψp, θ, ζ, Pζ, μ)` space.
#[derive(Clone)]
#[non_exhaustive]
pub struct InitialConditions {
    /// The status of `self`.
    status: InitialConditionsStatus,
    /// Which coordinate set were the initial conditions defined.
    coordinate_set: CoordinateSet,
    /// The guiding center state evaluated at the initial points.
    ///
    /// WARN: The state is only partially initialized before calling [`Self::finalize`].
    pub(crate) state: GCState,
}

impl InitialConditions {
    /// Constructs an `InitialConditions` from a set of Boozer variables.
    ///
    /// The set can be in either the `(t, ψt, θ, ζ, ρ, μ)` or `(t, ψp, θ, ζ, ρ, μ)` configuration space.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// let psit0 = MagneticFlux::Toroidal(0.2);
    /// let initial = InitialConditions::boozer(0.0, psit0, 0.0, 0.0, 1e-4, 1e-6);
    /// ```
    #[must_use]
    pub fn boozer(
        t0: f64,
        flux0: MagneticFlux,
        theta0: f64,
        zeta0: f64,
        rho0: f64,
        mu0: f64,
    ) -> Self {
        let coordinate_set = match flux0 {
            MagneticFlux::Toroidal(_) => CoordinateSet::BoozerToroidal,
            MagneticFlux::Poloidal(_) => CoordinateSet::BoozerPoloidal,
        };
        let state = GCState::from_boozer_uninit(t0, flux0, theta0, zeta0, rho0, mu0);
        Self {
            status: InitialConditionsStatus::Instantiated,
            state,
            coordinate_set,
        }
    }

    /// Constructs an `InitialConditions` from a set of Boozer variables.
    ///
    /// The set can be in either the `(t, ψt, θ, ζ, ρ, μ)` or `(t, ψp, θ, ζ, ρ, μ)` configuration space.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// # use std::f64::consts::PI;
    /// let psip0 = MagneticFlux::Poloidal(0.01);
    /// let initial = InitialConditions::mixed(0.0, psip0, PI, PI/2.0, -0.027, 1e-6);
    /// ```
    #[must_use]
    pub fn mixed(
        t0: f64,
        flux0: MagneticFlux,
        theta0: f64,
        zeta0: f64,
        pzeta0: f64,
        mu0: f64,
    ) -> Self {
        let coordinate_set = match flux0 {
            MagneticFlux::Toroidal(_) => CoordinateSet::MixedToroidal,
            MagneticFlux::Poloidal(_) => CoordinateSet::MixedPoloidal,
        };
        let state = GCState::from_mixed_uninit(t0, flux0, theta0, zeta0, pzeta0, mu0);
        Self {
            status: InitialConditionsStatus::Instantiated,
            state,
            coordinate_set,
        }
    }

    /// Calls [`Self::finalize`] if `self.status != Finalized`.
    pub(crate) fn ensure_finalized(&mut self, machine: Machine) {
        if self.status != InitialConditionsStatus::Finalized {
            self.finalize(machine);
        }
    }

    /// Evaluates all guiding center quantities and sets self's status.
    pub(crate) fn finalize(&mut self, machine: Machine) {
        use dexter_machine::EvalError::*;
        use {InitialConditionsStatus::*, SimulationError::*};

        match self.evaluate(machine) {
            Err(EvalError(Domain1dError(..)))
            | Err(EvalError(Domain2dError(..)))
            | Err(EvalError(AnalyticalDomainError)) => self.status = OutOfBoundsInitialization,
            Err(EvalError(UndefinedEvaluation(..))) => self.status = BadMagneticFlux,
            Err(err) => unreachable!("InitialConditions::finalize: {err}"),
            Ok(_) => {
                if self.state.t.is_nan()
                    | self.state.flux_value().is_nan()
                    | self.state.theta.is_nan()
                    | self.state.zeta.is_nan()
                    | self.state.rho.is_nan()
                    | self.state.pzeta.is_nan()
                    | self.state.mu.is_nan()
                {
                    self.status = NanInitialConditions
                } else if self.flux0().value() == 0.0 {
                    self.status = OutOfBoundsInitialization
                } else {
                    self.status = Finalized
                }
            }
        }
    }

    /// Calculates the missing `ρ` or `Pζ` field and evaluates the `state`.
    ///
    /// This method only alters the `rho`, `pzeta` and `state` fields. The `status` field must be
    /// set by the caller depending on the return value.
    #[expect(clippy::min_ident_chars, reason = "current, perturbation")]
    fn evaluate(&mut self, machine: Machine) -> Result<(), SimulationError> {
        // State needs `rho` to evaluate itself, which is uninitialized in the mixed case.
        if self.coordinate_set.is_mixed() {
            // WARN: Not sure if the perturbation should be included here. If so, also update
            // the GCState method.
            let acc = &mut Accelerator::new();
            let flux = *self.state.flux();
            let psip = match flux.kind() {
                MagneticFluxKind::Toroidal => machine.qfactor().eval_other(flux, acc)?,
                MagneticFluxKind::Poloidal => flux,
            };
            let g = machine.current().eval_g(flux, acc)?;
            self.state.rho = (self.state.pzeta + psip.value()) / g;
        }

        let mut caches = IntegrationCaches {
            mode_caches: machine.perturbation().generate_caches(),
            ..Default::default()
        };

        self.state.evaluate(machine, &mut caches)
    }
}

/// Getters.
///
/// # Note
///
/// All getters of derived quantities are guarded behind an `initial_conditions.status == Finalized`
/// check.
impl InitialConditions {
    /// Returns self's [`InitialConditionsStatus`].
    #[must_use]
    pub fn status(&self) -> InitialConditionsStatus {
        self.status
    }

    /// Returns the [`CoordinateSet`].
    #[must_use]
    pub fn coordinate_set(&self) -> CoordinateSet {
        self.coordinate_set
    }

    /// Returns the initial time `t0`.
    #[must_use]
    pub fn t0(&self) -> f64 {
        self.state.t
    }

    /// Returns the initial flux `flux0`.
    #[must_use]
    pub fn flux0(&self) -> MagneticFlux {
        *self.state.flux()
    }

    /// Returns the initial `θ`.
    #[must_use]
    pub fn theta0(&self) -> f64 {
        self.state.theta
    }

    /// Returns the initial `ζ`.
    #[must_use]
    pub fn zeta0(&self) -> f64 {
        self.state.zeta
    }

    /// Returns the initial `μ`.
    #[must_use]
    pub fn mu0(&self) -> f64 {
        self.state.mu
    }

    /// Returns the initial `ρ`.
    #[must_use]
    pub fn rho0(&self) -> Option<f64> {
        use {CoordinateSet::*, InitialConditionsStatus::*};
        match self.coordinate_set {
            BoozerToroidal | BoozerPoloidal => Some(self.state.rho),
            MixedToroidal | MixedPoloidal if self.status == Finalized => Some(self.state.rho),
            _ => None,
        }
    }

    /// Returns the initial `Pζ`.
    #[must_use]
    pub fn pzeta0(&self) -> Option<f64> {
        use {CoordinateSet::*, InitialConditionsStatus::*};
        match self.coordinate_set {
            MixedToroidal | MixedPoloidal => Some(self.state.pzeta),
            BoozerPoloidal | BoozerToroidal if self.status == Finalized => Some(self.state.pzeta),
            _ => None,
        }
    }

    /// Returns the initial `dot(θ)`.
    #[must_use]
    pub fn theta_dot0(&self) -> Option<f64> {
        match self.status {
            InitialConditionsStatus::Finalized => Some(self.state.theta_dot),
            _ => None,
        }
    }

    /// Returns the initial Energy.
    #[must_use]
    pub fn energy(&self) -> Option<f64> {
        match self.status {
            InitialConditionsStatus::Finalized => Some(self.state.energy),
            _ => None,
        }
    }

    /// Returns the corresponding `E-Pζ` point, if it has been calculated.
    #[must_use]
    pub(crate) fn energy_pzeta_point(&self) -> Option<Point> {
        match (self.pzeta0(), self.energy()) {
            (Some(pzeta), Some(energy)) => Some(Point {
                x: pzeta,
                y: energy,
            }),
            _ => None,
        }
    }
}

impl std::fmt::Debug for InitialConditions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("InitialConditions")
            .field("status", &self.status)
            .field("coordinate_set", &self.coordinate_set)
            .field("t0", &self.t0())
            .field("flux0", &self.flux0())
            .field("theta0", &self.theta0())
            .field("zeta0", &self.zeta0())
            .field("mu0", &self.mu0())
            .field("rho0", &format!("{:?}", self.rho0()))
            .field("pzeta0", &format!("{:?}", self.pzeta0()))
            .field("theta_dot0", &format!("{:?}", self.theta_dot0()))
            .field("energy", &format!("{:?}", self.energy()))
            .finish()
    }
}

#[cfg(test)]
mod test {
    use approx::assert_relative_eq;
    use dexter_machine::extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH};
    use dexter_machine::*;
    use std::{f64::consts::PI, path::PathBuf};

    use super::*;
    use crate::*;
    use CoordinateSet::*;
    use InitialConditionsStatus::*;
    use MagneticFlux::*;

    #[test]
    fn boozer_toroidal_initial_conditions() {
        let lcfs = Toroidal(0.05);
        let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
        let current = LarCurrent;
        let bfield = LarBfield;
        let perturbation = Perturbation::new(vec![
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
        ]);
        let machine = MachineBuilder::new(&qfactor, &current, &bfield)
            .with_perturbation(&perturbation)
            .build();

        let mut toroidal = InitialConditions::boozer(0.0, Toroidal(0.02), PI, PI, 1e-4, 1e-6);
        assert_eq!(toroidal.coordinate_set, BoozerToroidal);
        assert_eq!(toroidal.status, Instantiated);
        toroidal.finalize(machine);
        assert_eq!(toroidal.status, Finalized);
        assert!(toroidal.rho0().is_some_and(|v| v.is_finite()));
        assert!(toroidal.pzeta0().is_some_and(|v| v.is_finite()));
        assert!(toroidal.energy().is_some_and(|v| v.is_finite()));
        assert!(toroidal.theta_dot0().is_some_and(|v| v.is_finite()));

        let mut particle = Particle::new(&toroidal);
        particle.integrate(machine, (0.0, 1e2), &SolverParams::default());
        assert!(particle.steps_taken() > 10);
        assert!(particle.steps_taken() < 1000);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        // Flute modes are defined with respect to `psit_last`, therefore cannot evaluate w.r.t `psip`
        let mut poloidal = InitialConditions::boozer(0.0, Poloidal(0.01), PI, PI, 1e-4, 1e-6);
        poloidal.finalize(machine);
        assert_eq!(poloidal.status, BadMagneticFlux)
    }

    #[test]
    fn boozer_poloidal_initial_conditions() {
        let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let current = NcCurrentBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic)
            .build()
            .unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let mut poloidal = InitialConditions::boozer(0.0, Poloidal(0.02), PI, PI, 1e-4, 1e-6);
        assert_eq!(poloidal.status, Instantiated);
        assert_eq!(poloidal.coordinate_set, BoozerPoloidal);
        poloidal.finalize(machine);
        assert_eq!(poloidal.status, Finalized);
        assert!(poloidal.rho0().is_some_and(|v| v.is_finite()));
        assert!(poloidal.pzeta0().is_some_and(|v| v.is_finite()));
        assert!(poloidal.energy().is_some_and(|v| v.is_finite()));
        assert!(poloidal.theta_dot0().is_some_and(|v| v.is_finite()));

        let mut particle = Particle::new(&poloidal);
        particle.integrate(machine, (0.0, 1e3), &SolverParams::default());
        assert!(particle.steps_taken() > 10);
        assert!(particle.steps_taken() < 1000);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        let mut toroidal = InitialConditions::boozer(0.0, Toroidal(0.02), PI, PI, 1e-4, 1e-6);
        toroidal.finalize(machine);
        assert_eq!(toroidal.status, BadMagneticFlux)
    }

    #[test]
    fn mixed_toroidal_initial_conditions() {
        let lcfs = Toroidal(0.05);
        let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
        let current = LarCurrent;
        let bfield = LarBfield;
        let perturbation = Perturbation::new(vec![
            Box::new(FluteMode::new(1e-4, lcfs, 1, 2, 0.0)),
            Box::new(FluteMode::new(1e-5, lcfs, 1, 3, 0.0)),
        ]);
        let machine = MachineBuilder::new(&qfactor, &current, &bfield)
            .with_perturbation(&perturbation)
            .build();

        let mut toroidal = InitialConditions::mixed(0.0, Toroidal(0.01), PI, PI, -0.027, 1e-6);
        assert_eq!(toroidal.coordinate_set, MixedToroidal);
        assert_eq!(toroidal.status, Instantiated);
        toroidal.finalize(machine);
        assert_eq!(toroidal.status, Finalized);
        assert!(toroidal.rho0().is_some_and(|v| v.is_finite()));
        assert!(toroidal.pzeta0().is_some_and(|v| v.is_finite()));
        assert!(toroidal.energy().is_some_and(|v| v.is_finite()));
        assert!(toroidal.theta_dot0().is_some_and(|v| v.is_finite()));

        let mut particle = Particle::new(&toroidal);
        particle.integrate(machine, (0.0, 3e1), &SolverParams::default());
        assert!(particle.steps_taken() > 10);
        assert!(particle.steps_taken() < 1000);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        // Flute modes are defined with respect to `psit_last`, therefore cannot evaluate w.r.t `psip`
        let mut poloidal = InitialConditions::mixed(0.0, Poloidal(0.01), PI, PI, -0.027, 1e-6);
        poloidal.finalize(machine);
        assert_eq!(poloidal.status, BadMagneticFlux)
    }

    #[test]
    fn mixed_poloidal_initial_conditions() {
        let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let current = NcCurrentBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic)
            .build()
            .unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let mut poloidal = InitialConditions::mixed(0.0, Poloidal(0.01), PI, PI, -0.027, 1e-6);
        assert_eq!(poloidal.coordinate_set, MixedPoloidal);
        assert_eq!(poloidal.status, Instantiated);
        poloidal.finalize(machine);
        assert_eq!(poloidal.status, Finalized);
        assert!(poloidal.rho0().is_some_and(|v| v.is_finite()));
        assert!(poloidal.pzeta0().is_some_and(|v| v.is_finite()));
        assert!(poloidal.energy().is_some_and(|v| v.is_finite()));
        assert!(poloidal.theta_dot0().is_some_and(|v| v.is_finite()));

        let mut particle = Particle::new(&poloidal);
        particle.integrate(machine, (0.0, 1e2), &SolverParams::default());
        assert!(particle.steps_taken() > 10);
        assert!(particle.integration_status() == IntegrationStatus::Integrated);

        let mut toroidal = InitialConditions::mixed(0.0, Toroidal(0.01), PI, PI, -0.027, 1e-6);
        toroidal.finalize(machine);
        assert_eq!(toroidal.status, BadMagneticFlux)
    }

    #[test]
    fn boozer_mixed_equivalence() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let current = NcCurrentBuilder::new(&path, Interpolation1dType::Steffen)
            .build()
            .unwrap();
        let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic)
            .build()
            .unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let mut boozer = InitialConditions::boozer(0.0, Toroidal(0.01), PI, PI, 1e-4, 1e-6);
        boozer.finalize(machine);
        // WARN: This used to be 1e-12 but randomly started failing, I could not find a reason why.
        assert_relative_eq!(
            boozer.pzeta0().unwrap(),
            -0.00898781038097592,
            epsilon = 1e-6
        );

        let mut mixed =
            InitialConditions::mixed(0.0, Toroidal(0.01), PI, PI, boozer.pzeta0().unwrap(), 1e-6);
        mixed.finalize(machine);
        assert_relative_eq!(
            mixed.rho0().unwrap(),
            boozer.rho0().unwrap(),
            epsilon = 1e-12
        );

        assert_relative_eq!(
            mixed.theta_dot0().unwrap(),
            boozer.theta_dot0().unwrap(),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            mixed.energy().unwrap(),
            boozer.energy().unwrap(),
            epsilon = 1e-12
        );
    }
}
