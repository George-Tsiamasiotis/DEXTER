//! Representation of a particle

use std::time::Duration;

use derive_is_enum_variant::is_enum_variant as IsEnumVariant;
use equilibrium::{Bfield, Current, Perturbation, Qfactor};

use crate::routines::{Frequencies, close_theta_period, integrate, map_integrate};
use crate::state::Display;
use crate::{Evolution, IntegrationConfig, MappingParameters, SinglePeriodConfig, State};
use crate::{MappingConfig, ParticleError};

/// A set of a Particle's intial conditions.
#[derive(Clone, Debug)]
pub struct InitialConditions {
    /// The initial time.
    pub time0: f64,
    /// The initial `θ` angle.
    pub theta0: f64,
    /// The intial poloidal magnetic flux `ψp`.
    pub psip0: f64,
    /// The initial parallel gyro radius `ρ`.
    pub rho0: f64,
    /// The initial `ζ` angle.
    pub zeta0: f64,
    /// The magnetic moment `μ`.
    pub mu: f64,
}

/// The [`Particle`]'s integration status.
#[derive(Debug, Clone, Default, PartialEq, IsEnumVariant)]
pub enum IntegrationStatus {
    /// Initialized from [`InitialConditions`], not integrated.
    #[default]
    Initialized,
    /// Reached the end of the integration successfully.
    Integrated,
    /// Reached the end of the mapping successfully.
    Mapped,
    /// Integrated for 1 period of `θ-ψp`.
    SinglePeriodIntegrated,
    /// Escaped / Hit the wall.
    Escaped,
    /// NaN encountered during State evaluation.
    EvaluationNan,
    /// Timed out after a maximum number of steps.
    TimedOut(Duration),
    /// Intersections calculated from the mapping are invalid (The spacing between each
    /// intersection and its neighbors must be *exactly* 2π).
    InvalidIntersections,
    /// Integration/Mapping failed for unknown reasons.
    Failed(Box<str>),
}

/// Defines the Particle's orbit type from its θ-span.
#[derive(Debug, Default, Clone, IsEnumVariant)]
pub enum OrbitType {
    #[default]
    Undefined,
    Trapped,
    Passing,
}

/// Representation of a particle.
#[derive(Clone)]
pub struct Particle {
    /// The [`InitialConditions`] set of the particle.
    pub initial_conditions: InitialConditions,
    /// The initial [`State`] of the particle.
    pub(crate) initial_state: State,
    /// The final [`State`] of the particle.
    pub(crate) final_state: State,
    /// The time [`Evolution`] of the particle.
    pub evolution: Evolution,
    /// Status of the particle's integration.
    pub status: IntegrationStatus,
    /// The orbit type.
    pub orbit_type: OrbitType,
    /// The particle's `ωθ`, `ωζ` and `qkinetic`.
    pub frequencies: Frequencies,
}

impl Particle {
    /// Creates a new [`Particle`] from a set of [`InitialConditions`].
    ///
    /// Example
    ///
    /// ```
    /// # use std::path::PathBuf;
    /// # use particle::*;
    /// # use equilibrium::geometries::NcGeometryBuilder;
    /// #
    /// # let path = PathBuf::from("../equilibrium/lar_netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// let initial_conditions = InitialConditions {
    ///     time0: 0.0,
    ///     theta0: 0.0,
    ///     psip0: geometry.psip_wall() / 2.0,
    ///     rho0: 1e-4,
    ///     zeta0: 0.0,
    ///     mu: 0.0,
    /// };
    /// let mut particle = Particle::new(&initial_conditions);
    /// # Ok::<_, ParticleError>(())
    ///
    /// ```
    pub fn new(initial_conditions: &InitialConditions) -> Self {
        let initial_state = State::from_initial(initial_conditions);
        let mut evolution = Evolution::new();
        evolution.push_state(&initial_state);

        Self {
            initial_conditions: initial_conditions.to_owned(),
            evolution,
            initial_state,
            final_state: State::default(),
            status: IntegrationStatus::default(),
            orbit_type: OrbitType::default(),
            frequencies: Frequencies::default(),
        }
    }

    /// Integrates the particle, storing the calculated orbit in [`Particle::evolution`].
    ///
    /// Example
    ///
    /// ```
    /// # use std::path::PathBuf;
    /// # use particle::*;
    /// # use equilibrium::{geometries::*, qfactors::*, currents::*, bfields::*, harmonics::*, perturbations::*};
    /// #
    /// # let path = PathBuf::from("../equilibrium/lar_netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// # let perturbation = NcPerturbation::from_harmonics(&vec![
    /// #   NcHarmonicBuilder::new(&path, "steffen", 2, 1).build()?
    /// # ]);
    /// #
    /// # let initial_conditions = InitialConditions {
    /// #   time0: 0.0,
    /// #   theta0: 0.0,
    /// #   psip0: geometry.psip_wall() / 2.0,
    /// #   rho0: 1e-4,
    /// #   zeta0: 0.0,
    /// #   mu: 0.0,
    /// # };
    /// # let mut particle = Particle::new(&initial_conditions);
    /// #
    /// let config = IntegrationConfig {
    ///     method: SteppingMethod::EnergyAdaptiveStep,
    ///     max_steps: 10_000_000,
    ///     energy_rel_tol: 1e-12,
    ///     energy_abs_tol: 1e-13,
    ///     ..Default::default()
    /// };
    ///
    /// particle.integrate(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     (0.0, 1e6),
    ///     &config,
    /// );
    /// # Ok::<_, ParticleError>(())
    /// ```
    pub fn integrate(
        &mut self,
        qfactor: &impl Qfactor,
        current: &impl Current,
        bfield: &impl Bfield,
        perturbation: &impl Perturbation,
        t_eval: (f64, f64),
        config: &IntegrationConfig,
    ) {
        match integrate(self, qfactor, current, bfield, perturbation, t_eval, config) {
            Ok(()) => self.status = IntegrationStatus::Integrated,
            Err(error) => self.set_status_from_error(error),
        }
    }

    /// Integrates the particle, storing its intersections with the Poincare surface defined by
    /// [`MappingParameters`].
    ///
    /// Example
    ///
    /// ```
    /// # use std::path::PathBuf;
    /// # use particle::*;
    /// # use equilibrium::{geometries::*, qfactors::*, currents::*, bfields::*, harmonics::*, perturbations::*};
    /// #
    /// # let path = PathBuf::from("../equilibrium/lar_netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// # let perturbation = NcPerturbation::from_harmonics(&vec![
    /// #   NcHarmonicBuilder::new(&path, "steffen", 2, 1).build()?
    /// # ]);
    /// #
    /// # let initial_conditions = InitialConditions {
    /// #   time0: 0.0,
    /// #   theta0: 0.0,
    /// #   psip0: geometry.psip_wall() / 2.0,
    /// #   rho0: 1e-4,
    /// #   zeta0: 0.0,
    /// #   mu: 0.0,
    /// # };
    /// # let mut particle = Particle::new(&initial_conditions);
    /// #
    /// let config = MappingConfig {
    ///     method: SteppingMethod::ErrorAdaptiveStep,
    ///     max_steps: 100_000,
    ///     error_rel_tol: 1e-12,
    ///     error_abs_tol: 1e-13,
    ///     ..Default::default()
    /// };
    ///
    /// let params = MappingParameters::new(PoincareSection::ConstTheta, 0.0, 10);
    ///
    /// particle.map(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     &params,
    ///     &config,
    /// );
    /// # Ok::<_, ParticleError>(())
    /// ```
    pub fn map(
        &mut self,
        qfactor: &impl Qfactor,
        current: &impl Current,
        bfield: &impl Bfield,
        perturbation: &impl Perturbation,
        params: &MappingParameters,
        config: &MappingConfig,
    ) {
        match map_integrate(self, qfactor, current, bfield, perturbation, params, config) {
            Ok(()) => self.status = IntegrationStatus::Mapped,
            Err(error) => self.set_status_from_error(error),
        }
    }

    /// Integrates the particle for 1 `θ-ψp` period,  calculating its `ωθ`, `ωζ` and qkinetic.
    ///
    /// The orbit is stored in [`Particle::evolution`].
    ///
    /// Example
    ///
    /// ```
    /// # use std::path::PathBuf;
    /// # use particle::*;
    /// # use equilibrium::{geometries::*, qfactors::*, currents::*, bfields::*, harmonics::*, perturbations::*};
    /// #
    /// # let path = PathBuf::from("../equilibrium/lar_netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// # let perturbation = NcPerturbation::from_harmonics(&vec![
    /// #   NcHarmonicBuilder::new(&path, "steffen", 2, 1).build()?
    /// # ]);
    /// #
    /// # let initial_conditions = InitialConditions {
    /// #   time0: 0.0,
    /// #   theta0: 0.0,
    /// #   psip0: geometry.psip_wall() / 2.0,
    /// #   rho0: 1e-4,
    /// #   zeta0: 0.0,
    /// #   mu: 0.0,
    /// # };
    /// # let mut particle = Particle::new(&initial_conditions);
    /// #
    /// particle.single_period_integrate(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     &SinglePeriodConfig::default(),
    /// );
    /// # Ok::<_, ParticleError>(())
    /// ```
    pub fn single_period_integrate(
        &mut self,
        qfactor: &impl Qfactor,
        current: &impl Current,
        bfield: &impl Bfield,
        perturbation: &impl Perturbation,
        config: &SinglePeriodConfig,
    ) {
        match close_theta_period(self, qfactor, current, bfield, perturbation, config) {
            Ok(()) => self.status = IntegrationStatus::SinglePeriodIntegrated,
            Err(error) => self.set_status_from_error(error),
        }
    }

    /// Calculates the Particles OrbitType.
    pub(crate) fn calculate_orbit_type(&mut self) {
        // TODO: Decide how to setup up parameters
    }

    /// Sets the Particle's [`IntegrationStatus`] from a Result::Err() of an integration
    /// routine.
    pub(crate) fn set_status_from_error(&mut self, error: ParticleError) {
        use ParticleError::*;
        self.status = match error {
            EqError(..) => IntegrationStatus::Escaped,
            TimedOut(duration) => {
                self.evolution.duration = duration;
                IntegrationStatus::TimedOut(duration)
            }
            IntersectionError => IntegrationStatus::InvalidIntersections,
            EvaluationNaN => IntegrationStatus::EvaluationNan,
            // _ => IntegrationStatus::Failed(error.to_string().into_boxed_str()),
        }
    }
}

impl Particle {
    /// Returns the initial energy of the Particle, calculated from its initial State.
    pub fn initial_energy(&self) -> f64 {
        self.initial_state.energy()
    }

    /// Returns the final energy of the Particle, calculated from its final State.
    pub fn final_energy(&self) -> f64 {
        self.final_state.energy()
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("ψ-acc", &self.final_state.xacc)
            .field("θ-acc", &self.final_state.yacc)
            .field("hcache", &self.final_state.hcache.first().or(None))
            .field("Initial", &Display::from_state(&self.initial_state))
            .field(
                "Initial parallel energy",
                &self.initial_state.parallel_energy(),
            )
            .field(
                "Intial perpendicular energy",
                &self.initial_state.perpendicular_energy(),
            )
            .field("Initial energy", &self.initial_state.energy())
            .field("Final energy  ", &self.final_state.energy())
            .field("Status", &self.status)
            .field("Orbit type", &self.orbit_type)
            .field("Evolution", &self.evolution)
            .field("Frequencies", &self.frequencies)
            .finish()
    }
}
