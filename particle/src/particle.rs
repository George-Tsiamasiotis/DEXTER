//! Representation of a particle

use std::time::Duration;

use config::*;
use derive_is_enum_variant::is_enum_variant;
use equilibrium::{Bfield, Currents, Perturbation, Qfactor};

use crate::routines::{Frequencies, close_theta_period, integrate, map_integrate};
use crate::state::Display;
use crate::{Evolution, MappingParameters, State};

use crate::MagneticMoment;
use crate::{Flux, Length, ParticleError, Radians, Time};

/// A set of a Particle's intial conditions.
#[derive(Clone, Debug)]
pub struct InitialConditions {
    /// The initial time.
    pub time0: Time,
    /// The initial `θ` angle.
    pub theta0: Radians,
    /// The intial poloidal magnetic flux `ψp`.
    pub psip0: Flux,
    /// The initial parallel gyro radius `ρ`.
    pub rho0: Length,
    /// The initial `ζ` angle.
    pub zeta0: Flux,
    /// The magnetic moment `μ`.
    pub mu: MagneticMoment,
}

/// The [`Particle`]'s integration status.
#[derive(Debug, Clone, Default, PartialEq, is_enum_variant)]
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
    /// Timed out after [`config::MAX_STEPS`]
    TimedOut(Duration),
    /// Intersections calculated from the mapping are invalid (The spacing between each
    /// intersection and its neighbors must be *exactly* 2π).
    InvalidIntersections,
    /// Integration/Mapping failed for unknown reasons.
    Failed(Box<str>),
}

/// Defines the Particle's orbit type from its θ-span.
#[derive(Debug, Default, Clone, is_enum_variant)]
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
    pub initial_state: State,
    /// The final [`State`] of the particle.
    pub final_state: State,
    /// The [`Evolution`] time series of the particle.
    pub evolution: Evolution,
    /// Status of the particle's integration.
    pub status: IntegrationStatus,
    /// The orbit type.
    pub orbit_type: OrbitType,
    /// The particle's `ωθ`, `ωζ` and `qkinetic`.
    pub frequencies: Frequencies,
}

impl Particle {
    /// Creates a new [`Particle`] from the initial conditions.
    pub fn new(initial_conditions: &InitialConditions) -> Self {
        let initial_state = State::from_initial(initial_conditions);
        let mut evolution = Evolution::with_capacity(EVOLUTION_INIT_CAPACITY);
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

    /// Integrates the particle, storing the calculated orbit in [`Evolution`].
    pub fn integrate(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        currents: &Currents,
        perturbation: &Perturbation,
        t_eval: (Time, Time),
    ) {
        match integrate(self, qfactor, bfield, currents, perturbation, t_eval) {
            Ok(()) => self.status = IntegrationStatus::Integrated,
            Err(error) => self.set_status_from_error(error),
        }
    }

    /// Integrates the particle, storing its intersections with the Poincare surface defined by
    /// [`MappingParameters`].
    pub fn map(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        currents: &Currents,
        perturbation: &Perturbation,
        params: &MappingParameters,
    ) {
        match map_integrate(self, qfactor, bfield, currents, perturbation, params) {
            Ok(()) => self.status = IntegrationStatus::Mapped,
            Err(error) => self.set_status_from_error(error),
        }
    }

    /// Integrates the particle for 1 `θ-ψp` period,  calculating its `ωθ`, `ωζ` and qkinetic.
    pub fn calculate_frequencies(
        &mut self,
        qfactor: &Qfactor,
        bfield: &Bfield,
        currents: &Currents,
        perturbation: &Perturbation,
    ) {
        match close_theta_period(self, qfactor, bfield, currents, perturbation) {
            Ok(()) => self.status = IntegrationStatus::SinglePeriodIntegrated,
            Err(error) => self.set_status_from_error(error),
        }
    }

    /// Calculates the Particles OrbitType.
    pub(crate) fn calculate_orbit_type(&mut self) {
        let thetas = &self.evolution.theta;
        if thetas.is_empty() {
            self.orbit_type = OrbitType::Undefined;
            return;
        }
        let theta0 = thetas.first().expect("non empty");
        let thetaf = thetas.last().expect("non empty");

        use std::f64::consts::TAU;
        match (thetaf - theta0) > TAU - TRAPPED_THRESHOLD {
            true => self.orbit_type = OrbitType::Passing,
            false => self.orbit_type = OrbitType::Trapped,
        }
    }

    /// Sets the Particle's [`IntegrationStatus`] from a Result::Err() of an integration
    /// routine.
    pub(crate) fn set_status_from_error(&mut self, error: ParticleError) {
        use ParticleError::*;
        self.status = match error {
            EqError(..) => IntegrationStatus::Escaped, // FIXME:
            TimedOut(duration) => IntegrationStatus::TimedOut(duration),
            IntersectionError => IntegrationStatus::InvalidIntersections,
            _ => IntegrationStatus::Failed(error.to_string().into_boxed_str()),
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
