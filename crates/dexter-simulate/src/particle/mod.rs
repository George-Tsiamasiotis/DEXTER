//! Representation of a charged particle.

mod close;
mod evolution;
mod initial;
mod integrate;
mod intersect;

use crate::SolverParams;
pub use initial::{CoordinateSet, InitialConditions};
pub use intersect::{IntersectParams, Intersection};

// ===============================================================================================

use ndarray::Array1;
use rsl_interpolation::{Accelerator, Cache};
use std::time::Duration;

use dexter_common::export_array1D_getter_impl;
use dexter_equilibrium::{
    Bfield, Current, FluxCommute, Harmonic, HarmonicCache, Perturbation, Qfactor,
};
use evolution::Evolution;

/// Helper enum to define an [`InitialConditions`] set with respect to one of the flux
/// coordinates.
#[derive(Clone, Copy)]
pub enum InitialFlux {
    /// Initial flux `ψ0`.
    Toroidal(f64),
    /// Initial flux `ψp0`.
    Poloidal(f64),
}

impl InitialFlux {
    /// Returns the contained value, regardless of which variant.
    #[must_use]
    pub fn value(&self) -> f64 {
        match *self {
            Self::Toroidal(value) | Self::Poloidal(value) => value,
        }
    }
}

impl std::fmt::Debug for InitialFlux {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::Toroidal(psi0) => write!(f, "ψ0: {psi0}"),
            Self::Poloidal(psip0) => write!(f, "ψp0: {psip0}"),
        }
    }
}

// ===============================================================================================

/// Simple container for the equilibrium objects, to make passing them as parameters a bit easier.
pub(crate) struct EqObjects<'eq, Q, C, B, H>
where
    Q: Qfactor,
    C: Current,
    B: Bfield,
    H: Harmonic,
{
    /// The current configuration's [`Qfactor`].
    pub(crate) qfactor: &'eq Q,
    /// The current configuration's [`Current`].
    pub(crate) current: &'eq C,
    /// The current configuration's [`Bfield`].
    pub(crate) bfield: &'eq B,
    /// The current configuration's [`Perturbation`].
    pub(crate) perturbation: &'eq Perturbation<H>,
}

// ===============================================================================================

/// Container for the caching objects needed for the evaluations.
#[derive(Default, Debug)]
pub(crate) struct IntegrationCaches<C: HarmonicCache> {
    /// The `ψ` Accelerator. Only used when integrating with respect to `ψ`.
    pub(crate) psi_acc: Accelerator,
    /// The `ψp` Accelerator. Only used when integrating with respect to `ψp`.
    pub(crate) psip_acc: Accelerator,
    /// The `θ` Accelerator.
    pub(crate) theta_acc: Accelerator,
    /// The 2D Interpolation cache.
    pub(crate) spline_cache: Cache<f64>,
    /// The caches of the perturbation's harmonics.
    pub(crate) harmonic_caches: Vec<C>,
}

// ===============================================================================================

/// A [`Particle`]'s integration status.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntegrationStatus {
    /// Initialized by [`InitialConditions`], not integrated.
    Initialized,
    /// [`InitialConditions`] have not been fully calculated yet.
    PartlyInitialized,
    /// Invalid [`InitialConditions`]. May occur when using Mixed variables with objects that
    /// cannot define them.
    InvalidInitialConditions,
    /// [`InitialConditions`] were out of bounds.
    OutOfBoundsInitialization,
    /// Reached the end of the integration successfully.
    Integrated,
    /// Intersections calculation successful.
    Intersected,
    /// Integrated for a certain amount of `θ-ψ` periods.
    ClosedPeriods(usize),
    /// Escaped the last closed flux surface (LCFS).
    Escaped,
    /// Escaped when performing a step on the modified system.
    ///
    /// This indicates that something is wrong in Hénon's trick implementation.
    ModStateEscaped,
    /// Calculated some intersections correctly but also timed out.
    IntersectedTimedOut,
    /// Calculated invalid intersections.
    InvalidIntersections,
    /// Timed out after a maximum number of steps.
    TimedOut(Duration),
    /// Simulation failed for unknown reasons.
    Failed(Box<str>),
}

// ===============================================================================================

/// Accelerator/Cache stats of an integration routine.
#[derive(Debug, Default, Clone)]
pub struct ParticleCacheStats {
    /// The final state `ψ` Accelerator.
    pub psi_acc: Accelerator,
    /// The final state `ψp` Accelerator.
    pub psip_acc: Accelerator,
    /// The final state `θ` Accelerator.
    pub theta_acc: Accelerator,
    /// The sum of of the individual harmonic caches' hits.
    pub harmonic_cache_hits: usize,
    /// The sum of of the individual harmonic caches' misses.
    pub harmonic_cache_misses: usize,
}

// ===============================================================================================

/// A particle's calculated frequencies and `qkinetic`.
#[derive(Default, Clone)]
pub struct Frequencies {
    /// The particle's calculated `ωθ`.
    pub omega_theta: Option<f64>,
    /// The particle's calculated `ωζ`.
    pub omega_zeta: Option<f64>,
    /// The particle's calculated `qkinetic`.
    pub qkinetic: Option<f64>,
}

// ===============================================================================================
// ===============================================================================================

/// Representation of a charged particle.
#[derive(Clone)]
pub struct Particle {
    /// The [`InitialConditions`] set of the particle.
    initial_conditions: InitialConditions,
    /// Status of the particle's integration.
    integration_status: IntegrationStatus,
    /// The time evolution of the particle.
    evolution: Evolution,
    /// Stats about the particle's integration.
    stats: ParticleCacheStats,
    /// The particle's calculated `ωθ`, `ωζ` and `qkinetic`.
    frequencies: Frequencies,
    /// The particle's initial Energy in Normalized Units. The Energy depends both on the initial
    /// conditions and the equilibrium.
    initial_energy: Option<f64>,
    /// The particle's energy after the integration.
    final_energy: Option<f64>,
}

impl Particle {
    /// Creates a new [`Particle`] from a set of [`InitialConditions`].
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use dexter_simulate::*;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// # let psip_last = geometry.psip_last().unwrap();
    /// let psip0 = InitialFlux::Poloidal(0.5 * psip_last);
    /// let initial = InitialConditions::boozer(0.0, psip0, 0.0, 3.14, 1e-5, 1e-6);
    /// let mut particle = Particle::new(&initial);
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    #[must_use]
    pub fn new(initial_conditions: &InitialConditions) -> Self {
        use CoordinateSet::*;
        let integration_status = match initial_conditions.coordinate_set() {
            BoozerToroidal | BoozerPoloidal => IntegrationStatus::Initialized, // Always succeeds
            MixedToroidal | MixedPoloidal => IntegrationStatus::PartlyInitialized, // Needs `finalize()`
        };
        Self {
            initial_conditions: initial_conditions.to_owned(),
            integration_status,
            evolution: Evolution::default(),
            frequencies: Frequencies::default(),
            stats: ParticleCacheStats::default(),
            initial_energy: None,
            final_energy: None,
        }
    }
}

// Routines
impl Particle {
    /// Integrates the particle for a certain time interval.
    ///
    /// The time interval is in Normalized Units (inverse gyro-frequency).
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// # use std::path::PathBuf;
    /// #
    /// let path = PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// let lcfs = LastClosedFluxSurface::Toroidal(qfactor.psi_last().unwrap());
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, lcfs, 1, 1, 0.0),
    ///     CosHarmonic::new(1e-3, lcfs, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, lcfs, 1, 3, 0.0),
    /// ]);
    ///
    /// let psip0 = InitialFlux::Poloidal(0.015);
    /// let initial = InitialConditions::boozer(0.0, psip0, 0.0, 3.14, 1e-5, 1e-6);
    /// let mut particle = Particle::new(&initial);
    /// particle.integrate(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     (0.0, 1e2),
    ///     &SolverParams::default(),
    /// );
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    pub fn integrate<Q, C, B, H>(
        &mut self,
        qfactor: &Q,
        current: &C,
        bfield: &B,
        perturbation: &Perturbation<H>,
        teval: (f64, f64),
        solver_params: &SolverParams,
    ) where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let objects = EqObjects {
            qfactor,
            current,
            bfield,
            perturbation,
        };
        integrate::integrate(self, &objects, teval, solver_params);
    }

    /// Integrates the particle, calculating its intersections with a constant `θ` or `ζ` surface.
    ///
    /// The intersection surface, angle, and number of turns are configured with the helper struct
    /// [`IntersectParams`].
    ///
    /// Using the method described by [`Hénon`] we can force the solver to step exactly on the
    /// intersection surface.
    ///
    /// The differences between two consecutive values of the corresponding angle variable are
    /// guaranteed to be `2π +- ε`, where ε a number smaller than the solver's relative tolerance.
    ///
    /// [`Hénon`]: https://www.sciencedirect.com/science/article/abs/pii/0167278982900343
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// # use std::path::PathBuf;
    /// #
    /// let qfactor = UnityQfactor::new();
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::zero();
    ///
    /// let psi0 = InitialFlux::Toroidal(0.02);
    /// let initial = InitialConditions::boozer(0.0, psi0, 3.14, 0.0, 1e-4, 1e-6);
    /// let intersect_params = IntersectParams::new(Intersection::ConstTheta, 3.14, 10);
    ///
    /// let mut particle = Particle::new(&initial);
    /// particle.intersect(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     &intersect_params,
    ///     &SolverParams::default(),
    /// );
    ///
    /// assert_eq!(particle.steps_stored(), 10);
    /// # assert!(matches!(particle.integration_status(), IntegrationStatus::Intersected));
    /// # Ok::<_, SimulationError>(())
    ///
    /// ```
    pub fn intersect<Q, C, B, H>(
        &mut self,
        qfactor: &Q,
        current: &C,
        bfield: &B,
        perturbation: &Perturbation<H>,
        intersect_params: &IntersectParams,
        solver_params: &SolverParams,
    ) where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let objects = EqObjects {
            qfactor,
            current,
            bfield,
            perturbation,
        };
        intersect::intersect(self, &objects, intersect_params, solver_params);
    }

    /// Integrates the particle, for `periods` number of `θ-ψ` periods.
    ///
    /// # Example
    ///
    /// TODO:
    pub fn close<Q, C, B, H>(
        &mut self,
        qfactor: &Q,
        current: &C,
        bfield: &B,
        perturbation: &Perturbation<H>,
        periods: usize,
        solver_params: &SolverParams,
    ) where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let objects = EqObjects {
            qfactor,
            current,
            bfield,
            perturbation,
        };
        close::close(self, &objects, periods, solver_params);
    }
}

// Getters
impl Particle {
    /// Returns the Particle's [`InitialConditions`].
    #[must_use]
    pub fn initial_conditions(&self) -> InitialConditions {
        self.initial_conditions.clone()
    }

    /// Returns the Particle's [`IntegrationStatus`].
    #[must_use]
    pub fn integration_status(&self) -> IntegrationStatus {
        self.integration_status.clone()
    }

    /// Returns the total number of steps taken during the integration.
    ///
    /// This number is not necessarily the same as the number of steps stored.
    #[must_use]
    pub fn steps_taken(&self) -> usize {
        self.evolution.steps_taken
    }

    /// Returns the number of steps stored in the time series arrays.
    #[must_use]
    pub fn steps_stored(&self) -> usize {
        self.evolution.steps_stored()
    }

    /// Returns the duration of the integration routine.
    #[must_use]
    pub fn duration(&self) -> Duration {
        self.evolution.duration
    }

    /// Returns the particle's initial energy.
    #[must_use]
    pub fn initial_energy(&self) -> Option<f64> {
        self.initial_energy
    }

    /// Returns the particle's final energy.
    #[must_use]
    pub fn final_energy(&self) -> Option<f64> {
        self.final_energy
    }

    /// Returns the variance of the energy array.
    #[must_use]
    pub fn energy_var(&self) -> Option<f64> {
        self.evolution.energy_var()
    }

    /// Returns the particle's [`ParticleCacheStats`].
    #[must_use]
    pub fn cache_stats(&self) -> ParticleCacheStats {
        self.stats.clone()
    }

    /// Returns the particle's [`Frequencies`].
    #[must_use]
    pub fn frequencies(&self) -> Frequencies {
        self.frequencies.clone()
    }

    /// Returns the particle's calculated `ωθ`.
    #[must_use]
    pub fn omega_theta(&self) -> Option<f64> {
        self.frequencies.omega_theta
    }

    /// Returns the particle's calculated `ωζ`.
    #[must_use]
    pub fn omega_zeta(&self) -> Option<f64> {
        self.frequencies.omega_zeta
    }

    /// Returns the particle's calculated `qkinetic`.
    #[must_use]
    pub fn qkinetic(&self) -> Option<f64> {
        self.frequencies.qkinetic
    }

    /// Prints the particle's integration interpolation caches, [`Cache`] and
    /// [`Accelerator`].
    pub fn print_cache_stats(&self) {
        println!("{:#}", self.stats);
    }

    export_array1D_getter_impl!(t_array, evolution, t);
    export_array1D_getter_impl!(psi_array, evolution, psi_array);
    export_array1D_getter_impl!(psip_array, evolution, psip_array);
    export_array1D_getter_impl!(theta_array, evolution, theta_array);
    export_array1D_getter_impl!(zeta_array, evolution, zeta_array);
    export_array1D_getter_impl!(rho_array, evolution, rho_array);
    export_array1D_getter_impl!(mu_array, evolution, mu_array);
    export_array1D_getter_impl!(ptheta_array, evolution, ptheta_array);
    export_array1D_getter_impl!(pzeta_array, evolution, pzeta_array);
    export_array1D_getter_impl!(energy_array, evolution, energy_array);
}

impl std::fmt::Debug for Frequencies {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fn stringify(field: Option<f64>) -> String {
            if let Some(value) = field {
                format!("{value:.7}")
            } else {
                String::from("Not calculated")
            }
        }

        f.debug_struct("Frequencies")
            .field("omega_theta", &stringify(self.omega_theta))
            .field("omega_zeta", &stringify(self.omega_zeta))
            .field("qkinetic", &stringify(self.qkinetic))
            .finish()
    }
}

impl std::fmt::Debug for Particle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle")
            .field("initial conditions", &self.initial_conditions)
            .field("integration status", &self.integration_status)
            .field("evolution", &self.evolution)
            .field("frequencies", &self.frequencies)
            .field("initial energy", &self.initial_energy.unwrap_or(f64::NAN))
            .field("final energy  ", &self.final_energy.unwrap_or(f64::NAN))
            .field("energy variance", &self.energy_var().unwrap_or(f64::NAN))
            .finish()
    }
}

impl std::fmt::Display for ParticleCacheStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Particle Interpolation Caching stats")
            .field("ψ Accelerator", &self.psi_acc)
            .field("ψp Accelerator", &self.psip_acc)
            .field("θ Accelerator", &self.theta_acc)
            .field("total harmonic cache hits", &self.harmonic_cache_hits)
            .field("total harmonic cache misses", &self.harmonic_cache_misses)
            .finish()
    }
}
