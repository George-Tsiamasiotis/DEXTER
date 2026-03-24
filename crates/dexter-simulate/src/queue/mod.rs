//! Container type for batch handling Particles.

mod initials;
mod pbars;
mod stats;

pub use initials::QueueInitialConditions;
pub use initials::{poloidal_fluxes, toroidal_fluxes};

use pbars::{IntegratePbar, IntersectPbar};
use stats::QueueStats;

use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::ops::Index;
use std::slice::Iter;

use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, Perturbation, Qfactor};

use crate::{IntersectParams, Particle, SolverParams};

/// Indicates the routine Queue's particles executed.
#[derive(Default, Clone, Debug, PartialEq, Eq)]
pub enum Routine {
    /// [`Queue`] has not run any routine yet.
    #[default]
    None,
    /// [`Queue`] has run the [`Particle::integrate`] routine.
    Integrate,
    /// [`Queue`] has run the [`Particle::intersect`] routine.
    Intersect,
}

/// A collection of multiple [`Particles`](Particle), constructed from a
/// [`QueueInitialConditions`].
///
/// Offers the ability to batch [`Particle::integrate`] or [`Particle::intersect`] all the
/// contained particles, using multiple threads.
pub struct Queue {
    /// The sets of initial conditions.
    initial_conditions: QueueInitialConditions,
    /// The contained particles.
    particles: Vec<Particle>,
    /// Helper struct to keep track of a Queue's run statistics.
    stats: QueueStats,
    /// The routine Queue's particles executed.
    routine: Routine,
}

impl Queue {
    /// Creates a [`Queue`] from a [`QueueInitialConditions`].
    ///
    /// # Example
    ///
    /// Initial conditions in Boozer coordinates.
    /// ```
    /// # use dexter_simulate::*;
    /// use InitialFlux::*;
    /// let initial_conditions = QueueInitialConditions::boozer(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[1e-3, 2e-3],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// let queue = Queue::new(&initial_conditions);
    /// # Ok::<_, SimulationError>(())
    /// ```
    ///
    /// # Example
    ///
    /// Initial conditions in Mixed coordinates.
    /// ```
    /// # use dexter_simulate::*;
    /// use InitialFlux::*;
    /// let initial_conditions = QueueInitialConditions::mixed(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[-0.025, -0.027],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// let queue = Queue::new(&initial_conditions);
    /// # Ok::<_, SimulationError>(())
    /// ```
    #[must_use]
    pub fn new(initial_conditions: &QueueInitialConditions) -> Self {
        Self {
            initial_conditions: initial_conditions.clone(),
            particles: initial_conditions.to_particles(),
            stats: QueueStats::from_initial_conditions(initial_conditions),
            routine: Routine::None,
        }
    }

    /// Integrates all the contained particles for a specific time interval.
    ///
    /// # Example
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
    /// use InitialFlux::*;
    /// let initial_conditions = QueueInitialConditions::boozer(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[1e-4, 2e-4],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// let mut queue = Queue::new(&initial_conditions);
    /// queue.integrate(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     (0.0, 1e3),
    ///     &SolverParams::default(),
    /// );
    /// # Ok::<_, SimulationError>(())
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
        Q: Qfactor + FluxCommute + Send + Sync,
        C: Current + Send + Sync,
        B: Bfield + Send + Sync,
        H: Harmonic + Send + Sync,
    {
        let pbar = IntegratePbar::new(self);
        pbar.print_prelude();

        self.particles.par_iter_mut().for_each(|particle| {
            particle.integrate(qfactor, current, bfield, perturbation, teval, solver_params);
            pbar.inc(&particle.integration_status());
            pbar.print_stats();
        });
        pbar.finish();

        self.routine = Routine::Integrate;
        self.stats = QueueStats::from_completed_queue(self);
    }

    /// Integrates all the contained particle, calculating their intersections with a constant θ or ζ surface.
    ///
    /// Otherwise known as a `Poincare map`.
    ///
    /// The intersection surface, angle, and number of turns are configured with the helper struct [`IntersectParams`].
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// # use std::path::PathBuf;
    /// #
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.6);
    /// let qfactor = ParabolicQfactor::new(1.1, 4.2, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, lcfs, 1, 1, 0.0),
    ///     CosHarmonic::new(1e-3, lcfs, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, lcfs, 1, 3, 0.0),
    /// ]);
    ///
    /// use InitialFlux::*;
    /// let initial_conditions = QueueInitialConditions::boozer(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[1e-4, 2e-4],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// let mut queue = Queue::new(&initial_conditions);
    /// let intersect_params = IntersectParams::new(Intersection::ConstTheta, 0.0, 100);
    /// queue.intersect(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     &intersect_params,
    ///     &SolverParams::default(),
    /// );
    /// # Ok::<_, SimulationError>(())
    /// ```
    #[doc(alias = "poincare_map")]
    pub fn intersect<Q, C, B, H>(
        &mut self,
        qfactor: &Q,
        current: &C,
        bfield: &B,
        perturbation: &Perturbation<H>,
        intersect_params: &IntersectParams,
        solver_params: &SolverParams,
    ) where
        Q: Qfactor + FluxCommute + Send + Sync,
        C: Current + Send + Sync,
        B: Bfield + Send + Sync,
        H: Harmonic + Send + Sync,
    {
        let pbar = IntersectPbar::new(self, intersect_params);
        pbar.print_prelude();

        self.particles.par_iter_mut().for_each(|particle| {
            particle.intersect(
                qfactor,
                current,
                bfield,
                perturbation,
                intersect_params,
                solver_params,
            );
            pbar.inc(&particle.integration_status());
            pbar.print_stats();
        });
        pbar.finish();

        self.routine = Routine::Intersect;
        self.stats = QueueStats::from_completed_queue(self);
    }
}

/// Getters.
impl Queue {
    /// Returns the [`Queue`]'s [`QueueInitialConditions`].
    #[must_use]
    pub fn initial_conditions(&self) -> QueueInitialConditions {
        self.initial_conditions.clone()
    }

    /// Returns the [`Queue`]'s [`Routine`].
    #[must_use]
    pub fn routine(&self) -> Routine {
        self.routine.clone()
    }

    /// Returns the [`Queue`]'s number of contained [`Particle`]s.
    #[must_use]
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }

    /// Returns an iterator over the Queue's [`Particle`]s.
    ///
    /// The iterator yields all items from start to end.
    pub fn iter(&self) -> Iter<'_, Particle> {
        self.particles.iter()
    }

    /// Returns a [`Vec`] with all the contained particles.
    #[must_use]
    pub fn particles(&self) -> Vec<Particle> {
        self.particles.clone()
    }
}

impl Index<usize> for Queue {
    type Output = Particle;

    fn index(&self, index: usize) -> &Self::Output {
        self.particles.index(index)
    }
}

impl std::fmt::Debug for Queue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.stats.fmt(f)
    }
}
