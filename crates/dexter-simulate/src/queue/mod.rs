//! Container type for batch handling Particles.

mod initials;
mod pbars;
mod stats;

pub use initials::QueueInitialConditions;
pub use initials::{poloidal_fluxes, toroidal_fluxes};

use ndarray::Array1;
use ndarray_stats::QuantileExt;
use pbars::{ClosePbar, IntegratePbar, IntersectPbar};
use stats::QueueStats;

use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::ops::{Index, Range};
use std::slice::Iter;
use std::time::Duration;

use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, Perturbation, Qfactor};

use crate::queue::pbars::ClassifyPbar;
use crate::{EnergyPzetaPlane, IntersectParams, Particle, SolverParams};
use crate::{EnergyPzetaPosition, OrbitType};

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
    /// [`Queue`] has run the [`Particle::close`] routine.
    Close,
    /// [`Queue`] has run the [`Particle::classify`] routine.
    Classify,
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

    /// Creates a [`Queue`] from a slice of [`Particles`](Particle).
    ///
    /// The particles do not have to be initialized, but they must be defined on the same coordinate
    /// set (boozer/mixed).
    ///
    /// # Panics
    ///
    /// Panics if the slice is empty.
    ///
    /// # Example
    ///
    /// Initial conditions in Mixed coordinates.
    /// ```
    /// # use dexter_simulate::*;
    /// let psi0 = InitialFlux::Toroidal(0.02);
    /// let initial1 = InitialConditions::mixed(0.0, psi0, 1.0, 0.0, 0.0, 0.0);
    /// let initial2 = InitialConditions::mixed(0.0, psi0, 2.0, 0.0, 0.1, 0.0);
    /// let particle1 = Particle::new(&initial1);
    /// let particle2 = Particle::new(&initial2);
    ///
    /// let queue = Queue::from_particles(&[particle1, particle2]);
    ///
    /// assert_eq!(initial1.theta0(), queue[0].initial_conditions().theta0());
    /// assert_eq!(initial2.theta0(), queue[1].initial_conditions().theta0());
    /// # Ok::<_, SimulationError>(())
    /// ```
    #[must_use]
    pub fn from_particles(particles: &[Particle]) -> Self {
        assert!(!particles.is_empty(), "Particle slice cannot be empty");
        Self {
            initial_conditions: QueueInitialConditions::from_particles(particles),
            particles: particles.to_vec(),
            stats: QueueStats::default(),
            routine: Routine::None,
        }
    }
}

impl Queue {
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
    /// let lcfs = LastClosedFluxSurface::Toroidal(qfactor.psi_last());
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

    /// Integrates all the contained particles, calculating their intersections with a constant θ or ζ surface.
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

    /// Integrates all the contained particles for a specific number of periods and calculates
    /// their frequencies.
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
    /// let perturbation = Perturbation::zero();
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
    /// queue.close(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &perturbation,
    ///     1,
    ///     &SolverParams::default(),
    /// );
    /// # Ok::<_, SimulationError>(())
    /// ```
    pub fn close<Q, C, B, H>(
        &mut self,
        qfactor: &Q,
        current: &C,
        bfield: &B,
        perturbation: &Perturbation<H>,
        periods: usize,
        solver_params: &SolverParams,
    ) where
        Q: Qfactor + FluxCommute + Send + Sync,
        C: Current + Send + Sync,
        B: Bfield + Send + Sync,
        H: Harmonic + Send + Sync,
    {
        let pbar = ClosePbar::new(self);
        pbar.print_prelude();

        self.particles.par_iter_mut().for_each(|particle| {
            particle.close(
                qfactor,
                current,
                bfield,
                perturbation,
                periods,
                solver_params,
            );
            pbar.inc(&particle.integration_status());
            pbar.print_stats();
            particle.discard_arrays();
        });
        pbar.finish();

        self.routine = Routine::Close;
        self.stats = QueueStats::from_completed_queue(self);
    }

    /// Classifies all the contained particles' orbits. using their position on the (E, Pζ, μ=const)
    /// plane without integrating.
    ///
    /// # Note
    ///
    /// This method is experimental. It is exact for LAR equilibria and approximately correct for
    /// tokamak equilibria, depending on how shaped they are.
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
    /// queue.classify(&qfactor, &current, &bfield);
    /// # Ok::<_, SimulationError>(())
    /// ```
    pub fn classify<Q, C, B>(&mut self, qfactor: &Q, current: &C, bfield: &B)
    where
        Q: Qfactor + FluxCommute + Send + Sync,
        C: Current + Send + Sync,
        B: Bfield + Send + Sync,
    {
        let pbar = ClassifyPbar::new(self);
        pbar.print_prelude();

        self.particles.par_iter_mut().for_each(|particle| {
            particle.classify(qfactor, current, bfield);
            pbar.inc(&particle.orbit_type());
            pbar.print_stats();
        });
        pbar.finish();

        self.routine = Routine::Classify;
        self.stats = QueueStats::from_completed_queue(self);
    }

    /// Classifies all the contained particles' orbits. using their position on the (E, Pζ, μ=const)
    /// plane without integrating.
    ///
    /// This method is an optimization to [`Self::classify`]. Since the most common scenario is to
    /// classify particles with the same `μ`, we can generate the [`EnergyPzetaPlane`] only once
    /// and use it for all particles. This improves performance by 5-8 times.
    ///
    /// # Note
    ///
    /// This method is experimental. It is exact for LAR equilibria and approximately correct for
    /// tokamak equilibria, depending on how shaped they are.
    ///
    /// # Panics
    ///
    /// Panics if not every initial `μ` is equal with the others.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// # use std::path::PathBuf;
    /// #
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.6);
    /// let qfactor = ParabolicQfactor::new(1.1, 4.2, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    ///
    /// use InitialFlux::*;
    /// let initial_conditions = QueueInitialConditions::boozer(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[1e-4, 2e-4],
    ///     &[7e-6, 7e-6], // must all be equal
    /// )?;
    /// let mut queue = Queue::new(&initial_conditions);
    /// queue.classify_common_mu(&qfactor, &current, &bfield);
    /// # Ok::<_, SimulationError>(())
    /// ```
    #[expect(clippy::float_cmp, reason = "we need bit-to-bit equivalence")]
    pub fn classify_common_mu<Q, C, B>(&mut self, qfactor: &Q, current: &C, bfield: &B)
    where
        Q: Qfactor + FluxCommute + Send + Sync,
        C: Current + Send + Sync,
        B: Bfield + Send + Sync,
    {
        let mus = self.initial_conditions.mu_array_view();
        let mu = mus
            .first()
            .expect("Queue cannot be instantiated with empty arrays");
        assert!(
            mus.iter().all(|_mu| _mu == mu),
            "All initial `mu0` must be equal"
        );

        let plane = EnergyPzetaPlane::from_mu(qfactor, current, bfield, *mu);

        let pbar = ClassifyPbar::new(self);
        pbar.print_prelude();

        self.particles.par_iter_mut().for_each(|particle| {
            particle._classify(qfactor, current, bfield, Some(&plane));
            pbar.inc(&particle.orbit_type());
            pbar.print_stats();
        });
        pbar.finish();

        self.routine = Routine::Classify;
        self.stats = QueueStats::from_completed_queue(self);
    }
}

/// Slicing.
impl Queue {
    /// Iterates through `self`'s particles, keeping only the ones with an initial `Pζ` within the
    /// given `span`, discarding the rest.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// # use std::path::PathBuf;
    /// # use ndarray::Array1;
    /// #
    /// let num = 11;
    /// let psi0s = toroidal_fluxes(&Array1::linspace(0.0, 0.05, num).to_vec());
    /// let initial_conditions = QueueInitialConditions::mixed(
    ///     &Array1::zeros(num).to_vec(),
    ///     &psi0s.to_vec(),
    ///     &Array1::zeros(num).to_vec(),
    ///     &Array1::zeros(num).to_vec(),
    ///     &Array1::linspace(0.0, 1.0, num).to_vec(), // Pζ = 0.0, 0.1, ..., 0.9, 1.0
    ///     &Array1::zeros(num).to_vec(),
    /// )?;
    /// let mut queue = Queue::new(&initial_conditions);
    ///
    /// queue.retain_pzeta(0.45..0.64); // Only Pζ = 0.5 and Pζ = 0.6 remain
    /// assert_eq!(queue.particle_count(), 2);
    ///
    /// # Ok::<_, SimulationError>(())
    /// ```
    pub fn retain_pzeta(&mut self, span: Range<f64>) {
        self.particles.retain(|particle| {
            particle
                .initial_conditions()
                .pzeta0
                .is_some_and(|pzeta| span.contains(&pzeta))
        });
    }

    /// Iterates through `self`'s particles, keeping only the ones with an initial energy within the
    /// given `span`, discarding the rest.
    pub fn retain_energy(&mut self, span: Range<f64>) {
        self.particles.retain(|particle| {
            particle
                .initial_energy()
                .is_some_and(|energy| span.contains(&energy))
        });
    }

    /// Iterates through `self`'s particles, keeping only **one** particle in each `Pζ` *bin*.
    ///
    /// *Bins* are defined as the intervals:\
    /// `Pζmin <= Pζmin + ΔPζ <= Pζmin + 2ΔPζ <= ... <= Pζmin + (n-1)ΔPζ <= Pζmax`,\
    /// where `n=1..num_bins` and `ΔPζ = (Pζmax - Pζmin)/n`.
    ///
    /// When called after [`Queue::retain_energy`] with a small energy span, we can obtain a set of
    /// particles of approximately same energy and approximately equispaced `Pζ`s, depending on
    /// the `span` and `num_bins` parameters. This is useful since we essentially have no control
    /// over the particles' energies.
    ///
    /// # Note
    ///
    /// Particles are visited in order of instantiation, and the first particle that falls in an
    /// unfilled bin will be selected.
    pub fn bin_pzeta(&mut self, num_bins: usize) {
        // SAFETY: In the case that all Pζs are NaNs, which can happen with partly initialized
        // particles, all particles are simply discarded.
        let pzetas = Array1::from_iter(
            self.particles
                .iter()
                .map(|particle| particle.initial_conditions().pzeta0.unwrap_or(f64::NAN)),
        );
        let pzetamin = *pzetas.min_skipnan();
        let pzetamax = *pzetas.max_skipnan();

        let bins = Array1::linspace(pzetamin, pzetamax, num_bins + 1);
        let mut filled: Array1<bool> = Array1::from_elem(bins.len(), false);

        self.particles.retain(|particle| {
            let mut found = false;
            let pzeta = particle.initial_conditions().pzeta0.unwrap_or(f64::NAN);
            for n in 0..num_bins {
                if !filled[n] && (bins[n] <= pzeta) && (pzeta <= bins[n + 1]) {
                    filled[n] = true;
                    found = true;
                    break;
                }
            }
            found
        });
    }

    /// Iterates through `self`'s particles, keeping only the ones with an [`EnergyPzetaPosition`]
    /// that matches an element of `positions`.
    pub fn retain_energy_pzeta_positions(&mut self, positions: &[EnergyPzetaPosition]) {
        self.particles.retain(|particle| {
            let mut keep = false;
            for pos in positions {
                if particle.energy_pzeta_position() == *pos {
                    keep = true;
                    break;
                }
            }
            keep
        });
    }

    /// Iterates through `self`'s particles, keeping only the ones with an [`OrbitType`]
    /// that matches an element of `orbit_types`.
    pub fn retain_orbit_types(&mut self, orbit_types: &[OrbitType]) {
        self.particles.retain(|particle| {
            let mut keep = false;
            for typ in orbit_types {
                if particle.orbit_type() == *typ {
                    keep = true;
                    break;
                }
            }
            keep
        });
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

    /// Returns an [`Array1`] with all the particles' calculated **initial** energies.
    ///
    /// Particles are visited in order of instantiation.
    ///
    /// For particles whose energies have not been calculated, the corresponding element is
    /// replaced with NaN.
    #[must_use]
    pub fn energy_array(&self) -> Array1<f64> {
        Array1::from_iter(
            self.particles
                .iter()
                .map(|particle| particle.initial_energy().unwrap_or(f64::NAN)),
        )
    }

    /// Returns an [`Array1`] with the number of steps each particle has taken.
    ///
    /// Particles are visited in order of instantiation.
    #[must_use]
    pub fn steps_taken_array(&self) -> Array1<usize> {
        Array1::from_iter(self.particles.iter().map(Particle::steps_taken))
    }

    /// Returns an [`Array1`] with the number of steps each particle has stored.
    ///
    /// Particles are visited in order of instantiation.
    #[must_use]
    pub fn steps_stored_array(&self) -> Array1<usize> {
        Array1::from_iter(self.particles.iter().map(Particle::steps_stored))
    }

    /// Returns an [`Array1`] with all the particles' calculated `ωθ`.
    ///
    /// Particles are visited in order of instantiation.
    ///
    /// For particles where the [`Particle::close()`] failed, the corresponding element is
    /// replaced with `NaN`.
    #[must_use]
    pub fn omega_theta_array(&self) -> Array1<f64> {
        Array1::from_iter(
            self.particles
                .iter()
                .map(|particle| particle.omega_theta().unwrap_or(f64::NAN)),
        )
    }

    /// Returns an [`Array1`] with all the particles' calculated `ωζ`.
    ///
    /// Particles are visited in order of instantiation.
    ///
    /// For particles where the [`Particle::close()`] failed, the corresponding element is
    /// replaced with `NaN`.
    #[must_use]
    pub fn omega_zeta_array(&self) -> Array1<f64> {
        Array1::from_iter(
            self.particles
                .iter()
                .map(|particle| particle.omega_zeta().unwrap_or(f64::NAN)),
        )
    }

    /// Returns an [`Array1`] with all the particles' calculated `qkinetic`.
    ///
    /// Particles are visited in order of instantiation.
    ///
    /// For particles where the [`Particle::close()`] failed, the corresponding element is
    /// replaced with `NaN`.
    #[must_use]
    pub fn qkinetic_array(&self) -> Array1<f64> {
        Array1::from_iter(
            self.particles
                .iter()
                .map(|particle| particle.qkinetic().unwrap_or(f64::NAN)),
        )
    }

    /// Returns an [`Array1`] with all the particles' routine durations.
    ///
    /// Particles are visited in order of instantiation.
    #[must_use]
    pub fn durations(&self) -> Array1<Duration> {
        Array1::from_iter(self.particles.iter().map(Particle::duration))
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
