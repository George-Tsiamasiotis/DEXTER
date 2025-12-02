//! A collection of Particles.

use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use ndarray::{Array1, Array2, Axis};
use particle::{Flux, Radians};
use particle::{MappingParameters, Particle};
use utils::array2D_getter_impl;

use crate::progress_bars::{FrequenciesPbar, PoincarePbar};
use crate::{HeapInitialConditions, HeapStats, Result};

/// Describes the Routine by which the Heap's particle's where integrated.
#[non_exhaustive]
#[derive(Default, Clone, Debug)]
pub enum Routine {
    #[default]
    None,
    Integration,
    Poincare(MappingParameters),
    SinglePeriod,
}

/// A collections of mutliple [`Particle`]s, constructed from [`HeapInitialConditions`].
#[derive(Default)]
pub struct Heap {
    /// Initial conditions arrays.
    pub initials: HeapInitialConditions,
    /// Tracked [`Particle`]s.
    pub particles: Vec<Particle>,
    /// Describes the Routine by which the Heap's particle's where integrated.
    pub routine: Routine,
    /// Calculation results
    pub stats: HeapStats,
    /// The calculated θ angles, ignoring `PoincareSection`.
    pub thetas: Array2<Radians>,
    /// The calculated ζ angles, ignoring `PoincareSection`.
    pub zetas: Array2<Radians>,
    /// The calculated ψp flux, ignoring `PoincareSection`.
    pub psips: Array2<Flux>,
    /// The calculated ψ flux, ignoring `PoincareSection`.
    pub psis: Array2<Flux>,
}

impl Heap {
    /// Creates a [`Heap`], initializing a particle for each set of Initial Conditions.
    pub fn new(initials: &HeapInitialConditions) -> Self {
        let particles = initials.to_particles();
        Self {
            initials: initials.clone(),
            particles,
            stats: HeapStats::new(initials),
            ..Default::default()
        }
    }

    pub fn poincare(
        &mut self,
        qfactor: &Qfactor,
        currents: &Currents,
        bfield: &Bfield,
        perturbation: &Perturbation,
        params: &MappingParameters,
    ) -> Result<()> {
        let pbar = PoincarePbar::new(self, params);
        pbar.print_prelude();

        self.particles.par_iter_mut().try_for_each(|p| {
            let res = p.map(qfactor, bfield, currents, perturbation, params);
            pbar.inc(&p.status);
            pbar.print_stats();
            res
        })?;
        pbar.finish();

        self.routine = Routine::Poincare(*params);
        self.stats = HeapStats::from_heap(self);
        self.store_arrays(params)?;
        Ok(())
    }

    pub fn calculate_frequencies(
        &mut self,
        qfactor: &Qfactor,
        currents: &Currents,
        bfield: &Bfield,
        perturbation: &Perturbation,
    ) -> Result<()> {
        let pbar = FrequenciesPbar::new(self);
        pbar.print_prelude();

        self.particles.par_iter_mut().try_for_each(|p| {
            let res = p.calculate_frequencies(qfactor, bfield, currents, perturbation);
            pbar.inc(&p.status);
            pbar.print_stats();
            p.evolution.discard(); // We don't need the arrays anymore, avoid wasting memory
            res
        })?;
        pbar.finish();

        self.routine = Routine::SinglePeriod;
        self.stats = HeapStats::from_heap(self);
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.initials.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// Data extraction
impl Heap {
    array2D_getter_impl!(zetas, zetas, Radians);
    array2D_getter_impl!(psips, psips, Flux);
    array2D_getter_impl!(thetas, thetas, Radians);
    array2D_getter_impl!(psis, psis, Flux);

    /// Stores the Particle's time series and stacked into 2D arrays.
    fn store_arrays(&mut self, params: &MappingParameters) -> Result<()> {
        // We dont now how many particle's got completely integrated, so we push a new row for
        // every successful one.
        // We also include the initial point for now and drop it later, otherwise the code gets
        // ugly.
        let columns = params.intersections + 1;
        let shape = (0, columns);
        self.zetas = Array2::from_elem(shape, Radians::NAN);
        self.psips = Array2::from_elem(shape, Flux::NAN);
        self.thetas = Array2::from_elem(shape, Radians::NAN);
        self.psis = Array2::from_elem(shape, Flux::NAN);

        /// Copies the array of the calculated evolution `source` data into a new 1D array with length
        /// `columns` and pushes it to the 2D array `array`. If `len(source) < columns`, which will
        /// happen with escaped or timed out particles, the rest of the array is filled with NaN. This
        /// allows us to plot those particle's as well, while keeping all the data in the same 2D
        /// array.
        macro_rules! copy_and_fill_with_nan_and_push_row {
            ($particle:ident, $results_array:ident, $source:ident) => {
                assert!($particle.evolution.steps_stored() <= columns);
                self.$results_array.push_row(
                    Array1::from_shape_fn(columns, |i| {
                        $particle
                            .evolution
                            .$source()
                            .get(i)
                            .copied()
                            .unwrap_or(f64::NAN)
                    })
                    .view(),
                )?;
                // Still includes the initial point
                assert_eq!(self.$results_array.ncols(), params.intersections + 1);
            };
        }
        for p in self.particles.iter() {
            if should_be_plotted(p) {
                copy_and_fill_with_nan_and_push_row!(p, zetas, zeta);
                copy_and_fill_with_nan_and_push_row!(p, psips, psip);
                copy_and_fill_with_nan_and_push_row!(p, thetas, theta);
                copy_and_fill_with_nan_and_push_row!(p, psis, psi);
            }
        }

        // Remove intial points
        self.zetas.remove_index(Axis(1), 0);
        self.psips.remove_index(Axis(1), 0);
        self.thetas.remove_index(Axis(1), 0);
        self.psis.remove_index(Axis(1), 0);

        Ok(())
    }
}

/// Returns true if the particle should be plotted in the final Poincare map.
fn should_be_plotted(particle: &Particle) -> bool {
    let status_ok =
        particle.status.is_mapped() | particle.status.is_timed_out() | particle.status.is_escaped();
    // Drop particles that were initialized outside the wall
    let length_ok = particle.evolution.steps_stored() > 1;

    status_ok && length_ok
}

impl std::fmt::Debug for Heap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.stats.fmt(f)
    }
}
