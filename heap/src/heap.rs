//! A collection of Particles.

use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use particle::{MappingParameters, Particle, PoincareSection};

use crate::progress_bars::PoincarePbar;
use crate::{HeapInitialConditions, HeapStats, Result};

/// Describes the Routine by which the Heap's particle's where integrated.
#[non_exhaustive]
#[derive(Default, Clone, Debug)]
pub enum Routine {
    #[default]
    None,
    Integration,
    Poincare(PoincareSection),
    SinglePeriod,
}

/// A collections of mutliple [`Particle`]s, constructed from [`HeapInitialConditions`].
pub struct Heap {
    /// Initial conditions arrays.
    pub initials: HeapInitialConditions,
    /// Tracked [`Particle`]s.
    pub particles: Vec<Particle>,
    /// Describes the Routine by which the Heap's particle's where integrated.
    pub routine: Routine,
    /// Calculation results
    pub stats: HeapStats,
}

impl Heap {
    /// Creates a [`Heap`], initializing a particle for each set of Initial Conditions.
    pub fn new(initials: &HeapInitialConditions) -> Self {
        let particles = initials.to_particles();
        Self {
            initials: initials.clone(),
            particles,
            routine: Routine::default(),
            stats: HeapStats::new(initials),
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
            p.map(qfactor, bfield, currents, perturbation, params)
                .inspect(|()| {
                    pbar.inc(&p.status);
                    pbar.print_stats();
                })
        })?;
        pbar.finish();

        self.routine = Routine::Poincare(params.section);
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

impl std::fmt::Debug for Heap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.stats.fmt(f)
    }
}
