use std::time::Duration;

use particle::Particle;
use safe_unwrap::safe_unwrap;

use crate::{Heap, HeapInitialConditions, heap::Routine};

#[non_exhaustive]
#[derive(Default)]
pub struct HeapStats {
    routine: Routine,
    total_particles: usize,
    mapped: usize,
    escaped: usize,
    timedout: usize,
    failed: usize,
    /// Duration of the slowest particle.
    slowest: MapDuration,
    /// Duration of the fastest particle.
    fastest: MapDuration,
}
impl HeapStats {
    /// Creates a new [`HeapResults`], only containing information about the initial conditions.
    pub fn new(initial: &HeapInitialConditions) -> Self {
        Self {
            total_particles: initial.len(),
            ..Default::default()
        }
    }

    /// Creates a new Self from a Heap. This is needed since borrow checker exists
    pub fn from_heap(heap: &Heap) -> Self {
        let mut stat = Self::new(&heap.initials);
        stat.update(heap);
        stat
    }

    fn update(&mut self, heap: &Heap) {
        self.update_flags(heap);
        self.calculate_particle_nums(heap);
        self.calculate_durations(heap);
    }

    /// Updates various future flags.
    fn update_flags(&mut self, heap: &Heap) {
        self.routine = heap.routine.clone()
    }

    /// Counts the occurences of each [`IntegrationStatus`]'s variants.
    fn calculate_particle_nums(&mut self, heap: &Heap) {
        macro_rules! count_variants {
            ($is_enum:ident) => {
                heap.particles
                    .iter()
                    .filter(|p| p.status.$is_enum())
                    .count()
            };
        }
        self.mapped = count_variants!(is_mapped);
        self.escaped = count_variants!(is_escaped);
        self.timedout = count_variants!(is_timed_out);
        self.escaped = count_variants!(is_escaped);
        self.failed = count_variants!(is_failed);
        self.total_particles = heap.particles.len(); // Update just in case
    }

    /// Calculates and stores the fastest and slowest integrations.
    fn calculate_durations(&mut self, heap: &Heap) {
        let slowest = safe_unwrap!(
            "poincare.particles is non-empty",
            heap.particles.iter().max_by_key(|p| p.evolution.duration)
        );
        let fastest = safe_unwrap!(
            "poincare.particles is non-empty",
            heap.particles
                .iter()
                .filter(|p| p.evolution.steps_stored() > 0) // Drop invalid
                .min_by_key(|p| p.evolution.duration)
        );
        self.slowest = MapDuration::from(slowest);
        self.fastest = MapDuration::from(fastest);
    }
}

impl std::fmt::Debug for HeapStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HeapStats")
            .field("routine", &self.routine)
            .field("total_particles", &self.total_particles)
            .field("mapped", &self.mapped)
            .field("escaped", &self.escaped)
            .field("timedout", &self.timedout)
            .field("failed", &self.failed)
            .field("slowest", &self.slowest)
            .field("fastest", &self.fastest)
            .finish()
    }
}

// ===============================================================================================

/// Helper struct to display fastest and slowest particles
#[derive(Default)]
struct MapDuration {
    steps: usize,
    duration: Duration,
}

impl From<&Particle> for MapDuration {
    fn from(p: &Particle) -> Self {
        Self {
            steps: p.evolution.steps_taken(),
            duration: p.evolution.duration,
        }
    }
}

impl std::fmt::Debug for MapDuration {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "duration: {:?} ({} steps)", self.duration, self.steps)
    }
}
