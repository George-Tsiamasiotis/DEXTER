use std::time::Duration;

use particle::Particle;

use crate::{Heap, HeapInitialConditions, heap::Routine};

/// Keeps track of the number of Particles per `IntegrationStatus`.
#[derive(Default, Debug)]
struct IntegrationStatusNums {
    initialized: usize,
    integrated: usize,
    single_period: usize,
    mapped: usize,
    escaped: usize,
    evaluation_nan: usize,
    timedout: usize,
    invalid: usize,
    failed: usize,
}

/// Keeps track of the number of Particles per `OrbitType`.
#[derive(Default, Debug)]
struct OrbitTypeNums {
    undefined: usize,
    trapped: usize,
    passing: usize,
}

#[non_exhaustive]
#[derive(Default)]
pub struct HeapStats {
    routine: Routine,
    status_nums: IntegrationStatusNums,
    orbit_type_nums: OrbitTypeNums,
    total_particles: usize,
    /// Duration of the slowest particle.
    slowest: ParticleDuration,
    /// Duration of the fastest particle.
    fastest: ParticleDuration,
}

macro_rules! count_variants {
    ($heap:ident, $which:ident, $is_enum:ident) => {
        $heap
            .particles
            .iter()
            .filter(|p| p.$which.$is_enum())
            .count()
    };
}

impl HeapStats {
    /// Creates a new [`HeapStats`], only containing information about the initial conditions.
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
        self.calculate_orbit_types_nums(heap);
    }

    /// Updates various future flags.
    fn update_flags(&mut self, heap: &Heap) {
        self.routine = heap.routine.clone()
    }

    /// Counts the occurences of each [`IntegrationStatus`]'s variants.
    fn calculate_particle_nums(&mut self, heap: &Heap) {
        self.status_nums.initialized = count_variants!(heap, status, is_initialized);
        self.status_nums.integrated = count_variants!(heap, status, is_integrated);
        self.status_nums.mapped = count_variants!(heap, status, is_mapped);
        self.status_nums.escaped = count_variants!(heap, status, is_escaped);
        self.status_nums.evaluation_nan = count_variants!(heap, status, is_evaluation_nan);
        self.status_nums.timedout = count_variants!(heap, status, is_timed_out);
        self.status_nums.failed = count_variants!(heap, status, is_failed);
        self.status_nums.invalid = count_variants!(heap, status, is_invalid_intersections);
        self.status_nums.single_period = count_variants!(heap, status, is_single_period_integrated);
        self.total_particles = heap.particles.len(); // Update just in case
    }

    fn calculate_orbit_types_nums(&mut self, heap: &Heap) {
        self.orbit_type_nums.undefined = count_variants!(heap, orbit_type, is_undefined);
        self.orbit_type_nums.trapped = count_variants!(heap, orbit_type, is_trapped);
        self.orbit_type_nums.passing = count_variants!(heap, orbit_type, is_passing);
    }

    /// Calculates and stores the fastest and slowest integrations.
    fn calculate_durations(&mut self, heap: &Heap) {
        let slowest = heap.particles.iter().max_by_key(|p| p.evolution.duration);
        let fastest = heap
            .particles
            .iter()
            .filter(|p| p.evolution.steps_stored() > 0) // Drop invalid
            .min_by_key(|p| p.evolution.duration);
        self.slowest = ParticleDuration::from(slowest);
        self.fastest = ParticleDuration::from(fastest);
    }
}

impl std::fmt::Debug for HeapStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Heap") // Propagated to Heap
            .field("routine", &self.routine)
            .field("total_particles", &self.total_particles)
            .field("IntegrationStatusNums", &self.status_nums)
            .field("OrbitTypeNums", &self.orbit_type_nums)
            .field("slowest", &self.slowest)
            .field("fastest", &self.fastest)
            .finish()
    }
}

// ===============================================================================================

/// Helper struct to display fastest and slowest particles
#[derive(Default)]
struct ParticleDuration {
    steps: usize,
    duration: Duration,
}

impl From<Option<&Particle>> for ParticleDuration {
    fn from(p: Option<&Particle>) -> Self {
        let (duration, steps) = if let Some(p) = p {
            (p.evolution.duration, p.evolution.steps_taken())
        } else {
            (Duration::default(), 0)
        };
        Self { steps, duration }
    }
}

impl std::fmt::Debug for ParticleDuration {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "duration: {:?} ({} steps taken)",
            self.duration, self.steps
        )
    }
}
