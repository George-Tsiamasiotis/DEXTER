//! Helper struct to keep track of a Queue's run statistics.

#![expect(clippy::missing_docs_in_private_items, reason = "self-explanatory")]

use std::time::Duration;

use crate::{IntegrationStatus, Particle, Queue, QueueInitialConditions, queue::Routine};

/// Helper struct to keep track of [`IntegrationStatus`] counts in a [`Queue`].
#[derive(Default, Debug)]
struct IntegrationStatusCount {
    pub(crate) initialized: usize,
    pub(crate) out_of_bounds: usize,
    pub(crate) integrated: usize,
    pub(crate) escaped: usize,
    pub(crate) intersected: usize,
    pub(crate) intersected_timed_out: usize,
    pub(crate) invalid_intersections: usize,
    pub(crate) timed_out: usize,
    pub(crate) failed: usize,
}

/// Statistics container.
#[derive(Default)]
pub(crate) struct QueueStats {
    routine: Routine,
    integration_status_count: IntegrationStatusCount,
    particle_count: usize,
    slowest: ParticleDuration,
    fastest: ParticleDuration,
}

impl QueueStats {
    /// Instantiates a new [`QueueStats`]. Should be called at [`Queue`] initialization.
    pub(crate) fn from_initial_conditions(initial_conditions: &QueueInitialConditions) -> Self {
        Self {
            particle_count: initial_conditions.len(),
            ..Default::default()
        }
    }

    /// Fills the missing fields from a [`Queue`] that has completed its routine.
    ///
    /// Creating a new [`QueueStats`] is necessary since the previous one is already borrowed
    /// as mutable at the time that [`QueueStats`] must be updated.
    pub(crate) fn from_completed_queue(queue: &Queue) -> Self {
        Self::from_initial_conditions(&queue.initial_conditions()).update(queue)
    }

    /// Update self's fields.
    fn update(mut self, queue: &Queue) -> Self {
        self.update_flags(queue);
        self.count_integration_status(queue);
        self.calculate_durations(queue);
        self
    }

    /// Update any flag fields.
    pub(crate) fn update_flags(&mut self, queue: &Queue) {
        self.routine = queue.routine.clone()
    }

    /// Count how many times each [`IntegrationStatus`] variant appears in [`Queue`].
    pub(crate) fn count_integration_status(&mut self, queue: &Queue) {
        macro_rules! count_enum_variant {
            ($field:ident, $variant:ident) => {
                self.integration_status_count.$field = queue
                    .particles
                    .iter()
                    .filter(|particle| {
                        matches!(particle.integration_status(), IntegrationStatus::$variant)
                    })
                    .count();
            };
            ($field:ident, $variant:ident(..)) => {
                self.integration_status_count.$field = queue
                    .particles
                    .iter()
                    .filter(|particle| {
                        matches!(
                            particle.integration_status(),
                            IntegrationStatus::$variant(..)
                        )
                    })
                    .count();
            };
        }
        count_enum_variant!(initialized, Initialized);
        count_enum_variant!(out_of_bounds, OutOfBoundsInitialization);
        count_enum_variant!(integrated, Integrated);
        count_enum_variant!(escaped, Escaped);
        count_enum_variant!(intersected, Intersected);
        count_enum_variant!(intersected_timed_out, IntersectedTimedOut);
        count_enum_variant!(invalid_intersections, InvalidIntersections);
        count_enum_variant!(timed_out, TimedOut(..));
        count_enum_variant!(failed, Failed(..));
    }

    /// Calculates the fastest and slowest particles.
    pub(crate) fn calculate_durations(&mut self, queue: &Queue) {
        self.slowest = ParticleDuration::from(
            queue
                .particles
                .iter()
                .max_by_key(|particle| particle.duration()),
        );
        self.fastest = ParticleDuration::from(
            queue
                .particles
                .iter()
                .filter(|particle| particle.steps_stored() > 0)
                .min_by_key(|particle| particle.duration()),
        );
    }
}

impl std::fmt::Debug for QueueStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Queue") // Propagated to Queue's Debug
            .field("routine", &self.routine)
            .field("particle count", &self.particle_count)
            .field("integration status counts", &self.integration_status_count)
            .field("slowest", &self.slowest)
            .field("fastest", &self.fastest)
            .finish()
    }
}

// ===============================================================================================

/// Helper struct to calculate and store fastest and slowest durations.
#[derive(Default)]
struct ParticleDuration {
    /// The total number of steps taken.
    steps_taken: usize,
    /// The total duration.
    duration: Duration,
}

impl From<Option<&Particle>> for ParticleDuration {
    fn from(value: Option<&Particle>) -> Self {
        if let Some(particle) = value {
            Self {
                steps_taken: particle.steps_taken(),
                duration: particle.duration(),
            }
        } else {
            Self::default()
        }
    }
}

impl std::fmt::Debug for ParticleDuration {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "duration: {:?} ({} steps taken)",
            self.duration, self.steps_taken
        )
    }
}
