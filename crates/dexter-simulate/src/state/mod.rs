//! Definition of the state of a dynamical system.
//!
//! A `state` is a specific point on the current configuration space containing all the necessary
//! information for integrating the system, i.e dynamical coordinates, field values, time
//! derivatives etc.

mod energy;
mod guiding_center;

pub use energy::COMs;
pub(crate) use guiding_center::GCState;
