//! Numeric integration routines
//!
//! About Integration steps
//! -----------------------
//!
//! A. Setup
//!     1. Initialize return value `res`.
//!     2. Start `duration`
//!     3. Reset Evolution
//!     4. Evaluate `initial_state`
//!     5. Setup `state1` and `state2`
//!     6. Setup initial step
//! B. Main loop
//!     1. Time-out check (set status and break). Set `res` to the corresponding `ParticleError`
//!        variant.
//!     2. Routine success break. Set `res` to `Ok(())`.
//!     3. `state1` -> Step -> `state2`, steps_taken++
//!     4. Method specific checks
//!         - Push State
//!     5. `state1 = state2`
//! C. Finalization
//!     1. Routine specific final checks
//!     2. Evalulate `final_state`
//!     3. Finish() Evolution
//!     4. Calculate orbit type
//!     5. Save duration
//!     6. Return `res`

mod integrate;
mod map;
mod single_period_integrate;

mod henon;

pub(crate) use integrate::integrate;
pub(crate) use map::map_integrate;
pub(crate) use single_period_integrate::close_theta_period;

pub use map::{MappingParameters, PoincareSection};
pub use single_period_integrate::Frequencies;
