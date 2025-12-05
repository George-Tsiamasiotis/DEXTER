//! Numeric integration routines
//!
//! About Integration steps
//! -----------------------
//!
//! A. Setup
//!     1. Start `duration`
//!     2. Reset Evolution
//!     3. Evaluate `initial_state`
//!     4. Setup `state1` and `state2`
//!     5. Setup initial step
//! B. Main loop
//!     1. Time-out check (set status and break)
//!     2. Routine success break
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

mod frequencies;
mod integrate;
mod mapping;

mod henon;

pub(crate) use frequencies::close_theta_period;
pub(crate) use integrate::integrate;
pub(crate) use mapping::map_integrate;

pub use frequencies::Frequencies;
pub use mapping::{MappingParameters, PoincareSection};
