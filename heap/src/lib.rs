mod error;
mod heap;
mod initials;
mod progress_bars;
mod stats;

pub use error::HeapError;
pub use heap::Heap;
pub use initials::HeapInitialConditions;
pub use stats::HeapStats;

pub type Result<T> = std::result::Result<T, HeapError>;
