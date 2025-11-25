mod error;
mod heap;
mod initials;
mod progress_bars;

pub use error::HeapError;
pub use heap::Heap;
pub use initials::HeapInitialConditions;

pub type Result<T> = std::result::Result<T, HeapError>;
