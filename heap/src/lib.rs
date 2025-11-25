mod error;
mod initials;

pub use error::HeapError;
pub use initials::HeapInitialConditions;

pub type Result<T> = std::result::Result<T, HeapError>;
