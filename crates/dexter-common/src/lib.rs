//! Common utilities across the workspace.

mod macros;
mod threads;

pub use threads::{get_max_threads, set_num_threads};
