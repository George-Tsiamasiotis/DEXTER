//! Multithreading utility functions.

#![expect(clippy::unwrap_used, reason = "available_parallelism should not panic")]
#![expect(clippy::missing_panics_doc, reason = "should be fatal")]

use std::thread::available_parallelism;

/// Returns the device's number of availiale threads.
///
/// The number is obtained by calling [`std::thread::available_parallelism`].
#[must_use]
pub fn get_max_threads() -> usize {
    available_parallelism().unwrap().into()
}

/// Sets the global number of threads.
///
/// If `num` is either `0` or greater than the device's number of threads, it defaults to the
/// device's number of threads.
///
/// # Safety
///
/// The setting is passed to `rayon` as an environment variable through [`std::env::set_var`],
/// which is inherently unsafe. To avoid UB, this function should be called *once* at the *start of
/// a process*.
#[expect(unsafe_code, reason = "must set environment variable for rayon")]
pub fn set_num_threads(num: usize) {
    let max_threads: usize = available_parallelism().unwrap().into();

    let num_threads = if (1..=max_threads).contains(&num) {
        num
    } else {
        max_threads
    };

    // SAFETY: Is only called at the start of a process.
    unsafe {
        std::env::set_var("RAYON_NUM_THREADS", num_threads.to_string());
    }
}
