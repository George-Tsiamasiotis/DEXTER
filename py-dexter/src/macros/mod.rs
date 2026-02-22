//! Macros for exposing Rust methods to the Python API.
//!
//! It is necessary to create a new `#[pymethods]` impl block every time, since `#[pymethods]` does not
//! allow macros to be used inside it.

mod eval;
mod generics;
mod numpy;
mod pymethods;
