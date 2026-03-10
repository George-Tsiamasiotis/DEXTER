//! Python wrappers for rust error types.

use dexter::dexter_equilibrium::{EqError, EvalError};
use dexter::dexter_simulate::SimulationError;
use pyo3::PyErr;
use pyo3::exceptions::PyException;

/// Creates newtype wrappers around the foreign error types, to allow conversion to [`PyErr`] and
/// use of the `?` operator.
///
/// Source:
/// [`https://pyo3.rs/main/function/error-handling#foreign-rust-error-types`]
macro_rules! to_pyerr_impl {
    ($error_type:ident, $py_error_type: ident) => {
        #[derive(Debug)]
        pub struct $py_error_type(String);

        impl From<$py_error_type> for PyErr {
            fn from(error: $py_error_type) -> Self {
                PyException::new_err(error.0.to_string())
            }
        }

        impl From<$error_type> for $py_error_type {
            fn from(error: $error_type) -> Self {
                Self(error.to_string())
            }
        }

        impl From<PyErr> for $py_error_type {
            fn from(error: PyErr) -> Self {
                Self(error.to_string())
            }
        }
    };
}

to_pyerr_impl!(EqError, PyEqError);
to_pyerr_impl!(EvalError, PyEvalError);
to_pyerr_impl!(SimulationError, PySimulationError);
