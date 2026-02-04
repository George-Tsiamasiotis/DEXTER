//! Wrappers that expose the inner types' evaluation methods.
//!
//! We must create a new [`Accelerator`]/[`Cache`] in every call, but this is preferable than
//! messing up the newtype pattern by storing Accelerators inside the pywrappers, especially since
//! all file [`Qfactor`], [`Currents`], [`Bfield`], [`Harmonic`] and [`Perturbation`] require
//! a different number of Accelerators and [`Caches`].

/// Generates a 1D eval method from the wrapped Rust object.
#[macro_export]
macro_rules! py_eval1D {
    ($py_object:ident, $eval_method:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn $eval_method(&self, psip: f64) -> Result<f64, PyEqError> {
                let mut acc = Accelerator::new();
                Ok(self.0.$eval_method(psip, &mut acc)?)
            }
        }
    };
}

/// Generates a 2D eval method from the wrapped Rust object.
#[macro_export]
macro_rules! py_eval2D {
    ($py_object:ident, $eval_method:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn $eval_method(&self, psip: f64, theta: f64) -> Result<f64, PyEqError> {
                let mut xacc = Accelerator::new();
                let mut yacc = Accelerator::new();
                let mut cache = Cache::new();
                Ok(self
                    .0
                    .$eval_method(psip, theta, &mut xacc, &mut yacc, &mut cache)?)
            }
        }
    };
}
