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
            pub fn $eval_method(&self, flux: f64) -> Result<f64, PyEvalError> {
                let mut acc = Accelerator::new();
                Ok(self.0.$eval_method(flux, &mut acc)?)
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
            pub fn $eval_method(&self, flux: f64, theta: f64) -> Result<f64, PyEvalError> {
                let mut xacc = Accelerator::new();
                let mut yacc = Accelerator::new();
                let mut cache = Cache::new();
                Ok(self
                    .0
                    .$eval_method(flux, theta, &mut xacc, &mut yacc, &mut cache)?)
            }
        }
    };
}

/// Generates an eval method from the wrapped Harmonic object.
#[macro_export]
macro_rules! py_eval_harmonic {
    ($py_object:ident, $eval_method:ident) => {
        #[pymethods]
        impl $py_object {
            #[rustfmt::skip]
            pub fn $eval_method(&self, flux: f64, theta: f64, zeta: f64, t: f64) -> Result<f64, PyEvalError> {
                let mut cache = self.0.generate_cache();
                Ok(self
                    .0
                    .$eval_method(flux, theta, zeta, t, &mut cache)?)
            }
        }
    };
}

/// Generates an eval method from the wrapped Perturbation object.
#[macro_export]
macro_rules! py_eval_perturbation {
    ($py_object:ident, $eval_method:ident) => {
        #[pymethods]
        impl $py_object {
            #[rustfmt::skip]
            pub fn $eval_method(&self, flux: f64, theta: f64, zeta: f64, t: f64) -> Result<f64, PyEvalError> {
                // Unfortuneately we have to generate all the caches at each call
                let mut caches = self.0.generate_caches();
                Ok(self
                    .0
                    .$eval_method(flux, theta, zeta, t, &mut caches)?)
            }
        }
    };
}
