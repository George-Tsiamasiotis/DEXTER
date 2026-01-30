//! Wrappers that expose the inner types' [`ArrayD`] methods.

/// Generates a getter method that returns a 1D numpy array.
#[macro_export]
macro_rules! py_get_numpy1D {
    ($py_object:ident, $getter:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
                self.0.$getter().into_pyarray(py)
            }
        }
    };
}
