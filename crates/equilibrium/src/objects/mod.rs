pub mod bfields;
pub mod currents;
pub mod geometries;
pub mod harmonics;
pub mod perturbations;
pub mod qfactors;

/// Creates a getter method for extracting the flat Vec data as an Array2.
#[doc(hidden)]
#[macro_export]
macro_rules! fortran_vec_to_carray2d_impl {
    ($meth_name:ident, $field:ident, $var_name:ident) => {
        /// Returns the `R(ψp, θ)` data as a 2D array.
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "`(ψp, θ) data as an [`Array2<f64>`]"]
        pub fn $meth_name(&self) -> Array2<f64> {
            // Array is in Fortran order.
            let shape = (self.theta_data.len(), self.psip_data.len());
            Array2::from_shape_vec(shape, self.$field.to_vec())
                .expect("shape is correct by definition")
                .reversed_axes()
        }
    };
}
