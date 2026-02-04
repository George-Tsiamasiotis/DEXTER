/// Generates getters that return `[T]` fields to `Array1<T>`.
#[doc(hidden)]
#[macro_export]
macro_rules! vec_to_array1D_getter_impl {
    ($fun_name:ident, $($field:ident).+, $var_name:ident) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values as a 1D array." ]
        pub fn $fun_name(&self) -> Array1<f64> {
            Array1::from_vec(self.$($field).+.clone())
        }
    }
}

/// Creates a getter method for extracting the flat Vec data as an Array2.
/// The Vec is assumed to be in Fortran order, since it is intended for use by the splines.
#[doc(hidden)]
#[macro_export]
macro_rules! fortran_vec_to_carray2d_impl {
    ($meth_name:ident, $($field:ident).+, $var_name:ident) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values as a 2D array." ]
        pub fn $meth_name(&self) -> Array2<f64> {
            // Array is in Fortran order., so we must reverse the shape
            let actual_shape = self.shape();
            let shape = (actual_shape.1, actual_shape.0);
            Array2::from_shape_vec(shape, self.$($field).+.clone())
                .expect("Shape is correct by definition")
                .reversed_axes()
        }
    }
}

/// Generates getters for the fluxes' values at the wall
#[doc(hidden)]
#[macro_export]
macro_rules! fluxes_wall_value_getter_impl {
    () => {
        /// Returns the toroidal flux's value at the wall `ψ_wall`.
        pub fn psi_wall(&self) -> Option<f64> {
            self.psi.wall_value()
        }

        /// Returns the poloidal flux's value at the wall `ψp_wall`.
        pub fn psip_wall(&self) -> Option<f64> {
            self.psip.wall_value()
        }
    };
}

/// Generates getters for the fluxes' states.
#[doc(hidden)]
#[macro_export]
macro_rules! fluxes_state_getter_impl {
    () => {
        /// Returns the toroidal flux's state.
        pub fn psi_state(&self) -> NcFluxState {
            self.psi.state.clone()
        }

        /// Returns the poloidal flux's state.
        pub fn psip_state(&self) -> NcFluxState {
            self.psip.state.clone()
        }
    };
}
