/// Generates getters that return `[T]` fields to `Array1<T>`.
#[macro_export]
macro_rules! vec_to_array1D_getter_impl {
    ($fun_name:ident, $($field:ident).+) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($($field).+)]
        #[doc = "` values as a 1D array." ]
        pub fn $fun_name(&self) -> Array1<f64> {
            Array1::from_vec(self.$($field).+.clone())
        }
    }
}

/// Generates getters for the fluxes' values at the wall
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
