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

/// Generates getters that export a vec->Array1D getter of an inner field
#[doc(hidden)]
#[macro_export]
macro_rules! export_array1D_getter_impl {
    ($fun_name:ident, $($field:ident).+, $var_name:ident) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values as a 1D array." ]
        pub fn $fun_name(&self) -> Array1<f64> {
            self.$($field).+.$fun_name()
        }
    }
}
