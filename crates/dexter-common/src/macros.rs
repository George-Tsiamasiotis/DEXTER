//! Utility macros.

/// Generates a getter method for an `Array1<T>` field.
#[macro_export]
macro_rules! f64_slice_getter {
    ($fun_name:ident, $($field:ident).+, $var_name:ident) => {
        #[doc = "Returns a reference to the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values." ]
        pub fn $fun_name(&self) -> &[f64] {
            &self.$($field).+
        }
    };
    ($fun_name:ident, $($field:ident).+(), $var_name:ident) => {
        #[doc = "Returns a reference to the `"]
        #[doc = stringify!($var_name)]
        #[doc = "` values." ]
        pub fn $fun_name(&self) -> &[f64] {
            self.$($field).+()
        }
    }
}
