/// Generates getters that return a primitive field.
#[macro_export]
macro_rules! primitive_getter_impl {
    ($field:ident, $typ:ident, $doc:literal) => {
        #[doc = $doc]
        pub fn $field(&self) -> $typ {
            self.$field.clone()
        }
    };
}

/// Generates getters that return an Option<`primitive`>.
#[macro_export]
macro_rules! fallible_primitive_getter_impl {
    ($field:ident, $typ:ident, $doc:literal) => {
        #[doc = $doc]
        pub fn $field(&self) -> Option<$typ> {
            self.$field.clone()
        }
    };
}

/// Generates getters that return `[T]` fields to `Array1<T>`.
#[macro_export]
macro_rules! array1D_getter_impl {
    ($fun_name:ident, $($field:ident).+) => {
        #[doc = "Returns the `"]
        #[doc = stringify!($($field).+)]
        #[doc = "` array." ]
        pub fn $fun_name(&self) -> Array1<f64> {
            Array1::from(self.$($field).+.clone())
        }
    }
}
