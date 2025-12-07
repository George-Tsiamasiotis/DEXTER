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
