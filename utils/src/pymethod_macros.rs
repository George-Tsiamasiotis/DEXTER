//! Miscellaneous getter methods
//!
//! It is necessary to create a new `#[pymethods]` impl block every time, since `#[pymethods]` does not
//! allow macros to be used inside it.

/// Generates a impl Debug block, using the inner type's Debug representation.
#[macro_export]
macro_rules! py_debug_impl {
    ($py_object:ident) => {
        impl std::fmt::Debug for $py_object {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                self.0.fmt(f)
            }
        }
    };
}

/// Generates a `__repr__` method, corresponding to the inner type's `Debug` representation.
#[macro_export]
macro_rules! py_repr_impl {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            pub fn __repr__(&self) -> String {
                format!("{:#?}", self.0)
            }
        }
    };
}

/// Generates a python getter method for a primitive type, which is expected to be a pub field
#[macro_export]
macro_rules! py_get_primitive_field {
    ($py_object:ident, $field:ident, $typ:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $field(&self) -> $typ {
                self.0.$field as $typ
            }
        }
    };
}

/// Exports a getter method on the exposed python wrapper
#[macro_export]
macro_rules! py_export_getter {
    ($py_object:ident, $getter:ident, $typ:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter(&self) -> $typ {
                self.0.$getter() as $typ
            }
        }
    };
}

/// Generates a getter method for the NetCDF's path field.
#[macro_export]
macro_rules! py_get_path {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn path(&self) -> String {
                String::from(safe_unwrap!(
                    "file already opened",
                    self.0.path.clone().into_os_string().into_string()
                ))
            }
        }
    };
}

/// Generates a getter method for the object's interpolation type.
#[macro_export]
macro_rules! py_get_typ {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn typ(&self) -> String {
                String::from(self.0.typ.clone())
            }
        }
    };
}
