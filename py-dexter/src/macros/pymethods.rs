//! Miscellaneous getter methods

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

/// Generates a getter method for the NetCDF's path field.
#[macro_export]
macro_rules! py_get_path {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn path(&self) -> String {
                String::from(
                    self.0
                        .path()
                        .clone()
                        .into_os_string()
                        .into_string()
                        .unwrap(), // Safe: path exists since object exists.
                )
            }
        }
    };
}

/// Exports a getter method of a public field on the exposed python wrapper
#[macro_export]
macro_rules! py_export_pub_field {
    ($py_object:ident, $getter:ident, $T:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter(&self) -> $T {
                self.0.$getter as $T
            }
        }
    };
    // In python, this translates to `Optional[T]`, which is equivalent to
    // `None | T`
    ($py_object:ident, $getter:ident, Option<$T:ident>) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter(&self) -> Option<$T> {
                self.0.$getter as Option<$T>
            }
        }
    };
}

/// Exports a getter method on the exposed python wrapper
#[macro_export]
macro_rules! py_export_getter {
    ($py_object:ident, $getter:ident, $T:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter(&self) -> $T {
                self.0.$getter() as $T
            }
        }
    };
    // In python, this translates to `Optional[T]`, which is equivalent to
    // `None | T`
    ($py_object:ident, $getter:ident, Option<$T:ident>) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter(&self) -> Option<$T> {
                self.0.$getter() as Option<$T>
            }
        }
    };
}

/// Generates a python getter that returns the debug String ("{:?}") of an enum
#[macro_export]
macro_rules! py_get_enum_string {
    ($py_object:ident, $getter:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn $getter(&self) -> String {
                format!("{:?}", self.0.$getter())
            }
        }
    };
}

/// Generates a python getter that returns the netCDF convention version as a String
#[macro_export]
macro_rules! py_get_netcdf_version {
    ($py_object:ident) => {
        #[pymethods]
        impl $py_object {
            #[getter]
            pub fn get_netcdf_version(&self) -> String {
                self.0.netcdf_version().to_string()
            }
        }
    };
}
