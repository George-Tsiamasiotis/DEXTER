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

/// Generates getters for the fluxes' values.
#[doc(hidden)]
#[macro_export]
macro_rules! fluxes_values_array_getter_impl {
    () => {
        /// Returns the toroidal flux's values as a 1D array, if they exist.
        pub fn psi_array(&self) -> Option<Array1<f64>> {
            self.psi
                .values()
                .map(|values| Array1::from(Vec::from(values)))
        }

        /// Returns the poloidal flux's values as a 1D array, if they exist.
        pub fn psip_array(&self) -> Option<Array1<f64>> {
            self.psip
                .values()
                .map(|values| Array1::from(Vec::from(values)))
        }
    };
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

/// Generates a getter for the object's equilibrium type.
#[doc(hidden)]
#[macro_export]
macro_rules! equilibrium_type_getter_impl {
    () => {
        /// Returns the object's [`EquilibriumType`].
        pub fn equilibrium_type(&self) -> EquilibriumType {
            self.equilibrium_type.clone()
        }
    };
}

/// Generates a getter for the object's path to netCDF file.
#[doc(hidden)]
#[macro_export]
macro_rules! netcdf_path_getter_impl {
    () => {
        /// Returns the netCDF file's path.
        pub fn path(&self) -> PathBuf {
            self.path.clone()
        }
    };
}

/// Generates a getter for the object's Interpolation types
#[doc(hidden)]
#[macro_export]
macro_rules! interp_type_getter_impl {
    // 1D Interpolation
    (1) => {
        /// Returns the interpolation type.
        pub fn interp_type(&self) -> String {
            self.interp_type.clone()
        }
    };
    // 2D Interpolation
    (2) => {
        /// Returns the 1D interpolation type.
        pub fn interp1d_type(&self) -> String {
            self.interp1d_type.clone()
        }

        /// Returns the 2D interpolation type.
        pub fn interp2d_type(&self) -> String {
            self.interp2d_type.clone()
        }
    };
}

/// Generates a getter for the object's netCDF convention version.
#[doc(hidden)]
#[macro_export]
macro_rules! netcdf_version_getter_impl {
    () => {
        /// Returns the object's [`Version`](semver::Version).
        pub fn netcdf_version(&self) -> semver::Version {
            self.netcdf_version.clone()
        }
    };
}

/// Generates a getter of the shape of 2D numerical objects.
#[doc(hidden)]
#[macro_export]
macro_rules! shape2d_getter_impl {
    () => {
        /// Returns the the (ψ/ψp, θ) shape of the 2D arrays, depending on the state of each
        /// flux coordinate. If both coordinates are "good", they are guaranteed to be of the same
        /// length.
        pub fn shape(&self) -> (usize, usize) {
            // One of the 2 is guaranteed to be non-zero.
            let psi_len = match self.psi.state {
                NcFluxState::None => 0,
                _ => self.psi.uvalues().len(),
            };
            let psip_len = match self.psip.state {
                NcFluxState::None => 0,
                _ => self.psip.uvalues().len(),
            };
            // If they both exist, they are guaranteed to have the same length.
            let xlen = psi_len.max(psip_len);
            (xlen, self.theta_values.len())
        }
    };
}
