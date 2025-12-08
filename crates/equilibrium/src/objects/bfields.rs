//! Representation of an equilibrium's magnetic field.

use common::array1D_getter_impl;
use ndarray::{Array1, Array2};
use rsl_interpolation::{Accelerator, Cache, DynInterpolation2d, Interp2dType, make_interp2d_type};
use std::f64::consts::TAU;
use std::path::PathBuf;

use crate::Bfield;
use crate::Result;
use crate::fortran_vec_to_carray2d_impl;
use crate::{Flux, Radians};

/// Used to create a [`NcBfield`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
pub struct NcBfieldBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 2D [`DynInterpolation2d`], in case-insensitive string format.
    typ: String,
}

impl NcBfieldBuilder {
    /// Creates a new [`NcBfieldBuilder`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::bfields;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = bfields::NcBfieldBuilder::new(&path, "bicubic");
    /// ```
    pub fn new(path: &PathBuf, typ: &str) -> Self {
        Self {
            path: path.clone(),
            typ: typ.into(),
        }
    }

    /// Creates a new [`NcBfield`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::bfields;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let bfield = bfields::NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// # Ok::<_, equilibrium::EqError>(())
    /// ```
    pub fn build(self) -> Result<NcBfield> {
        NcBfield::build(self)
    }
}

// ===============================================================================================

/// Magnetic field reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcBfieldBuilder`].
#[non_exhaustive]
pub struct NcBfield {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 2D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.Interp2dType.html#implementors
    typ: String,

    /// The `ψp` data array.
    psip_data: Vec<Flux>,
    /// The `θ` data array.
    theta_data: Vec<Radians>,
    /// The `B` data array, flattened in C order.
    b_data: Vec<f64>,

    /// Interpolator over the `b` values, as a function of ψp, θ.
    b_interp: DynInterpolation2d<f64>,
}

/// Creation
impl NcBfield {
    /// Constructs an [`NcBfield`] from [`NcBfieldBuilder`].
    pub(crate) fn build(builder: NcBfieldBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP_NORM)?.to_vec();
        let theta_data = extract_1d_array(&f, THETA)?.to_vec();
        let b_data = extract_2d_array(&f, B_NORM)?;

        // `Interp.za` must be in Fortran order.
        let b_data_fortran_flat = b_data
            .flatten_with_order(ndarray::Order::ColumnMajor)
            .to_vec();

        let b_interp = make_interp2d_type(&builder.typ)?.build(
            &psip_data,
            &theta_data,
            &b_data_fortran_flat,
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: builder.typ,
            psip_data,
            theta_data,
            b_data: b_data_fortran_flat,
            b_interp,
        })
    }
}

/// Interpolation
#[rustfmt::skip] // pretty!
impl Bfield for NcBfield {
    fn b(&self, psip: Flux, theta: Radians, xacc: &mut Accelerator, yacc: &mut Accelerator, cache: &mut Cache<f64>) -> Result<f64> {
        Ok(self.b_interp.eval(&self.psip_data,&self.theta_data,&self.b_data, psip, theta.rem_euclid(TAU), xacc, yacc, cache)?)
    }

    fn db_dpsip(&self, psip: Flux, theta: Radians, xacc: &mut Accelerator, yacc: &mut Accelerator, cache: &mut Cache<f64>) -> Result<f64> {
        Ok(self.b_interp.eval_deriv_x(&self.psip_data, &self.theta_data, &self.b_data, psip, theta.rem_euclid(TAU), xacc, yacc, cache)?)
    }

    fn db_dtheta(&self, psip: Flux, theta: Radians, xacc: &mut Accelerator, yacc: &mut Accelerator, cache: &mut Cache<f64>) -> Result<f64> {
        Ok(self.b_interp.eval_deriv_y(&self.psip_data, &self.theta_data, &self.b_data, psip, theta.rem_euclid(TAU), xacc, yacc, cache,)?)
    }
}

/// Getters
impl NcBfield {
    /// Returns the netCDF file's path.
    pub fn path(&self) -> PathBuf {
        self.path.clone()
    }

    /// Returns the interpolation type.
    pub fn typ(&self) -> String {
        self.typ.clone()
    }

    /// Returns the shape of the `b` array.
    pub fn shape(&self) -> (usize, usize) {
        (self.psip_data.len(), self.theta_data.len())
    }

    array1D_getter_impl!(psip_data, psip_data);
    array1D_getter_impl!(theta_data, theta_data);
    fortran_vec_to_carray2d_impl!(b_data, b_data, B);
}

impl std::fmt::Debug for NcBfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcBfield")
            .field("path", &self.path())
            .field("typ", &self.typ())
            .field("shape", &self.shape())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::extract::STUB_TEST_NETCDF_PATH;

    /// Uninitialized Cache bug,
    ///
    /// This caused NaNs appearing when evaluating with ψp < psip_data[1] AND θ < theta_data[1]
    /// with a newly created Cache.
    ///
    /// Fixed in rsl-interpolation 0.1.16-pre.3
    #[test]
    fn test_nan_near_first_flux_surface_theta0() {
        let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
        let typ = "bicubic";
        let builder = NcBfieldBuilder::new(&path, typ);
        let bfield = builder.build().unwrap();
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();

        let psip = 1e-10;
        let theta = 1e-10;

        let val = bfield
            .b(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        assert!(val.is_finite());
    }
}
