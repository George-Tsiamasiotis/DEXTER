//! Representation of an equilibrium's plasma current.

use common::array1D_getter_impl;
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::path::{Path, PathBuf};

use crate::Current;
use crate::Result;

/// Used to create a [`NcCurrent`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
pub struct NcCurrentBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`], in case-insensitive string format.
    typ: String,
}

impl NcCurrentBuilder {
    /// Creates a new [`NcCurrentBuilder`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcCurrentBuilder::new(&path, "cubic");
    /// ```
    pub fn new(path: &Path, typ: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            typ: typ.into(),
        }
    }

    /// Creates a new [`NcCurrent`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let current = NcCurrentBuilder::new(&path, "cubic").build()?;
    /// # Ok::<_, equilibrium::EqError>(())
    /// ```
    pub fn build(self) -> Result<NcCurrent> {
        NcCurrent::build(self)
    }
}

// ===============================================================================================

/// Plasma current reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcCurrentBuilder`].
#[non_exhaustive]
pub struct NcCurrent {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,

    /// The `ψp` data array.
    psip_data: Vec<f64>,
    /// The `g` data array.
    g_data: Vec<f64>,
    /// The `I` data array.
    i_data: Vec<f64>,

    /// Interpolator over the `g` values, as a function of ψp.
    g_interp: DynInterpolation<f64>,
    /// Interpolator over the `I` values, as a function of ψp.
    i_interp: DynInterpolation<f64>,
}

/// Creation
impl NcCurrent {
    /// Constructs an [`NcCurrent`] from [`NcCurrentBuilder`].
    pub(crate) fn build(builder: NcCurrentBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP_NORM)?.to_vec();
        let g_data = extract_1d_array(&f, G_NORM)?.to_vec();
        let i_data = extract_1d_array(&f, I_NORM)?.to_vec();

        let g_interp = make_interp_type(&builder.typ)?.build(&psip_data, &g_data)?;
        let i_interp = make_interp_type(&builder.typ)?.build(&psip_data, &i_data)?;

        Ok(Self {
            path: path.to_owned(),
            typ: builder.typ,
            psip_data,
            g_data,
            i_data,
            g_interp,
            i_interp,
        })
    }
}

/// Interpolation
impl Current for NcCurrent {
    fn g(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .g_interp
            .eval(&self.psip_data, &self.g_data, psip, acc)?)
    }

    fn i(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .i_interp
            .eval(&self.psip_data, &self.i_data, psip, acc)?)
    }

    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .g_interp
            .eval_deriv(&self.psip_data, &self.g_data, psip, acc)?)
    }

    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .i_interp
            .eval_deriv(&self.psip_data, &self.i_data, psip, acc)?)
    }
}

/// Getters
impl NcCurrent {
    /// Returns the netCDF file's path.
    pub fn path(&self) -> PathBuf {
        self.path.clone()
    }

    /// Returns the interpolation type.
    pub fn typ(&self) -> String {
        self.typ.clone()
    }

    /// Returns the number of data points.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.psip_data.len()
    }

    array1D_getter_impl!(psip_data, psip_data);
    array1D_getter_impl!(g_data, g_data);
    array1D_getter_impl!(i_data, i_data);
}

impl std::fmt::Debug for NcCurrent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcCurrent")
            .field("path", &self.path())
            .field("typ", &self.typ())
            .field("len", &self.len())
            .finish()
    }
}

// ===============================================================================================

/// Large Aspect Ratio approximation, with `g=1` and `I=0`.
///
/// # Example
/// ```
/// # use equilibrium::*;
/// let current = LarCurrent;
/// ```
#[derive(Debug)]
pub struct LarCurrent;

impl Current for LarCurrent {
    /// Always returns `1.0`.
    #[allow(unused_variables)]
    fn g(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }

    /// Always returns `0.0`.
    #[allow(unused_variables)]
    fn i(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }

    /// Always returns `0.0`.
    #[allow(unused_variables)]
    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }

    /// Always returns `0.0`.
    #[allow(unused_variables)]
    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }
}
