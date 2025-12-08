//! Representation of an equilibrium's plasma current.

use std::path::PathBuf;

use common::array1D_getter_impl;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};

use ndarray::Array1;

use crate::Current;
use crate::Flux;
use crate::Result;

/// Used to create a [`NcCurrent`].
pub struct NcCurrentBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,
}

impl NcCurrentBuilder {
    /// Creates a new [`NcCurrentBuilder`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let builder = NcCurrentBuilder::new(&path, "cubic");
    /// ```
    pub fn new(path: &PathBuf, typ: &str) -> Self {
        Self {
            path: path.clone(),
            typ: typ.into(),
        }
    }

    /// Creates a new [`NcQfactor`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let current = NcCurrentBuilder::new(&path, "cubic").build()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn build(self) -> Result<NcCurrent> {
        NcCurrent::build(self)
    }
}

// ===============================================================================================

/// Plasma current reconstructed from a netCDF file.
pub struct NcCurrent {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,

    /// The `ψp` data array.
    psip_data: Vec<Flux>,
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

        let psip_data = extract_1d_array(&f, PSIP_NORM)?
            .as_standard_layout()
            .to_vec();
        let g_data = extract_1d_array(&f, G_NORM)?.as_standard_layout().to_vec();
        let i_data = extract_1d_array(&f, I_NORM)?.as_standard_layout().to_vec();

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
    fn g(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .g_interp
            .eval(&self.psip_data, &self.g_data, psip, acc)?)
    }

    fn i(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .i_interp
            .eval(&self.psip_data, &self.i_data, psip, acc)?)
    }

    fn dg_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .g_interp
            .eval_deriv(&self.psip_data, &self.g_data, psip, acc)?)
    }

    fn di_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .i_interp
            .eval_deriv(&self.psip_data, &self.i_data, psip, acc)?)
    }
}

// Getters
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

    array1D_getter_impl!(psip_data, psip_data, Flux);
    array1D_getter_impl!(g_data, g_data, f64);
    array1D_getter_impl!(i_data, i_data, f64);
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::extract::STUB_NETCDF_PATH;

    fn create_nc_current() -> NcCurrent {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        let typ = "steffen";
        NcCurrentBuilder::new(&path, typ).build().unwrap()
    }

    #[test]
    fn test_current_creation() {
        let c = create_nc_current();
        let _ = format!("{c:?}");
    }

    #[test]
    fn test_getters() {
        let c = create_nc_current();
        c.path();
        c.typ();
        c.len();

        assert_eq!(c.psip_data().ndim(), 1);
        assert_eq!(c.g_data().ndim(), 1);
        assert_eq!(c.i_data().ndim(), 1);
    }

    #[test]
    fn test_extraction_methods() {}

    #[test]
    fn test_spline_evaluation() {
        let c = create_nc_current();
        let mut acc = Accelerator::new();

        let psip = 0.015;
        c.g(psip, &mut acc).unwrap();
        c.i(psip, &mut acc).unwrap();
        c.di_dpsip(psip, &mut acc).unwrap();
        c.dg_dpsip(psip, &mut acc).unwrap();
    }
}
