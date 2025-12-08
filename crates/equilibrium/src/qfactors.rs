//! Representation of an equilibrium's q-factor.

use std::path::PathBuf;

use common::array1D_getter_impl;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};

use ndarray::Array1;

use crate::Flux;
use crate::Qfactor;
use crate::Result;

/// Used to create a [`NcQfactor`].
pub struct NcQfactorBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,
}

impl NcQfactorBuilder {
    /// Creates a new [`NcQfactorBuilder`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let builder = NcQfactorBuilder::new(&path, "cubic");
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
    /// let qfactor = NcQfactorBuilder::new(&path, "cubic").build()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn build(self) -> Result<NcQfactor> {
        NcQfactor::build(self)
    }
}

// ===============================================================================================

/// q-factor reconstructed from a netCDF file.
pub struct NcQfactor {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,

    /// The `ψp` data array.
    psip_data: Vec<Flux>,
    /// The `q` data array.
    q_data: Vec<f64>,
    /// The `ψ` data array.
    psi_data: Vec<Flux>,

    /// Interpolator over the `q` values, as a function of ψp.
    q_interp: DynInterpolation<f64>,
    /// Interpolator over the `ψ` values, as a function of ψp.
    psi_interp: DynInterpolation<f64>,
}

/// Creation
impl NcQfactor {
    /// Constructs an [`NcQfactor`] from [`NcQfactorBuilder`].
    pub(crate) fn build(builder: NcQfactorBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP_NORM)?
            .as_standard_layout()
            .to_vec();
        let psi_data = extract_1d_array(&f, PSI_NORM)?
            .as_standard_layout()
            .to_vec();
        let q_data = extract_1d_array(&f, Q)?.as_standard_layout().to_vec();

        let q_interp = make_interp_type(&builder.typ)?.build(&psip_data, &q_data)?;
        let psi_interp = make_interp_type(&builder.typ)?.build(&psip_data, &psi_data)?;

        Ok(Self {
            path: path.to_owned(),
            typ: builder.typ,
            psip_data,
            q_data,
            psi_data,
            q_interp,
            psi_interp,
        })
    }
}

/// Interpolation
impl Qfactor for NcQfactor {
    fn q(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .q_interp
            .eval(&self.psip_data, &self.q_data, psip, acc)?)
    }

    fn psi(&self, psip: Flux, acc: &mut Accelerator) -> Result<Flux> {
        Ok(self
            .psi_interp
            .eval(&self.psip_data, &self.psi_data, psip, acc)?)
    }

    fn dpsi_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .psi_interp
            .eval_deriv(&self.psip_data, &self.psi_data, psip, acc)?)
    }
}

/// Getters
impl NcQfactor {
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
    array1D_getter_impl!(psi_data, psi_data, Flux);
    array1D_getter_impl!(q_data, q_data, f64);
}

impl std::fmt::Debug for NcQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcQfactor")
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

    fn create_nc_qfactor() -> NcQfactor {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        let typ = "steffen";
        NcQfactorBuilder::new(&path, typ).build().unwrap()
    }

    #[test]
    fn test_qfactor_creation() {
        let q = create_nc_qfactor();
        let _ = format!("{q:?}");
    }

    #[test]
    fn test_getters() {
        let q = create_nc_qfactor();
        q.path();
        q.typ();
        q.len();

        assert_eq!(q.psip_data().ndim(), 1);
        assert_eq!(q.q_data().ndim(), 1);
        assert_eq!(q.psi_data().ndim(), 1);
    }

    #[test]
    fn test_spline_evaluation() {
        let q = create_nc_qfactor();
        let mut acc = Accelerator::new();

        let psip = 0.015;
        q.q(psip, &mut acc).unwrap();
        q.psi(psip, &mut acc).unwrap();
        q.dpsi_dpsip(psip, &mut acc).unwrap();
    }
}
