//! Representation of an equilibrium's q-factor.

use common::array1D_getter_impl;
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::path::{Path, PathBuf};

use crate::Qfactor;
use crate::Result;

/// Used to create an [`NcQfactor`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
pub struct NcQfactorBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`], in case-insensitive string format.
    typ: String,
}

impl NcQfactorBuilder {
    /// Creates a new [`NcQfactorBuilder`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcQfactorBuilder::new(&path, "cubic");
    /// ```
    pub fn new(path: &Path, typ: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            typ: typ.into(),
        }
    }

    /// Creates a new [`NcQfactor`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, "cubic").build()?;
    /// # Ok::<_, equilibrium::EqError>(())
    /// ```
    pub fn build(self) -> Result<NcQfactor> {
        NcQfactor::build(self)
    }
}

// ===============================================================================================

/// q-factor reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcQfactorBuilder`].
#[non_exhaustive]
pub struct NcQfactor {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`], in case-insensitive string format.
    typ: String,

    /// The `ψp` data array.
    psip_data: Vec<f64>,
    /// The `q` data array.
    q_data: Vec<f64>,
    /// The `ψ` data array.
    psi_data: Vec<f64>,

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

        let psip_data = extract_1d_array(&f, PSIP_NORM)?.to_vec();
        let psi_data = extract_1d_array(&f, PSI_NORM)?.to_vec();
        let q_data = extract_1d_array(&f, Q)?.to_vec();

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
    fn q(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .q_interp
            .eval(&self.psip_data, &self.q_data, psip, acc)?)
    }

    fn psi(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .psi_interp
            .eval(&self.psip_data, &self.psi_data, psip, acc)?)
    }

    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
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

    array1D_getter_impl!(psip_data, psip_data);
    array1D_getter_impl!(psi_data, psi_data);
    array1D_getter_impl!(q_data, q_data);
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

impl Clone for NcQfactor {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            typ: self.typ.clone(),
            psip_data: self.psip_data.clone(),
            q_data: self.q_data.clone(),
            psi_data: self.psi_data.clone(),
            q_interp: make_interp_type(&self.typ())
                .expect("already succeeded")
                .build(&self.psip_data().to_vec(), &self.q_data().to_vec())
                .expect("already succeeded"),
            psi_interp: make_interp_type(&self.typ())
                .expect("already succeeded")
                .build(&self.psip_data().to_vec(), &self.psi_data().to_vec())
                .expect("already succeeded"),
        }
    }
}

// ===============================================================================================

/// Analytical q-factor profile with `q=1`.
///
/// # Example
/// ```
/// # use equilibrium::*;
/// let qfactor = UnityQfactor;
/// ```
#[derive(Debug, Clone)]
pub struct UnityQfactor;

impl Qfactor for UnityQfactor {
    /// Always returns `1.0`.
    #[allow(unused_variables)]
    fn q(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }

    /// Always returns `psip`.
    #[allow(unused_variables)]
    fn psi(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(psip)
    }

    /// Always returns `1.0`.
    #[allow(unused_variables)]
    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }
}
