//! Representation of an equilibrium's plasma current.

use crate::{fluxes_wall_value_getter_impl, vec_to_array1D_getter_impl};
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::path::{Path, PathBuf};

use crate::{Current, NcFlux, NcFluxState, Result};

/// Used to create a [`NcCurrent`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
pub struct NcCurrentBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`] type (case-insensitive).
    typ: String,
}

impl NcCurrentBuilder {
    /// Creates a new [`NcCurrentBuilder`] from a netCDF file at `path`, with `typ` interpolation
    /// type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
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
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let current = NcCurrentBuilder::new(&path, "cubic").build()?;
    /// Ok::<_, EqError>(())
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

    psi: NcFlux,
    psip: NcFlux,

    g_values: Vec<f64>,
    g_of_psi_interp: Option<DynInterpolation<f64>>,
    g_of_psip_interp: Option<DynInterpolation<f64>>,

    i_values: Vec<f64>,
    i_of_psi_interp: Option<DynInterpolation<f64>>,
    i_of_psip_interp: Option<DynInterpolation<f64>>,
}

/// Creation
impl NcCurrent {
    /// Constructs an [`NcCurrent`] from [`NcCurrentBuilder`].
    fn build(builder: NcCurrentBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;

        let g_values = extract_1d_array(&f, G_NORM)?.to_vec();
        let i_values = extract_1d_array(&f, I_NORM)?.to_vec();
        let psi = NcFlux::toroidal(&f);
        let psip = NcFlux::poloidal(&f);

        // Create interpolators, if possible
        use NcFluxState::Good;
        let g_of_psi_interp = match psi.state {
            Good => Some(make_interp_type(&builder.typ)?.build(&psi.values, &g_values)?),
            _ => None,
        };
        let i_of_psi_interp = match psi.state {
            Good => Some(make_interp_type(&builder.typ)?.build(&psi.values, &i_values)?),
            _ => None,
        };

        let g_of_psip_interp = match psip.state {
            Good => Some(make_interp_type(&builder.typ)?.build(&psip.values, &g_values)?),
            _ => None,
        };
        let i_of_psip_interp = match psip.state {
            Good => Some(make_interp_type(&builder.typ)?.build(&psip.values, &i_values)?),
            _ => None,
        };

        Ok(Self {
            path,
            typ: builder.typ,
            psi,
            psip,
            g_values,
            i_values,
            g_of_psi_interp,
            g_of_psip_interp,
            i_of_psi_interp,
            i_of_psip_interp,
        })
    }
}

/// Evaluations
impl Current for NcCurrent {
    fn g_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psi.state == NcFluxState::Good);
        match &self.g_of_psi_interp {
            Some(interp) => Ok(interp.eval(&self.psi.values, &self.g_values, psi, acc)?),
            None => panic!("Flux coordinate ψ is not well defined."),
        }
    }

    fn g_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psip.state == NcFluxState::Good);
        match &self.g_of_psip_interp {
            Some(interp) => Ok(interp.eval(&self.psip.values, &self.g_values, psip, acc)?),
            None => panic!("Flux coordinate ψp is not well defined."),
        }
    }

    fn i_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psi.state == NcFluxState::Good);
        match &self.i_of_psi_interp {
            Some(interp) => Ok(interp.eval(&self.psi.values, &self.i_values, psi, acc)?),
            None => panic!("Flux coordinate ψ is not well defined."),
        }
    }

    fn i_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psip.state == NcFluxState::Good);
        match &self.i_of_psip_interp {
            Some(interp) => Ok(interp.eval(&self.psip.values, &self.i_values, psip, acc)?),
            None => panic!("Flux coordinate ψp is not well defined."),
        }
    }

    fn dg_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psi.state == NcFluxState::Good);
        match &self.g_of_psi_interp {
            Some(interp) => Ok(interp.eval_deriv(&self.psi.values, &self.g_values, psi, acc)?),
            None => panic!("Flux coordinate ψ is not well defined."),
        }
    }

    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psip.state == NcFluxState::Good);
        match &self.g_of_psip_interp {
            Some(interp) => Ok(interp.eval_deriv(&self.psip.values, &self.g_values, psip, acc)?),
            None => panic!("Flux coordinate ψp is not well defined."),
        }
    }

    fn di_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psi.state == NcFluxState::Good);
        match &self.i_of_psi_interp {
            Some(interp) => Ok(interp.eval_deriv(&self.psi.values, &self.i_values, psi, acc)?),
            None => panic!("Flux coordinate ψ is not well defined."),
        }
    }

    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        debug_assert!(self.psip.state == NcFluxState::Good);
        match &self.i_of_psip_interp {
            Some(interp) => Ok(interp.eval_deriv(&self.psip.values, &self.i_values, psip, acc)?),
            None => panic!("Flux coordinate ψp is not well defined."),
        }
    }
}

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
        self.g_values.len()
    }

    fluxes_wall_value_getter_impl!();
    vec_to_array1D_getter_impl!(psi_array, psi.values);
    vec_to_array1D_getter_impl!(psip_array, psip.values);
    vec_to_array1D_getter_impl!(g_array, g_values);
    vec_to_array1D_getter_impl!(i_array, i_values);
}

impl std::fmt::Debug for NcCurrent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcCurrent")
            .field("path", &self.path())
            .field("typ", &self.typ())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .finish()
    }
}
