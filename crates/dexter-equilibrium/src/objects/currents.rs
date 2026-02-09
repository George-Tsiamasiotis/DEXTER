//! Representation of an equilibrium's plasma current.

use crate::{
    equilibrium_type_getter_impl, fluxes_state_getter_impl, fluxes_values_array_getter_impl,
    fluxes_wall_value_getter_impl, interp_type_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl, vec_to_array1D_getter_impl,
};
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::path::{Path, PathBuf};

use crate::flux::{NcFlux, NcFluxState};
use crate::{Current, EqError, EquilibriumType, Result};

// ===============================================================================================

/// Analytical Large Aspect Ratio Current with g=1 and I=0.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
pub struct LarCurrent {
    equilibrium_type: EquilibriumType,
}

impl LarCurrent {
    /// Creates a new `LarCurrent`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let current = LarCurrent::new();
    /// ```
    #[allow(clippy::new_without_default, reason = "Just confuses things")]
    pub fn new() -> Self {
        Self {
            equilibrium_type: EquilibriumType::Analytical,
        }
    }

    equilibrium_type_getter_impl!();
}

#[allow(unused_variables)]
impl Current for LarCurrent {
    fn g_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }

    fn g_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }

    fn i_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }

    fn i_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }

    fn dg_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }

    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }

    fn di_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }

    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(0.0)
    }
}

impl std::fmt::Debug for LarCurrent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Large Aspect Ratio Current with g=1 and I=0")
            .finish()
    }
}

// ===============================================================================================

/// Used to create a [`NcCurrent`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
pub struct NcCurrentBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`] type (case-insensitive).
    interp_type: String,
}

impl NcCurrentBuilder {
    /// Creates a new [`NcCurrentBuilder`] from a netCDF file at `path`, with `interp_type`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcCurrentBuilder::new(&path, "cubic");
    /// ```
    pub fn new(path: &Path, interp_type: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type: interp_type.into(),
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

/// Numerical plasma current reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcCurrentBuilder`].
#[non_exhaustive]
pub struct NcCurrent {
    /// Path to the netCDF file.
    path: PathBuf,
    netcdf_version: semver::Version,

    equilibrium_type: EquilibriumType,
    interp_type: String,

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
        let netcdf_version = extract_version(&f)?;

        let psi = NcFlux::toroidal(&f);
        let psip = NcFlux::poloidal(&f);
        let g_values = extract_1d_array(&f, G_NORM)?.to_vec();
        let i_values = extract_1d_array(&f, I_NORM)?.to_vec();

        // Create interpolators, if possible
        use NcFluxState::Good;
        #[rustfmt::skip]
        let g_of_psi_interp = match psi.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psi.uvalues(), &g_values)?),
            _ => None,
        };
        #[rustfmt::skip]
        let i_of_psi_interp = match psi.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psi.uvalues(), &i_values)?),
            _ => None,
        };

        #[rustfmt::skip]
        let g_of_psip_interp = match psip.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psip.uvalues(), &g_values)?),
            _ => None,
        };
        #[rustfmt::skip]
        let i_of_psip_interp = match psip.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psip.uvalues(), &i_values)?),
            _ => None,
        };

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
            netcdf_version,
            path,
            interp_type: builder.interp_type,
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

impl Current for NcCurrent {
    fn g_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.g_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psi.uvalues(), &self.g_values, psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("g(ψ)".into())),
        }
    }

    fn g_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.g_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psip.uvalues(), &self.g_values, psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("g(ψp)".into())),
        }
    }

    fn i_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.i_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psi.uvalues(), &self.i_values, psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("I(ψ)".into())),
        }
    }

    fn i_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.i_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psip.uvalues(), &self.i_values, psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("I(ψp)".into())),
        }
    }

    fn dg_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.g_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv(self.psi.uvalues(), &self.g_values, psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("dg(ψ)/dψ".into())),
        }
    }

    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.g_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv(self.psip.uvalues(), &self.g_values, psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("dg(ψp)/dψp".into())),
        }
    }

    fn di_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.i_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv(self.psi.uvalues(), &self.i_values, psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("dI(ψ)/dψ".into())),
        }
    }

    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.i_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv(self.psip.uvalues(), &self.i_values, psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("dI(ψp)/dψp".into())),
        }
    }
}

/// Getters
impl NcCurrent {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(1);
    fluxes_wall_value_getter_impl!();
    fluxes_state_getter_impl!();
    fluxes_values_array_getter_impl!();
    vec_to_array1D_getter_impl!(g_array, g_values, g);
    vec_to_array1D_getter_impl!(i_array, i_values, I);
}

impl std::fmt::Debug for NcCurrent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcCurrent")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("equilibrium type", &self.equilibrium_type())
            .field("interpolation type", &self.interp_type())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_current(path_str: &str) -> NcCurrent {
        let path = PathBuf::from(path_str);
        let builder = NcCurrentBuilder::new(&path, "steffen");
        builder.build().unwrap()
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::TOROIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let current = create_nc_current(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(current.psi_state(), NcFluxState::Good);
        assert_eq!(current.psip_state(), NcFluxState::Bad);
        assert!(current.g_of_psi_interp.is_some());
        assert!(current.i_of_psi_interp.is_some());
        assert!(current.g_of_psip_interp.is_none());
        assert!(current.i_of_psip_interp.is_none());

        assert!(current.psi_array().is_some());
        assert!(current.psip_array().is_some());
    }

    #[test]
    fn good_psi_evals() {
        let current = create_nc_current(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        current.g_of_psi(0.01, &mut acc).unwrap();
        current.i_of_psi(0.01, &mut acc).unwrap();
        current.dg_dpsi(0.01, &mut acc).unwrap();
        current.di_dpsi(0.01, &mut acc).unwrap();
    }

    #[test]
    fn bad_psip_evals() {
        let current = create_nc_current(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EqError::UndefinedEvaluation as err;
        matches!(current.g_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(current.i_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(current.dg_dpsip(0.01, &mut acc), Err(err(..)));
        matches!(current.di_dpsip(0.01, &mut acc), Err(err(..)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let current = create_nc_current(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(current.psi_state(), NcFluxState::Bad);
        assert_eq!(current.psip_state(), NcFluxState::Good);
        assert!(current.g_of_psi_interp.is_none());
        assert!(current.i_of_psi_interp.is_none());
        assert!(current.g_of_psip_interp.is_some());
        assert!(current.i_of_psip_interp.is_some());

        assert!(current.psi_array().is_some());
        assert!(current.psip_array().is_some());
    }

    #[test]
    fn good_psip_evals() {
        let current = create_nc_current(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        current.g_of_psip(0.01, &mut acc).unwrap();
        current.i_of_psip(0.01, &mut acc).unwrap();
        current.dg_dpsip(0.01, &mut acc).unwrap();
        current.di_dpsip(0.01, &mut acc).unwrap();
    }

    #[test]
    fn bad_psi_evals() {
        let current = create_nc_current(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EqError::UndefinedEvaluation as err;
        matches!(current.g_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(current.i_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(current.dg_dpsi(0.01, &mut acc), Err(err(..)));
        matches!(current.di_dpsi(0.01, &mut acc), Err(err(..)));
    }
}
