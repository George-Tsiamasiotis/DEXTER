//! Representation of an equilibrium's magnetic_field.

use crate::{
    equilibrium_type_getter_impl, fluxes_state_getter_impl, fluxes_values_array_getter_impl,
    fluxes_wall_value_getter_impl, fortran_vec_to_carray2d_impl, interp_type_getter_impl,
    netcdf_path_getter_impl, netcdf_version_getter_impl, shape2d_getter_impl,
    vec_to_array1D_getter_impl,
};
use ndarray::{Array1, Array2, Order::ColumnMajor};
use rsl_interpolation::{Accelerator, Cache, DynInterpolation2d, Interp2dType, make_interp2d_type};
use std::path::{Path, PathBuf};

use crate::flux::{NcFlux, NcFluxState};
use crate::{Bfield, EqError, EquilibriumType, Result};

// ===============================================================================================

/// Analytical Large Aspect Ratio magnetic field with B(ψ, θ) = 1 - sqrt(2ψ)cos(θ).
///
/// The LAR magnetic field can only be expressed with respect to the toroidal flux ψ.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
pub struct LarBfield {
    equilibrium_type: EquilibriumType,
}

impl LarBfield {
    /// Creates a new `LarBfield`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let bfield = LarBfield::new();
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
impl Bfield for LarBfield {
    fn b_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(1.0 - (2.0 * psi).sqrt() * theta.cos())
    }

    fn b_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Err(EqError::UndefinedEvaluation("B(ψp, θ)".into()))
    }

    fn db_dpsi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(-theta.cos() / (2.0 * psi).sqrt())
    }

    fn db_dpsip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Err(EqError::UndefinedEvaluation("dB(ψp, θ)/dψp".into()))
    }

    fn db_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok((2.0 * psi).sqrt() * theta.sin())
    }

    fn db_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Err(EqError::UndefinedEvaluation("dB(ψp, θ)/dθ".into()))
    }
}

impl std::fmt::Debug for LarBfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Large Aspect Ratio Bfield with B(ψ, θ) = 1 - sqrt(2ψ)cos(θ).")
            .finish()
    }
}

// ===============================================================================================

/// Used to create a [`NcBfield`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
pub struct NcBfieldBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 2D [`DynInterpolation2d`], in case-insensitive string format.
    interp_type: String,
}

impl NcBfieldBuilder {
    /// Creates a new [`NcBfieldBuilder`] from a netCDF file at `path`, with 2D interpolation
    /// type `interp_type`.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcBfieldBuilder::new(&path, "bicubic");
    /// ```
    pub fn new(path: &Path, interp_type: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type: interp_type.into(),
        }
    }

    /// Creates a new [`NcBfield`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// # Ok::<_, EqError>(())
    /// ```
    pub fn build(self) -> Result<NcBfield> {
        NcBfield::build(self)
    }
}

// ===============================================================================================

/// Numerical magnetic field profile reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcBfieldBuilder`].
#[non_exhaustive]
pub struct NcBfield {
    /// Path to the netCDF file.
    path: PathBuf,
    netcdf_version: semver::Version,

    equilibrium_type: EquilibriumType,
    interp_type: String,

    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,

    /// The boozer toroidal angle `θ` in [rads].
    theta_values: Vec<f64>,
    psi: NcFlux,
    psip: NcFlux,

    /// The `B` values, flattened in F order.
    b_values_fortran_flat: Vec<f64>,
    /// B(ψ, θ) interpolatior.
    b_of_psi_interp: Option<DynInterpolation2d<f64>>,
    /// B(ψp, θ) interpolatior.
    b_of_psip_interp: Option<DynInterpolation2d<f64>>,
}

/// Creation
impl NcBfield {
    /// Constructs an [`NcBfield`] from an [`NcBfieldBuilder`].
    pub(crate) fn build(builder: NcBfieldBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;
        let netcdf_version = extract_version(&f)?;

        let theta_values = extract_1d_array(&f, THETA)?.to_vec();
        let psi = NcFlux::toroidal(&f);
        let psip = NcFlux::poloidal(&f);

        let baxis = extract_scalar(&f, BAXIS)?;
        let b_values_fortran_flat = extract_2d_array(&f, B_NORM)?
            .flatten_with_order(ColumnMajor)
            .to_vec();

        // Create interpolators, if possible
        use NcFluxState::Good;
        let b_of_psi_interp = match psi.state {
            Good => Some(make_interp2d_type(&builder.interp_type)?.build(
                psi.uvalues(),
                &theta_values,
                &b_values_fortran_flat,
            )?),
            _ => None,
        };
        let b_of_psip_interp = match psip.state {
            Good => Some(make_interp2d_type(&builder.interp_type)?.build(
                psip.uvalues(),
                &theta_values,
                &b_values_fortran_flat,
            )?),
            _ => None,
        };

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            theta_values,
            psi,
            psip,
            baxis,
            b_values_fortran_flat,
            b_of_psi_interp,
            b_of_psip_interp,
        })
    }
}

impl Bfield for NcBfield {
    fn b_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.b_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(
                self.psi.uvalues(),
                &self.theta_values,
                &self.b_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("B(ψ, θ)".into())),
        }
    }

    fn b_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.b_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(
                self.psip.uvalues(),
                &self.theta_values,
                &self.b_values_fortran_flat,
                psip,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("B(ψp, θ)".into())),
        }
    }

    fn db_dpsi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.b_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv_x(
                self.psi.uvalues(),
                &self.theta_values,
                &self.b_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("dB(ψ, θ)/dψ".into())),
        }
    }

    fn db_dpsip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.b_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv_x(
                self.psip.uvalues(),
                &self.theta_values,
                &self.b_values_fortran_flat,
                psip,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("dB(ψp, θ)/dψp".into())),
        }
    }

    fn db_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.b_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv_y(
                self.psi.uvalues(),
                &self.theta_values,
                &self.b_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("dB(ψ, θ)/dθ".into())),
        }
    }

    fn db_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.b_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv_y(
                self.psip.uvalues(),
                &self.theta_values,
                &self.b_values_fortran_flat,
                psip,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("dB(ψp, θ)/dθ".into())),
        }
    }
}

/// Getters
impl NcBfield {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(1);

    /// Returns the magnetic field strength on the axis `B0` **in \[T\]**.
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    shape2d_getter_impl!();
    fluxes_wall_value_getter_impl!();
    fluxes_state_getter_impl!();
    fluxes_values_array_getter_impl!();
    vec_to_array1D_getter_impl!(theta_array, theta_values, θ);
    fortran_vec_to_carray2d_impl!(b_array, b_values_fortran_flat, B);
}

impl std::fmt::Debug for NcBfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcBfield")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("equilibrium type", &self.equilibrium_type())
            .field("interpolation type", &self.interp_type())
            .field("baxis [T]", &self.baxis)
            .field("shape (ψ/ψp, θ)", &self.shape())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_bfield(path_str: &str) -> NcBfield {
        let path = PathBuf::from(&path_str);
        let builder = NcBfieldBuilder::new(&path, "bicubic");
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
        let bfield = create_nc_bfield(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(bfield.psi_state(), NcFluxState::Good);
        assert_eq!(bfield.psip_state(), NcFluxState::Bad);
        assert!(bfield.b_of_psi_interp.is_some());
        assert!(bfield.b_of_psip_interp.is_none());

        assert!(bfield.psi_array().is_some());
        assert!(bfield.psip_array().is_some());
    }

    #[test]
    fn good_psi_evals() {
        let b = create_nc_bfield(TOROIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        b.b_of_psi(p, t, &mut a1, &mut a2, &mut c).unwrap();
        b.db_dpsi(p, t, &mut a1, &mut a2, &mut c).unwrap();
        b.db_of_psi_dtheta(p, t, &mut a1, &mut a2, &mut c).unwrap();
    }

    #[test]
    fn bad_psip_evals() {
        let b = create_nc_bfield(TOROIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        use EqError::UndefinedEvaluation as err;
        matches!(b.b_of_psip(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(b.db_dpsip(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(
            b.db_of_psip_dtheta(p, t, &mut a1, &mut a2, &mut c),
            Err(err(..))
        );
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let bfield = create_nc_bfield(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(bfield.psi_state(), NcFluxState::Bad);
        assert_eq!(bfield.psip_state(), NcFluxState::Good);
        assert!(bfield.b_of_psi_interp.is_none());
        assert!(bfield.b_of_psip_interp.is_some());

        assert!(bfield.psi_array().is_some());
        assert!(bfield.psip_array().is_some());
    }

    #[test]
    fn good_psip_evals() {
        let b = create_nc_bfield(POLOIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        b.b_of_psip(p, t, &mut a1, &mut a2, &mut c).unwrap();
        b.db_dpsip(p, t, &mut a1, &mut a2, &mut c).unwrap();
        b.db_of_psip_dtheta(p, t, &mut a1, &mut a2, &mut c).unwrap();
    }

    #[test]
    fn bad_psi_evals() {
        let b = create_nc_bfield(POLOIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        use EqError::UndefinedEvaluation as err;
        matches!(b.b_of_psi(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(b.db_dpsi(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(
            b.db_of_psi_dtheta(p, t, &mut a1, &mut a2, &mut c),
            Err(err(..))
        );
    }
}
