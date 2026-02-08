//! Representation of an equilibrium's q-factor profile.

use crate::{
    equilibrium_type_getter_impl, fluxes_state_getter_impl, fluxes_values_array_getter_impl,
    fluxes_wall_value_getter_impl, interp_type_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl, vec_to_array1D_getter_impl,
};
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::path::{Path, PathBuf};

use crate::flux::{NcFlux, NcFluxState};
use crate::{EqError, EquilibriumType, FluxCommute, Qfactor, Result};

// ===============================================================================================

/// Analytical q-factor profile of q = 1 and ψ=ψp.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
pub struct UnityQfactor {
    equilibrium_type: EquilibriumType,
}

impl UnityQfactor {
    /// Creates a new `UnityQfactor`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let qfactor = UnityQfactor::new();
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
impl FluxCommute for UnityQfactor {
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(psi)
    }

    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(psip)
    }
}

#[allow(unused_variables)]
impl Qfactor for UnityQfactor {
    fn q_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }

    fn q_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }

    fn dpsip_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }

    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(1.0)
    }
}

impl std::fmt::Debug for UnityQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("q-factor of q = 1 and ψ=ψp").finish()
    }
}

// ===============================================================================================

/// Analytical q-factor of parabolic q(ψ) profile.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
pub struct ParabolicQfactor {
    equilibrium_type: EquilibriumType,
    /// The `q` value on the magnetic axis.
    qaxis: f64,
    /// The `q` value at the wall.
    qwall: f64,
    /// The value of the toroidal flux `ψ` at the wall, in Normalized units.
    psi_wall: f64,
}

impl ParabolicQfactor {
    /// Creates a new `ParabolicQfactor`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let qfactor = ParabolicQfactor::new(1.1, 3.8, 0.45);
    /// ```
    #[allow(clippy::new_without_default, reason = "Just confuses things")]
    pub fn new(qaxis: f64, qwall: f64, psi_wall: f64) -> Self {
        assert!(psi_wall.signum() == 1.0, "`ψwall must be positive");
        Self {
            equilibrium_type: EquilibriumType::Analytical,
            qaxis,
            qwall,
            psi_wall,
        }
    }

    equilibrium_type_getter_impl!();

    /// Returns the value of `q` on the magnetic axis.
    pub fn qaxis(&self) -> f64 {
        self.qaxis
    }

    /// Returns the value of `q` on the wall.
    pub fn qwall(&self) -> f64 {
        self.qwall
    }

    /// Returns the toroidal flux's value at the wall `ψ_wall`.
    pub fn psi_wall(&self) -> f64 {
        self.psi_wall
    }

    /// Returns the poloidal flux's value at the wall `ψp_wall`.
    pub fn psip_wall(&self) -> f64 {
        self.psip_of_psi(self.psi_wall(), &mut Accelerator::new())
            .expect("If this fails, something is wrong with the formula")
    }
}

#[allow(unused_variables)]
impl FluxCommute for ParabolicQfactor {
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        let atan_arg = psi * (self.qwall - self.qaxis).sqrt() / (self.psi_wall * self.qaxis.sqrt());
        let coef = self.psi_wall / (self.qaxis * (self.qwall - self.qaxis)).sqrt();
        Ok(coef * atan_arg.atan())
    }

    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        let tan_arg = (self.qaxis * (self.qwall - self.qaxis)).sqrt() * psip / self.psi_wall;
        let coef = self.psi_wall * self.qaxis.sqrt() / (self.qwall - self.qaxis).sqrt();
        Ok(coef * tan_arg.tan())
    }
}

#[allow(unused_variables)]
// TODO: Cache reoccurring values when sure the formulas are correct.
impl Qfactor for ParabolicQfactor {
    fn q_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.qaxis + (self.qwall - self.qaxis) * (psi / self.psi_wall).powi(2))
    }

    fn q_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        let psi = self.psi_of_psip(psip, acc)?;
        // Create a new Accelerator for `psi` else the `psip` Accelerator looses its state.
        self.q_of_psi(psi, &mut Accelerator::new())
    }

    fn dpsip_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        let denom = self.qaxis * self.psi_wall.powi(2) + (self.qwall - self.qaxis) * psi.powi(2);
        Ok(self.psi_wall.powi(2) / denom)
    }

    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        let cos_arg = (self.qaxis * (self.qwall - self.qaxis)).sqrt() * psip / self.psi_wall;
        Ok(self.qaxis / cos_arg.cos().powi(2))
    }
}

impl std::fmt::Debug for ParabolicQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ParabolicQfactor: q-factor of parabolic q(ψ) profile.")
            .field("equilibrium_type", &self.equilibrium_type)
            .field("qaxis", &self.qaxis)
            .field("qwall", &self.qwall)
            .field("psi_wall", &self.psi_wall)
            .finish()
    }
}

// ===============================================================================================

/// Used to create a [`NcQfactor`].
///
/// Exists for future configuration flexibility.
pub struct NcQfactorBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`] type (case-insensitive).
    interp_type: String,
}

impl NcQfactorBuilder {
    /// Creates a new [`NcQfactorBuilder`] from a netCDF file at `path`, with `typ` interpolation
    /// type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcQfactorBuilder::new(&path, "cubic");
    /// ```
    pub fn new(path: &Path, interp_type: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type: interp_type.into(),
        }
    }

    /// Creates a new [`NcQfactor`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, "cubic").build()?;
    /// Ok::<_, EqError>(())
    /// ```
    pub fn build(self) -> Result<NcQfactor> {
        NcQfactor::build(self)
    }
}

/// Numerical q-factor profile reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcQfactorBuilder`].
pub struct NcQfactor {
    /// Path to the netCDF file.
    path: PathBuf,
    netcdf_version: semver::Version,

    equilibrium_type: EquilibriumType,
    interp_type: String,

    psi: NcFlux,
    psip: NcFlux,

    psip_of_psi_interp: Option<DynInterpolation<f64>>,
    psi_of_psip_interp: Option<DynInterpolation<f64>>,

    q_values: Vec<f64>,
    q_of_psi_interp: Option<DynInterpolation<f64>>,
    q_of_psip_interp: Option<DynInterpolation<f64>>,
}

/// Creation
impl NcQfactor {
    fn build(builder: NcQfactorBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;
        let netcdf_version = extract_version(&f)?;

        let psi = NcFlux::toroidal(&f);
        let psip = NcFlux::poloidal(&f);
        let q_values = extract_1d_array(&f, Q)?.to_vec();

        // TODO: Integrate q/ι in case one of the fluxes is missing.

        // Create interpolators, if possible
        use NcFluxState::Good;
        #[rustfmt::skip]
        let psip_of_psi_interp = match psi.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psi.uvalues(), psip.uvalues())?),
            _ => None,
        };
        #[rustfmt::skip]
        let psi_of_psip_interp = match psip.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psip.uvalues(), psi.uvalues())?),
            _ => None,
        };

        #[rustfmt::skip]
        let q_of_psi_interp = match psi.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psi.uvalues(), &q_values)?),
            _ => None,
        };
        #[rustfmt::skip]
        let q_of_psip_interp = match psip.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psip.uvalues(), &q_values)?),
            _ => None,
        };

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            psi,
            psip,
            psip_of_psi_interp,
            psi_of_psip_interp,
            q_values,
            q_of_psi_interp,
            q_of_psip_interp,
        })
    }
}

impl FluxCommute for NcQfactor {
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psip_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psi.uvalues(), self.psip.uvalues(), psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("ψp(ψ)".into())),
        }
    }

    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psi_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psip.uvalues(), self.psi.uvalues(), psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("ψ(ψp)".into())),
        }
    }
}

impl Qfactor for NcQfactor {
    fn q_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.q_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psi.uvalues(), &self.q_values, psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("q(ψ)".into())),
        }
    }

    fn q_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.q_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(self.psip.uvalues(), &self.q_values, psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("q(ψp)".into())),
        }
    }

    fn dpsip_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psip_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv(self.psi.uvalues(), self.psip.uvalues(), psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("dψp(ψ)/dψ".into())),
        }
    }

    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psi_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval_deriv(self.psip.uvalues(), self.psi.uvalues(), psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("dψ(ψp)/dψp".into())),
        }
    }
}

/// Getters
impl NcQfactor {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(1);

    /// Returns the value of q on the magnetic axis.
    pub fn qaxis(&self) -> f64 {
        self.q_values.first().copied().expect("Always exists")
    }

    /// Returns the value of q on the wall.
    pub fn qwall(&self) -> f64 {
        self.q_values.last().copied().expect("Always exists")
    }

    fluxes_wall_value_getter_impl!();
    fluxes_state_getter_impl!();
    fluxes_values_array_getter_impl!();
    vec_to_array1D_getter_impl!(q_array, q_values, q);
}

impl std::fmt::Debug for NcQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcQfactor")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("equilibrium type", &self.equilibrium_type())
            .field("interpolation type", &self.interp_type())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .field("qaxis", &self.qaxis())
            .field("qwall", &self.qwall())
            .field("ψwall", &self.psi_wall())
            .field("ψpwall", &self.psip_wall())
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_qfactor(path_str: &str) -> NcQfactor {
        let path = PathBuf::from(&path_str);
        let builder = NcQfactorBuilder::new(&path, "steffen");
        builder.build().unwrap()
    }

    /// Make sure that dψ(ψp)/dψp and q(ψ) are close enough.
    pub(super) fn test_dpsi_dpsip_q_closeness(qfactor: &impl Qfactor, psip_wall: f64, qwall: f64) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psips = Array1::linspace(0.02 * psip_wall, 0.98 * psip_wall, 100);

        let mut acc = Accelerator::new();
        for psip in psips.iter().copied() {
            assert!(
                (qfactor.q_of_psip(psip, &mut acc).unwrap()
                    - qfactor.dpsi_dpsip(psip, &mut acc).unwrap())
                .abs()
                    < qwall * 1e-4
            )
        }
    }

    /// Make sure that dψp(ψ)/dψ and i(ψ) are close enough.
    pub(super) fn test_dpsip_dpsi_iota_closeness(
        qfactor: &impl Qfactor,
        psi_wall: f64,
        qwall: f64,
    ) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psis = Array1::linspace(0.02 * psi_wall, 0.98 * psi_wall, 100);

        let mut acc = Accelerator::new();
        for psi in psis.iter().copied() {
            assert!(
                (qfactor.iota_of_psi(psi, &mut acc).unwrap()
                    - qfactor.dpsip_dpsi(psi, &mut acc).unwrap())
                .abs()
                    < qwall * 1e-4
            )
        }
    }
}

#[cfg(test)]
mod test_derivatives_closeness {
    use crate::extract::TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn unity_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = UnityQfactor::new();
        test_dpsi_dpsip_q_closeness(&qfactor, 0.45, 3.9);
    }

    #[test]
    fn unity_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = UnityQfactor::new();
        test_dpsip_dpsi_iota_closeness(&qfactor, 0.45, 3.9);
    }

    #[test]
    fn parabolic_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, 0.45);
        let psip_wall = qfactor.psip_wall();
        let qwall = qfactor.qwall();
        test_dpsi_dpsip_q_closeness(&qfactor, psip_wall, qwall);
    }

    #[test]
    fn parabolic_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, 0.45);
        let psip_wall = qfactor.psip_wall();
        let qwall = qfactor.qwall();
        test_dpsip_dpsi_iota_closeness(&qfactor, psip_wall, qwall);
    }

    #[test]
    fn nc_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        let psip_wall = qfactor.psip_wall().unwrap();
        let qwall = qfactor.qwall();
        test_dpsi_dpsip_q_closeness(&qfactor, psip_wall, qwall);
    }

    #[test]
    fn nc_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        let psi_wall = qfactor.psi_wall().unwrap();
        let qwall = qfactor.qwall();
        test_dpsip_dpsi_iota_closeness(&qfactor, psi_wall, qwall);
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::TOROIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(qfactor.psi_state(), NcFluxState::Good);
        assert_eq!(qfactor.psip_state(), NcFluxState::Bad);
        assert!(qfactor.q_of_psi_interp.is_some());
        assert!(qfactor.q_of_psip_interp.is_none());
        assert!(qfactor.psip_of_psi_interp.is_some());
        assert!(qfactor.psi_of_psip_interp.is_none());

        assert!(qfactor.psi_array().is_some());
        assert!(qfactor.psip_array().is_some());
    }

    #[test]
    fn good_psi_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        qfactor.q_of_psi(0.01, &mut acc).unwrap();
        qfactor.iota_of_psi(0.01, &mut acc).unwrap();
        qfactor.psip_of_psi(0.01, &mut acc).unwrap();
        qfactor.dpsip_dpsi(0.01, &mut acc).unwrap();
    }

    #[test]
    fn bad_psip_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EqError::UndefinedEvaluation as err;
        matches!(qfactor.q_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.iota_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.psi_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.dpsi_dpsip(0.01, &mut acc), Err(err(..)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(qfactor.psi_state(), NcFluxState::Bad);
        assert_eq!(qfactor.psip_state(), NcFluxState::Good);
        assert!(qfactor.q_of_psi_interp.is_none());
        assert!(qfactor.q_of_psip_interp.is_some());
        assert!(qfactor.psip_of_psi_interp.is_none());
        assert!(qfactor.psi_of_psip_interp.is_some());

        assert!(qfactor.psi_array().is_some());
        assert!(qfactor.psip_array().is_some());
    }

    #[test]
    fn good_psip_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        qfactor.q_of_psip(0.01, &mut acc).unwrap();
        qfactor.iota_of_psip(0.01, &mut acc).unwrap();
        qfactor.psi_of_psip(0.01, &mut acc).unwrap();
        qfactor.dpsi_dpsip(0.01, &mut acc).unwrap();
    }

    #[test]
    fn bad_psi_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EqError::UndefinedEvaluation as err;
        matches!(qfactor.q_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.iota_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.psip_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.dpsip_dpsi(0.01, &mut acc), Err(err(..)));
    }
}
