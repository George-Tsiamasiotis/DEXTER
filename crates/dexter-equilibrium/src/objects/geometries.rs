//! Representation of an equilibrium's general geometry.

use crate::{
    equilibrium_type_getter_impl, fluxes_state_getter_impl, fluxes_wall_value_getter_impl,
    fortran_vec_to_carray2d_impl, interp_type_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl, vec_to_array1D_getter_impl,
};
use ndarray::{Array1, Array2, Order::ColumnMajor};
use rsl_interpolation::{
    Accelerator, Cache, DynInterpolation, DynInterpolation2d, Interp2dType, InterpType,
    make_interp_type, make_interp2d_type,
};
use std::path::{Path, PathBuf};

use crate::flux::{NcFlux, NcFluxState};
use crate::{EqError, EquilibriumType, Geometry, Result};

/// Used to create an [`NcGeometry`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
pub struct NcGeometryBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`], in case-insensitive string format.
    interp1d_type: String,
    /// 2D [`DynInterpolation2d`], in case-insensitive string format.
    interp2d_type: String,
}

impl NcGeometryBuilder {
    /// Creates a new [`NcGeometryBuilder`] from a netCDF file at `path`, with 1D interpolation
    /// type `typ1d` and 2D interpolation type `typ2d`.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcGeometryBuilder::new(&path, "akima", "bicubic");
    /// ```
    pub fn new(path: &Path, interp1d_type: &str, interp2d_type: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            interp1d_type: interp1d_type.into(),
            interp2d_type: interp2d_type.into(),
        }
    }

    /// Creates a new [`NcGeometry`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let geometry = NcGeometryBuilder::new(&path, "akima", "bicubic").build()?;
    /// # Ok::<_, EqError>(())
    /// ```
    pub fn build(self) -> Result<NcGeometry> {
        NcGeometry::build(self)
    }
}

// ===============================================================================================

/// Describes the general geometry of the equilibrium.
///
/// Stores fluxes, angles and lab variables' data, and provides interpolation methods between them.
///
/// Should be created with an [`NcGeometryBuilder`].
#[non_exhaustive]
pub struct NcGeometry {
    /// Path to the netCDF file.
    path: PathBuf,
    netcdf_version: semver::Version,

    equilibrium_type: EquilibriumType,
    /// Interpolation type of the 1D quantities
    interp1d_type: String,
    /// Interpolation type of the 2D quantities
    interp2d_type: String,

    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,
    /// The horizontal position of the magnetic axis `R0` in [m].
    raxis: f64,
    /// The vertical position of the magnetic axis in [m].
    zaxis: f64,
    /// The geometrical axis (device major radius) in [m].
    rgeo: f64,

    /// The boozer toroidal angle `θ` in [rads].
    theta_values: Vec<f64>,
    psi: NcFlux,
    psip: NcFlux,

    psip_of_psi_interp: Option<DynInterpolation<f64>>,
    psi_of_psip_interp: Option<DynInterpolation<f64>>,

    /// The radial coordinate r in [m].
    r_values: Vec<f64>,
    r_of_psi_interp: Option<DynInterpolation<f64>>,
    r_of_psip_interp: Option<DynInterpolation<f64>>,
    psi_of_r_interp: Option<DynInterpolation<f64>>,
    psip_of_r_interp: Option<DynInterpolation<f64>>,

    /// The `R` coordinate in [m], flattened in F order.
    rlab_values_fortran_flat: Vec<f64>,
    /// R(ψ, θ) interpolatior.
    rlab_of_psi_interp: Option<DynInterpolation2d<f64>>,
    /// R(ψp, θ) interpolatior.
    rlab_of_psip_interp: Option<DynInterpolation2d<f64>>,

    /// The `Z` coordinate in [m], flattened in F order.
    zlab_values_fortran_flat: Vec<f64>,
    /// Z(ψ, θ) interpolatior.
    zlab_of_psi_interp: Option<DynInterpolation2d<f64>>,
    /// Z(ψp, θ) interpolatior.
    zlab_of_psip_interp: Option<DynInterpolation2d<f64>>,

    /// The VMEC output to Boozer Jacobian in [m/T], flattened in F order.
    jacobian_values_fortran_flat: Vec<f64>,
    /// J(ψ, θ) interpolatior.
    jacobian_of_psi_interp: Option<DynInterpolation2d<f64>>,
    /// J(ψp, θ) interpolatior.
    jacobian_of_psip_interp: Option<DynInterpolation2d<f64>>,
}

/// Creation
impl NcGeometry {
    /// Constructs an [`NcGeometry`] from an [`NcGeometryBuilder`]
    pub(crate) fn build(builder: NcGeometryBuilder) -> Result<Self> {
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
        let raxis = extract_scalar(&f, RAXIS)?;
        let zaxis = extract_scalar(&f, ZAXIS)?;
        let rgeo = extract_scalar(&f, RGEO)?;
        let r_values = extract_1d_array(&f, R)?.to_vec();
        let rlab_values_fortran_flat = extract_2d_array(&f, RLAB)?
            .flatten_with_order(ColumnMajor)
            .to_vec();
        let zlab_values_fortran_flat = extract_2d_array(&f, ZLAB)?
            .flatten_with_order(ColumnMajor)
            .to_vec();
        let jacobian_values_fortran_flat = extract_2d_array(&f, JACOBIAN)?
            .flatten_with_order(ColumnMajor)
            .to_vec();

        // Create interpolators, if possible
        use NcFluxState::Good;
        let psip_of_psi_interp = match psi.state {
            Good => {
                Some(make_interp_type(&builder.interp1d_type)?.build(&psi.values, &psip.values)?)
            }
            _ => None,
        };
        let psi_of_psip_interp = match psip.state {
            Good => {
                Some(make_interp_type(&builder.interp1d_type)?.build(&psip.values, &psi.values)?)
            }
            _ => None,
        };

        let r_of_psi_interp = match psi.state {
            Good => Some(make_interp_type(&builder.interp1d_type)?.build(&psi.values, &r_values)?),
            _ => None,
        };
        let r_of_psip_interp = match psip.state {
            Good => Some(make_interp_type(&builder.interp1d_type)?.build(&psip.values, &r_values)?),
            _ => None,
        };

        // r-values are guaranteed to be in inreasing order, however the corresponding flux values
        // might no exist.
        let psi_of_r_interp = match psi.state {
            NcFluxState::None => None,
            _ => Some(make_interp_type(&builder.interp1d_type)?.build(&r_values, &psi.values)?),
        };
        let psip_of_r_interp = match psip.state {
            NcFluxState::None => None,
            _ => Some(make_interp_type(&builder.interp1d_type)?.build(&r_values, &psip.values)?),
        };

        let rlab_of_psi_interp = match psi.state {
            Good => Some(make_interp2d_type(&builder.interp2d_type)?.build(
                &psi.values,
                &theta_values,
                &rlab_values_fortran_flat,
            )?),
            _ => None,
        };
        let rlab_of_psip_interp = match psip.state {
            Good => Some(make_interp2d_type(&builder.interp2d_type)?.build(
                &psip.values,
                &theta_values,
                &rlab_values_fortran_flat,
            )?),
            _ => None,
        };
        let zlab_of_psi_interp = match psi.state {
            Good => Some(make_interp2d_type(&builder.interp2d_type)?.build(
                &psi.values,
                &theta_values,
                &zlab_values_fortran_flat,
            )?),
            _ => None,
        };
        let zlab_of_psip_interp = match psip.state {
            Good => Some(make_interp2d_type(&builder.interp2d_type)?.build(
                &psip.values,
                &theta_values,
                &zlab_values_fortran_flat,
            )?),
            _ => None,
        };

        let jacobian_of_psi_interp = match psi.state {
            Good => Some(make_interp2d_type(&builder.interp2d_type)?.build(
                &psi.values,
                &theta_values,
                &jacobian_values_fortran_flat,
            )?),
            _ => None,
        };
        let jacobian_of_psip_interp = match psip.state {
            Good => Some(make_interp2d_type(&builder.interp2d_type)?.build(
                &psip.values,
                &theta_values,
                &jacobian_values_fortran_flat,
            )?),
            _ => None,
        };

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
            netcdf_version,
            path,
            interp1d_type: builder.interp1d_type,
            interp2d_type: builder.interp2d_type,
            theta_values,
            psi,
            psip,
            psi_of_psip_interp,
            psip_of_psi_interp,
            baxis,
            raxis,
            zaxis,
            rgeo,
            r_values,
            r_of_psi_interp,
            r_of_psip_interp,
            psi_of_r_interp,
            psip_of_r_interp,
            rlab_values_fortran_flat,
            rlab_of_psi_interp,
            rlab_of_psip_interp,
            zlab_values_fortran_flat,
            zlab_of_psi_interp,
            zlab_of_psip_interp,
            jacobian_values_fortran_flat,
            jacobian_of_psi_interp,
            jacobian_of_psip_interp,
        })
    }
}

/// Evaluations
impl Geometry for NcGeometry {
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psip_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(&self.psi.values, &self.psip.values, psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("ψp(ψ)".into())),
        }
    }

    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psi_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(&self.psip.values, &self.psi.values, psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("ψ(ψp)".into())),
        }
    }

    fn r_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.r_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(&self.psi.values, &self.r_values, psi, acc)?),
            None => Err(EqError::UndefinedEvaluation("r(ψ)".into())),
        }
    }

    fn r_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.r_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(&self.psip.values, &self.r_values, psip, acc)?),
            None => Err(EqError::UndefinedEvaluation("r(ψp)".into())),
        }
    }

    fn psi_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psi_of_r_interp.as_ref() {
            Some(i) => Ok(i.eval(&self.r_values, &self.psi.values, r, acc)?),
            None => Err(EqError::UndefinedEvaluation("ψ(r)".into())),
        }
    }

    fn psip_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64> {
        match self.psip_of_r_interp.as_ref() {
            Some(i) => Ok(i.eval(&self.r_values, &self.psip.values, r, acc)?),
            None => Err(EqError::UndefinedEvaluation("ψp(r)".into())),
        }
    }

    fn rlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.rlab_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(
                &self.psi.values,
                &self.theta_values,
                &self.rlab_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("R(ψ, θ)".into())),
        }
    }

    fn rlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.rlab_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(
                &self.psip.values,
                &self.theta_values,
                &self.rlab_values_fortran_flat,
                psip,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("R(ψp, θ)".into())),
        }
    }

    fn zlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.zlab_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(
                &self.psi.values,
                &self.theta_values,
                &self.zlab_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("Z(ψ, θ)".into())),
        }
    }

    fn zlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.zlab_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(
                &self.psip.values,
                &self.theta_values,
                &self.zlab_values_fortran_flat,
                psip,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("Z(ψp, θ)".into())),
        }
    }

    fn jacobian_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.jacobian_of_psi_interp.as_ref() {
            Some(i) => Ok(i.eval(
                &self.psi.values,
                &self.theta_values,
                &self.jacobian_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("J(ψ, θ)".into())),
        }
    }

    fn jacobian_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        match self.jacobian_of_psip_interp.as_ref() {
            Some(i) => Ok(i.eval(
                &self.psip.values,
                &self.theta_values,
                &self.jacobian_values_fortran_flat,
                psip,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?),
            None => Err(EqError::UndefinedEvaluation("J(ψp, θ)".into())),
        }
    }
}

/// Getters
impl NcGeometry {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(2);

    /// Returns the magnetic field strength on the axis `B0` **in \[T\]**.
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    /// Returns the horizontal position of the magnetic axis `R0` **in \[m\]**.
    pub fn raxis(&self) -> f64 {
        self.raxis
    }

    /// Returns the vertical position of the magnetic axis **in \[m\]**.
    pub fn zaxis(&self) -> f64 {
        self.zaxis
    }

    /// Returns the geometrical axis (device major radius) **in \[m\]**.
    pub fn rgeo(&self) -> f64 {
        self.rgeo
    }

    /// Returns the `r` coordinate's value at the wall **in \[m\]**, if it exists.
    pub fn rwall(&self) -> Option<f64> {
        self.r_values.last().copied()
    }

    /// Returns the the (ψ/ψp, θ) shape of the 2D arrays, depending on the state of each
    /// flux coordinate. If both coordinates are "good", they are guaranteed to be of the same
    /// length.
    pub fn shape(&self) -> (usize, usize) {
        // One of the 2 is guaranteed to be non-zero.
        let psi_len = match self.psi.state {
            NcFluxState::None => 0,
            _ => self.psi.values.len(),
        };
        let psip_len = match self.psip.state {
            NcFluxState::None => 0,
            _ => self.psip.values.len(),
        };
        // If they both exist, they are guaranteed to have the same length.
        let xlen = psi_len.max(psip_len);
        (xlen, self.theta_values.len())
    }

    fluxes_wall_value_getter_impl!();
    fluxes_state_getter_impl!();
    vec_to_array1D_getter_impl!(psi_array, psi.values, ψ);
    vec_to_array1D_getter_impl!(psip_array, psip.values, ψp);
    vec_to_array1D_getter_impl!(theta_array, theta_values, θ);
    vec_to_array1D_getter_impl!(r_array, r_values, r);
    fortran_vec_to_carray2d_impl!(rlab_array, rlab_values_fortran_flat, R);
    fortran_vec_to_carray2d_impl!(zlab_array, zlab_values_fortran_flat, Z);
    fortran_vec_to_carray2d_impl!(jacobian_array, jacobian_values_fortran_flat, J);
}

impl std::fmt::Debug for NcGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcGeometry")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("equilibrium type", &self.equilibrium_type())
            .field("1D interpolation type", &self.interp1d_type())
            .field("2D interpolation type", &self.interp2d_type())
            .field("baxis [m]", &self.baxis)
            .field("raxis [m]", &self.raxis)
            .field("zaxis [m]", &self.zaxis)
            .field("rgeo [m]", &self.rgeo)
            .field("rwall [m]", &self.rwall())
            .field("shape (ψ/ψp, θ)", &self.shape())
            .field("psi_wall", &self.psi_wall())
            .field("psip_wall", &self.psip_wall())
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_geometry(path_str: &str) -> NcGeometry {
        let path = PathBuf::from(&path_str);
        let builder = NcGeometryBuilder::new(&path, "steffen", "bicubic");
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
        let geometry = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(geometry.psi_state(), NcFluxState::Good);
        assert_eq!(geometry.psip_state(), NcFluxState::Bad);
        assert!(geometry.psip_of_psi_interp.is_some());
        assert!(geometry.r_of_psi_interp.is_some());
        assert!(geometry.rlab_of_psi_interp.is_some());
        assert!(geometry.zlab_of_psi_interp.is_some());
        assert!(geometry.jacobian_of_psi_interp.is_some());

        assert!(geometry.psi_of_psip_interp.is_none());
        assert!(geometry.r_of_psip_interp.is_none());
        assert!(geometry.rlab_of_psip_interp.is_none());
        assert!(geometry.zlab_of_psip_interp.is_none());
        assert!(geometry.jacobian_of_psip_interp.is_none());

        assert!(geometry.psi_of_r_interp.is_some());
        assert!(geometry.psip_of_r_interp.is_some());
    }

    #[test]
    fn good_psi_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        g.psip_of_psi(p, &mut a1).unwrap();
        g.r_of_psi(p, &mut a1).unwrap();
        g.rlab_of_psi(p, t, &mut a1, &mut a2, &mut c).unwrap();
        g.zlab_of_psi(p, t, &mut a1, &mut a2, &mut c).unwrap();
        g.jacobian_of_psi(p, t, &mut a1, &mut a2, &mut c).unwrap();
    }

    #[test]
    fn bad_psip_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        use EqError::UndefinedEvaluation as err;
        matches!(g.psi_of_psip(p, &mut a1), Err(err(..)));
        matches!(g.r_of_psip(p, &mut a1), Err(err(..)));
        matches!(g.rlab_of_psip(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(g.zlab_of_psip(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(
            g.jacobian_of_psip(p, t, &mut a1, &mut a2, &mut c),
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
        let geometry = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(geometry.psi_state(), NcFluxState::Bad);
        assert_eq!(geometry.psip_state(), NcFluxState::Good);
        assert!(geometry.psip_of_psi_interp.is_none());
        assert!(geometry.r_of_psi_interp.is_none());
        assert!(geometry.rlab_of_psi_interp.is_none());
        assert!(geometry.zlab_of_psi_interp.is_none());
        assert!(geometry.jacobian_of_psi_interp.is_none());

        assert!(geometry.psi_of_psip_interp.is_some());
        assert!(geometry.r_of_psip_interp.is_some());
        assert!(geometry.rlab_of_psip_interp.is_some());
        assert!(geometry.zlab_of_psip_interp.is_some());
        assert!(geometry.jacobian_of_psip_interp.is_some());

        assert!(geometry.psi_of_r_interp.is_some());
        assert!(geometry.psip_of_r_interp.is_some());
    }

    #[test]
    fn good_psip_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        g.psi_of_psip(p, &mut a1).unwrap();
        g.r_of_psip(p, &mut a1).unwrap();
        g.rlab_of_psip(p, t, &mut a1, &mut a2, &mut c).unwrap();
        g.zlab_of_psip(p, t, &mut a1, &mut a2, &mut c).unwrap();
        g.jacobian_of_psip(p, t, &mut a1, &mut a2, &mut c).unwrap();
    }

    #[test]
    fn bad_psi_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        use EqError::UndefinedEvaluation as err;
        matches!(g.psip_of_psi(p, &mut a1), Err(err(..)));
        matches!(g.r_of_psi(p, &mut a1), Err(err(..)));
        matches!(g.rlab_of_psi(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(g.zlab_of_psi(p, t, &mut a1, &mut a2, &mut c), Err(err(..)));
        matches!(
            g.jacobian_of_psi(p, t, &mut a1, &mut a2, &mut c),
            Err(err(..))
        );
    }
}
