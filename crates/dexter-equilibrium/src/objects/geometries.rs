//! Representation of an equilibrium's general geometry.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    debug_assert_non_negative_r, equilibrium_type_getter_impl, fluxes_state_getter_impl,
    fluxes_values_array_getter_impl, fluxes_wall_value_getter_impl, fortran_vec_to_carray2d_impl,
    interp_type_getter_impl, netcdf_path_getter_impl, netcdf_version_getter_impl,
    shape2d_getter_impl,
};
use dexter_common::vec_to_array1D_getter_impl;
use ndarray::{Array1, Array2, Order::ColumnMajor};
use rsl_interpolation::{
    Accelerator, Cache, DynInterpolation, DynInterpolation2d, Interp2dType, InterpType,
    make_interp_type, make_interp2d_type,
};
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::objects::nc_flux::{NcFlux, NcFluxState};
use crate::{EqError, EvalError};
use crate::{EquilibriumType, FluxCommute, Geometry};

// ===============================================================================================

/// Analytical Large Aspect Ratio Geometry of a circular device.
///
/// The definitions are not very strict at the moment.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
#[derive(Debug)]
pub struct LarGeometry {
    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,
    /// The horizontal position of the magnetic axis `R0` in [m].
    raxis: f64,
    /// The minor radius `rwall` in [m].
    rwall: f64,
    /// The toroidal flux's value at the wall `ψ_wall` in Normalized units.
    psi_wall: f64,
}

impl LarGeometry {
    /// Creates a new `LarCurrent`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let geometry = LarGeometry::new(2.0, 1.75, 0.5);
    /// ```
    #[must_use]
    pub fn new(baxis: f64, raxis: f64, rwall: f64) -> Self {
        let psi_wall_si = baxis * rwall.powi(2) / 2.0;
        let psi_wall = psi_wall_si / (baxis * raxis.powi(2));
        Self {
            equilibrium_type: EquilibriumType::Analytical,
            baxis,
            raxis,
            rwall,
            psi_wall,
        }
    }
}

impl Geometry for LarGeometry {
    fn r_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!((2.0 * psi).sqrt()))
    }

    fn r_of_psip(&self, _: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "r(ψp) (defined through q)".into(),
        ))
    }

    fn psi_of_r(&self, r: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_r!(r);
        Ok(debug_assert_is_finite!(r.powi(2) / 2.0))
    }

    fn psip_of_r(&self, _: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "ψp(r) (defined through q)".into(),
        ))
    }

    fn rlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        _: &mut Accelerator,
        _: &mut Accelerator,
        _: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(
            self.raxis + self.raxis * (2.0 * psi).sqrt() * theta.cos()
        ))
    }

    fn rlab_of_psip(
        &self,
        _: f64,
        _: f64,
        _: &mut Accelerator,
        _: &mut Accelerator,
        _: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "R(ψp, θ) (defined through q)".into(),
        ))
    }

    fn zlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        _: &mut Accelerator,
        _: &mut Accelerator,
        _: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(
            self.raxis * (2.0 * psi).sqrt() * theta.sin()
        ))
    }

    fn zlab_of_psip(
        &self,
        _: f64,
        _: f64,
        _: &mut Accelerator,
        _: &mut Accelerator,
        _: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "Z(ψp, θ) (defined through q)".into(),
        ))
    }

    fn jacobian_of_psi(
        &self,
        _: f64,
        _: f64,
        _: &mut Accelerator,
        _: &mut Accelerator,
        _: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "J(ψ, θ) (defined through q, g, I and B)".into(),
        ))
    }

    fn jacobian_of_psip(
        &self,
        _: f64,
        _: f64,
        _: &mut Accelerator,
        _: &mut Accelerator,
        _: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "J(ψp, θ) (defined through q, g, I and B)".into(),
        ))
    }

    fn rlab_wall(&self) -> Array1<f64> {
        use core::f64::consts::PI;
        let arr = Array1::linspace(0.0, 2.0 * PI, 1000);
        let mut acc2 = Accelerator::new();
        let mut acc1 = Accelerator::new();
        let mut cache = Cache::new();
        arr.mapv(|theta| {
            match self.rlab_of_psi(self.psi_wall, theta, &mut acc2, &mut acc1, &mut cache) {
                Ok(rlab_wall_value) => rlab_wall_value,
                Err(_) => unreachable!("Expression is analytical, cannot fail"),
            }
        })
    }

    fn zlab_wall(&self) -> Array1<f64> {
        use core::f64::consts::PI;
        let arr = Array1::linspace(0.0, 2.0 * PI, 1000);
        let mut acc1 = Accelerator::new();
        let mut acc2 = Accelerator::new();
        let mut cache = Cache::new();
        arr.mapv(|theta| {
            match self.zlab_of_psi(self.psi_wall, theta, &mut acc1, &mut acc2, &mut cache) {
                Ok(zlab_wall_value) => zlab_wall_value,
                Err(_) => unreachable!("Expression is analytical, cannot fail"),
            }
        })
    }
}

impl LarGeometry {
    equilibrium_type_getter_impl!();

    /// Returns the magnetic field strength on the axis `B0` in **\[T\]**.
    #[must_use]
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    /// Returns the horizontal position of the magnetic axis `R0` in **\[m\]**.
    #[must_use]
    pub fn raxis(&self) -> f64 {
        self.raxis
    }

    /// Returns the vertical position of the magnetic axis **in \[m\]**.
    #[must_use]
    pub fn zaxis(&self) -> f64 {
        0.0
    }

    /// Returns the geometrical axis (device major radius) **in \[m\]**.
    #[must_use]
    pub fn rgeo(&self) -> f64 {
        self.raxis
    }

    /// Returns the `r` coordinate's value at the wall in **\[m\]**.
    #[must_use]
    pub fn rwall(&self) -> f64 {
        self.rwall
    }

    /// Returns the todoidal flux's `ψ` value at the wall.
    #[must_use]
    pub fn psi_wall(&self) -> f64 {
        self.psi_wall
    }
}

// ===============================================================================================

/// Used to create an [`NcGeometry`].
///
/// Exists for future configuration flexibility.
#[non_exhaustive]
#[derive(Debug)]
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
    /// type `interp1d_type` and 2D interpolation type `interp2d_type`.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcGeometryBuilder::new(&path, "akima", "bicubic");
    /// ```
    #[must_use]
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
    ///
    /// # Errors
    ///
    /// Returns an [`EqError`] if it fails to build the [`NcGeometry`].
    pub fn build(self) -> Result<NcGeometry, EqError> {
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
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,

    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// Interpolation type of the 1D quantities.
    interp1d_type: String,
    /// Interpolation type of the 2D quantities.
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
    /// The toroidal flux coordinate.
    psi: NcFlux,
    /// The poloidal flux coordinate.
    psip: NcFlux,

    /// The `ψp(ψ)` interpolator.
    psip_of_psi_interp: Option<DynInterpolation<f64>>,
    /// The `ψ(ψp)` interpolator.
    psi_of_psip_interp: Option<DynInterpolation<f64>>,

    /// The radial coordinate r in [m].
    r_values: Vec<f64>,
    /// The `r(ψ)` interpolator.
    r_of_psi_interp: Option<DynInterpolation<f64>>,
    /// The `r(ψp)` interpolator.
    r_of_psip_interp: Option<DynInterpolation<f64>>,
    /// The `ψ(r)` interpolator.
    psi_of_r_interp: Option<DynInterpolation<f64>>,
    /// The `ψp(r)` interpolator.
    psip_of_r_interp: Option<DynInterpolation<f64>>,

    /// The `R` coordinate in [m], flattened in F order.
    rlab_values_fortran_flat: Vec<f64>,
    /// `R(ψ, θ)` interpolator.
    rlab_of_psi_interp: Option<DynInterpolation2d<f64>>,
    /// `R(ψp, θ)` interpolator.
    rlab_of_psip_interp: Option<DynInterpolation2d<f64>>,

    /// The `Z` coordinate in [m], flattened in F order.
    zlab_values_fortran_flat: Vec<f64>,
    /// `Z(ψ, θ)` interpolator.
    zlab_of_psi_interp: Option<DynInterpolation2d<f64>>,
    /// `Z(ψp, θ)` interpolator.
    zlab_of_psip_interp: Option<DynInterpolation2d<f64>>,

    /// The VMEC output to Boozer Jacobian in [m/T], flattened in F order.
    jacobian_values_fortran_flat: Vec<f64>,
    /// `J(ψ, θ)` interpolator.
    jacobian_of_psi_interp: Option<DynInterpolation2d<f64>>,
    /// `J(ψp, θ)` interpolator.
    jacobian_of_psip_interp: Option<DynInterpolation2d<f64>>,
}

/// Creation.
impl NcGeometry {
    /// Constructs an [`NcGeometry`] from an [`NcGeometryBuilder`].
    pub(crate) fn build(builder: NcGeometryBuilder) -> Result<Self, EqError> {
        use crate::extract;
        use crate::extract::netcdf_fields::{NC_BAXIS, NC_RAXIS, NC_RGEO, NC_ZAXIS};
        use crate::extract::netcdf_fields::{NC_JACOBIAN, NC_R, NC_RLAB, NC_THETA, NC_ZLAB};

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let theta_values = extract::array_1d(&file, NC_THETA)?.to_vec();
        let psi = NcFlux::toroidal(&file);
        let psip = NcFlux::poloidal(&file);

        let baxis: f64 = extract::scalar(&file, NC_BAXIS)?;
        let raxis: f64 = extract::scalar(&file, NC_RAXIS)?;
        let zaxis: f64 = extract::scalar(&file, NC_ZAXIS)?;
        let rgeo: f64 = extract::scalar(&file, NC_RGEO)?;
        let r_values = extract::array_1d(&file, NC_R)?.to_vec();
        let rlab_values_fortran_flat = extract::array_2d(&file, NC_RLAB)?
            .flatten_with_order(ColumnMajor)
            .to_vec();
        let zlab_values_fortran_flat = extract::array_2d(&file, NC_ZLAB)?
            .flatten_with_order(ColumnMajor)
            .to_vec();
        let jacobian_values_fortran_flat = extract::array_2d(&file, NC_JACOBIAN)?
            .flatten_with_order(ColumnMajor)
            .to_vec();

        debug_assert!(baxis.is_finite(), "NaN 'baxis' field encountered");
        debug_assert!(raxis.is_finite(), "NaN 'raxis' field encountered");
        debug_assert!(zaxis.is_finite(), "NaN 'zaxis' field encountered");
        debug_assert!(rgeo.is_finite(), "NaN 'rgeo' field encountered");
        debug_assert_all_finite_values(&theta_values);
        debug_assert_all_finite_values(&r_values);
        debug_assert_all_finite_values(&rlab_values_fortran_flat);
        debug_assert_all_finite_values(&zlab_values_fortran_flat);
        debug_assert_all_finite_values(&jacobian_values_fortran_flat);

        // Create interpolators, if possible
        use NcFluxState::Good;
        let psip_of_psi_interp = if (psi.state() == Good) & (psip.state() != NcFluxState::None) {
            Some(make_interp_type(&builder.interp1d_type)?.build(psi.uvalues(), psip.uvalues())?)
        } else {
            None
        };
        let psi_of_psip_interp = if (psip.state() == Good) & (psi.state() != NcFluxState::None) {
            Some(make_interp_type(&builder.interp1d_type)?.build(psip.uvalues(), psi.uvalues())?)
        } else {
            None
        };

        let r_of_psi_interp = if psi.state() == Good {
            Some(make_interp_type(&builder.interp1d_type)?.build(psi.uvalues(), &r_values)?)
        } else {
            None
        };
        let r_of_psip_interp = if psip.state() == Good {
            Some(make_interp_type(&builder.interp1d_type)?.build(psip.uvalues(), &r_values)?)
        } else {
            None
        };

        // Neither the fluxes or `r` is guaranteed to exist.
        // If `r` exists, then it is guaranteed it's in increasing order.
        let psi_of_r_interp = match psi.state() {
            NcFluxState::None => None,
            _ => make_interp_type(&builder.interp1d_type)?
                .build(&r_values, psi.uvalues())
                .ok(),
        };
        let psip_of_r_interp = match psip.state() {
            NcFluxState::None => None,
            _ => make_interp_type(&builder.interp1d_type)?
                .build(&r_values, psip.uvalues())
                .ok(),
        };

        let rlab_of_psi_interp = if psi.state() == Good {
            Some(make_interp2d_type(&builder.interp2d_type)?.build(
                psi.uvalues(),
                &theta_values,
                &rlab_values_fortran_flat,
            )?)
        } else {
            None
        };
        let rlab_of_psip_interp = if psip.state() == Good {
            Some(make_interp2d_type(&builder.interp2d_type)?.build(
                psip.uvalues(),
                &theta_values,
                &rlab_values_fortran_flat,
            )?)
        } else {
            None
        };

        let zlab_of_psi_interp = if psi.state() == Good {
            Some(make_interp2d_type(&builder.interp2d_type)?.build(
                psi.uvalues(),
                &theta_values,
                &zlab_values_fortran_flat,
            )?)
        } else {
            None
        };
        let zlab_of_psip_interp = if psip.state() == Good {
            Some(make_interp2d_type(&builder.interp2d_type)?.build(
                psip.uvalues(),
                &theta_values,
                &zlab_values_fortran_flat,
            )?)
        } else {
            None
        };

        let jacobian_of_psi_interp = if psi.state() == Good {
            Some(make_interp2d_type(&builder.interp2d_type)?.build(
                psi.uvalues(),
                &theta_values,
                &jacobian_values_fortran_flat,
            )?)
        } else {
            None
        };
        let jacobian_of_psip_interp = if psip.state() == Good {
            Some(make_interp2d_type(&builder.interp2d_type)?.build(
                psip.uvalues(),
                &theta_values,
                &jacobian_values_fortran_flat,
            )?)
        } else {
            None
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

impl FluxCommute for NcGeometry {
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.psip_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                self.psip.uvalues(),
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψp(ψ)".into())),
        }
    }

    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.psi_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                self.psi.uvalues(),
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψ(ψp)".into())),
        }
    }
}

impl Geometry for NcGeometry {
    fn r_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.r_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.r_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("r(ψ)".into())),
        }
    }

    fn r_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.r_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.r_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("r(ψp)".into())),
        }
    }

    fn psi_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_r!(r);
        match self.psi_of_r_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                &self.r_values,
                self.psi.uvalues(),
                r,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψ(r)".into())),
        }
    }

    fn psip_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_r!(r);
        match self.psip_of_r_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                &self.r_values,
                self.psip.uvalues(),
                r,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψp(r)".into())),
        }
    }

    fn rlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.rlab_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.theta_values,
                &self.rlab_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("R(ψ, θ)".into())),
        }
    }

    fn rlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.rlab_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.theta_values,
                &self.rlab_values_fortran_flat,
                psip,
                theta,
                psip_acc,
                theta_acc,
                cache,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("R(ψp, θ)".into())),
        }
    }

    fn zlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.zlab_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.theta_values,
                &self.zlab_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("Z(ψ, θ)".into())),
        }
    }

    fn zlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.zlab_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.theta_values,
                &self.zlab_values_fortran_flat,
                psip,
                theta,
                psip_acc,
                theta_acc,
                cache,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("Z(ψp, θ)".into())),
        }
    }

    fn jacobian_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.jacobian_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.theta_values,
                &self.jacobian_values_fortran_flat,
                psi,
                theta,
                psi_acc,
                theta_acc,
                cache,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("J(ψ, θ)".into())),
        }
    }

    fn jacobian_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.jacobian_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.theta_values,
                &self.jacobian_values_fortran_flat,
                psip,
                theta,
                psip_acc,
                theta_acc,
                cache,
            )?)),
            None => Err(EvalError::UndefinedEvaluation("J(ψp, θ)".into())),
        }
    }

    fn rlab_wall(&self) -> Array1<f64> {
        // Get the last row of the C-ordered `rlab_array`
        self.rlab_array().row(self.shape().0 - 1).to_owned()
    }

    fn zlab_wall(&self) -> Array1<f64> {
        // Get the last row of the C-ordered `zlab_array`
        self.zlab_array().row(self.shape().0 - 1).to_owned()
    }
}

/// Getters.
impl NcGeometry {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(2);

    /// Returns the magnetic field strength on the axis `B0` **in \[T\]**.
    #[must_use]
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    /// Returns the horizontal position of the magnetic axis `R0` **in \[m\]**.
    #[must_use]
    pub fn raxis(&self) -> f64 {
        self.raxis
    }

    /// Returns the vertical position of the magnetic axis **in \[m\]**.
    #[must_use]
    pub fn zaxis(&self) -> f64 {
        self.zaxis
    }

    /// Returns the geometrical axis (device major radius) **in \[m\]**.
    #[must_use]
    pub fn rgeo(&self) -> f64 {
        self.rgeo
    }

    /// Returns the `r` coordinate's value at the wall **in \[m\]**.
    #[must_use]
    pub fn rwall(&self) -> f64 {
        match self.r_values.last().copied() {
            Some(rwall) => rwall,
            None => unreachable!("NcGeometry cannot be created if `r_values` dont exist"),
        }
    }

    shape2d_getter_impl!();
    fluxes_wall_value_getter_impl!();
    fluxes_state_getter_impl!();
    fluxes_values_array_getter_impl!();
    vec_to_array1D_getter_impl!(theta_array, theta_values, theta);
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
            .field("baxis [T]", &self.baxis)
            .field("raxis [m]", &self.raxis)
            .field("zaxis [m]", &self.zaxis)
            .field("rgeo [m]", &self.rgeo)
            .field("rwall [m]", &self.rwall())
            .field("shape (ψ/ψp, θ)", &self.shape())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .finish()
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
        assert_eq!(geometry.psi.state(), NcFluxState::Good);
        assert_eq!(geometry.psip.state(), NcFluxState::Bad);
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

        assert!(geometry.psi_array().is_some());
        assert!(geometry.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psi_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        assert!(g.psip_of_psi(p, &mut a1).unwrap().is_finite());
        assert!(g.r_of_psi(p, &mut a1).unwrap().is_finite());
        assert!(g.rlab_of_psi(p, t, &mut a1, &mut a2, &mut c).unwrap().is_finite());
        assert!(g.zlab_of_psi(p, t, &mut a1, &mut a2, &mut c).unwrap().is_finite());
        assert!(g.jacobian_of_psi(p, t, &mut a1, &mut a2, &mut c).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        use EvalError::UndefinedEvaluation as err;
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
        assert_eq!(geometry.psi.state(), NcFluxState::Bad);
        assert_eq!(geometry.psip.state(), NcFluxState::Good);
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

        assert!(geometry.psi_array().is_some());
        assert!(geometry.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psip_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        assert!(g.psi_of_psip(p, &mut a1).unwrap().is_finite());
        assert!(g.r_of_psip(p, &mut a1).unwrap().is_finite());
        assert!(g.rlab_of_psip(p, t, &mut a1, &mut a2, &mut c).unwrap().is_finite());
        assert!(g.zlab_of_psip(p, t, &mut a1, &mut a2, &mut c).unwrap().is_finite());
        assert!(g.jacobian_of_psip(p, t, &mut a1, &mut a2, &mut c).unwrap().is_finite());
    }

    #[test]
    fn bad_psi_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let mut a1 = Accelerator::new();
        let mut a2 = Accelerator::new();
        let mut c = Cache::new();
        let p = 0.01;
        let t = 3.14;
        use EvalError::UndefinedEvaluation as err;
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
