//! Representation of a machine's general geometry.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_flux, debug_assert_non_negative_r,
    fortran_vec_to_carrayview2d_impl, nc_fluxes_values_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl,
};
use core::f64::consts::PI;
use ndarray::{Array1, ArrayView2, Order::ColumnMajor};
use rsl_interpolation::{Accelerator, Accelerator2d};
use std::hint::cold_path;
use std::path::{Path, PathBuf};

use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{EvalError, MachineError};
use crate::{Geometry, MachineObject, MachineType};
use crate::{Interpolation1dType, Interpolation2dType};
use crate::{MagneticFlux, MagneticFlux::*};
use dexter_common::{DynInterpolator, DynInterpolator2d, make_interp, make_interp2d};

// ===============================================================================================

/// Analytical Large Aspect Ratio Geometry of a circular device.
///
/// The definitions are not very strict at the moment.
///
/// # Note
///
/// No ψ bounds checks are performed in evaluations.
#[non_exhaustive]
#[derive(Debug)]
pub struct LarGeometry {
    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,
    /// The horizontal position of the magnetic axis `R0` in [m].
    raxis: f64,
    /// The `r` coordinate's value at the last closed flux surface `rlast` in [m].
    ///
    /// In LAR configuration, `rlast` coincides with the device's minor radius.
    rlast: f64,
    /// The value of the last closed toroidal flux surface `ψt_last` in Normalized units.
    ///
    /// In LAR configuration, the wall coincides with the last closed flux surface.
    psit_last: f64,
}

impl LarGeometry {
    /// Creates a new `LarGeometry`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// let geometry = LarGeometry::new(2.0, 1.75, 0.5);
    /// ```
    #[must_use]
    pub fn new(baxis: f64, raxis: f64, rlast: f64) -> Self {
        let psit_last_si = baxis * rlast.powi(2) / 2.0;
        let psit_last = psit_last_si / (baxis * raxis.powi(2));
        Self {
            baxis,
            raxis,
            rlast,
            psit_last,
        }
    }

    /// Returns the last closed toroidal flux surface.
    #[must_use]
    pub fn psit_last(&self) -> MagneticFlux {
        Toroidal(self.psit_last)
    }
}

impl MachineObject for LarGeometry {
    fn machine_type(&self) -> MachineType {
        MachineType::Analytical
    }

    fn psit_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Bad
    }
}

impl Geometry for LarGeometry {
    fn baxis(&self) -> f64 {
        self.baxis
    }

    fn raxis(&self) -> f64 {
        self.raxis
    }

    fn zaxis(&self) -> f64 {
        0.0
    }

    fn rlast(&self) -> f64 {
        self.rlast
    }

    fn rgeo(&self) -> f64 {
        self.raxis
    }

    fn eval_r(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        match flux {
            // Must normalize `r`
            Toroidal(psit) => Ok(debug_assert_is_finite!(self.raxis * (2.0 * psit).sqrt())),
            Poloidal(_) => Err(EvalError::UndefinedEvaluation(
                "r(ψp) (defined through q)".into(),
            )),
        }
    }

    fn eval_psit_of_r(&self, r: f64, _: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        debug_assert_non_negative_r!(r);
        // Must normalize `r` first
        Ok(Toroidal(debug_assert_is_finite!(
            (r / self.raxis).powi(2) / 2.0
        )))
    }

    fn eval_psip_of_r(&self, _: f64, _: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        Err(EvalError::UndefinedEvaluation(
            "ψp(r) (defined through q)".into(),
        ))
    }

    fn eval_rlab(
        &self,
        flux: MagneticFlux,
        theta: f64,
        _: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        match flux {
            Toroidal(psit) => Ok(debug_assert_is_finite!(
                self.raxis + self.raxis * (2.0 * psit).sqrt() * theta.cos()
            )),
            Poloidal(_) => Err(EvalError::UndefinedEvaluation(
                "R(ψp, θ) (defined through q)".into(),
            )),
        }
    }

    fn eval_zlab(
        &self,
        flux: MagneticFlux,
        theta: f64,
        _: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        match flux {
            Toroidal(psit) => Ok(debug_assert_is_finite!(
                self.raxis * (2.0 * psit).sqrt() * theta.sin()
            )),
            Poloidal(_) => Err(EvalError::UndefinedEvaluation(
                "Z(ψp, θ) (defined through q)".into(),
            )),
        }
    }

    fn eval_jacobian(
        &self,
        flux: MagneticFlux,
        _: f64,
        _: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let msg = format!("J({}, θ) (defined through q, g, I and B)", flux.symbol());
        Err(EvalError::UndefinedEvaluation(msg))
    }

    fn rlab_last(&self) -> Array1<f64> {
        let arr = Array1::linspace(0.0, 2.0 * PI, 1000);
        let acc = &mut Accelerator2d::new();
        arr.mapv(
            |theta| match self.eval_rlab(Toroidal(self.psit_last), theta, acc) {
                Ok(rlab_lcfs_value) => rlab_lcfs_value,
                Err(_) => unreachable!("Expression is analytical, cannot fail"),
            },
        )
    }

    fn zlab_last(&self) -> Array1<f64> {
        let arr = Array1::linspace(0.0, 2.0 * PI, 1000);
        let acc = &mut Accelerator2d::new();
        arr.mapv(
            |theta| match self.eval_zlab(Toroidal(self.psit_last), theta, acc) {
                Ok(zlab_lcfs_value) => zlab_lcfs_value,
                Err(_) => unreachable!("Expression is analytical, cannot fail"),
            },
        )
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
    /// The 1D interpolation type.
    interp1d_type: Interpolation1dType,
    /// The 2D interpolation type.
    interp2d_type: Interpolation2dType,
}

impl NcGeometryBuilder {
    /// Creates a new `NcGeometryBuilder` from a netCDF file at `path`, with 1D interpolation
    /// type `interp1d_type` and 2D interpolation type `interp2d_type`.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let interp1d_type = Interpolation1dType::Akima;
    /// let interp2d_type = Interpolation2dType::Bicubic;
    /// let builder = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type);
    /// ```
    #[must_use]
    pub fn new(
        path: &Path,
        interp1d_type: Interpolation1dType,
        interp2d_type: Interpolation2dType,
    ) -> Self {
        Self {
            path: path.to_path_buf(),
            interp1d_type,
            interp2d_type,
        }
    }

    /// Creates a new [`NcGeometry`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let interp1d_type = Interpolation1dType::Akima;
    /// let interp2d_type = Interpolation2dType::Bicubic;
    /// let builder = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type);
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`MachineError`] if it fails to build the [`NcGeometry`].
    pub fn build(self) -> Result<NcGeometry, MachineError> {
        NcGeometry::build(self)
    }
}

// ===============================================================================================

/// Geometry of a realistic configuration.
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

    /// The 1D interpolation type.
    interp1d_type: Interpolation1dType,
    /// The 2D interpolation type.
    interp2d_type: Interpolation2dType,

    /// Magnetic field strength on the axis `B0` in [T].
    baxis: f64,
    /// The horizontal position of the magnetic axis `R0` in [m].
    raxis: f64,
    /// The vertical position of the magnetic axis in [m].
    zaxis: f64,
    /// The geometrical axis (device major radius) in [m].
    rgeo: f64,

    /// The boozer toroidal angle `θ` in [rads].
    theta_values: Box<[f64]>,
    /// The toroidal flux coordinate.
    psit: NcFlux,
    /// The poloidal flux coordinate.
    psip: NcFlux,

    /// The radial coordinate r in [m].
    r_values: Box<[f64]>,
    /// The `r(ψt)` interpolator.
    r_of_psit_interp: Option<DynInterpolator>,
    /// The `r(ψp)` interpolator.
    r_of_psip_interp: Option<DynInterpolator>,
    /// The `ψt(r)` interpolator.
    psit_of_r_interp: Option<DynInterpolator>,
    /// The `ψp(r)` interpolator.
    psip_of_r_interp: Option<DynInterpolator>,

    /// The `R` coordinate in [m], flattened in F order.
    rlab_values_fortran_flat: Box<[f64]>,
    /// `R(ψt, θ)` interpolator.
    rlab_of_psit_interp: Option<DynInterpolator2d>,
    /// `R(ψp, θ)` interpolator.
    rlab_of_psip_interp: Option<DynInterpolator2d>,

    /// The `Z` coordinate in [m], flattened in F order.
    zlab_values_fortran_flat: Box<[f64]>,
    /// `Z(ψt, θ)` interpolator.
    zlab_of_psit_interp: Option<DynInterpolator2d>,
    /// `Z(ψp, θ)` interpolator.
    zlab_of_psip_interp: Option<DynInterpolator2d>,

    /// The VMEC output to Boozer Jacobian in [m/T], flattened in F order.
    jacobian_values_fortran_flat: Option<Box<[f64]>>,
    /// `J(ψt, θ)` interpolator.
    jacobian_of_psit_interp: Option<DynInterpolator2d>,
    /// `J(ψp, θ)` interpolator.
    jacobian_of_psip_interp: Option<DynInterpolator2d>,
}

/// Creation.
impl NcGeometry {
    /// Constructs an `NcGeometry` from an [`NcGeometryBuilder`].
    pub(crate) fn build(builder: NcGeometryBuilder) -> Result<Self, MachineError> {
        use crate::extract;
        use crate::extract::netcdf_fields::{NC_BAXIS, NC_RAXIS, NC_RGEO, NC_ZAXIS};
        use crate::extract::netcdf_fields::{NC_JACOBIAN, NC_R, NC_RLAB, NC_THETA, NC_ZLAB};

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let theta_values = extract::array_1d(&file, NC_THETA)?.to_vec();
        let psit = NcFlux::toroidal(&file);
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
        let jacobian_values_fortran_flat = extract::array_2d(&file, NC_JACOBIAN)
            .map(|array| array.flatten_with_order(ColumnMajor).to_vec())
            .ok();

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let r_of_psit_interp = match psit.state() {
            Good => make_interp(builder.interp1d_type, psit.uvalues(), &r_values).ok(),
            _ => None,
        };
        let r_of_psip_interp = match psip.state() {
            Good => make_interp(builder.interp1d_type, psip.uvalues(), &r_values).ok(),
            _ => None,
        };

        // Neither the fluxes or `r` is guaranteed to exist.
        // If `r` exists, then it is guaranteed it's in increasing order.
        let psit_of_r_interp = match psit.state() {
            FluxCoordinateState::NoValues => None,
            _ => make_interp(builder.interp1d_type, &r_values, psit.uvalues()).ok(),
        };
        let psip_of_r_interp = match psip.state() {
            FluxCoordinateState::NoValues => None,
            _ => make_interp(builder.interp1d_type, &r_values, psip.uvalues()).ok(),
        };

        let rlab_of_psit_interp = match psit.state() {
            Good => make_interp2d(
                builder.interp2d_type,
                psit.uvalues(),
                &theta_values,
                &rlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };
        let rlab_of_psip_interp = match psip.state() {
            Good => make_interp2d(
                builder.interp2d_type,
                psip.uvalues(),
                &theta_values,
                &rlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };

        let zlab_of_psit_interp = match psit.state() {
            Good => make_interp2d(
                builder.interp2d_type,
                psit.uvalues(),
                &theta_values,
                &zlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };
        let zlab_of_psip_interp = match psip.state() {
            Good => make_interp2d(
                builder.interp2d_type,
                psip.uvalues(),
                &theta_values,
                &zlab_values_fortran_flat,
            )
            .ok(),
            _ => None,
        };

        let jacobian_of_psit_interp = match psit.state() {
            Good => jacobian_values_fortran_flat.as_ref().and_then(|values| {
                make_interp2d(builder.interp2d_type, psit.uvalues(), &theta_values, values).ok()
            }),
            _ => None,
        };
        let jacobian_of_psip_interp = match psip.state() {
            Good => jacobian_values_fortran_flat.as_ref().and_then(|values| {
                make_interp2d(builder.interp2d_type, psip.uvalues(), &theta_values, values).ok()
            }),
            _ => None,
        };

        Ok(Self {
            netcdf_version,
            path,
            interp1d_type: builder.interp1d_type,
            interp2d_type: builder.interp2d_type,
            theta_values: theta_values.into_boxed_slice(),
            psit,
            psip,
            baxis,
            raxis,
            zaxis,
            rgeo,
            r_values: r_values.into_boxed_slice(),
            r_of_psit_interp,
            r_of_psip_interp,
            psit_of_r_interp,
            psip_of_r_interp,
            rlab_values_fortran_flat: rlab_values_fortran_flat.into_boxed_slice(),
            rlab_of_psit_interp,
            rlab_of_psip_interp,
            zlab_values_fortran_flat: zlab_values_fortran_flat.into_boxed_slice(),
            zlab_of_psit_interp,
            zlab_of_psip_interp,
            jacobian_values_fortran_flat: jacobian_values_fortran_flat.map(Vec::into_boxed_slice),
            jacobian_of_psit_interp,
            jacobian_of_psip_interp,
        })
    }
}

impl MachineObject for NcGeometry {
    fn machine_type(&self) -> MachineType {
        MachineType::Numerical
    }

    fn psit_state(&self) -> FluxCoordinateState {
        self.psit.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip.state()
    }
}

impl Geometry for NcGeometry {
    fn baxis(&self) -> f64 {
        self.baxis
    }

    fn raxis(&self) -> f64 {
        self.raxis
    }

    fn zaxis(&self) -> f64 {
        self.zaxis
    }

    fn rgeo(&self) -> f64 {
        self.rgeo
    }

    fn rlast(&self) -> f64 {
        match self.r_values.last().copied() {
            Some(rlast) => rlast,
            None => unreachable!("NcGeometry cannot be created if `r_values` dont exist"),
        }
    }

    fn eval_r(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let interp_opt = match flux {
            Toroidal(_) => self.r_of_psit_interp.as_ref(),
            Poloidal(_) => self.r_of_psip_interp.as_ref(),
        };
        let Some(interp) = interp_opt else {
            cold_path();
            let msg = format!("r({})", flux.symbol());
            return Err(EvalError::UndefinedEvaluation(msg));
        };
        // This cannot panic. If `interp` is `Some` then the corresponding values exist.
        let (val, xa) = match flux {
            Toroidal(val) => (val, self.psit.uvalues()),
            Poloidal(val) => (val, self.psip.uvalues()),
        };
        let ya = &self.theta_values;
        Ok(debug_assert_is_finite!(interp.eval(xa, ya, val, acc)?))
    }

    fn eval_psit_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        debug_assert_non_negative_r!(r);
        match self.psit_of_r_interp.as_ref() {
            Some(interp) => Ok(Toroidal(debug_assert_is_finite!(interp.eval(
                &self.r_values,
                self.psit.uvalues(),
                r,
                acc
            )?))),
            None => Err(EvalError::UndefinedEvaluation("ψ(r)".into())),
        }
    }

    fn eval_psip_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        debug_assert_non_negative_r!(r);
        match self.psip_of_r_interp.as_ref() {
            Some(interp) => Ok(Poloidal(debug_assert_is_finite!(interp.eval(
                &self.r_values,
                self.psip.uvalues(),
                r,
                acc
            )?))),
            None => Err(EvalError::UndefinedEvaluation("ψp(r)".into())),
        }
    }

    fn eval_rlab(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let interp_opt = match flux {
            Toroidal(_) => self.rlab_of_psit_interp.as_ref(),
            Poloidal(_) => self.rlab_of_psip_interp.as_ref(),
        };
        let Some(interp) = interp_opt else {
            cold_path();
            let msg = format!("R({}, θ)", flux.symbol());
            return Err(EvalError::UndefinedEvaluation(msg));
        };
        // This cannot panic. If `interp` is `Some` then the corresponding values exist.
        let (val, xa) = match flux {
            Toroidal(val) => (val, self.psit.uvalues()),
            Poloidal(val) => (val, self.psip.uvalues()),
        };
        let ya = &self.theta_values;
        let za = &self.rlab_values_fortran_flat;
        Ok(debug_assert_is_finite!(
            interp.eval(xa, ya, za, val, theta, acc)?
        ))
    }

    fn eval_zlab(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let interp_opt = match flux {
            Toroidal(_) => self.zlab_of_psit_interp.as_ref(),
            Poloidal(_) => self.zlab_of_psip_interp.as_ref(),
        };
        let Some(interp) = interp_opt else {
            cold_path();
            let msg = format!("Z({}, θ)", flux.symbol());
            return Err(EvalError::UndefinedEvaluation(msg));
        };
        // This cannot panic. If `interp` is `Some` then the corresponding values exist.
        let (val, xa) = match flux {
            Toroidal(val) => (val, self.psit.uvalues()),
            Poloidal(val) => (val, self.psip.uvalues()),
        };
        let ya = &self.theta_values;
        let za = &self.zlab_values_fortran_flat;
        Ok(debug_assert_is_finite!(
            interp.eval(xa, ya, za, val, theta, acc)?
        ))
    }

    #[expect(clippy::manual_let_else, reason = "clarity")]
    fn eval_jacobian(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let interp_opt = match flux {
            Toroidal(_) => self.jacobian_of_psit_interp.as_ref(),
            Poloidal(_) => self.jacobian_of_psip_interp.as_ref(),
        };
        let Some(interp) = interp_opt else {
            cold_path();
            let msg = format!("J({}, θ)", flux.symbol());
            return Err(EvalError::UndefinedEvaluation(msg));
        };
        // This cannot panic. If `interp` is `Some` then the corresponding values exist.
        let (val, xa) = match flux {
            Toroidal(val) => (val, self.psit.uvalues()),
            Poloidal(val) => (val, self.psip.uvalues()),
        };
        let ya = &self.theta_values;
        let za = match self.jacobian_values_fortran_flat.as_ref() {
            Some(values) => values,
            None => unreachable!("interp is Some"),
        };
        Ok(debug_assert_is_finite!(
            interp.eval(xa, ya, za, val, theta, acc)?
        ))
    }

    fn rlab_last(&self) -> Array1<f64> {
        // Get the last row of the C-ordered `rlab_array`
        self.rlab_array().row(self.shape().0 - 1).to_owned()
    }

    fn zlab_last(&self) -> Array1<f64> {
        // Get the last row of the C-ordered `zlab_array`
        self.zlab_array().row(self.shape().0 - 1).to_owned()
    }
}

/// Getters.
impl NcGeometry {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();

    /// Returns the 1D interpolation type.
    #[must_use]
    pub fn interp1d_type(&self) -> Interpolation1dType {
        self.interp1d_type
    }

    /// Returns the 2D interpolation type.
    #[must_use]
    pub fn interp2d_type(&self) -> Interpolation2dType {
        self.interp2d_type
    }

    /// Returns the (ψ, θ) shape of the 2D arrays.
    #[must_use]
    pub fn shape(&self) -> (usize, usize) {
        let psit_len = match self.psit_state() {
            FluxCoordinateState::NoValues => 0,
            _ => self.psit.uvalues().len(),
        };
        let psip_len = match self.psip_state() {
            FluxCoordinateState::NoValues => 0,
            _ => self.psip.uvalues().len(),
        };
        let xlen = psit_len.max(psip_len);
        (xlen, self.theta_values.len())
    }

    /// Returns the last closed toroidal flux surface.
    #[must_use]
    pub fn psit_last(&self) -> Option<MagneticFlux> {
        self.psit.last_value().map(Toroidal)
    }

    /// Returns the last closed poloidal flux surface.
    #[must_use]
    pub fn psip_last(&self) -> Option<MagneticFlux> {
        self.psip.last_value().map(Poloidal)
    }

    /// Returns the Jacobian values as a 2D array view.
    #[must_use]
    pub fn jacobian_array(&self) -> Option<ArrayView2<'_, f64>> {
        self.jacobian_values_fortran_flat.as_ref().map(|values| {
            let actual_shape = self.shape();
            let shape = (actual_shape.1, actual_shape.0);
            match ArrayView2::from_shape(shape, values) {
                Ok(view) => view.reversed_axes(),
                Err(_) => unreachable!("Shape is correct by definition"),
            }
        })
    }

    /// Returns a reference to the `θ` values.
    #[must_use]
    pub fn theta_values(&self) -> &[f64] {
        &self.theta_values
    }

    /// Returns a reference to the `r` values.
    #[must_use]
    pub fn r_values(&self) -> &[f64] {
        &self.r_values
    }

    nc_fluxes_values_getter_impl!();
    fortran_vec_to_carrayview2d_impl!(rlab_array, rlab_values_fortran_flat, R);
    fortran_vec_to_carrayview2d_impl!(zlab_array, zlab_values_fortran_flat, Z);
}

impl std::fmt::Debug for NcGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcGeometry")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("1D interpolation type", &self.interp1d_type())
            .field("2D interpolation type", &self.interp2d_type())
            .field("baxis [T]", &self.baxis)
            .field("raxis [m]", &self.raxis)
            .field("zaxis [m]", &self.zaxis)
            .field("rgeo [m]", &self.rgeo)
            .field("rlast [m]", &self.rlast())
            .field("shape (ψ, θ)", &self.shape())
            .field("ψt", &self.psit)
            .field("ψp", &self.psip)
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_geometry(path_str: &str) -> NcGeometry {
        let path = PathBuf::from(&path_str);
        let builder = NcGeometryBuilder::new(
            &path,
            Interpolation1dType::Steffen,
            Interpolation2dType::Bicubic,
        );
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
        assert_eq!(geometry.psit_state(), FluxCoordinateState::Good);
        assert_eq!(geometry.psip_state(), FluxCoordinateState::Bad);

        assert_eq!(geometry.psit.state(), FluxCoordinateState::Good);
        assert_eq!(geometry.psip.state(), FluxCoordinateState::Bad);
        assert!(geometry.r_of_psit_interp.is_some());
        assert!(geometry.rlab_of_psit_interp.is_some());
        assert!(geometry.zlab_of_psit_interp.is_some());
        assert!(geometry.jacobian_of_psit_interp.is_some());

        assert!(geometry.r_of_psip_interp.is_none());
        assert!(geometry.rlab_of_psip_interp.is_none());
        assert!(geometry.zlab_of_psip_interp.is_none());
        assert!(geometry.jacobian_of_psip_interp.is_none());

        assert!(geometry.psit_of_r_interp.is_some());
        assert!(geometry.psip_of_r_interp.is_some());

        assert!(geometry.psit_values().is_some());
        assert!(geometry.psip_values().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psit_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let th = 3.14;
        let flux = Toroidal(0.01);
        assert!(g.eval_r(flux, acc1).unwrap().is_finite());
        assert!(g.eval_rlab(flux, th, acc2).unwrap().is_finite());
        assert!(g.eval_zlab(flux, th, acc2).unwrap().is_finite());
        assert!(g.eval_jacobian(flux, th, acc2).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let g = create_nc_geometry(TOROIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let th = 3.14;
        let flux = Poloidal(0.01);
        use EvalError::UndefinedEvaluation as err;
        matches!(g.eval_r(flux, acc1), Err(err(..)));
        matches!(g.eval_rlab(flux, th, acc2), Err(err(..)));
        matches!(g.eval_zlab(flux, th, acc2), Err(err(..)));
        matches!(g.eval_jacobian(flux, th, acc2), Err(err(..)));
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
        assert_eq!(geometry.psit_state(), FluxCoordinateState::Bad);
        assert_eq!(geometry.psip_state(), FluxCoordinateState::Good);

        assert_eq!(geometry.psit.state(), FluxCoordinateState::Bad);
        assert_eq!(geometry.psip.state(), FluxCoordinateState::Good);
        assert!(geometry.r_of_psit_interp.is_none());
        assert!(geometry.rlab_of_psit_interp.is_none());
        assert!(geometry.zlab_of_psit_interp.is_none());
        assert!(geometry.jacobian_of_psit_interp.is_none());

        assert!(geometry.r_of_psip_interp.is_some());
        assert!(geometry.rlab_of_psip_interp.is_some());
        assert!(geometry.zlab_of_psip_interp.is_some());
        assert!(geometry.jacobian_of_psip_interp.is_some());

        assert!(geometry.psit_of_r_interp.is_some());
        assert!(geometry.psip_of_r_interp.is_some());

        assert!(geometry.psit_values().is_some());
        assert!(geometry.psip_values().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psip_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let th = 3.14;
        let flux = Poloidal(0.01);
        assert!(g.eval_r(flux, acc1).unwrap().is_finite());
        assert!(g.eval_rlab(flux, th, acc2).unwrap().is_finite());
        assert!(g.eval_zlab(flux, th, acc2).unwrap().is_finite());
        assert!(g.eval_jacobian(flux, th, acc2).unwrap().is_finite());
    }

    #[test]
    fn bad_psit_evals() {
        let g = create_nc_geometry(POLOIDAL_TEST_NETCDF_PATH);
        let acc1 = &mut Accelerator::new();
        let acc2 = &mut Accelerator2d::new();
        let th = 3.14;
        let flux = Toroidal(0.01);
        use EvalError::UndefinedEvaluation as err;
        matches!(g.eval_r(flux, acc1), Err(err(..)));
        matches!(g.eval_rlab(flux, th, acc2), Err(err(..)));
        matches!(g.eval_zlab(flux, th, acc2), Err(err(..)));
        matches!(g.eval_jacobian(flux, th, acc2), Err(err(..)));
    }
}
