//! Object for conversion from normalized to lab quantities
//!
//! We do not really care about performance here, creating accelerators in every eval call
//! is much simpler.

use std::f64::consts::TAU;
use std::path::PathBuf;

use common::array1D_getter_impl;
use rsl_interpolation::{
    Accelerator, Cache, DynInterpolation, DynInterpolation2d, Interp2dType, InterpType,
    make_interp_type, make_interp2d_type,
};

use ndarray::{Array1, Array2};

use crate::fortran_vec_to_carray2d_impl;
use crate::{Flux, Length, Radians};
use crate::{Geometry, Result};

/// Used to create an [`NcGeometry`].
pub struct NcGeometryBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ1d: String,
    /// 2D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.Interp2dType.html#implementors
    typ2d: String,
}

impl NcGeometryBuilder {
    /// Creates a new [`NcGeometryBuilder`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let builder = NcGeometryBuilder::new(&path, "akima", "bicubic");
    /// ```
    pub fn new(path: &PathBuf, typ1d: &str, typ2d: &str) -> Self {
        Self {
            path: path.clone(),
            typ1d: typ1d.into(),
            typ2d: typ2d.into(),
        }
    }

    /// Creates a new [`NcGeometry`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geometry = NcGeometryBuilder::new(&path, "akima", "bicubic").build()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn build(self) -> Result<NcGeometry> {
        NcGeometry::build(self)
    }
}

// ===============================================================================================

/// Describes the general geometry of the equilibrium.
///
/// Stores fluxes, angles and lab variables' data, and provides interpolation methods between them.
pub struct NcGeometry {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ1d: String,
    /// 2D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.Interp2dType.html#implementors
    typ2d: String,

    /// Magnetic field strength on the axis `B0` **in \[T\]**.
    baxis: f64,
    /// The horizontal position of the magnetic axis `R0` **in \[m\]**.
    raxis: Length,
    /// The vertical position of the magnetic axis **in \[m\]**.
    zaxis: Length,
    /// The geometrical axis (device major radius) **in \[m\]**.
    rgeo: Length,

    /// The boozer toroidal angle `θ` **in \[rads\]**.
    theta_data: Vec<Radians>,
    /// The poloidal flux `ψp` **in Normalized Units**.
    psip_data: Vec<Flux>,
    /// The toroidal flux `ψ` **in Normalized Units**.
    psi_data: Vec<Flux>,
    /// The radial coordinate r **in \[m\]**.
    r_data: Vec<Length>,

    /// R(ψp, θ): The `R` coordinate with respect to boozer coordinates **in \[m\]**, flattened
    /// in F order.
    rlab_data_fortran_flat: Vec<Length>,
    /// Z(ψp, θ): The `Z` coordinate with respect to boozer coordinates **in \[m\]**, flattened
    /// in F order.
    zlab_data_fortran_flat: Vec<Length>,
    /// J(ψp, θ): The VMEC output to Boozer Jacobian in **\[ m/T \]**, flattened in F order.
    jacobian_data_fortran_flat: Vec<f64>,

    /// Interpolator of `ψp(r)` **in \[m\]**.
    psip_of_r_interp: DynInterpolation<f64>,
    /// Interpolator of `r(ψp)` **in \[m\]**.
    r_of_psip_interp: DynInterpolation<f64>,

    /// Interpolator over the R coordinate, as a function of ψp, θ.
    rlab_interp: DynInterpolation2d<f64>,
    /// Interpolator over the Z coordinate, as a function of ψp, θ.
    zlab_interp: DynInterpolation2d<f64>,
    /// Interpolator over the Jacobian, as a function of ψp, θ.
    jacobian_interp: DynInterpolation2d<f64>,
}

/// Creation
impl NcGeometry {
    /// Constructs an [`NcGeometry`] from [`NcGeometryBuilder`].
    pub(crate) fn build(builder: NcGeometryBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;

        let baxis = extract_scalar(&f, BAXIS)?;
        let raxis = extract_scalar(&f, RAXIS)?;
        let zaxis = extract_scalar(&f, ZAXIS)?;
        let rgeo = extract_scalar(&f, RGEO)?;
        let psip_data = extract_1d_array(&f, PSIP_NORM)?
            .as_standard_layout()
            .to_vec();
        let psi_data = extract_1d_array(&f, PSI_NORM)?
            .as_standard_layout()
            .to_vec();
        let r_data = extract_1d_array(&f, R)?.as_standard_layout().to_vec();
        let theta_data = extract_1d_array(&f, THETA)?.as_standard_layout().to_vec();
        let rlab_data = extract_2d_array(&f, RLAB)?;
        let zlab_data = extract_2d_array(&f, ZLAB)?;
        let jacobian_data = extract_2d_array(&f, JACOBIAN)?;

        // Interpolator's `za` input must be in Fortran order.
        let order = ndarray::Order::ColumnMajor;
        let rlab_data_fortran_flat = rlab_data.flatten_with_order(order).to_vec();
        let zlab_data_fortran_flat = zlab_data.flatten_with_order(order).to_vec();
        let jacobian_data_fortran_flat = jacobian_data.flatten_with_order(order).to_vec();

        let r_of_psip_interp = make_interp_type(&builder.typ1d)?.build(&psip_data, &r_data)?;

        let psip_of_r_interp = make_interp_type(&builder.typ1d)?.build(&r_data, &psip_data)?;

        let rlab_interp = make_interp2d_type(&builder.typ2d)?.build(
            &psip_data,
            &theta_data,
            &rlab_data_fortran_flat,
        )?;

        let zlab_interp = make_interp2d_type(&builder.typ2d)?.build(
            &psip_data,
            &theta_data,
            &zlab_data_fortran_flat,
        )?;

        let jacobian_interp = make_interp2d_type(&builder.typ2d)?.build(
            &psip_data,
            &theta_data,
            &jacobian_data_fortran_flat,
        )?;

        Ok(Self {
            path,
            typ1d: builder.typ1d,
            typ2d: builder.typ2d,
            baxis,
            raxis,
            zaxis,
            rgeo,
            psip_data,
            psi_data,
            theta_data,
            r_data,
            rlab_data_fortran_flat,
            zlab_data_fortran_flat,
            jacobian_data_fortran_flat,
            psip_of_r_interp,
            r_of_psip_interp,
            rlab_interp,
            zlab_interp,
            jacobian_interp,
        })
    }
}

/// Interpolation
impl Geometry for NcGeometry {
    fn r(&self, psip: Flux) -> Result<Length> {
        let mut acc = Accelerator::new();
        Ok(self
            .r_of_psip_interp
            .eval(&self.psip_data, &self.r_data, psip, &mut acc)?)
    }

    fn psip(&self, r: Length) -> Result<Flux> {
        let mut acc = Accelerator::new();
        Ok(self
            .psip_of_r_interp
            .eval(&self.r_data, &self.psip_data, r, &mut acc)?)
    }

    fn rlab(&self, psip: Flux, theta: Radians) -> Result<f64> {
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        Ok(self.rlab_interp.eval(
            &self.psip_data,
            &self.theta_data,
            &self.rlab_data_fortran_flat,
            psip,
            theta.rem_euclid(TAU),
            &mut xacc,
            &mut yacc,
            &mut cache,
        )?)
    }

    fn zlab(&self, psip: Flux, theta: Radians) -> Result<f64> {
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        Ok(self.zlab_interp.eval(
            &self.psip_data,
            &self.theta_data,
            &self.zlab_data_fortran_flat,
            psip,
            theta.rem_euclid(TAU),
            &mut xacc,
            &mut yacc,
            &mut cache,
        )?)
    }

    fn jacobian(&self, psip: Flux, theta: Radians) -> Result<f64> {
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        Ok(self.jacobian_interp.eval(
            &self.psip_data,
            &self.theta_data,
            &self.jacobian_data_fortran_flat,
            psip,
            theta.rem_euclid(TAU),
            &mut xacc,
            &mut yacc,
            &mut cache,
        )?)
    }
}

/// Getters
impl NcGeometry {
    /// Returns the netCDF file's path.
    pub fn path(&self) -> PathBuf {
        self.path.clone()
    }

    /// Returns the 1D interpolation type.
    pub fn typ1d(&self) -> String {
        self.typ1d.clone()
    }

    /// Returns the 2D interpolation type.
    pub fn typ2d(&self) -> String {
        self.typ2d.clone()
    }

    /// Returns the shape of the `b` array.
    pub fn shape(&self) -> (usize, usize) {
        (self.psip_data.len(), self.theta_data.len())
    }

    /// Returns the magnetic field strength on the axis `B0` **in \[T\]**.
    pub fn baxis(&self) -> f64 {
        self.baxis
    }

    /// Retruns the horizontal position of the magnetic axis `R0` **in \[m\]**.
    pub fn raxis(&self) -> f64 {
        self.raxis
    }

    /// Retruns the vertical position of the magnetic axis **in \[m\]**.
    pub fn zaxis(&self) -> f64 {
        self.zaxis
    }

    /// Returns the geometrical axis (device major radius) **in \[m\]**.
    pub fn rgeo(&self) -> f64 {
        self.rgeo
    }

    /// Returns the tokamak's minor radius `r_wall` **in \[m\]**.
    pub fn r_wall(&self) -> f64 {
        // `r_data` is always non-empty, otherwise `NcGeometry` cannot be constructed
        self.r_data.last().copied().expect("array non-empty")
    }

    /// Returns the poloidal flux's value at the wall `ψp_wall` **in Normalized Units**.
    pub fn psip_wall(&self) -> f64 {
        // `psip_data` is always non-empty, otherwise `NcGeometry` cannot be constructed
        self.psip_data.last().copied().expect("array non-empty")
    }

    /// Returns the toroidal flux's value at the wall `ψ_wall` **in Normalized Units**.
    pub fn psi_wall(&self) -> f64 {
        // `psi_data` is always non-empty, otherwise `NcGeometry` cannot be constructed
        self.psi_data.last().copied().expect("array non-empty")
    }

    array1D_getter_impl!(theta_data, theta_data, Radians);
    array1D_getter_impl!(psip_data, psip_data, Flux);
    array1D_getter_impl!(psi_data, psi_data, Flux);
    array1D_getter_impl!(r_data, r_data, Length);

    fortran_vec_to_carray2d_impl!(rlab_data, rlab_data_fortran_flat, R);
    fortran_vec_to_carray2d_impl!(zlab_data, zlab_data_fortran_flat, Z);
    fortran_vec_to_carray2d_impl!(jacobian_data, jacobian_data_fortran_flat, J);
}

impl std::fmt::Debug for NcGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcGeometry")
            .field("path", &self.path())
            .field("typ 1D", &self.typ1d())
            .field("typ 2D", &self.typ2d())
            .field("Baxis [T]", &format!("{:.7}", self.baxis()))
            .field("Raxis [m]", &format!("{:.7}", self.raxis()))
            .field("Zaxis [m]", &format!("{:.7}", self.zaxis()))
            .field("Rgeo [m]", &format!("{:.7}", self.rgeo()))
            .field("ψp_wall", &format!("{:.7}", self.psip_wall()))
            .field("ψ_wall", &format!("{:.7}", self.psi_wall()))
            .field("r_wall", &format!("{:.7}", self.r_wall()))
            .field("shape", &self.shape())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::extract::STUB_NETCDF_PATH;

    fn create_nc_geometry() -> NcGeometry {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        let typ1d = "steffen";
        let typ2d = "bicubic";
        NcGeometryBuilder::new(&path, typ1d, typ2d).build().unwrap()
    }

    #[test]
    fn test_geometry_creation() {
        let g = create_nc_geometry();
        let _ = format!("{g:?}");
    }

    #[test]
    fn test_getters() {
        let g = create_nc_geometry();
        g.path();
        g.typ1d();
        g.typ2d();
        g.baxis();
        g.raxis();
        g.zaxis();
        g.rgeo();
        g.psip_wall();
        g.psi_wall();
        g.r_wall();
        g.shape();

        assert_eq!(g.psip_data().ndim(), 1);
        assert_eq!(g.psi_data().ndim(), 1);
        assert_eq!(g.r_data().ndim(), 1);
        assert_eq!(g.theta_data().ndim(), 1);
        assert_eq!(g.rlab_data().ndim(), 2);
        assert_eq!(g.zlab_data().ndim(), 2);
        assert_eq!(g.jacobian_data().ndim(), 2);
    }

    #[test]
    fn test_spline_evaluation() {
        let g = create_nc_geometry();

        let r = 0.1;
        let psip = 0.015;
        let theta = 0.0;
        g.r(psip).unwrap();
        g.psip(r).unwrap();
        g.rlab(psip, theta).unwrap();
        g.zlab(psip, theta).unwrap();
        g.jacobian(psip, theta).unwrap();
    }
}
