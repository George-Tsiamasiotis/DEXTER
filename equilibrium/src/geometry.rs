//! Object for conversion from normalized to lab quantities
//!
//! Use the lower level interpolators here, to avoid storing the same data many times.
//!
//! We do not really care about performance here, creating accelerators in every eval call
//! is much simpler.

use std::f64::consts::TAU;
use std::path::PathBuf;

use ndarray::{Array1, Array2};
use rsl_interpolation::{
    Accelerator, Cache, DynInterpolation, DynInterpolation2d, Interp2dType, InterpType,
};
use rsl_interpolation::{make_interp_type, make_interp2d_type};
use utils::array1D_getter_impl;

use crate::Result;
use crate::{Flux, Length, Radians};

/// Stores fluxes, angles and lab variables' data, and provides interpolation methods between them.
pub struct Geometry {
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

    /// R(ψp, θ): The `R` coordinate with respect to boozer coordinates **in \[m\]**.
    rlab_flat_data: Vec<Length>,
    /// Z(ψp, θ): The `Z` coordinate with respect to boozer coordinates **in \[m\]**.
    zlab_flat_data: Vec<Length>,
    /// J(ψp, θ): The VMEC output to Boozer Jacobian in **\[ m/T \]**.
    jacobian_flat_data: Vec<f64>,

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

// Creation
impl Geometry {
    /// Constructs a [`Geometry`] from a netCDF file, with `typ1d` 1D interpolation type, and
    /// `typ2d` 2D interpolation type.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ1d: &str, typ2d: &str) -> Result<Self> {
        use crate::extract::*;
        use config::netcdf_fields::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(path)?;
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
        let rlab_data = extract_2d_array(&f, RLAB)?.to_owned();
        let zlab_data = extract_2d_array(&f, ZLAB)?.to_owned();
        let jacobian_data = extract_2d_array(&f, JACOBIAN)?.to_owned();

        // Interpolator's `za` input must be in Fortran order.
        let order = ndarray::Order::ColumnMajor;
        let rlab_flat_data = rlab_data.flatten_with_order(order).to_vec();
        let zlab_flat_data = zlab_data.flatten_with_order(order).to_vec();
        let jacobian_flat_data = jacobian_data.flatten_with_order(order).to_vec();

        let r_of_psip_interp = make_interp_type(typ1d)?.build(&psip_data, &r_data)?;

        let psip_of_r_interp = make_interp_type(typ1d)?.build(&r_data, &psip_data)?;

        let rlab_interp =
            make_interp2d_type(typ2d)?.build(&psip_data, &theta_data, &rlab_flat_data)?;

        let zlab_interp =
            make_interp2d_type(typ2d)?.build(&psip_data, &theta_data, &zlab_flat_data)?;

        let jacobian_interp =
            make_interp2d_type(typ2d)?.build(&psip_data, &theta_data, &jacobian_flat_data)?;

        Ok(Self {
            path,
            typ1d: typ1d.into(),
            typ2d: typ2d.into(),
            baxis,
            raxis,
            zaxis,
            rgeo,
            psip_data,
            psi_data,
            theta_data,
            r_data,
            rlab_flat_data,
            zlab_flat_data,
            jacobian_flat_data,
            psip_of_r_interp,
            r_of_psip_interp,
            rlab_interp,
            zlab_interp,
            jacobian_interp,
        })
    }
}

/// Interpolation
impl Geometry {
    /// Calculates the radial coordinate `r(ψp)` **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let r = geom.r(0.015)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn r(&self, psip: Flux) -> Result<Length> {
        let mut acc = Accelerator::new();
        Ok(self
            .r_of_psip_interp
            .eval(&self.psip_data, &self.r_data, psip, &mut acc)?)
    }

    /// Calculates the poloidal flux `ψp(r)`, where r is **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let psip = geom.psip(0.45)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn psip(&self, r: Length) -> Result<Flux> {
        let mut acc = Accelerator::new();
        Ok(self
            .psip_of_r_interp
            .eval(&self.r_data, &self.psip_data, r, &mut acc)?)
    }

    /// Calculates `R(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let rlab = geom.rlab(0.015, 2.0*PI)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn rlab(&self, psip: Flux, theta: Radians) -> Result<f64> {
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        Ok(self.rlab_interp.eval(
            &self.psip_data,
            &self.theta_data,
            &self.rlab_flat_data,
            psip,
            theta.rem_euclid(TAU),
            &mut xacc,
            &mut yacc,
            &mut cache,
        )?)
    }

    /// Calculates `Z(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let zlab = geom.zlab(0.015, 2.0*PI)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn zlab(&self, psip: Flux, theta: Radians) -> Result<f64> {
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        Ok(self.zlab_interp.eval(
            &self.psip_data,
            &self.theta_data,
            &self.zlab_flat_data,
            psip,
            theta.rem_euclid(TAU),
            &mut xacc,
            &mut yacc,
            &mut cache,
        )?)
    }

    /// Calculates the Jacobian `J(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let j = geom.jacobian(0.015, 2.0*PI)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn jacobian(&self, psip: Flux, theta: Radians) -> Result<f64> {
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        Ok(self.jacobian_interp.eval(
            &self.psip_data,
            &self.theta_data,
            &self.jacobian_flat_data,
            psip,
            theta.rem_euclid(TAU),
            &mut xacc,
            &mut yacc,
            &mut cache,
        )?)
    }
}

/// Getters
impl Geometry {
    /// Returns the `R(ψp, θ)` data as a 2D array, in C order.
    pub fn rlab_data(&self) -> Array2<f64> {
        // Array is in Fortran order.
        let shape = (self.theta_data.len(), self.psip_data.len());
        Array2::from_shape_vec(shape, self.rlab_flat_data.clone())
            .expect("shape is correct by definition")
            .reversed_axes()
    }

    /// Returns the `Z(ψp, θ)` data as a 2D array, in C order.
    pub fn zlab_data(&self) -> Array2<f64> {
        // Array is in Fortran order.
        let shape = (self.theta_data.len(), self.psip_data.len());
        Array2::from_shape_vec(shape, self.zlab_flat_data.clone())
            .expect("shape is correct by definition")
            .reversed_axes()
    }

    /// Returns the `J(ψp, θ)` data as a 2D array, in C order.
    pub fn jacobian_data(&self) -> Array2<f64> {
        // Array is in Fortran order.
        let shape = (self.theta_data.len(), self.psip_data.len());
        Array2::from_shape_vec(shape, self.jacobian_flat_data.clone())
            .expect("shape is correct by definition")
            .reversed_axes()
    }

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
        // `r_data` is always non-empty, otherwise `Geometry` cannot be constructed
        self.r_data.last().copied().expect("array non-empty")
    }

    /// Returns the poloidal flux's value at the wall `ψp_wall` **in Normalized Units**.
    pub fn psip_wall(&self) -> f64 {
        // `psip_data` is always non-empty, otherwise `Geometry` cannot be constructed
        self.psip_data.last().copied().expect("array non-empty")
    }

    /// Returns the toroidal flux's value at the wall `ψ_wall` **in Normalized Units**.
    pub fn psi_wall(&self) -> f64 {
        // `psi_data` is always non-empty, otherwise `Geometry` cannot be constructed
        self.psi_data.last().copied().expect("array non-empty")
    }

    array1D_getter_impl!(theta_data, theta_data, Radians);
    array1D_getter_impl!(psip_data, psip_data, Flux);
    array1D_getter_impl!(psi_data, psi_data, Flux);
    array1D_getter_impl!(r_data, r_data, Length);
}

impl std::fmt::Debug for Geometry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Geometry")
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
    use config::STUB_NETCDF_PATH;

    fn create_geometry() -> Geometry {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        Geometry::from_dataset(&path, "steffen", "bicubic").unwrap()
    }

    #[test]
    fn test_geometry_creation() {
        let g = create_geometry();
        let _ = format!("{g:?}");
    }

    #[test]
    fn test_getters() {
        let g = create_geometry();
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
        let g = create_geometry();

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
