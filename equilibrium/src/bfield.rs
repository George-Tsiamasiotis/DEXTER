use std::path::PathBuf;

use rsl_interpolation::{Accelerator, Cache, DynSpline2d, make_spline2d};
use utils::array1D_getter_impl;

use crate::Result;
use crate::{Flux, Radians};

use ndarray::{Array1, Array2, meshgrid};

/// Magnetic field reconstructed from a netCDF file.
pub struct Bfield {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 2D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.Interp2dType.html#implementors
    typ: String,
    /// Spline over the magnetic field strength data, as a function of ψp, θ.
    b_spline: DynSpline2d<f64>,
}

// Creation
impl Bfield {
    /// Constructs a [`Bfield`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str) -> Result<Self> {
        use crate::extract::*;
        use config::netcdf_fields::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP_NORM)?
            .as_standard_layout()
            .to_owned();
        let theta_data = extract_1d_array(&f, THETA)?.as_standard_layout().to_owned();

        let b_data = extract_2d_array(&f, B_NORM)?;

        // `Spline.za` is in Fortran order.
        let order = ndarray::Order::ColumnMajor;
        let b_data_flat = b_data.flatten_with_order(order).to_owned();

        // `extract_array()` has already checked if the arrays are empty
        let b_spline = make_spline2d(
            typ,
            psip_data.as_slice().expect("array is non-empty"),
            theta_data.as_slice().expect("array is non-empty"),
            b_data_flat.as_slice().expect("array is non-empty"),
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            b_spline,
        })
    }
}

// Interpolation
impl Bfield {
    /// Calculates `B(ψp, θ)`,
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
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let b =  bfield.b(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn b(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self.b_spline.eval(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕B(ψp, θ) /𝜕𝜃`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let db_dtheta = bfield.db_dtheta(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn db_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_y(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕B(ψp, θ) /𝜕ψp`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let db_dpsip = bfield.db_dpsip(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn db_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_x(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕²B(ψp, θ) /𝜕𝜓p²`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dpsip2(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn d2b_dpsip2(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_xx(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕²B(ψp, θ) /𝜕θ²`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dtheta2(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn d2b_dtheta2(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_yy(psip, mod2pi(theta), xacc, yacc, cache)?)
    }

    /// Calculates `𝜕²B(ψp, θ) /𝜕ψp𝜕θ`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
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
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dpsip_dtheta(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn d2b_dpsip_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64> {
        Ok(self
            .b_spline
            .eval_deriv_xy(psip, mod2pi(theta), xacc, yacc, cache)?)
    }
}

/// Getters
impl Bfield {
    /// Returns the `R(ψp, θ)` data as a 2D array.
    pub fn b_data(&self) -> Array2<f64> {
        // Array is in Fortran order.
        let shape = (self.b_spline.ya.len(), self.b_spline.xa.len());
        Array2::from_shape_vec(shape, self.b_spline.za.to_vec())
            .expect("shape is correct by definition")
            .reversed_axes()
    }

    /// Returns the `𝜕B(ψp, θ) /𝜕ψp` data as a 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    pub fn db_dpsip_data(&self) -> Array2<f64> {
        let psip_data = self.psip_data();
        let theta_data = self.theta_data();
        let grid = meshgrid((&psip_data, &theta_data), ndarray::MeshIndex::IJ);

        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        // If this fails, there is something wrong with the data
        Array2::from_shape_fn(self.shape(), |(i, j)| {
            (self.db_dpsip(
                grid.0[[i, j]],
                grid.1[[i, j]],
                &mut xacc,
                &mut yacc,
                &mut cache,
            ))
            .expect("something is wrong with the data")
        })
    }

    /// Returns the `𝜕B(ψp, θ) /𝜕𝜃` data as a 2D array.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    pub fn db_dtheta_data(&self) -> Array2<f64> {
        let psip_data = self.psip_data();
        let theta_data = self.theta_data();
        let grid = meshgrid((&psip_data, &theta_data), ndarray::MeshIndex::IJ);

        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();
        // If this fails, there is something wrong with the data
        Array2::from_shape_fn(self.shape(), |(i, j)| {
            (self.db_dtheta(
                grid.0[[i, j]],
                grid.1[[i, j]],
                &mut xacc,
                &mut yacc,
                &mut cache,
            ))
            .expect("something is wrong with the data")
        })
    }

    /// Returns the netCDF file's path.
    pub fn path(&self) -> PathBuf {
        self.path.clone()
    }

    /// Returns the interpolation type.
    pub fn typ(&self) -> String {
        self.typ.clone()
    }

    /// Returns the shape of the `b` array.
    pub fn shape(&self) -> (usize, usize) {
        (self.b_spline.xa.len(), self.b_spline.ya.len())
    }

    array1D_getter_impl!(psip_data, b_spline.xa, Flux);
    array1D_getter_impl!(theta_data, b_spline.ya, Flux);
}

/// Returns θ % 2π.
fn mod2pi(theta: f64) -> f64 {
    use std::f64::consts::TAU;
    theta.rem_euclid(TAU)
}

impl std::fmt::Debug for Bfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Bfield")
            .field("path", &self.path())
            .field("typ", &self.typ())
            .field("shape", &self.shape())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use config::STUB_NETCDF_PATH;

    fn create_bfield() -> Bfield {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        Bfield::from_dataset(&path, "bicubic").unwrap()
    }

    #[test]
    fn test_bfield_creation() {
        let b = create_bfield();
        let _ = format!("{b:?}");
    }

    #[test]
    fn test_getters() {
        let b = create_bfield();
        b.path();
        b.typ();
        b.shape();

        assert_eq!(b.psip_data().ndim(), 1);
        assert_eq!(b.theta_data().ndim(), 1);
        assert_eq!(b.b_data().ndim(), 2);
        assert_eq!(b.db_dpsip_data().ndim(), 2);
        assert_eq!(b.db_dtheta_data().ndim(), 2);
    }

    #[test]
    fn test_spline_evaluation() {
        let b = create_bfield();
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();

        let psip = 0.015;
        let theta = 0.0;
        b.b(psip, theta, &mut xacc, &mut yacc, &mut cache).unwrap();
        b.db_dpsip(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.db_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.d2b_dpsip2(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.d2b_dtheta2(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        b.d2b_dpsip_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
    }
}
