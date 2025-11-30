use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline, make_spline};
use utils::array1D_getter_impl;

use crate::Flux;
use crate::Result;

use ndarray::Array1;
use safe_unwrap::safe_unwrap;

/// q-factor reconstructed from a netCDF file.
pub struct Qfactor {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,
    /// Spline over the q-factor data, as a function of ψp.
    q_spline: DynSpline<f64>,
    /// Spline over the toroidal flux data, as a function of ψp.
    psi_spline: DynSpline<f64>,
}

// Creation
impl Qfactor {
    /// Constructs a [`Qfactor`] from a netCDF file at `path`, with spline of `typ` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
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
        let psi_data = extract_1d_array(&f, PSI_NORM)?
            .as_standard_layout()
            .to_owned();
        // R_NORM isn't really useful at the moment
        let q_data = extract_1d_array(&f, Q)?.as_standard_layout().to_owned();

        let q_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", q_data.as_slice()),
        )?;
        let psi_spline = make_spline(
            typ,
            safe_unwrap!("array is non-empty", psip_data.as_slice()),
            safe_unwrap!("array is non-empty", psi_data.as_slice()),
        )?;

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            q_spline,
            psi_spline,
        })
    }
}

// Interpolation
impl Qfactor {
    /// Calculates the q-factor `q(ψp)`.
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
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let q = qfactor.q(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn q(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.q_spline.eval(psip, acc)?)
    }

    /// Calculates the toroidal flux `ψ(ψp)`.
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
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let psi = qfactor.psi(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn psi(&self, psip: Flux, acc: &mut Accelerator) -> Result<Flux> {
        Ok(self.psi_spline.eval(psip, acc)?)
    }

    /// Calculates the derivative `dψ/dψp`.
    ///
    /// This value should always equal `q(ψp)`.
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
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dpsi_dpsip = qfactor.dpsi_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn dpsi_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.psi_spline.eval_deriv(psip, acc)?)
    }
}

// Getters
impl Qfactor {
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
        self.q_spline.xa.len()
    }

    /// Returns the poloidal flux's value at the wall `ψp_wall` **in Normalized Units**.
    pub fn psip_wall(&self) -> f64 {
        safe_unwrap!("array is non-empty", self.q_spline.xa.last().copied())
    }

    /// Returns the toroidal flux's value at the wall `ψ_wall` **in Normalized Units**.
    pub fn psi_wall(&self) -> f64 {
        safe_unwrap!("array is non-empty", self.psi_spline.ya.last().copied())
    }

    array1D_getter_impl!(psip_data, q_spline.xa, Flux);
    array1D_getter_impl!(psi_data, psi_spline.ya, Flux);
    array1D_getter_impl!(q_data, q_spline.ya, f64);
}

impl std::fmt::Debug for Qfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Qfactor")
            .field("path", &self.path())
            .field("typ", &self.typ())
            .field("ψp_wall", &format!("{:.7}", self.psip_wall()))
            .field("ψ_wall", &format!("{:.7}", self.psi_wall()))
            .field("len", &self.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use config::STUB_NETCDF_PATH;

    fn create_qfactor() -> Qfactor {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        Qfactor::from_dataset(&path, "akima").unwrap()
    }

    #[test]
    fn test_qfactor_creation() {
        let q = create_qfactor();
        let _ = format!("{q:?}");
    }

    #[test]
    fn test_getters() {
        let q = create_qfactor();
        q.path();
        q.typ();
        q.psip_wall();
        q.psi_wall();
        q.len();

        assert_eq!(q.psip_data().ndim(), 1);
        assert_eq!(q.q_data().ndim(), 1);
        assert_eq!(q.psi_data().ndim(), 1);
    }

    #[test]
    fn test_spline_evaluation() {
        let q = create_qfactor();
        let mut acc = Accelerator::new();

        let psip = 0.015;
        q.q(psip, &mut acc).unwrap();
        q.psi(psip, &mut acc).unwrap();
    }
}
