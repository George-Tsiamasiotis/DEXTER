//! Definition of the `FluxCoordinateArray` object.

use dexter_common::make_interp;
use ndarray::{ArrayView1, Axis};

use crate::extract::netcdf_fields::{NC_PSIP_NORM, NC_PSIT_NORM};
use crate::{Accelerator, Interpolation1dType, MachineError, MagneticFluxArray};
use crate::{MagneticFluxKind, extract};

/// Defines whether or not a Flux coordinate can be used for evaluations.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[expect(clippy::exhaustive_enums, reason = "No other possible states")]
pub enum FluxCoordinateState {
    /// Good coordinate: values exist and are increasing. Can be used as an x coordinate for
    /// evaluations.
    Good,
    /// Bad coordinate: values exist but are not monotonic. Can only be used as y-data.
    Bad,
    /// Does not exist or is empty (fallback to when an extraction error occurs).
    NoValues,
}

impl FluxCoordinateState {
    /// Returns `Good` if the values are strictly monotonic, and `Bad` otherwise.
    fn from_values(values: &[f64]) -> Self {
        let arr = ArrayView1::from(values);
        // Use abs() to handle negative q-factors; we only care for monotonicity
        if arr
            .abs()
            .diff(1, Axis(0))
            .iter()
            .all(|diff| diff.is_sign_positive())
        {
            Self::Good
        } else {
            Self::Bad
        }
    }
}

/// Representation of a Flux Coordinate created from a netCDF file. Contains the Flux's values
/// and state.
#[derive(Clone)]
pub(crate) struct FluxCoordinateArray {
    /// The extracted flux values, if the extraction was successful.
    values: Option<Box<[f64]>>,
    /// The state of the coordinate.
    state: FluxCoordinateState,
    /// The kind of the stored flux values.
    kind: MagneticFluxKind,
}

impl FluxCoordinateArray {
    /// Creates the toroidal flux `ψt` representation.
    pub(crate) fn toroidal(file: &netcdf::File) -> Self {
        Self::build(file, NC_PSIT_NORM, MagneticFluxKind::Toroidal)
    }

    /// Creates the poroidal flux `ψp` representation.
    pub(crate) fn poloidal(file: &netcdf::File) -> Self {
        Self::build(file, NC_PSIP_NORM, MagneticFluxKind::Poloidal)
    }

    /// Creates the Flux Coordinate from raw values.
    pub(crate) fn from_raw_values(values: &[f64], kind: MagneticFluxKind) -> Self {
        Self {
            values: Some(values.to_vec().into_boxed_slice()),
            state: FluxCoordinateState::from_values(values),
            kind,
        }
    }

    /// Extracts the values and sets the state.
    ///
    /// It is not guaranteed that both fluxes (good or bad) exist in a dataset, so we dont want to
    /// return an error in case the corresponding field is not found in the netCDF file.
    fn build(file: &netcdf::File, netcdf_field: &str, kind: MagneticFluxKind) -> Self {
        if let Ok(array) = extract::array_1d(file, netcdf_field)
            && !array.is_empty()
        {
            let values = array.to_vec().into_boxed_slice();
            let state = FluxCoordinateState::from_values(&values);
            Self {
                values: Some(values),
                state,
                kind,
            }
        } else {
            Self {
                values: None,
                state: FluxCoordinateState::NoValues,
                kind,
            }
        }
    }

    /// Returns the state of the flux coordinate.
    pub(crate) fn state(&self) -> FluxCoordinateState {
        self.state
    }

    /// Returns the flux values, if the exist.
    pub(crate) fn values(&self) -> Option<&[f64]> {
        self.values.as_deref()
    }

    /// Returns a slice to the wrapped values, panicking if they do not exist.
    ///
    /// Since the values are referenced by the interpolators all the time, we do the unwrapping
    /// here to avoid messing up the interpolation calls, which results in very ugly code.
    ///
    /// This method should only be used in places where we know the values exists, for example
    /// under an `state != FluxCoordinateState::NoValues` guard.
    pub(crate) fn uvalues(&self) -> &[f64] {
        match self.values.as_ref() {
            Some(values) => values,
            None => unreachable!("Should only be called when values are guaranteed to exist"),
        }
    }

    /// Returns the flux value at the last closed flux surface.
    pub(crate) fn last_value(&self) -> Option<f64> {
        self.values
            .as_ref()
            .and_then(|values| values.last().copied())
    }
}

impl From<MagneticFluxArray> for FluxCoordinateArray {
    fn from(array: MagneticFluxArray) -> Self {
        let values = array.values().to_vec().into_boxed_slice();
        let state = FluxCoordinateState::from_values(&values);
        Self {
            values: Some(values),
            state,
            kind: array.kind(),
        }
    }
}

impl From<&MagneticFluxArray> for FluxCoordinateArray {
    fn from(array: &MagneticFluxArray) -> Self {
        let values = array.values().to_vec().into_boxed_slice();
        let state = FluxCoordinateState::from_values(&values);
        Self {
            values: Some(values),
            state,
            kind: array.kind(),
        }
    }
}

/// Integrators.
impl FluxCoordinateArray {
    /// Checks if the input `ψ` and `q` arrays are valid.
    fn check_input_data(flux_values: &[f64], q_values: &[f64]) -> Result<(), MachineError> {
        if flux_values.is_empty() || q_values.is_empty() {
            return Err(MachineError::DataError(
                "Input arrays cannot be empty".into(),
            ));
        }
        if flux_values[0] != 0.0 {
            return Err(MachineError::DataError(
                "The first value of 'flux_array' must be 0.0".into(),
            ));
        }
        if flux_values.len() != q_values.len() {
            return Err(MachineError::DataError(
                "Input arrays must have the same size".into(),
            ));
        }
        if !ArrayView1::from(flux_values)
            .diff(1, Axis(0))
            .iter()
            .all(|diff| diff.is_sign_positive())
        {
            return Err(MachineError::DataError(
                "Input magnetic flux array must be monotonic".into(),
            ));
        }
        Ok(())
    }

    /// Creates a `FluxCoordinateArray` of kind `Poloidal` by integrating `ι(ψt)` as defined by the
    /// `psit_values` and `i_values` slices.
    ///
    /// # Errors
    ///
    /// This method returns a [`MachineError`] if [`make_interp`] fails.
    pub(crate) fn psip_from_psit(
        psit_values: &[f64],
        q_values: &[f64],
        interp_type: Interpolation1dType,
    ) -> Result<Self, MachineError> {
        Self::check_input_data(psit_values, q_values)?;

        let i_values = &q_values.iter().map(|q| q.recip()).collect::<Box<[f64]>>();
        let i_of_psit_interp = make_interp(interp_type, psit_values, i_values)?;
        let acc = &mut Accelerator::new();

        let psip_values: Box<[f64]> = psit_values
            .iter()
            .map(|psit| {
                i_of_psit_interp
                    .eval_integ(psit_values, i_values, 0.0, *psit, acc)
                    .expect("Cannot fail by definition of the interpolator")
            })
            .collect();

        Ok(Self::from_raw_values(
            &psip_values,
            MagneticFluxKind::Poloidal,
        ))
    }

    /// Creates a `FluxCoordinateArray` of kind `Toroidal` by integrating `q(ψp)` as defined by the
    /// `psip_values` and `q_values` slices.
    ///
    /// # Errors
    ///
    /// This method returns a [`MachineError`] if [`make_interp`] fails.
    pub(crate) fn psit_from_psip(
        psip_values: &[f64],
        q_values: &[f64],
        interp_type: Interpolation1dType,
    ) -> Result<Self, MachineError> {
        Self::check_input_data(psip_values, q_values)?;

        let q_of_psip_interp = make_interp(interp_type, psip_values, q_values)?;
        let acc = &mut Accelerator::new();

        let psit_values: Box<[f64]> = psip_values
            .iter()
            .map(|psit| {
                q_of_psip_interp
                    .eval_integ(psip_values, q_values, 0.0, *psit, acc)
                    .expect("Cannot fail by definition of the interpolator")
            })
            .collect();
        Ok(Self::from_raw_values(
            &psit_values,
            MagneticFluxKind::Toroidal,
        ))
    }
}

impl std::fmt::Debug for FluxCoordinateArray {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let len = match self.values.as_ref() {
            Some(values) => format!("{}", values.len()),
            None => String::from("No values"),
        };
        let last = match self.last_value() {
            Some(value) => format!("{value:.6e}"),
            None => String::from("No value"),
        };
        f.debug_struct("FluxCoordinateArray")
            .field("state", &self.state)
            .field("kind", &self.kind)
            .field("len", &len)
            .field("last value", &last)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::extract::{TEST_NETCDF_PATH, open};

    use super::*;

    fn create_empty_flux() -> FluxCoordinateArray {
        FluxCoordinateArray {
            values: None,
            state: FluxCoordinateState::NoValues,
            kind: MagneticFluxKind::Toroidal,
        }
    }

    #[test]
    #[should_panic]
    fn get_uvalues_empty() {
        let flux = create_empty_flux();
        let _ = format!("{flux:?}");
        let _ = flux.uvalues();
    }

    #[test]
    fn get_values_empty() {
        let flux = create_empty_flux();
        let _ = format!("{flux:?}");
        assert!(flux.values().is_none());
    }

    #[test]
    fn non_existing_flux_variable() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let file = open(&path).unwrap();
        let noflux = FluxCoordinateArray::build(&file, "NOT_A_FIELD", MagneticFluxKind::Toroidal);
        assert_eq!(noflux.state, FluxCoordinateState::NoValues);
        assert!(noflux.values.is_none());
    }
}
