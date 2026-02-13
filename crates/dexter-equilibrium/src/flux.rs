//! Definition of the `NcFlux` object.

use ndarray::{ArrayView1, Axis};

use crate::extract::extract_1d_array;
use crate::extract::netcdf_fields::*;

/// Defines whether or not a Flux coordinate can be used for evaluations.
#[derive(Debug, PartialEq, Clone)]
pub enum NcFluxState {
    /// Good coordinate: values exist and are increasing. Can be used as an x coordinate for
    /// evaluations.
    Good,
    /// Bad coordinate: values exist but are not monotonic. Can only be used as y-data.
    Bad,
    /// Does not exist or is empty (fallback to when an extraction error occurs).
    None,
}

impl NcFluxState {
    /// Returns Self::Good if the values are strictly increasing, and Self::Bad otherwise.
    fn from_array(arr: ArrayView1<f64>) -> Self {
        // Use abs() to handle negative q-factors
        match arr.abs().diff(1, Axis(0)).iter().all(|d| d.signum() == 1.0) {
            true => Self::Good,
            false => Self::Bad,
        }
    }
}

/// Contains the Flux's values and state.
#[derive(Clone)]
pub struct NcFlux {
    values: Option<Vec<f64>>,
    pub(crate) state: NcFluxState,
}

impl NcFlux {
    /// Creates the toroidal flux `ψ` representation.
    pub(crate) fn toroidal(file: &netcdf::File) -> Self {
        Self::build(file, PSI_NORM)
    }

    /// Creates the poroidal flux `ψp` representation.
    pub(crate) fn poloidal(file: &netcdf::File) -> Self {
        Self::build(file, PSIP_NORM)
    }

    /// Extracts the values and sets the state.
    fn build(file: &netcdf::File, netcdf_field: &str) -> Self {
        match extract_1d_array(file, netcdf_field) {
            Ok(arr) => Self {
                values: Some(arr.to_vec()),
                state: NcFluxState::from_array(arr.view()),
            },
            Err(_) => Self {
                values: None,
                state: NcFluxState::None,
            },
        }
    }

    /// Returns the flux values, if the exist.
    pub(crate) fn values(&self) -> Option<&[f64]> {
        self.values.as_deref()
    }

    /// Returns a slice to the wrapped values, panicking if they do not exist.
    ///
    /// Since the values are referenced by the interpolators all the time, we do the unwrapping
    /// here to avoid messing up the interpolation calls
    ///
    /// This method should only be used in places where we know the values exists, for example
    /// under an `state != NcFluxState::None` guard.
    pub(crate) fn uvalues(&self) -> &[f64] {
        match self.values.as_ref() {
            Some(values) => values,
            None => unreachable!("Should only be called when values are guaranteed to exist"),
        }
    }

    /// Returns the flux's value at the wall, if it exists.
    pub(crate) fn wall_value(&self) -> Option<f64> {
        self.values
            .as_ref()
            .and_then(|values| values.last().copied())
    }
}

impl std::fmt::Debug for NcFlux {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let len = match self.values.as_ref() {
            Some(values) => format!("{}", values.len()),
            None => String::from("No values"),
        };
        let wall = match self.wall_value() {
            Some(value) => format!("{}", value),
            None => String::from("No value"),
        };
        f.debug_struct("Flux")
            .field("state", &self.state)
            .field("len", &len)
            .field("wall value", &wall)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::extract::{TEST_NETCDF_PATH, open};

    use super::*;

    fn create_empty_flux() -> NcFlux {
        NcFlux {
            values: None,
            state: NcFluxState::None,
        }
    }

    #[test]
    #[should_panic]
    fn get_uvalues_empty() {
        let flux = create_empty_flux();
        let _ = format!("{flux:?}");
        flux.uvalues();
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
        let noflux = NcFlux::build(&file, "NOT_A_FIELD");
        assert_eq!(noflux.state, NcFluxState::None);
        assert!(noflux.values.is_none());
    }
}
