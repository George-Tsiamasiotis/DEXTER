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
        match arr.diff(1, Axis(0)).iter().all(|d| d.signum() == 1.0) {
            true => Self::Good,
            false => Self::Bad,
        }
    }
}

/// Contains the Flux's values and state.
///
/// Even if the state is None (error trying to extract the values), we set the values as an empty
/// vec instead of defining it as an Option<Vec<f64>>, since this greatly simplifies the code.
/// Moreover, the logic on the evaluation sites checks the state anyway, so an extra `is_some()`
/// would be redundant. Finally, even if an evaluation is called, the empty vec would cause a
/// panic, which is what we want.
#[derive(Clone)]
pub struct NcFlux {
    pub(crate) values: Vec<f64>,
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
                values: arr.to_vec(),
                state: NcFluxState::from_array(arr.view()),
            },
            Err(_) => Self {
                values: Vec::new(),
                state: NcFluxState::None,
            },
        }
    }

    /// Returns the flux's value at the wall, if it exists.
    pub(crate) fn wall_value(&self) -> Option<f64> {
        self.values.last().copied()
    }
}

impl std::fmt::Debug for NcFlux {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Flux")
            .field("state", &self.state)
            .field("len", &self.values.len())
            .field("wall value", &self.wall_value())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use crate::extract::{TEST_NETCDF_PATH, open};

    use super::*;

    #[test]
    fn non_existing_flux_variable() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let file = open(&path).unwrap();
        let noflux = NcFlux::build(&file, "NOT_A_FIELD");
        assert_eq!(noflux.state, NcFluxState::None);
        assert!(noflux.values.is_empty());
    }
}
