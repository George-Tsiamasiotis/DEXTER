//! Machine representation objects.

pub(crate) mod flux_coordinate_array;
pub(crate) mod getters;

pub(crate) mod bfields;
pub(crate) mod currents;
pub(crate) mod efields;
pub(crate) mod flute_mode;
pub(crate) mod geometries;
pub(crate) mod nc_flute_mode;
pub(crate) mod perturbation;
pub(crate) mod qfactors;

pub(crate) use flux_coordinate_array::{FluxCoordinateArray, FluxCoordinateState};

/// Describes the type of machine.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MachineType {
    /// Machine reconstructed from numerical data.
    ///
    /// Evaluations are calculated by interpolation.
    Numerical,
    /// Machine described by analytical formulas.
    ///
    /// Evaluations are calculated by simply evaluating the formulas.
    Analytical,
}

/// Debug-asserts that `r` is non-negative.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_non_negative_r {
    ($r: expr) => {
        debug_assert!($r >= 0.0, "Encountered negative r")
    };
}

/// Debug-asserts that `flux` is non-negative.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_non_negative_flux {
    ($flux: expr) => {
        debug_assert!($flux.value() >= 0.0, "Encountered negative magnetic flux")
    };
}

/// Debug-asserts that `psit` is non-negative.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_non_negative_psit {
    ($psit: expr) => {
        debug_assert!($psit >= 0.0, "Encountered negative ψt")
    };
}

/// Debug-asserts that `psip` is non-negative.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_non_negative_psip {
    ($psip: expr) => {
        debug_assert!($psip >= 0.0, "Encountered negative ψp")
    };
}

/// Debug-asserts that `theta` is 2π-modulo.
///
/// This is needed since the B array padding extends the interpolation domain.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_is_2pi_modulo {
    ($theta: expr) => {
        debug_assert!((0.0..=TAU).contains(&$theta), "Encountered non-2π-modulo θ")
    };
}

/// Debug-asserts that the wrapped expression `is_finite()` and returns it.
///
/// Use `std::debug_assert` and `std::dbg` syntax.
#[doc(hidden)]
#[macro_export]
macro_rules! debug_assert_is_finite {
    ($arg: expr) => {
        if std::cfg!(debug_assertions) {
            if $arg.is_finite() {
                $arg
            } else {
                panic!(
                    "[{}:{}:{}]: NaN value encountered",
                    std::file!(),
                    std::line!(),
                    std::column!(),
                )
            }
        } else {
            $arg
        }
    };
}
