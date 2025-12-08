// ==================== NetCDF

/// Stub NetCDF file path to be used for unit tests.
///
/// Must be *relative* to the crate's `CARGO_MANIFEST_DIR`. Can be created with
/// `equilibrium/scripts/npz_to_netcdf.py` from an `.npz` file.
pub const STUB_NETCDF_PATH: &str = "../../data/stub_netcdf.nc";

// ==================== Solver

/// The maximum amount of steps a particle can make before terminating its integration.
pub const MAX_STEPS: usize = 1_000_000;

/// The maximum amount of steps a particle can make when integrating for a single period.
pub const SINGLE_PERIOD_MAX_STEPS: usize = 100_000;

/// The initial time step for the RKF45 adaptive step method. Should be small enough to account
/// for fast particles. The value is empirical.
pub const RKF45_FIRST_STEP: f64 = 1e-1;

/// The safety factor of the solver. Should be less than 1.0
pub const SAFETY_FACTOR: f64 = 0.9;

/// The relative tolerance of the energy error in every step.
pub const ENERGY_REL_TOL: f64 = 1e-10;

/// The relative tolerance of the stepping error in every step.
pub const STEP_REL_TOL: f64 = 1e-11;

/// The θ-span "relaxation" factor. When integrating for a single period, we can treat particles
/// where |θ_initial - θ_final| ≈ 2π (approx) as passing. The critirion is
/// `|θ[-1] - θ[0]| > TAU - TRAPPED_THRESHOLD`
pub const TRAPPED_THRESHOLD: f64 = 1e-7;

// ==================== Mapping

/// The maximum allowed absolute difference between two consecutive intersections.
pub const MAP_THRESHOLD: f64 = 1e-9;

// ==================== Misc

/// The starting capacity of the Evolution time series vectors.
pub const EVOLUTION_INIT_CAPACITY: usize = 2000;
