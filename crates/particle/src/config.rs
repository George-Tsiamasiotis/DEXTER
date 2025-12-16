#[allow(unused_imports)] // doc
use crate::Particle;

#[derive(Debug, Clone)]
/// The method used to calculate the next optimal step.
pub enum SteppingMethod {
    /// Forces the step size to be small enough so that the Energy difference from step to step is
    /// under a certain threshold. The tolerances can be adjusted with the `energy_rel_tol` and
    /// `energy_abs_tol` fields.
    EnergyAdaptiveStep,
    /// Classic RK error estimation : Adjust the step size to minimize the local truncation error.
    ErrorAdaptiveStep,
}

// ===============================================================================================

/// Ensures that all methods' configurations have the fields required by the Stepper.
pub(crate) trait StepperConfig {
    fn method(&self) -> &SteppingMethod;
    fn energy_rel_tol(&self) -> f64;
    fn energy_abs_tol(&self) -> f64;
    fn error_rel_tol(&self) -> f64;
    fn error_abs_tol(&self) -> f64;
    fn safety_factor(&self) -> f64;
}

#[rustfmt::skip]
macro_rules! stepper_config_impl {
    ($object:ident) => {
        impl StepperConfig for $object {
            #[inline(always)]
            fn method(&self) -> &SteppingMethod { &self.method }
            #[inline(always)]
            fn energy_rel_tol(&self) -> f64 { self.energy_rel_tol }
            #[inline(always)]
            fn energy_abs_tol(&self) -> f64 { self.energy_abs_tol }
            #[inline(always)]
            fn error_rel_tol(&self) -> f64 { self.error_rel_tol }
            #[inline(always)]
            fn error_abs_tol(&self) -> f64 { self.error_abs_tol }
            #[inline(always)]
            fn safety_factor(&self) -> f64 { self.safety_factor }
        }
    };
}

stepper_config_impl!(IntegrationConfig);
stepper_config_impl!(MappingConfig);
stepper_config_impl!(SinglePeriodConfig);

// ===============================================================================================

/// Defines the parameters of the [`Particle::integrate`] routine.
///
/// See [`IntegrationConfig::default`] for the default values.
#[derive(Debug, Clone)]
pub struct IntegrationConfig {
    /// The optimal step calculation method.
    pub method: SteppingMethod,
    /// The maximum amount of steps a particle can make before terminating its integration.
    pub max_steps: usize,
    /// The initial time step for the RKF45 adaptive step method. The value is empirical.
    pub first_step: f64,
    /// The safety factor of the solver. Should be less than 1.0
    pub safety_factor: f64,
    /// The relative tolerance of the energy difference in every step.
    pub energy_rel_tol: f64,
    /// The absolute tolerance of the energy difference in every step.
    pub energy_abs_tol: f64,
    /// The relative tolerance of the local truncation error in every step.
    pub error_rel_tol: f64,
    /// The absolute tolerance of the local truncation error in every step.
    pub error_abs_tol: f64,
}

impl Default for IntegrationConfig {
    fn default() -> Self {
        Self {
            method: SteppingMethod::EnergyAdaptiveStep,
            max_steps: 1_000_000,
            first_step: 1e-1,
            safety_factor: 0.9,
            energy_rel_tol: 1e-10,
            energy_abs_tol: 1e-12,
            error_rel_tol: 1e-12,
            error_abs_tol: 1e-14,
        }
    }
}

// ===============================================================================================

/// Defines the parameters of the [`Particle::map`] routine.
///
/// See [`MappingConfig::default`] for the default values.
#[derive(Debug, Clone)]
pub struct MappingConfig {
    /// The optimal step calculation method.
    pub method: SteppingMethod,
    /// The maximum amount of steps a particle can make before terminating its integration.
    pub max_steps: usize,
    /// The initial time step for the RKF45 adaptive step method. The value is empirical.
    pub first_step: f64,
    /// The safety factor of the solver. Should be less than 1.0
    pub safety_factor: f64,
    /// The relative tolerance of the energy difference in every step.
    pub energy_rel_tol: f64,
    /// The absolute tolerance of the energy difference in every step.
    pub energy_abs_tol: f64,
    /// The relative tolerance of the local truncation error in every step.
    pub error_rel_tol: f64,
    /// The absolute tolerance of the local truncation error in every step.
    pub error_abs_tol: f64,
    /// The maximum allowed absolute difference between the difference of two consecutive
    /// intersections and 2π.
    pub map_threshold: f64,
}

impl Default for MappingConfig {
    fn default() -> Self {
        Self {
            method: SteppingMethod::EnergyAdaptiveStep,
            max_steps: 1_000_000,
            first_step: 1e-1,
            safety_factor: 0.9,
            energy_rel_tol: 1e-10,
            energy_abs_tol: 1e-12,
            error_rel_tol: 1e-12,
            error_abs_tol: 1e-14,
            map_threshold: 1e-9,
        }
    }
}

// ===============================================================================================

/// Defines the parameters of the [`Particle::single_period_integrate`] routine.
///
/// See [`SinglePeriodConfig::default`] for the default values.
#[derive(Debug, Clone)]
pub struct SinglePeriodConfig {
    /// The optimal step calculation method.
    pub method: SteppingMethod,
    /// The maximum amount of steps a particle can make before terminating its integration.
    pub max_steps: usize,
    /// The initial time step for the RKF45 adaptive step method. The value is empirical.
    pub first_step: f64,
    /// The safety factor of the solver. Should be less than 1.0
    pub safety_factor: f64,
    /// The relative tolerance of the energy difference in every step.
    pub energy_rel_tol: f64,
    /// The absolute tolerance of the energy difference in every step.
    pub energy_abs_tol: f64,
    /// The relative tolerance of the local truncation error in every step.
    pub error_rel_tol: f64,
    /// The absolute tolerance of the local truncation error in every step.
    pub error_abs_tol: f64,
}

impl Default for SinglePeriodConfig {
    fn default() -> Self {
        Self {
            method: SteppingMethod::EnergyAdaptiveStep,
            max_steps: 1_000_000,
            first_step: 1e-1,
            safety_factor: 0.9,
            energy_rel_tol: 1e-10,
            energy_abs_tol: 1e-12,
            error_rel_tol: 1e-12,
            error_abs_tol: 1e-14,
        }
    }
}
