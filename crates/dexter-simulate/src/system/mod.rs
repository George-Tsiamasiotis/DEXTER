mod rkf45;

use crate::SteppingMethod;
pub(crate) use rkf45::Stepper;

/// Ensures that all methods' configurations have the fields required by the Stepper.
pub(crate) trait StepperParams {
    fn method(&self) -> &SteppingMethod;
    fn energy_rel_tol(&self) -> f64;
    fn energy_abs_tol(&self) -> f64;
    fn error_rel_tol(&self) -> f64;
    fn error_abs_tol(&self) -> f64;
    fn safety_factor(&self) -> f64;
}

#[doc(hidden)]
#[macro_export]
macro_rules! stepper_params_trait_impl {
    (&$object:ident) => {
#[rustfmt::skip]
        impl StepperParams for &$object {
            fn method(&self) -> &SteppingMethod { &self.method }
            fn energy_rel_tol(&self) -> f64 { self.energy_rel_tol }
            fn energy_abs_tol(&self) -> f64 { self.energy_abs_tol }
            fn error_rel_tol(&self) -> f64 { self.error_rel_tol }
            fn error_abs_tol(&self) -> f64 { self.error_abs_tol }
            fn safety_factor(&self) -> f64 { self.safety_factor }
        }
    };
}
