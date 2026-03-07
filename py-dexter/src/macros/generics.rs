/// Manually monomorphizes Partlcle's `integrate()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_particle_integrate_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty, $P:ty) => {
        #[pymethods]
        impl PyParticle {
            #[allow(non_snake_case)]
            pub fn $fun_name<'py>(
                &mut self,
                qfactor: &$Q,
                current: &$C,
                bfield: &$B,
                per: &$P,
                teval: (f64, f64),
                // SolverParams
                method: Option<Bound<'py, PyAny>>,
                max_steps: Option<usize>,
                first_step: Option<f64>,
                safety_factor: Option<f64>,
                energy_rel_tol: Option<f64>,
                energy_abs_tol: Option<f64>,
                error_rel_tol: Option<f64>,
                error_abs_tol: Option<f64>,
            ) -> PyResult<()> {
                // Create default and replace with passed values
                let mut solver_params = SolverParams::default();
                if let Some(method) = method {
                    solver_params.method = Self::resolve_stepping_method(method)?
                }
                max_steps.inspect(|v| solver_params.max_steps = *v);
                first_step.inspect(|v| solver_params.first_step = *v);
                safety_factor.inspect(|v| solver_params.safety_factor = *v);
                energy_rel_tol.inspect(|v| solver_params.energy_rel_tol = *v);
                energy_abs_tol.inspect(|v| solver_params.energy_abs_tol = *v);
                error_rel_tol.inspect(|v| solver_params.error_rel_tol = *v);
                error_abs_tol.inspect(|v| solver_params.error_abs_tol = *v);
                self.0.integrate(
                    &qfactor.0,
                    &current.0,
                    &bfield.0,
                    &per.0,
                    teval,
                    &solver_params,
                );
                Ok(())
            }
        }
    };
}

/// Manually monomorphizes Partlcle's `intersect()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_particle_intersect_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty, $P:ty) => {
        #[pymethods]
        impl PyParticle {
            #[allow(non_snake_case)]
            pub fn $fun_name<'py>(
                &mut self,
                qfactor: &$Q,
                current: &$C,
                bfield: &$B,
                per: &$P,
                intersect_params: &PyIntersectParams,
                // SolverParams
                method: Option<Bound<'py, PyAny>>,
                max_steps: Option<usize>,
                first_step: Option<f64>,
                safety_factor: Option<f64>,
                energy_rel_tol: Option<f64>,
                energy_abs_tol: Option<f64>,
                error_rel_tol: Option<f64>,
                error_abs_tol: Option<f64>,
            ) -> PyResult<()> {
                // Create default and replace with passed values
                let mut solver_params = SolverParams::default();
                if let Some(method) = method {
                    solver_params.method = Self::resolve_stepping_method(method)?
                }
                max_steps.inspect(|v| solver_params.max_steps = *v);
                first_step.inspect(|v| solver_params.first_step = *v);
                safety_factor.inspect(|v| solver_params.safety_factor = *v);
                energy_rel_tol.inspect(|v| solver_params.energy_rel_tol = *v);
                energy_abs_tol.inspect(|v| solver_params.energy_abs_tol = *v);
                error_rel_tol.inspect(|v| solver_params.error_rel_tol = *v);
                error_abs_tol.inspect(|v| solver_params.error_abs_tol = *v);
                self.0.intersect(
                    &qfactor.0,
                    &current.0,
                    &bfield.0,
                    &per.0,
                    &intersect_params.0,
                    &solver_params,
                );
                Ok(())
            }
        }
    };
}
