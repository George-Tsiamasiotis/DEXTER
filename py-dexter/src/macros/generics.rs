/// Manually monomorphizes Partlcle's `integrate()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_particle_integration_impl {
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
                // IntegrationParams
                method: Option<Bound<'py, PyAny>>,
                max_steps: Option<usize>,
                first_step: Option<f64>,
                safety_factor: Option<f64>,
                energy_rel_tol: Option<f64>,
                energy_abs_tol: Option<f64>,
                error_rel_tol: Option<f64>,
                error_abs_tol: Option<f64>,
            ) -> PyResult<()> {
                let mut params = IntegrationParams::default();
                if let Some(method) = method {
                    params.method = Self::resolve_stepping_method(method)?
                }
                max_steps.inspect(|v| params.max_steps = *v);
                first_step.inspect(|v| params.first_step = *v);
                safety_factor.inspect(|v| params.safety_factor = *v);
                energy_rel_tol.inspect(|v| params.energy_rel_tol = *v);
                energy_abs_tol.inspect(|v| params.energy_abs_tol = *v);
                error_rel_tol.inspect(|v| params.error_rel_tol = *v);
                error_abs_tol.inspect(|v| params.error_abs_tol = *v);
                self.0
                    .integrate(&qfactor.0, &current.0, &bfield.0, &per.0, teval, &params);
                Ok(())
            }
        }
    };
}
