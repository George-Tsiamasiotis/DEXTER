/// Manually monomorphizes Partlcle's `integrate()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_particle_integrate_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty, $P:ty) => {
        #[pymethods]
        impl PyParticle {
            #[allow(non_snake_case)]
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
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
                    solver_params.method = resolve_stepping_method(method)?
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
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
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
                    solver_params.method = resolve_stepping_method(method)?
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

/// Manually monomorphizes Partlcle's `close()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_particle_close_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty, $P:ty) => {
        #[pymethods]
        impl PyParticle {
            #[allow(non_snake_case)]
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
            pub fn $fun_name<'py>(
                &mut self,
                qfactor: &$Q,
                current: &$C,
                bfield: &$B,
                per: &$P,
                periods: usize,
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
                    solver_params.method = resolve_stepping_method(method)?
                }
                max_steps.inspect(|v| solver_params.max_steps = *v);
                first_step.inspect(|v| solver_params.first_step = *v);
                safety_factor.inspect(|v| solver_params.safety_factor = *v);
                energy_rel_tol.inspect(|v| solver_params.energy_rel_tol = *v);
                energy_abs_tol.inspect(|v| solver_params.energy_abs_tol = *v);
                error_rel_tol.inspect(|v| solver_params.error_rel_tol = *v);
                error_abs_tol.inspect(|v| solver_params.error_abs_tol = *v);
                self.0.close(
                    &qfactor.0,
                    &current.0,
                    &bfield.0,
                    &per.0,
                    periods,
                    &solver_params,
                );
                Ok(())
            }
        }
    };
}

// ===============================================================================================

/// Manually monomorphizes Queue's `integrate()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_queue_integrate_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty, $P:ty) => {
        #[pymethods]
        impl PyQueue {
            #[allow(non_snake_case)]
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
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
                    solver_params.method = resolve_stepping_method(method)?
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

/// Manually monomorphizes Queue's `intersect()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_queue_intersect_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty, $P:ty) => {
        #[pymethods]
        impl PyQueue {
            #[allow(non_snake_case)]
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
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
                    solver_params.method = resolve_stepping_method(method)?
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

/// Manually monomorphizes Queue's `close()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_queue_close_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty, $P:ty) => {
        #[pymethods]
        impl PyQueue {
            #[allow(non_snake_case)]
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
            pub fn $fun_name<'py>(
                &mut self,
                qfactor: &$Q,
                current: &$C,
                bfield: &$B,
                per: &$P,
                periods: usize,
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
                    solver_params.method = resolve_stepping_method(method)?
                }
                max_steps.inspect(|v| solver_params.max_steps = *v);
                first_step.inspect(|v| solver_params.first_step = *v);
                safety_factor.inspect(|v| solver_params.safety_factor = *v);
                energy_rel_tol.inspect(|v| solver_params.energy_rel_tol = *v);
                energy_abs_tol.inspect(|v| solver_params.energy_abs_tol = *v);
                error_rel_tol.inspect(|v| solver_params.error_rel_tol = *v);
                error_abs_tol.inspect(|v| solver_params.error_abs_tol = *v);
                self.0.close(
                    &qfactor.0,
                    &current.0,
                    &bfield.0,
                    &per.0,
                    periods,
                    &solver_params,
                );
                Ok(())
            }
        }
    };
}

/// ==============================================================================================
/// ==============================================================================================

/// Manually monomorphizes COMs' `energy_of_psi_grid()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_energy_of_psi_grid_impl {
    ($fun_name:ident, $Q:ty, $C:ty, $B:ty) => {
        #[pymethods]
        impl PyCOMs {
            #[allow(non_snake_case)]
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
            pub fn $fun_name<'py>(
                &self,
                py: Python<'py>,
                qfactor: &$Q,
                current: &$C,
                bfield: &$B,
                psi_values: Vec<f64>,
                theta_values: Vec<f64>,
            ) -> PyResult<Bound<'py, PyArray2<f64>>> {
                Ok(self
                    .0
                    .energy_of_psi_grid(
                        &qfactor.0,
                        &current.0,
                        &bfield.0,
                        &Array1::from_vec(theta_values),
                        &Array1::from_vec(psi_values),
                    )
                    .or_else(|err| Err(PyErr::new::<PyValueError, _>(err.to_string())))?
                    .into_pyarray(py))
            }
        }
    };
}

/// Manually monomorphizes COMs' `energy_of_psip_grid()` method for the specific equilibrium objects.
#[macro_export]
macro_rules! generic_energy_of_psip_grid_impl {
    ($fun_name:ident, $C:ty, $B:ty) => {
        #[pymethods]
        impl PyCOMs {
            #[allow(non_snake_case)]
            #[expect(clippy::too_many_arguments, reason = "python kwargs")]
            pub fn $fun_name<'py>(
                &self,
                py: Python<'py>,
                current: &$C,
                bfield: &$B,
                psip_values: Vec<f64>,
                theta_values: Vec<f64>,
            ) -> PyResult<Bound<'py, PyArray2<f64>>> {
                Ok(self
                    .0
                    .energy_of_psip_grid(
                        &current.0,
                        &bfield.0,
                        &Array1::from_vec(theta_values),
                        &Array1::from_vec(psip_values),
                    )
                    .or_else(|err| Err(PyErr::new::<PyValueError, _>(err.to_string())))?
                    .into_pyarray(py))
            }
        }
    };
}
