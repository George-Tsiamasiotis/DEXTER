mod macros;
mod pyerror;

mod pycommon;
mod pylibrium;
mod pysimulate;

use pyo3::prelude::*;

use dexter::dexter_equilibrium::constants::DEFAULT_THETA_PADDING_WIDTH;
use dexter::dexter_simulate::constants::{
    DEFAULT_ENERGY_ABS_TOL, DEFAULT_ENERGY_REL_TOL, DEFAULT_ERROR_ABS_TOL, DEFAULT_ERROR_REL_TOL,
    DEFAULT_FIRST_STEP, DEFAULT_MAX_STEPS, DEFAULT_SAFETY_FACTOR, DEFAULT_STEPPING_METHOD,
};

#[pymodule]
#[rustfmt::skip]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Constants
    m.add("DEFAULT_THETA_PADDING_WIDTH", DEFAULT_THETA_PADDING_WIDTH)?;
    m.add("DEFAULT_STEPPING_METHOD", DEFAULT_STEPPING_METHOD.to_string())?;
    m.add("DEFAULT_MAX_STEPS", DEFAULT_MAX_STEPS)?;
    m.add("DEFAULT_FIRST_STEP", DEFAULT_FIRST_STEP)?;
    m.add("DEFAULT_SAFETY_FACTOR", DEFAULT_SAFETY_FACTOR)?;
    m.add("DEFAULT_ENERGY_REL_TOL", DEFAULT_ENERGY_REL_TOL)?;
    m.add("DEFAULT_ENERGY_ABS_TOL", DEFAULT_ENERGY_ABS_TOL)?;
    m.add("DEFAULT_ERROR_REL_TOL", DEFAULT_ERROR_REL_TOL)?;
    m.add("DEFAULT_ERROR_ABS_TOL", DEFAULT_ERROR_ABS_TOL)?;
    // Constants
    m.add_function(wrap_pyfunction!(pycommon::get_max_threads, m)?)?; // Equilibrium
    m.add_function(wrap_pyfunction!(pycommon::set_num_threads, m)?)?; // Equilibrium
    m.add_class::<pylibrium::PyLastClosedFluxSurface>()?;
    m.add_class::<pylibrium::PyLarGeometry>()?;
    m.add_class::<pylibrium::PyNcGeometry>()?;
    m.add_class::<pylibrium::PyNcQfactor>()?;
    m.add_class::<pylibrium::PyUnityQfactor>()?;
    m.add_class::<pylibrium::PyParabolicQfactor>()?;
    m.add_class::<pylibrium::PyLarCurrent>()?;
    m.add_class::<pylibrium::PyNcCurrent>()?;
    m.add_class::<pylibrium::PyLarBfield>()?;
    m.add_class::<pylibrium::PyNcBfield>()?;
    m.add_class::<pylibrium::PyCosHarmonic>()?;
    m.add_class::<pylibrium::PyNcHarmonic>()?;
    m.add_class::<pylibrium::PyCosPerturbation>()?;
    m.add_class::<pylibrium::PyNcPerturbation>()?;
    // Simulation
    m.add_class::<pysimulate::PyCOMs>()?;
    m.add_class::<pysimulate::PyInitialFlux>()?;
    m.add_class::<pysimulate::PyInitialConditions>()?;
    m.add_class::<pysimulate::PyIntersectParams>()?;
    m.add_class::<pysimulate::PyParticle>()?;
    m.add_class::<pysimulate::PyInitialFluxArray>()?;
    m.add_class::<pysimulate::PyQueueInitialConditions>()?;
    m.add_class::<pysimulate::PyQueue>()?;
    Ok(())
}
