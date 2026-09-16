mod args;
mod common;
mod error;
mod machine;
mod macros;
mod simulate;

pub use args::*;
pub use common::*;
pub use error::*;
pub use machine::*;
pub use simulate::*;

pub type Result<T> = std::result::Result<T, DexterError>;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(common::get_max_threads, m)?)?;
    m.add_function(wrap_pyfunction!(common::set_num_threads, m)?)?;
    m.add_class::<machine::PyMagneticFlux>()?;
    m.add_class::<machine::PyGeometry>()?;
    m.add_class::<machine::PyQfactor>()?;
    m.add_class::<machine::PyCurrent>()?;
    m.add_class::<machine::PyBfield>()?;
    m.add_class::<machine::PyMode>()?;
    m.add_class::<machine::PyPerturbation>()?;
    m.add_class::<simulate::PyInitialConditions>()?;
    m.add_class::<simulate::PySolverParams>()?;
    m.add_class::<simulate::PyIntersectParams>()?;
    m.add_class::<simulate::PyParticle>()?;
    m.add_class::<simulate::PyMagneticFluxArray>()?;
    m.add_class::<simulate::PyQueueInitialConditions>()?;
    m.add_class::<simulate::PyQueue>()?;
    Ok(())
}
