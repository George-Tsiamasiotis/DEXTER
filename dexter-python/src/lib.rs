//! Module and sub-module final exports
//!
//! # Notes
//!
//! - Use _<module-name>() to avoid ambiguity and overwrite the name with
//!   `#[pymodule(name = "<module-name>")]`.

mod pylibrium;

mod pyerrors;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Equilibrium
    m.add_class::<pylibrium::PyQfactor>()?;
    m.add_class::<pylibrium::PyQfactor>()?;
    m.add_class::<pylibrium::PyCurrents>()?;
    m.add_class::<pylibrium::PyBfield>()?;
    m.add_class::<pylibrium::PyHarmonic>()?;
    m.add_class::<pylibrium::PyPerturbation>()?;
    Ok(())
}
