mod macros;
mod pyerror;

mod pylibrium;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Equilibrium
    m.add_class::<pylibrium::PyNcCurrent>()?;
    Ok(())
}
