//! Module and sub-module final exports

mod macros;
mod pylibrium;

mod pyerrors;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Equilibrium
    m.add_class::<pylibrium::PyNcGeometry>()?;
    m.add_class::<pylibrium::PyNcQfactor>()?;
    m.add_class::<pylibrium::PyUnityQfactor>()?;
    m.add_class::<pylibrium::PyNcCurrent>()?;
    m.add_class::<pylibrium::PyLarCurrent>()?;
    m.add_class::<pylibrium::PyNcBfield>()?;
    m.add_class::<pylibrium::PyNcHarmonic>()?;
    m.add_class::<pylibrium::PyNcPerturbation>()?;
    Ok(())
}
