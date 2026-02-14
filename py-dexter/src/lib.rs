mod macros;
mod pyerror;

mod pylibrium;

use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Equilibrium
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
    Ok(())
}
