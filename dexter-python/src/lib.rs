//! Module and sub-module final exports

mod pyheap;
mod pylibrium;
mod pyparticle;

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
    // Particle
    m.add_class::<pyparticle::PyInitialConditions>()?;
    m.add_class::<pyparticle::PyMappingParameters>()?;
    m.add_class::<pyparticle::PyEvolution>()?;
    m.add_class::<pyparticle::PyFrequencies>()?;
    m.add_class::<pyparticle::PyParticle>()?;
    // Heap
    m.add_class::<pyheap::PyHeapInitialConditions>()?;
    Ok(())
}
