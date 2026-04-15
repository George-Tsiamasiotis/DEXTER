mod macros;
mod pyerror;

mod pycommon;

mod pybfields;
mod pycurrents;
mod pygeometries;
mod pyharmonics;
mod pyperturbation;
mod pyqfactors;

mod pylibrium_misc;

mod pycoms;
mod pyparticle;
mod pyqueue;

use pyo3::prelude::*;

#[pymodule]
#[rustfmt::skip]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Equilibrium
    m.add_function(wrap_pyfunction!(pycommon::get_max_threads, m)?)?;
    m.add_function(wrap_pyfunction!(pycommon::set_num_threads, m)?)?;
    m.add_class::<pylibrium_misc::PyLastClosedFluxSurface>()?;
    m.add_class::<pygeometries::PyLarGeometry>()?;
    m.add_class::<pygeometries::PyNcGeometry>()?;
    m.add_class::<pyqfactors::PyNcQfactor>()?;
    m.add_class::<pyqfactors::PyUnityQfactor>()?;
    m.add_class::<pyqfactors::PyParabolicQfactor>()?;
    m.add_class::<pycurrents::PyLarCurrent>()?;
    m.add_class::<pycurrents::PyNcCurrent>()?;
    m.add_class::<pybfields::PyLarBfield>()?;
    m.add_class::<pybfields::PyNcBfield>()?;
    m.add_class::<pyharmonics::PyCosHarmonic>()?;
    m.add_class::<pyharmonics::PyNcHarmonic>()?;
    m.add_class::<pyperturbation::PyCosPerturbation>()?;
    m.add_class::<pyperturbation::PyNcPerturbation>()?;
    // Simulation
    m.add_class::<pycoms::PyParabola>()?;
    m.add_class::<pycoms::PyEnergyPzetaPlane>()?;
    m.add_class::<pycoms::PyCOMs>()?;
    m.add_class::<pyparticle::PyInitialFlux>()?;
    m.add_class::<pyparticle::PyInitialConditions>()?;
    m.add_class::<pyparticle::PyIntersectParams>()?;
    m.add_class::<pyparticle::PyParticle>()?;
    m.add_class::<pyqueue::PyInitialFluxArray>()?;
    m.add_class::<pyqueue::PyQueueInitialConditions>()?;
    m.add_class::<pyqueue::PyQueue>()?;
    Ok(())
}
