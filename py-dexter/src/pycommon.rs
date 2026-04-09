use pyo3::prelude::*;

#[pyfunction(name = "_py_set_num_threads")]
pub fn set_num_threads(num: usize) {
    dexter::dexter_simulate::set_num_threads(num);
}

#[pyfunction(name = "_py_get_max_threads")]
pub fn get_max_threads() -> usize {
    dexter::dexter_simulate::get_max_threads()
}
