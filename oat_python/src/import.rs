//! Import data from Python

use itertools::Itertools;
use ordered_float::OrderedFloat;
use pyo3::pyfunction;
use pyo3::prelude::*;
use pyo3::types::PyType;

use sprs::CsMatBase;




pub fn import_sparse_matrix(py: Python, scipy_csr: &PyAny )     
    -> PyResult< CsMatBase<OrderedFloat<f64>, usize, Vec<usize>, Vec<usize>, Vec<OrderedFloat<f64>>> >    
{
    // Check if the object is an instance of csr_matrix

    let shape: (usize,usize) = scipy_csr.getattr("shape").ok().unwrap().extract().ok().unwrap();
    let indptr: Vec<usize> = scipy_csr.getattr("indptr").ok().unwrap().extract().ok().unwrap();
    let indices: Vec<usize> = scipy_csr.getattr("indices").ok().unwrap().extract().ok().unwrap();
    let data: Vec< f64 > = scipy_csr.getattr("data").ok().unwrap().extract().ok().unwrap();
    let data = data.into_iter().map(|v| OrderedFloat(v)).collect_vec();

    return Ok( CsMatBase::new(
        shape, // shape: 
        indptr,
        indices,
        data,
    ) )
}

