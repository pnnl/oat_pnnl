//! # oat_python: python wrappers for the oat_rust repository

extern crate serde;
extern crate serde_json;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3_log;



/// Return the transpose of a list of lists
/// 
/// We regard the input as a sparse 0-1 matrix in vector-of-rowvectors format
#[pyfunction]
fn transpose_listlist( vecvec: Vec<Vec<usize>>) -> PyResult<Vec<Vec<usize>>> {
    // precompute the number of columns of the untransposed matrix
    // note we have to add 1, since arrays are 0-indexed and we
    // want to be able to use the max value as an index
    let ncol = vecvec.iter().flatten().max().unwrap_or(&0).clone() + 1; 
    let mut transposed = vec![vec![]; ncol];

    for (rowind, row) in vecvec.iter().enumerate() {
        for colind in row {
            transposed[*colind].push(rowind)
        }
    }
    Ok(transposed)
}

/// Return the transpose of a list of lists (SUBROUTINE)
/// 
/// We regard the input as a sparse 0-1 matrix in vector-of-rowvectors format
pub fn unique_row_indices_helper( vecvec:& Vec<Vec<usize>>) -> Vec<usize> {
    let mut uindices = Vec::new();
    let mut include;
    for (rowind, row) in vecvec.iter().enumerate() {
        include = true;
        for priorind in uindices.iter() {            
            if row == &vecvec[*priorind] { include = false; break }
        }
        if include { uindices.push(rowind) };
    }
    uindices
}

/// Return the transpose of a list of lists
/// 
/// We regard the input as a sparse 0-1 matrix in vector-of-rowvectors format
#[pyfunction]
fn unique_row_indices( vecvec: Vec<Vec<usize>>) -> PyResult<Vec<usize>> {
    Ok(unique_row_indices_helper( & vecvec))
}

/// Return the transpose of a list of lists
/// 
/// We regard the input as a sparse 0-1 matrix in vector-of-rowvectors format
#[pyfunction]
fn unique_rows( vecvec: Vec<Vec<usize>>) -> PyResult<Vec<Vec<usize>>> {
    let uindices = unique_row_indices_helper(&vecvec);
    let urows = uindices.iter().map(|x| vecvec[*x].clone() );
    Ok(urows.collect())
}


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn oat_python(_py: Python, m: &PyModule) -> PyResult<()> {
    pyo3_log::init(); 
    m.add_function(wrap_pyfunction!(unique_row_indices, m)?)?;
    m.add_function(wrap_pyfunction!(unique_rows, m)?)?;
    m.add_function(wrap_pyfunction!(transpose_listlist, m)?)?;

    Ok(())
}
