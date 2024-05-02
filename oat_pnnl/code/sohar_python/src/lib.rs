extern crate serde;
extern crate serde_json;
use ordered_float::OrderedFloat;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3_log;
use sohar_rust::simplicial::from::graph_weighted::data_structures::CliqueBoundaryMatrix;
use sohar_rust::simplicial::from::graph_weighted::user_interface::CliqueParams;
use sohar_rust::simplicial::from::graph_weighted::user_interface::cliques_in_order;
use sohar_rust::simplicial::from::relation::UseClearing;
use sohar_rust::simplicial::simplices::Simplex;
use sohar_rust::simplicial::simplices::SimplexFiltered;
use oat_rust::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
use oat_rust::matrices::display::{print_indexed_minor_views, print_indexed_major_views};
use std::io;
use std::time::Instant;

use num::rational::Ratio;
use itertools::Itertools;
use oat_rust;
use oat_rust::matrices::matrix_oracle_traits::OracleMajorAscend;
use oat_rust::matrices::matrix_oracle_traits::OracleMinorDescend;
use oat_rust::matrices::operations::umatch::row_major::new_umatchrowmajor_with_clearing;
use oat_rust::rings::operator_structs::ring_native::field_rational_i64;
// use oat_rust::utilities::cw_complexes::simplices_weighted::clique::CliqueBoundaryMatrix;
// use oat_rust::utilities::homology::clique::CliqueParams;
// use oat_rust::utilities::homology::clique::cliques_in_order;
// use oat_rust::utilities::homology::dowker::UseClearing;
// use oat_rust::utilities::cw_complexes::simplices_weighted::clique::Simplex;
// use oat_rust::utilities::cw_complexes::simplices_weighted::clique::SimplexIter;
use oat_rust::utilities::partial_order::OrderComparatorAutoGt;
use oat_rust::utilities::partial_order::OrderComparatorAutoLt;
use oat_rust::utilities::partial_order::OrderComparatorLtByKey;
use oat_rust::utilities::partial_order::is_sorted_strictly;

// use numpy::ndarray::{Array1, ArrayD, ArrayView1, ArrayViewD, ArrayViewMutD, Zip};
// use numpy::{
//     datetime::{units, Timedelta},
//     Complex64, IntoPyArray, PyArray1, PyArrayDyn, PyReadonlyArray1, PyReadonlyArrayDyn,
//     PyReadwriteArray1, PyReadwriteArrayDyn,
// };


// #[pyfunction]
// fn axpy_py<'py>(
//     py: Python<'py>,
//     x: PyReadonlyArrayDyn<'_, f64>,
// ) -> f64 {
//     let x = x.as_array();
//     // let y = y.as_array();
//     // let z = axpy(a, x, y);
//     // z.into_pyarray(py)
//     x[[0,0]]
// }



/// Compute basis of cycle representatives, over the 2-element field.
/// 
/// Input is a relation formatted vec-of-rowvec matrix.
#[pyfunction]
fn homology_basis_from_dowker( 
            dowker_simplices_vec_format: Vec<Vec<usize>>, 
            maxdim: usize
        ) 
    -> PyResult<Vec<Vec<Vec<(Vec<usize>,(i64,i64))>>>> {
    // precompute the number of columns of the untransposed matrix
    // note we have to add 1, since arrays are 0-indexed and we
    // want to be able to use the max value as an index
    let basis = 
        sohar_rust::simplicial::from::relation::homology_basis_from_dowker(
            & dowker_simplices_vec_format, 
            maxdim,
            UseClearing::Yes,
        );
    let convert = |x: Ratio<i64> | (x.numer().clone(), x.denom().clone());
    let basis_new = basis.iter().map(
                |x|
                x.iter().map(
                    |y|
                    y.iter().map(
                        |z|
                        (z.0.clone(), convert(z.1))
                    ).collect_vec()
                ).collect_vec()
            ).collect_vec();

    Ok(basis_new)
}


/// Compute basis of cycle representatives, over the rationals
/// 
/// Input is a relation formatted as a vec-of-rowvec matrix.
#[pyfunction]
fn persistent_homology_basis_from_clique( 
            dissimilarity_matrix: Vec<Vec<f64>>, 
            maxdim: usize,
            maxdis: Option<f64>,  
        ) 
    -> PyResult<
            (
                Vec<Vec<(f64, f64)>>,
                Vec<Vec<Vec<(Vec<usize>,(i64,i64))>>>,
                Vec<Vec<Vec<(Vec<usize>,(i64,i64))>>>,
            )>
    {

    // let maxdis = match use_enclosing_radius{
    //     true => Some(maxdis),
    //     false => None,
    // };
    for i in 0 .. dissimilarity_matrix.len() {
        for j in i .. dissimilarity_matrix.len() {
            assert_eq!( dissimilarity_matrix[i][j], dissimilarity_matrix[j][i] );
        }
    }
 

    let ring_operator = field_rational_i64();
    let params = CliqueParams::new(
                & dissimilarity_matrix, 
                maxdim, 
                maxdis, 
                ring_operator.clone(),
            );
    println!("CliqueParams: {:?}", &params.maxdis() );
    let boundary_matrix = CliqueBoundaryMatrix::new( & params );
    let boundary_matrix_ref = & boundary_matrix;    
    let keymaj_vec = cliques_in_order( & params );
    // let keymin_vec = cliques_in_order( & CliqueParams::new(
    //             dissimilarity_matrix, maxdim+1, maxdis, ring_operator.clone()
    //         ) ) ;


    // A GOOD CHECK BUT IT CAN TAKE TIME
    // verify_viewmajorascend_compatible_with_viewminordescend(
    //         boundary_matrix_ref,
    //         keymin_vec.iter().cloned(),
    //         keymaj_vec.iter().cloned(),
    //     );       

    // println!("(oat_rust_python) keymaj_vec = "); 
    // for entry in keymaj_vec.iter() { println!("{:?}", entry) }


    // println!("press any key to continue");
    // let mut guess = String::new();
    // io::stdin()
    //     .read_line(&mut guess)
    //     .expect("Failed to read line");    

    let iter_keymaj = keymaj_vec.iter().cloned();    

    // A GOOD CHECK BUT IT CAN TAKE TIME
    // println!("check that oracle has strictly sorted rows");
    // // print_indexed_major_views( & boundary_matrix_ref, iter_keymaj.clone() );  // print the major views       
    // for keymaj in iter_keymaj.clone() {        
    //     assert!( is_sorted_strictly( 
    //                                     &boundary_matrix_ref.view_major_ascend(keymaj.clone()).collect_vec() , 
    //                                     &OrderComparatorLtByKey::new( OrderComparatorAutoLt::new() ) 
    //                                 ) );
    // }

    // println!("press enter to continue");
    // let mut guess = String::new();
    // io::stdin()
    //     .read_line(&mut guess)
    //     .expect("Failed to read line");        

    // A GOOD CHECK BUT IT CAN TAKE TIME
    // println!("check that oracle has strictly sorted columns");
    // // print_indexed_minor_views( & boundary_matrix_ref, iter_keymaj.clone() );  // print the major views        
    // for keymaj in iter_keymaj.clone() {
    //     assert!( is_sorted_strictly(    &boundary_matrix_ref.view_minor_descend(keymaj).iter().cloned().collect_vec() , 
    //                                     &OrderComparatorLtByKey::new( OrderComparatorAutoGt::new() )  // NOTE THAT HERE WE USE GT 
    //                                 ) );
    // }    

    // println!("press enter to continue");
    // let mut guess = String::new();
    // io::stdin()
    //     .read_line(&mut guess)
    //     .expect("Failed to read line");       

    // println!("starting umatch");
    let umatch = new_umatchrowmajor_with_clearing(
            boundary_matrix_ref, 
            iter_keymaj.clone(), 
            ring_operator.clone(), 
            OrderComparatorAutoLt::new(), 
            OrderComparatorAutoLt::new(), 
        );      

    // println!("start build bd matrix");
    // let boundary_matrix = oat_rust::utilities::homology::clique::get_clique_boundary_matrix(
    //     dissimilarity_matrix,
    //     maxdis,
    //     field_rational_i64()        
    // );

    // println!("start umatch");    
    // let umatch = oat_rust::utilities::homology::clique::umatch_from_clique(
    //     & boundary_matrix,
    //     maxdim,
    // );

       
    // let iter_keymaj = 
    //     (0..maxdim+1).map(
    //         |dim|
    //         {
    //             let mut vec = SimplexIter::new( 
    //                         dim,
    //                         & umatch.array_mapping_ref().dismat,
    //                         umatch.array_mapping_ref().maxdis,                
    //                     )
    //                     .collect_vec();
    //             vec.sort_by(|x,y| y.filvalue.cmp(&x.filvalue));
    //             vec
    //         }
    //     )
    //     .flatten();    


    // println!("setting up to unpack");  
    let dim_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.dim();
    let fil_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.fil().into_inner();    
    let barcode = umatch.barcode( iter_keymaj, dim_fn, fil_fn, true , true);

    // println!("getting intervals, reps, bounding chains");
    let (intervals, representatives, bounding_chains ) = barcode.unwrap();

    let convert_coefficients = |x: Ratio<i64> | (x.numer().clone(), x.denom().clone());
    let convert_simplex = |x: SimplexFiltered< OrderedFloat<f64> > | x.vertices.iter().map(|y| y.clone() as usize ).collect_vec();
    
    // println!("start reformat reps");   

    let now = Instant::now();        
    let represntatives_new = representatives.unwrap().iter().map(
                |x| // each x is a list of cycle reps
                x.iter().map(
                    |y| // each y is a cycle rep
                    y.iter().map(
                        |z| // each z is an entry in a cycle rep
                        ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
                    ).collect_vec()
                ).collect_vec()
            ).collect_vec();
    println!("time compute reps {}", now.elapsed().as_millis());            

    // println!("start reformat bounding chains");                 
    let now = Instant::now();            
    let bounding_chains_new = bounding_chains.unwrap().iter().map(
                |x| // each x is a list of cycle reps
                x.iter().map(
                    |y| // each y is a cycle rep
                    y.iter().map(
                        |z| // each z is an entry in a cycle rep
                        ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
                    ).collect_vec()
                ).collect_vec()
            ).collect_vec();           
    println!("time compute bounding chains {}", now.elapsed().as_millis());                

    Ok( ( intervals, represntatives_new, bounding_chains_new ) )
}




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
fn sohar_python(_py: Python, m: &PyModule) -> PyResult<()> {
    pyo3_log::init();
    m.add_function(wrap_pyfunction!(persistent_homology_basis_from_clique, m)?)?;        
    m.add_function(wrap_pyfunction!(homology_basis_from_dowker, m)?)?;    
    m.add_function(wrap_pyfunction!(unique_row_indices, m)?)?;
    m.add_function(wrap_pyfunction!(unique_rows, m)?)?;
    m.add_function(wrap_pyfunction!(transpose_listlist, m)?)?;

    Ok(())
}
