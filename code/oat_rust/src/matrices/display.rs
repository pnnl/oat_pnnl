//! Display data from a matrix oracle.

use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, IndicesAndCoefficients};
use itertools::Itertools;
use std::fmt::Debug;

use super::matrix_oracle_traits::OracleMinorDescend;



/// Print one major view of the matrix (represented as a vector of entries) for each item in the iterator.
/// 
/// ```
/// // import the relevant crates
/// use oat_rust::matrices::matrix_types::vec_of_vec::{VecOfVecSimple}; // a particular type of sparse matrix
/// use oat_rust::matrices::matrix_oracle_traits::OracleMajorAscend; // the trait that defines the command for ascending major views
/// use oat_rust::matrices::display::print_indexed_major_views;
///         
/// // define a row-major sparse matrix that represents the follow matrix
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVecSimple::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
/// 
/// // define an iterator that runs over the major keys (i.e. the row indices)
/// let iter_keymaj = 0..2;  
/// 
/// // print the major views specified by the iterator.  this should show the following:
/// // $ major_view 0: [(1, 5), (2, 5)]
/// // $ major_view 1: [(2, 7)]
/// print_indexed_major_views( &&matrix, iter_keymaj ); // we pass `&&matrix` because `&matrix` implements the oracle trait and we have to pass a *reference* to something that implements the oracle trait
/// ```
pub fn print_indexed_major_views
            < Matrix, IterKeyMaj >
            ( matrix: & Matrix, iter_keymaj: IterKeyMaj ) 
    where   
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMaj:                 Clone + Debug,        
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   Debug,
        IterKeyMaj:                     IntoIterator< Item = Matrix::KeyMaj >,        
{
    for keymaj in iter_keymaj {        
        println!("major_view for key {:?}", keymaj.clone());
        for entry in matrix.view_major_ascend( keymaj ).into_iter()  { 
            println!("{:?}", entry) 
        }
    }
}

/// Print one major view of the matrix (represented as a vector of entries) for each item in the iterator.
pub fn print_indexed_minor_views
            < Matrix, IterKeyMin >
            ( matrix: & Matrix, iter_keymin: IterKeyMin ) 
    where   
        Matrix:                         OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMin:                 Clone + Debug,        
        Matrix::ViewMinorDescend:       IntoIterator,
        Matrix::ViewMinorDescendEntry:  Debug,
        IterKeyMin:                     IntoIterator< Item = Matrix::KeyMin >,        
{
    for keymin in iter_keymin {
        println!("minor_view {:?}: {:?}", keymin.clone(), matrix.view_minor_descend( keymin ).into_iter().collect_vec()  );
    }
}









#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_print_indexed_major_views() {
        use crate::matrices::matrix_types::vec_of_vec::{VecOfVecSimple}; // a particular type of sparse matrix
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend; // the trait that defines the command for ascending major views
        use crate::matrices::display::print_indexed_major_views;

        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVecSimple::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
        let iter_keymaj = 0..2;  // iterates over the major keys, i.e. over the row indices

        // print the major views specified by the iterator
        print_indexed_major_views( &&matrix, iter_keymaj ); // we pass `&&matrix` because `&matrix` implements the oracle trait and we have to pass a *reference* to something that implements the oracle trait
        // this should show the following:
    }


}
