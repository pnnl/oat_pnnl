//! Internal vectors are stored in unsorted order
//! 
//! See the [parent module](super) for details on what vector-of-vectors format means.
//! 
//! # Warning: Under construction
//! 
//! This module is under construction following some refactoring changes.
//! If you think you'll need the functionality that this module is intended to provide,
//! there are two options: 
//! 
//! (1) Check out the source code.  You might be able to implement the missing functionality yourself!
//!     If you do, we'd love for you to share your code with the community.  Reach out to us for help
//!     doing that!   **Some of the source code has been commented out, so check out the source file
//!     for bits and pieces that might be useful.**
//! 
//! (2) Write to us and let us know. Knowing what people care about helps us prioritize the
//!     code we work on!

use crate::algebra::vectors::entries::{KeyValGet, KeyValTypes};
use crate::algebra::matrices::query::{   ViewRow,
                                        ViewRowAscend,
                                        ViewRowDescend, ViewColDescend, IndicesAndCoefficients, ViewColAscend, ViewCol, MatrixEntry, MatrixOracle,
                                        // ViewCol, 
                                        // ViewColAscend,
                                        // ViewColDescend,
                                        // WhichMajor,
                                        // MajorDimension
                                    };

                                    use crate::utilities::binary_search::{find_sorted_binary_oracle};
                                    use crate::utilities::order::{JudgePartialOrder, is_sorted_strictly, OrderOperatorByKey, };
                                    use crate::utilities::statistics::histogram;


use rand::Rng;                                          // we use this module to generate random elements
use rand::distributions::{Bernoulli, Distribution};     // we use this module to generate random elements

use std::iter::{Rev, Cloned};
use std::marker::PhantomData;
use std::slice::Iter;
use itertools::Itertools;

use super::super::bimajor::{MatrixBimajor, MatrixBimajorData};

//  Unsorted vector of vectors
//  -------------------------------------------------------------------------------------------------------

/// An vector of unsorted vectors, representing a sparse matrix.  
/// 
/// See the [parent module](super) for details on what vector-of-vectors format means.
/// 
/// Each of the internal vectors should have entries sorted in asecneding order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::unsorted::VecOfVecUnsorted;
/// use oat_rust::algebra::matrices::query::*;
/// 
/// // Create a new vec-of-vec matrix.
/// let matrix  =   VecOfVecUnsorted::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                 );
/// ```
pub struct VecOfVecUnsorted
                < IndexCoeffPair >
    where
        IndexCoeffPair:             KeyValTypes,
    // where   IndexCoeffPair:     KeyValGet < ColIndex, Coefficient >,
{
    pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
    // phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
    // phantom_minkey: PhantomData< ColIndex >,
    // phantom_snzval: PhantomData< Coefficient >,
}


impl    < 'a, IndexCoeffPair, >
        
        VecOfVecUnsorted 
            < IndexCoeffPair, >

    where
        IndexCoeffPair:             KeyValTypes,            
        
        // where   IndexCoeffPair:    KeyValGet < ColIndex, Coefficient >        

{
    // Make new `VecOfVecUnsorted` from an existing vector of vectors. 
    pub fn new( vecvec: Vec < Vec < IndexCoeffPair > > ) -> Self
    {        
        VecOfVecUnsorted{   
                    vec_of_vec:     vecvec,                    
                    // phantom_kvpair: PhantomData,
                    // phantom_minkey: PhantomData,
                    // phantom_snzval: PhantomData,                    
                }
    }
}

// // IndicesAndCoefficients
// impl < 'a, IndexCoeffPair, >

//     IndicesAndCoefficients for
//     &'a VecOfVecUnsorted< IndexCoeffPair, >

//     where
//         IndexCoeffPair:             KeyValTypes,    

// {   
//     type ColIndex = IndexCoeffPair::Key; 
//     type RowIndex = usize; 
//     type Coefficient = IndexCoeffPair::Val; 
//     type EntryMajor = &'a IndexCoeffPair;
//     type EntryMinor = ( usize, IndexCoeffPair::Val );
// }    


// // ViewRow
// impl < 'a, IndexCoeffPair, >
    
//     ViewRow  for 
    
//     &'a VecOfVecUnsorted 
//         < IndexCoeffPair, >

//     where
//         IndexCoeffPair:             KeyValTypes,   

//     // where   IndexCoeffPair:     KeyValGet < ColIndex, Coefficient >,
// {       
//     type ViewMajor          =   &'a[IndexCoeffPair];
//     type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;

//     fn view_major( & self, index: usize ) -> &'a[IndexCoeffPair] {
//         return self.vec_of_vec[index].as_slice() //.cloned()
//     } 
// }




























#[cfg(test)]
mod tests {    

    use crate::algebra::matrices::display::print_indexed_major_views;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_vec_of_vec_construction() {

        let _matrix  =   VecOfVecUnsorted::new(
                                                vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ]
                                                    );                                      
    }
}
