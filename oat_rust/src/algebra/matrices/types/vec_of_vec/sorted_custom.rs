//!  Internal vectors are stored in a *customized* sorted order
//! 
//! See the [parent module](super) for details on what vector-of-vectors format means.
//! 
//! 
//! # Warning: Under construction
//! 
//! This module is under construction following some refactoring changes.
//! If you think you'll need the functionality that this module is intended to provide,
//! there are two options: 
//! 
//! (1) Check out the source code.  You might be able to implement the missing functionality yourself!
//!     If you do, we'd love for you to share your code with the community.  Reach out to us for help
//!     doing that!  **Some of the source code has been commented out, so check out the source file
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

                                    use crate::utilities::order::{JudgePartialOrder, is_sorted_strictly, OrderOperatorByKey, };

use std::iter::{Rev, Cloned};




/// An vector of strictly sorted vectors, representing a sparse matrix.  
/// 
/// By "strictly sorted" we mean that (i) the entries in each vector appear in ascending order,
/// according to index, and (ii) within each vector, consequtive entries have distinct indices.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted_custom::*;
/// use oat_rust::algebra::matrices::query::*;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// 
/// // Standard method to create a row-major vec-of-vec matrix (second function 
/// // argument specifies the order on entries).
/// let matrix  =   VecOfVec::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                     OrderOperatorByKey::new(),
///                 );
/// 
/// // Streamlined method to create a row-major vec-of-vec matrix (order of entires 
/// // is automatically inferred from the partial order on indices).
/// let matrix  =   vecvec_with_defualt_order(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                 );
/// 
/// ```
#[derive(Debug, Clone)]
pub struct VecOfVec
                < IndexCoeffPair, OrderOperator, >

    // < ColIndex, Coefficient, IndexCoeffPair, OrderOperator, >

    // where   IndexCoeffPair:     KeyValGet < ColIndex, Coefficient, >,
    //         OrderOperator:    JudgePartialOrder<  IndexCoeffPair >,


    where   OrderOperator:    JudgePartialOrder<  IndexCoeffPair >,   
            IndexCoeffPair:     KeyValTypes,         
{
    vec_of_vec:         Vec< Vec< IndexCoeffPair > >,
    order_operator:   OrderOperator,
    // pub phantom_minkey: PhantomData< ColIndex >,
    // pub phantom_snzval: PhantomData< Coefficient >,
}


impl    < IndexCoeffPair, OrderOperator, >
        
        VecOfVec 
            < IndexCoeffPair, OrderOperator, >
        
        where   // IndexCoeffPair:     KeyValGet < ColIndex, Coefficient >,
                OrderOperator:    JudgePartialOrder<  IndexCoeffPair >,
                IndexCoeffPair:     KeyValTypes,
                // ColIndex:             ,
                // Coefficient:             ,
                // Self:               'a,        

{
    /// Make a new `VecOfVec`.  Throws an error if one of the internal vectors is not sorted
    /// in strictly asecending order, by index.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted_custom::VecOfVec;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // This line should not panic.
    /// let matrix = VecOfVec::new( vec![ vec![ (0,5), (1,6) ] ], OrderOperatorAuto::new() ); // no panic
    /// 
    /// // This line should panic.
    /// // let matrix = VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ], OrderOperatorAuto::new() );
    /// // And here is the test that confirms it.
    /// let result = std::panic::catch_unwind(|| {VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ], OrderOperatorAuto::new() )} );
    /// assert!(result.is_err()); 
    /// ```
    pub fn new( vecvec: Vec < Vec < IndexCoeffPair > >, mut order_operator: OrderOperator ) -> Self
    {
        for vec in vecvec.iter() {
            if ! is_sorted_strictly( vec, &mut order_operator ) {
                panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
            }
        }
        
        VecOfVec{   
                    vec_of_vec:         vecvec,                    
                    order_operator,
                    // phantom_minkey:     PhantomData,
                    // phantom_snzval:     PhantomData,                    
                }
    }

    /// Return the internally stored `Vec<Vec<EntryType>>` struct and the explicit order comparator.  Consumes the wrapper `VecOfVec` struct.
    pub fn decompose( self ) -> ( Vec < Vec < IndexCoeffPair > >, OrderOperator ) {
        ( self.vec_of_vec, self.order_operator, )
    }

    /// Append a vector to the internally stored `Vec<Vec<EntryType>>` struct.  The new vector must be sorted in strictly ascending order; the function panics, otherwise.
    pub fn push_vec( &mut self, new_sorted_vec: Vec < IndexCoeffPair > ) {
        if ! is_sorted_strictly( & new_sorted_vec, &mut self.order_operator ) {
            panic!("Attempt to append a non-strictly-sorted vector to `VecOfVec`.  To proceed, the new vector must be sorted in *strictly* ascending order, by index.") 
        }        
        self.vec_of_vec.push( new_sorted_vec );
    }

    /// Returns an immutable reference to the `vec_number`th internal vector.
    pub fn peek_at( & self, vec_number: usize ) -> &Vec< IndexCoeffPair > {        
        & self.vec_of_vec[ vec_number ]
    }

}

// THIS IS NOT AN ASSOCIATED FUNCTION!
pub fn vecvec_with_defualt_order 
        < IndexCoeffPair, >
        ( vecvec: Vec < Vec < IndexCoeffPair > > )
        ->
        VecOfVec 
            < IndexCoeffPair, OrderOperatorByKey< IndexCoeffPair::Key, IndexCoeffPair::Val, IndexCoeffPair >, >

        where   IndexCoeffPair:         KeyValGet < IndexCoeffPair::Key, IndexCoeffPair::Val >,
                IndexCoeffPair:         KeyValTypes,
                IndexCoeffPair::Key:    PartialOrd,
    
{

    for vec in vecvec.iter() {
        if ! is_sorted_strictly( vec, &mut OrderOperatorByKey::new() ) {
            panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
        }
    }

    VecOfVec::new(   
                vecvec,                    
                OrderOperatorByKey::new() 
        )
}

impl < 'a, IndexCoeffPair, OrderOperator, >

    IndicesAndCoefficients for 
    &'a VecOfVec< IndexCoeffPair, OrderOperator, >

    where   
        OrderOperator:    JudgePartialOrder<  IndexCoeffPair >,
        IndexCoeffPair:     KeyValTypes,             

{ 
    type EntryMajor = &'a IndexCoeffPair;
    type EntryMinor = (usize, IndexCoeffPair::Val);       
    type ColIndex = IndexCoeffPair::Key; 
    type RowIndex = usize; 
    type Coefficient = IndexCoeffPair::Val; 
}  



// impl < 'a, IndexCoeffPair, OrderOperator, >
    
//     ViewRow  for 
    
//     &'a VecOfVec
//         < IndexCoeffPair, OrderOperator, >

//     where   //IndexCoeffPair:     KeyValGet < ColIndex, Coefficient >,
//             OrderOperator:    JudgePartialOrder<  IndexCoeffPair >,    
//             IndexCoeffPair:     KeyValTypes,        
// {
//     type ViewMajor          =   &'a[IndexCoeffPair];
//     type ViewMajorIntoIter  =   std::slice::Iter< 'a, IndexCoeffPair >;

//     fn view_major( & self, index: usize ) -> &'a[IndexCoeffPair] {
//         return self.vec_of_vec[index].as_slice()
//     } 
// }


// impl < 'a, IndexCoeffPair, OrderOperator, >
    
//     ViewRowAscend  for 
    
//     &'a VecOfVec 
//         < IndexCoeffPair, OrderOperator, >

//     where   // IndexCoeffPair:     KeyValGet < ColIndex, Coefficient >,
//             OrderOperator:    JudgePartialOrder<  IndexCoeffPair >, 
//             IndexCoeffPair:     KeyValTypes,                 
// {
//     type ViewMajorAscend            =   &'a[IndexCoeffPair];
//     type ViewMajorAscendIntoIter    =   std::slice::Iter< 'a, IndexCoeffPair >;
        
//     /// Assumes that entries in each vector are sorted in ascending order.
//     // fn view_major_ascend( & self, index: usize ) -> & Vec< IndexCoeffPair > {
//     //     return self.view_major( index )
//     // } 
//     fn view_major_ascend( & self, index: usize ) -> &'a[IndexCoeffPair] {
//         self.view_major( index )
//     }     
// }

// impl < 'a, IndexCoeffPair, OrderOperator, > 
    
//     ViewRowDescend  for 
    
//     &'a VecOfVec 
//         < IndexCoeffPair, OrderOperator, >

//     where   //IndexCoeffPair:     KeyValGet < ColIndex, Coefficient >,
//             OrderOperator:    JudgePartialOrder<  IndexCoeffPair >,
//             IndexCoeffPair:     KeyValTypes,              
// {
//     type ViewMajorDescend           =   Rev<std::slice::Iter<'a, IndexCoeffPair>>;
//     type ViewMajorDescendIntoIter   =   < Self::ViewMajorDescend as IntoIterator >::IntoIter;
        
//     /// Assumes that entries in each vector are sorted in ascending order.    
//     fn view_major_descend( & self, index: usize ) -> Rev<std::slice::Iter<'a, IndexCoeffPair >> {
//         return self.vec_of_vec[index].iter().rev()
//     } 
// }
















#[cfg(test)]
mod tests {    

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

}
