//! Vec-of-Vec sparse matrices.
//! 
//! A "vec-of-vec" is a vector of vectors representing the nonzero entries of a sparse matrix.
//! For example, if `A` is the following matrix
//! 
//! ```ignore
//! 5 6
//! 7 0
//! ```
//! 
//! Then the row-major vec-of-vec representation of `A` would be 
//! 
//! ```ignore
//! vec![
//!     vec![ (0, 5), (1, 6) ],
//!     vec![ (0, 7) ],
//! ]
//! ```
//! 
//! and the column-major vec-of-vec representation of `A` would be
//! 
//! ```ignore
//! vec![
//!     vec![ (0, 5), (1, 7) ],
//!     vec![ (0, 6) ],
//! ]
//! ```
//! 
//! This module contains several variants of the vec-of-vec data structure for use
//! with the [matrix oracle traits](crate::matrices::matrix_oracle_traits).
//! 
//! # Design elements
//! 
//! There are three separate vec-of-vec data structures, each motivated by a different use.
//! 
//! 
//! 
//! - [VecOfVec](VecOfVec)
//!   - Raison D'Etre
//!     - *Lifetimes are required by the matrix oracle traits.*  
//!       To the best of their knowledge (please let them know if information is found to the contrary), the developers believe that a `Vec<Vec<EntryType>>` object does not have an explicit
//!       lifetime parameter.  Wrapping a `Vec<Vec<EntryType>>` in a `VecOfVec<..>`  struct
//!       provides a lifetime parameter which can be used in implementing the [matrix oracle traits](crate::matrices::matrix_oracle_traits).
//!     - *Strict order is required for safe, efficient implementation of ascending and descending oracles.*  
//!                 To facilitate the implementation of the [OracleMajorAscend](crate::matrices::matrix_oracle_traits::OracleMajorAscend) and
//!                  [OracleMajorDescend](crate::matrices::matrix_oracle_traits::OracleMajorDescend)
//!                  traits in an efficient manner, entries in the inner vectors should appear in
//!                  in strictly sorted order according to index 
//!                  (either ascending or descending; arbitrarily, we chose ascending).
//!                  Thus, if we want a vec-of-vec data structure that can efficiently implement matrix oracle
//!                  traits, then we need some way to ensure that the inner vectors are sorted, for *safety*.
//!                  The `VecOfVec` accomplishes this by restricting user access to the `Vec<Vec<EntryType>>` which it stores internally.
//!                  One can only generate new instances of the `VecOfVec` struct by 
//!                  using the constructor `new` function, and this function panics if it receives a vec-of-vec
//!                  where one of the inner vectors is not sorted.
//!     - *Explicit order comparators make order of entries unambiguous, and afford greater flexibility.*  
//!        Since we allow matrix indices to have essentially any type, it is necessary to be explit about 
//!                  how order is defined.  This is the reason why `VecOfVec` objects store an internal 
//!                  `OrderComparator` object, which determines when one index is strictly less than another.
//!                  Typically, an order comapartor requires zero memory.  Moreover, it's necessary to provide
//!                  an order comaparator when constructig a `VecOfVec` (so that the `new` function can
//!                  check that entries are sorted), so requiring an order comparator as an argument in the 
//!                  `new` constructor adds no additional burden to the user.
//!
//!   - *Methods and challenges of editing data stored in a `VecOfVec`.*  
//!     - Due to the requirements of object
//!         safety, the user is highly limited in their ability to modify that data stored internally by a `VecOfVec`.
//!         The most general route is to retrieve the internally stored `Vec<Vec<EntryType>>`, modify it, then 
//!         wrap it in a new `VecOfVec` struct.  Note, however, that this will incur the cost of (i) re-sorting
//!         each internal vector, or (ii) verifying that each internal vector is sorted, if it is already sorted.
//!   - *Alternatives*
//!     - If the restrictions imposed by safety requirements on the `VecOfVec` struct are overly onerous, 
//!       consider using [VecOfVecUnsorted](VecOfVecUnsorted).  This does not implement the [OracleMajorAscend](crate::matrices::matrix_oracle_traits::OracleMajorAscend) or
//!       [OracleMajorDescend](crate::matrices::matrix_oracle_traits::OracleMajorDescend) methods, but it is much easier to modifiy.
//! 
//! - [VecOfVecUnsorted](VecOfVecUnsorted)
//!   - Raison D'Etre
//!     - The requirement that internal vectors be strictly sorted may place an undue borden on the
//!                  user in some cases, for example
//!       - when one wishes to make frequent updates to the matrix, without re-sorting and re-checking each internal vector each time
//!       - when defining an order comparator is difficult or onerous
//!     Nevertheless, one may still wish to implement `OracleMajor` on 
//!                  Unlike a struct of type `Vec< Vec< EntryType > >`, a struct of type [VecOfVecUnsorted](VecOfVecUnsorted) has a lifetime
//!                  parameter.  We use that parameter in the implementation of the 
//!                  [matrix oracle traits](crate::matrices::matrix_oracle_traits).
//!     - **Issue** Could the same be achieved with a reference, `&'a Vec<Vec<EntryType>>`?
//!                  
//! - [VecOfVecSimple](VecOfVecSimple) 
//!   - Raison D'Etre  
//!     
//!     - The Vec-of-vec format is used extensively throughout the documentation and unit tests in this library
//!       because it is one of the most human readable sparse vector formats available.  However, this struct
//!       has a number of type parameters; more type parameters make code less readable and examples/tests harder to write.
//!                   The [VecOfVecSimple](VecOfVecSimple) struct has fewer type parameters; it is easer to read, construct,
//!                   and analyze.  (The price you pay for this simplicity is flexibility/generality, but in unit tests
//!         this matters less).
//! 
//! 
//! 

use serde_json::value::Index;

use crate::matrices::matrix_oracle_traits::{   OracleMajor,
                                        OracleMajorAscend,
                                        OracleMajorDescend, OracleMinorDescend, IndicesAndCoefficients, OracleMinorAscend, OracleMinor,
                                        // OracleMinor, 
                                        // OracleMinorAscend,
                                        // OracleMinorDescend,
                                        // WhichMajor,
                                        // MajorDimension
                                    };
use crate::entries::{KeyValGet, KeyValAssociatedTypes};
use crate::utilities::partial_order::{StrictlyLess, is_sorted_strictly, OrderComparatorAutoLtByKey, OrderComparatorAutoLtByFirstTupleEntry};
use std::iter::{Rev, Cloned};
use std::marker::PhantomData;
use std::slice::Iter;
use itertools::Itertools;




//  Unsorted vector of vectors
//  -------------------------------------------------------------------------------------------------------

/// An vector of unsorted vectors, representing a sparse matrix.  
/// 
/// Each of the internal vectors should have entries sorted in asecneding order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::*;
/// use oat_rust::matrices::matrix_oracle_traits::*;
/// 
/// // Create a new vec-of-vec matrix.
/// let matrix  =   VecOfVecUnsorted::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                 );
/// ```
pub struct VecOfVecUnsorted
                < IndexCoeffPair >
    where
        IndexCoeffPair:             KeyValAssociatedTypes,
    // where   IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
{
    pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
    // phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
    // phantom_minkey: PhantomData< KeyMin >,
    // phantom_snzval: PhantomData< SnzVal >,
}


impl    < 'a, IndexCoeffPair, >
        
        VecOfVecUnsorted 
            < IndexCoeffPair, >

    where
        IndexCoeffPair:             KeyValAssociatedTypes,            
        
        // where   IndexCoeffPair:    KeyValGet < KeyMin, SnzVal >        

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

// IndicesAndCoefficients
impl < 'a, IndexCoeffPair, >

    IndicesAndCoefficients for
    &'a VecOfVecUnsorted< IndexCoeffPair, >

    where
        IndexCoeffPair:             KeyValAssociatedTypes,    

{   type KeyMin = IndexCoeffPair::Key; type KeyMaj = usize; type SnzVal = IndexCoeffPair::Val;  }    


// OracleMajor
impl < 'a, IndexCoeffPair, >
    
    OracleMajor  for 
    
    &'a VecOfVecUnsorted 
        < IndexCoeffPair, >

    where
        IndexCoeffPair:             KeyValAssociatedTypes,   

    // where   IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
{       
    type ViewMajor          =   &'a[IndexCoeffPair];
    type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajorEntry     =   &'a IndexCoeffPair;

    fn view_major( & self, index: usize ) -> &'a[IndexCoeffPair] {
        return self.vec_of_vec[index].as_slice() //.cloned()
    } 
}




//  Vector of vectors (STRICTLY SORTED)
//  -------------------------------------------------------------------------------------------------------




/// An vector of strictly sorted vectors, representing a sparse matrix.  
/// 
/// By "strictly sorted" we mean that (i) the entries in each vector appear in ascending order,
/// according to index, and (ii) within each vector, consequtive entries have distinct indices.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::*;
/// use oat_rust::matrices::matrix_oracle_traits::*;
/// use oat_rust::utilities::partial_order::OrderComparatorAutoLtByKey;
/// 
/// // Standard method to create a row-major vec-of-vec matrix (second function 
/// // argument specifies the order on entries).
/// let matrix  =   VecOfVec::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                     OrderComparatorAutoLtByKey::new(),
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
pub struct  VecOfVec
                < IndexCoeffPair, OrderComparator, >

    // < KeyMin, SnzVal, IndexCoeffPair, OrderComparator, >

    // where   IndexCoeffPair:     KeyValGet < KeyMin, SnzVal, >,
    //         OrderComparator:    StrictlyLess<  IndexCoeffPair >,


    where   OrderComparator:    StrictlyLess<  IndexCoeffPair >,   
            IndexCoeffPair:     KeyValAssociatedTypes,         
{
    vec_of_vec:         Vec< Vec< IndexCoeffPair > >,
    order_comparator:   OrderComparator,
    // pub phantom_minkey: PhantomData< KeyMin >,
    // pub phantom_snzval: PhantomData< SnzVal >,
}


impl    < IndexCoeffPair, OrderComparator, >
        
        VecOfVec 
            < IndexCoeffPair, OrderComparator, >
        
        where   // IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
                OrderComparator:    StrictlyLess<  IndexCoeffPair >,
                IndexCoeffPair:     KeyValAssociatedTypes,
                // KeyMin:             ,
                // SnzVal:             ,
                // Self:               'a,        

{
    /// Make a new `VecOfVec`.  Throws an error if one of the internal vectors is not sorted
    /// in strictly asecending order, by index.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVec;
    /// use oat_rust::utilities::partial_order::OrderComparatorAutoAnyType;
    /// 
    /// // This line should not panic.
    /// let matrix = VecOfVec::new( vec![ vec![ (0,5), (1,6) ] ], OrderComparatorAutoAnyType::new() ); // no panic
    /// 
    /// // This line should panic.
    /// // let matrix = VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ], OrderComparatorAutoAnyType::new() );
    /// // And here is the test that confirms it.
    /// let result = std::panic::catch_unwind(|| {VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ], OrderComparatorAutoAnyType::new() )} );
    /// assert!(result.is_err()); 
    /// ```
    pub fn new( vecvec: Vec < Vec < IndexCoeffPair > >, mut order_comparator: OrderComparator ) -> Self
    {
        for vec in vecvec.iter() {
            if ! is_sorted_strictly( & vec, &mut order_comparator ) {
                panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
            }
        }
        
        VecOfVec{   
                    vec_of_vec:         vecvec,                    
                    order_comparator:   order_comparator,
                    // phantom_minkey:     PhantomData,
                    // phantom_snzval:     PhantomData,                    
                }
    }

    /// Return the internally stored `Vec<Vec<EntryType>>` struct and the explicit order comparator.  Consumes the wrapper `VecOfVec` struct.
    pub fn decompose( self ) -> ( Vec < Vec < IndexCoeffPair > >, OrderComparator ) {
        ( self.vec_of_vec, self.order_comparator, )
    }

    /// Append a vector to the internally stored `Vec<Vec<EntryType>>` struct.  The new vector must be sorted in strictly ascending order; the function panics, otherwise.
    pub fn push_vec( &mut self, new_sorted_vec: Vec < IndexCoeffPair > ) {
        if ! is_sorted_strictly( & new_sorted_vec, &mut self.order_comparator ) {
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
            < IndexCoeffPair, OrderComparatorAutoLtByKey< IndexCoeffPair::Key, IndexCoeffPair::Val, IndexCoeffPair >, >

        where   IndexCoeffPair:         KeyValGet < IndexCoeffPair::Key, IndexCoeffPair::Val >,
                IndexCoeffPair:         KeyValAssociatedTypes,
                IndexCoeffPair::Key:    PartialOrd,
    
{

    for vec in vecvec.iter() {
        if ! is_sorted_strictly( & vec, &mut OrderComparatorAutoLtByKey::new() ) {
            panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
        }
    }

    VecOfVec::new(   
                vecvec,                    
                OrderComparatorAutoLtByKey::new() 
        )
}

impl < 'a, IndexCoeffPair, OrderComparator, >

    IndicesAndCoefficients for 
    &'a VecOfVec< IndexCoeffPair, OrderComparator, >

    where   
        OrderComparator:    StrictlyLess<  IndexCoeffPair >,
        IndexCoeffPair:     KeyValAssociatedTypes,             

{ type KeyMin = IndexCoeffPair::Key; type KeyMaj = usize; type SnzVal = IndexCoeffPair::Val; }  



impl < 'a, IndexCoeffPair, OrderComparator, >
    
    OracleMajor  for 
    
    &'a VecOfVec
        < IndexCoeffPair, OrderComparator, >

    where   //IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
            OrderComparator:    StrictlyLess<  IndexCoeffPair >,    
            IndexCoeffPair:     KeyValAssociatedTypes,        
{
    type ViewMajor          =   &'a[IndexCoeffPair];
    type ViewMajorIntoIter  =   std::slice::Iter< 'a, IndexCoeffPair >;
    type ViewMajorEntry     =   &'a IndexCoeffPair;

    fn view_major( & self, index: usize ) -> &'a[IndexCoeffPair] {
        return self.vec_of_vec[index].as_slice()
    } 
}


impl < 'a, IndexCoeffPair, OrderComparator, >
    
    OracleMajorAscend  for 
    
    &'a VecOfVec 
        < IndexCoeffPair, OrderComparator, >

    where   // IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
            OrderComparator:    StrictlyLess<  IndexCoeffPair >, 
            IndexCoeffPair:     KeyValAssociatedTypes,                 
{
    type ViewMajorAscend            =   &'a[IndexCoeffPair];
    type ViewMajorAscendIntoIter    =   std::slice::Iter< 'a, IndexCoeffPair >;
    type ViewMajorAscendEntry       =   &'a IndexCoeffPair;
        
    /// Assumes that entries in each vector are sorted in ascending order.
    // fn view_major_ascend( & self, index: usize ) -> & Vec< IndexCoeffPair > {
    //     return self.view_major( index )
    // } 
    fn view_major_ascend( & self, index: usize ) -> &'a[IndexCoeffPair] {
        return self.view_major( index )
    }     
}

impl < 'a, IndexCoeffPair, OrderComparator, > 
    
    OracleMajorDescend  for 
    
    &'a VecOfVec 
        < IndexCoeffPair, OrderComparator, >

    where   //IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
            OrderComparator:    StrictlyLess<  IndexCoeffPair >,
            IndexCoeffPair:     KeyValAssociatedTypes,              
{
    type ViewMajorDescend           =   Rev<std::slice::Iter<'a, IndexCoeffPair>>;
    type ViewMajorDescendIntoIter   =   < Self::ViewMajorDescend as IntoIterator >::IntoIter;
    type ViewMajorDescendEntry      =   &'a IndexCoeffPair;
        
    /// Assumes that entries in each vector are sorted in ascending order.    
    fn view_major_descend( & self, index: usize ) -> Rev<std::slice::Iter<'a, IndexCoeffPair >> {
        return self.vec_of_vec[index].iter().rev()
    } 
}








//  Vector of vectors FROM BOROW (STRICTLY SORTED)
//  -------------------------------------------------------------------------------------------------------




/// An vector of strictly sorted vectors, representing a sparse matrix.  
/// 
/// By "strictly sorted" we mean that (i) the entries in each vector appear in ascending order,
/// according to index, and (ii) within each vector, consequtive entries have distinct indices.
/// 
/// Unlike a [`VecOfVec`], this struct holds only a *reference* to a vector of vectors,
/// internally (hence the `Borrow` in `VecOfVecFromBorrow`).  
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::{VecOfVecFromBorrow, vecvec_with_defualt_order};
/// use oat_rust::matrices::matrix_oracle_traits::*;
/// use oat_rust::utilities::partial_order::OrderComparatorAutoLtByKey;
/// 
/// // Standard method to create a row-major vec-of-vec matrix (second function 
/// // argument specifies the order on entries).
/// let matrix  =   VecOfVecFromBorrow::new(
///                     & vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                     OrderComparatorAutoLtByKey::new(),
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
pub struct VecOfVecFromBorrow

    < 'a, IndexCoeffPair, OrderComparator, >

    where   OrderComparator:    StrictlyLess<  IndexCoeffPair >,  
            IndexCoeffPair:     KeyValAssociatedTypes,
{
    vec_of_vec:         &'a Vec< Vec< IndexCoeffPair > >,
    order_comparator:   OrderComparator,
    // phantom_keymin:     PhantomData< KeyMin >,
    // phantom_snzval:     PhantomData< SnzVal >,
}


impl    < 'a, IndexCoeffPair, OrderComparator, >
        
        VecOfVecFromBorrow
        < 'a, IndexCoeffPair, OrderComparator, >
        
        where   OrderComparator:    StrictlyLess<  IndexCoeffPair >,
                IndexCoeffPair:     KeyValAssociatedTypes,        

{
    /// Make a new `VecOfVecFromBorrow`.  Throws an error if one of the internal vectors is not sorted
    /// in strictly asecending order, by index.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::matrices::matrix_types::vec_of_vec::{VecOfVec, VecOfVecFromBorrow};
    /// use oat_rust::utilities::partial_order::OrderComparatorAutoAnyType;
    /// 
    /// // This line should not panic.
    /// let matrix = VecOfVecFromBorrow::new( & vec![ vec![ (0,5), (1,6) ] ], OrderComparatorAutoAnyType::new() ); // no panic
    /// 
    /// // This line should panic.
    /// // let matrix = VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ], OrderComparatorAutoAnyType::new() );
    /// // And here is the test that confirms it.
    /// let result = std::panic::catch_unwind(|| {VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ], OrderComparatorAutoAnyType::new() )} );
    /// assert!(result.is_err()); 
    /// ```
    pub fn new( vecvec: &'a Vec < Vec < IndexCoeffPair > >, mut order_comparator: OrderComparator ) -> Self
    {
        for vec in vecvec.iter() {
            if ! is_sorted_strictly( & vec, &mut order_comparator ) {
                panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
            }
        }
        
        VecOfVecFromBorrow{   
                    vec_of_vec:         vecvec,                    
                    order_comparator:   order_comparator, 
                    // phantom_keymin:     PhantomData,
                    // phantom_snzval:     PhantomData,                 
                }
    }

    /// Constructs a `VecOfVecFromBorrow` out of a reference to a `VecOfVecSimple` **without checking order of entries in each internal vector**.
    /// 
    /// This operation is safe because the data stored in a `VecOfVecSimple` obeys the same ordering rules that a `VecOfVec` obeys.
    pub fn from_vecofvecsimple< KeyMin, SnzVal> ( vecvec_simple: &'a VecOfVecSimple< KeyMin, SnzVal > ) 
        -> 
        VecOfVecFromBorrow< 'a, (KeyMin, SnzVal), OrderComparatorAutoLtByFirstTupleEntry > 
        where
            KeyMin:     PartialOrd
        {
            VecOfVecFromBorrow{ vec_of_vec: vecvec_simple.vec_of_vec_borrowed(), order_comparator: OrderComparatorAutoLtByFirstTupleEntry::new(), }
        }

    /// Return the internally stored reference `&'a Vec<Vec<EntryType>>` and the explicit order comparator.  Consumes the wrapper `VecOfVec` struct.
    pub fn decompose( self ) -> ( &'a Vec < Vec < IndexCoeffPair > >, OrderComparator ) {
        ( self.vec_of_vec, self.order_comparator, )
    }

    /// Returns an immutable reference to the `vec_number`th internal vector.
    pub fn peek_at( & self, vec_number: usize ) -> &Vec< IndexCoeffPair > {        
        & self.vec_of_vec[ vec_number ]
    }

}

// THIS IS NOT AN ASSOCIATED FUNCTION!
pub fn vecvec_with_defualt_order_from_ref
        < 'a, KeyMin, SnzVal, IndexCoeffPair, >
        ( vecvec: &'a Vec < Vec < IndexCoeffPair > > )
        ->
        VecOfVecFromBorrow
        < 'a, IndexCoeffPair, OrderComparatorAutoLtByKey< KeyMin, SnzVal, IndexCoeffPair >, >

        where   IndexCoeffPair:     KeyValAssociatedTypes< Key = KeyMin, Val = SnzVal >,
                IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
                KeyMin:             PartialOrd,
                SnzVal:             ,
    
{

    for vec in vecvec.iter() {
        if ! is_sorted_strictly( & vec, &mut OrderComparatorAutoLtByKey::new() ) {
            panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
        }
    }

    VecOfVecFromBorrow::new(   
                vecvec,                    
                OrderComparatorAutoLtByKey::new() 
        )
}


//  IndicesAndCoefficients
impl < 'a, IndexCoeffPair, OrderComparator, >
    
    IndicesAndCoefficients  for 
    
    VecOfVecFromBorrow
        < 'a, IndexCoeffPair, OrderComparator, >

    where   //IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
            OrderComparator:    StrictlyLess<  IndexCoeffPair >,
            IndexCoeffPair:     KeyValAssociatedTypes,

{   type KeyMin = IndexCoeffPair::Key; type KeyMaj = usize; type SnzVal = IndexCoeffPair::Val; }            


//  OracleMajor
impl < 'a, IndexCoeffPair, OrderComparator, >
    
    OracleMajor  for 
    
    VecOfVecFromBorrow
    < 'a, IndexCoeffPair, OrderComparator, >

    where   //IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
            OrderComparator:    StrictlyLess<  IndexCoeffPair >, 
            IndexCoeffPair:     KeyValAssociatedTypes,               
{
    type ViewMajor          =   &'a[IndexCoeffPair];
    type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajorEntry     =   &'a IndexCoeffPair;    
        
    fn view_major( & self, index: usize ) -> &'a[IndexCoeffPair] {
        return self.vec_of_vec[index].as_slice()
    } 
}


impl < 'a, IndexCoeffPair, OrderComparator, >
    
    OracleMajorAscend  for 
    
    VecOfVecFromBorrow
    < 'a, IndexCoeffPair, OrderComparator, >

    where   // IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
            OrderComparator:    StrictlyLess<  IndexCoeffPair >,   
            IndexCoeffPair:     KeyValAssociatedTypes,             
{
    type ViewMajorAscend            =   &'a[IndexCoeffPair];
    type ViewMajorAscendIntoIter    =   std::slice::Iter< 'a, IndexCoeffPair >;
    type ViewMajorAscendEntry       =   &'a IndexCoeffPair;
    
    /// Assumes that entries in each vector are sorted in ascending order.
    // fn view_major_ascend( & self, index: usize ) -> & Vec< IndexCoeffPair > {
    //     return self.view_major( index )
    // } 
    fn view_major_ascend( & self, index: usize ) -> &'a[IndexCoeffPair] {
        return self.view_major( index )
    }     
}

impl < 'a, IndexCoeffPair, OrderComparator, >
    
    OracleMajorDescend  for 
    
    VecOfVecFromBorrow
    < 'a, IndexCoeffPair, OrderComparator, >

    where   //IndexCoeffPair:     KeyValGet < KeyMin, SnzVal >,
            OrderComparator:    StrictlyLess<  IndexCoeffPair >,
            IndexCoeffPair:     KeyValAssociatedTypes,            
{
    type ViewMajorDescend           =   Rev<std::slice::Iter<'a, IndexCoeffPair>>;
    type ViewMajorDescendIntoIter   =   Rev<std::slice::Iter<'a, IndexCoeffPair>>;
    type ViewMajorDescendEntry      =   &'a IndexCoeffPair;    
        
    /// Assumes that entries in each vector are sorted in ascending order.    
    fn view_major_descend( & self, index: usize ) -> Rev<std::slice::Iter<'a, IndexCoeffPair >> {
        return self.vec_of_vec[index].iter().rev()
    } 
}










//  VecOfVecSimple: A simplified VecOfVec convenient for small examples / human readability
//  -------------------------------------------------------------------------------------------------------

/// Essentially the same as [`VecOfVec`], but with the added requriements that (i) index-coefficient pairs be stored
/// as tuples, and (ii) order on indices is the default implemented by the `PartialOrd` trait.  
/// These restrictions sometimes allow one to omit a large number of type parameters when working with [`VecOfVecSimple`]
/// versus [`VecOfVec`].
/// 
/// Each of the internal vectors should have entries sorted in *strictly* asecneding order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::*;
/// use oat_rust::matrices::matrix_oracle_traits::*;
/// use std::marker::PhantomData;
/// 
/// // Create a new matrix.
/// let matrix  =   VecOfVecSimple::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                 );
/// 
/// 
/// 
/// ```
#[derive(Debug, Clone)]
pub struct  VecOfVecSimple
            < KeyMin, SnzVal >

{
    vec_of_vec:         Vec< Vec< ( KeyMin, SnzVal ) > >,
    order_comparator:   OrderComparatorAutoLtByKey < KeyMin, SnzVal, ( KeyMin, SnzVal ) >,
    // phantom_lifetime:   PhantomData< &'a KeyMin >,
}

impl < KeyMin, SnzVal >
        
        VecOfVecSimple
            < KeyMin, SnzVal > 

{
    /// Make a new `VecOfVecSimple`.  Throws an error if one of the internal vectors is not sorted
    /// in strictly asecending order, by index.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
    /// 
    /// // This line should not panic.
    /// let matrix = VecOfVecSimple::new( vec![ vec![ (0,5), (1,6) ] ] ); // no panic
    /// 
    /// // This line should panic.
    /// // let matrix = VecOfVecSimple::new( vec![ vec![ (1,6), (0,5) ] ] );
    /// // And here is the test that confirms it.
    /// let result = std::panic::catch_unwind(|| {VecOfVecSimple::new( vec![ vec![ (1,6), (0,5) ] ] )} );
    /// assert!(result.is_err()); 
    /// ```
    pub fn new( vecvec: Vec < Vec < ( KeyMin, SnzVal ) > > ) -> Self  
        where   KeyMin: PartialOrd,
                ( KeyMin, SnzVal ):     KeyValGet< KeyMin, SnzVal >, // if we comment this out then  we get an error sayint that `OrderComparatorAutoLtByKey` doesn't implement the `StrictlyLess` trait; probably this has something to do with the fact that `OrderComparatorAutoLtByKey` autoimplements `OrderComparatorAutoLtByKey` for structs that implement `KeyValGet`
    {
        for vec in vecvec.iter() {
            if ! is_sorted_strictly( & vec, &mut OrderComparatorAutoLtByKey::new() ) {
                panic!("Attempt to construct `VecOfVecSimple` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
            }
        }        

        VecOfVecSimple{   
                    vec_of_vec:         vecvec,   
                    order_comparator:   OrderComparatorAutoLtByKey::new(),                 
                    // phantom_lifetime:   PhantomData,                  
                }
    }

    /// Returns a clone of the internally stored order comparator.
    pub fn clone_order_comparator( &self ) -> OrderComparatorAutoLtByKey < KeyMin, SnzVal, ( KeyMin, SnzVal ) > 
        where OrderComparatorAutoLtByKey < KeyMin, SnzVal, ( KeyMin, SnzVal ) > :   Clone
    { self.order_comparator.clone() }

    /// Returns an immutable view of the `Vec< Vec< (usize, SnzVal) > >` that stores the entries of of the matrix, internally.
    pub fn vec_of_vec_borrowed( &self ) -> & Vec< Vec< (KeyMin, SnzVal) > > { & self.vec_of_vec }


    /// Append a vector to the internally stored `Vec<Vec<EntryType>>` struct.  The new vector must be sorted in strictly ascending order; the function panics, otherwise.
    pub fn push_vec( &mut self, new_sorted_vec: Vec < (KeyMin, SnzVal) > ) 
        where 
            OrderComparatorAutoLtByKey<KeyMin, SnzVal, (KeyMin, SnzVal)>: StrictlyLess< (KeyMin, SnzVal)>
    {
        if ! is_sorted_strictly( & new_sorted_vec, &mut self.order_comparator ) {
            panic!("Attempt to append a non-strictly-sorted vector to `VecOfVecSimple`.  To proceed, the new vector must be sorted in *strictly* ascending order, by index.") 
        }        
        self.vec_of_vec.push( new_sorted_vec );
    }    

    /// Creates a [`VecOfVecSimple`] from an iterable that runs over iterables that run over tuples.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
    /// use oat_rust::matrices::matrix_oracle_traits::OracleMajorAscend;
    /// 
    /// let iter = (0..2).map( |x| vec![(x,x)] );
    /// let vec_of_vec = VecOfVecSimple::from_iterable( iter );
    /// itertools::assert_equal( (& vec_of_vec).view_major_ascend(0), vec![(0,0)]);
    /// itertools::assert_equal( (& vec_of_vec).view_major_ascend(1), vec![(1,1)]) 
    /// ```
    pub fn from_iterable< I >( iter: I ) -> VecOfVecSimple< KeyMin, SnzVal > 
        where
            I:          IntoIterator,
            I::Item:    IntoIterator< Item = (KeyMin, SnzVal) >,
            KeyMin:     Clone + PartialOrd,
            SnzVal:     Clone,
    {
        let vecvec =    iter.into_iter().map( |x| x.into_iter().collect_vec() ).collect_vec();
        VecOfVecSimple::new( vecvec )
    }


}


//  ORACLE IMPLEMENTATIONS FOR &'a VecOfVecSimple
//  ---------------------------------------------


//  IndicesAndCoefficients
impl < 'a, KeyMin, SnzVal > 

    IndicesAndCoefficients for 
    &'a VecOfVecSimple < KeyMin, SnzVal >

{ type KeyMin = KeyMin; type KeyMaj = usize; type SnzVal = SnzVal; }   


// OracleMajor
impl < 'a, KeyMin, SnzVal > 
    
    OracleMajor for 
    
    &'a VecOfVecSimple 
        < KeyMin, SnzVal > 

    where   KeyMin:     Clone,    
            SnzVal:     Clone,
{      
    type ViewMajor          =   Cloned< Iter< 'a, (KeyMin, SnzVal) > >;
    type ViewMajorIntoIter  =   Cloned< Iter< 'a, (KeyMin, SnzVal) > >;
    type ViewMajorEntry     =   (KeyMin, SnzVal);

    fn view_major( & self, index: usize ) -> Cloned< Iter< 'a, (KeyMin, SnzVal) > > {
        return self.vec_of_vec[index].iter().cloned()
    } 
}

// impl < 'a, KeyMin, SnzVal > OracleRefInherit for &'a VecOfVecSimple < KeyMin, SnzVal > {}

impl < 'a, KeyMin, SnzVal > 
    
    OracleMajorAscend  for 
    
    &'a VecOfVecSimple < KeyMin, SnzVal > 

    where   KeyMin:     Clone,    
            SnzVal:     Clone,    
{
    type ViewMajorAscend            =   Cloned< Iter< 'a, (KeyMin, SnzVal) > >;
    type ViewMajorAscendIntoIter    =   Cloned< Iter< 'a, (KeyMin, SnzVal) > >;
    type ViewMajorAscendEntry       =   (KeyMin, SnzVal);

    /// Assumes that entries in each vector are sorted in ascending order.
    fn view_major_ascend( & self, index: usize ) -> Cloned< Iter< 'a, (KeyMin, SnzVal) > > {
        return self.view_major( index )
    } 
}

impl < 'a, KeyMin, SnzVal > 
    
    OracleMajorDescend  for 
    
    &'a VecOfVecSimple < KeyMin, SnzVal > 

    where   KeyMin:     Clone,    
            SnzVal:     Clone,  
{
    type ViewMajorDescend           =   Cloned<Rev<std::slice::Iter<'a, (KeyMin, SnzVal)>>>;
    type ViewMajorDescendIntoIter   =   Cloned<Rev<std::slice::Iter<'a, (KeyMin, SnzVal)>>>;
    type ViewMajorDescendEntry      =   (KeyMin, SnzVal);

    /// Assumes that entries in each vector are sorted in ascending order.    
    fn view_major_descend( & self, index: usize ) -> Cloned<Rev<std::slice::Iter<'a, (KeyMin, SnzVal)>>> {
        return self.vec_of_vec[index].iter().rev().cloned()
    } 
}

/// Represents a minor view of a `VecOfVecSimple`, with entries appearing in descending order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::{VecOfVecSimple, VecOfVecSimpleViewMinorDescend};
/// use oat_rust::matrices::matrix_oracle_traits::OracleMinorDescend;
///         
/// // define a row-major sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVecSimple::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_descend( 0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_descend( 1 ) );
/// itertools::assert_equal( vec![ (1,7), (0,5) ], (& matrix).view_minor_descend( 2 ) );  
/// ```
pub struct VecOfVecSimpleViewMinorDescend< 'a, KeyMin, SnzVal > 
    where   KeyMin:     Clone,    
            SnzVal:     Clone,  
{
    vec_of_vec:             &'a VecOfVecSimple< KeyMin, SnzVal >,
    keymaj:                 usize,
    keymin:                 KeyMin,
    phantom_keymin:         PhantomData< KeyMin >,
}

// Iterator
impl < 'a, KeyMin, SnzVal > 

        Iterator for 

        VecOfVecSimpleViewMinorDescend< 'a, KeyMin, SnzVal >

    where   KeyMin:     Clone + PartialEq,    
            SnzVal:     Clone,          
{
    type Item = (usize, SnzVal);

    fn next( &mut self ) -> Option< Self::Item > {
        // extract the underlying vector of vectors
        let vecvec_data =  self.vec_of_vec.vec_of_vec_borrowed();                
        // (assuming the matrix is row-major) scan the rows of the matrix from bottom to top
        while self.keymaj > 0 {
            // drop the row number by 1
            self.keymaj -= 1;            
            // get the row
            let view_major = & vecvec_data[ self.keymaj ];
            // scan the row to see if it contains an entry of form ( my_keymin, snzval )
            for ( keymin, snzval ) in view_major {
                // if it does, then return ( row_number, snzval )
                if keymin != & self.keymin { continue }
                else { return Some( ( self.keymaj, snzval.clone() ) ) }
            }
        }
        return None
    }
}        

/// Represents a minor view of a `VecOfVecSimple`, with entries appearing in ascending order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::{VecOfVecSimple, VecOfVecSimpleViewMinorAscend};
/// use oat_rust::matrices::matrix_oracle_traits::OracleMinorAscend;
///         
/// // define a row-major sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVecSimple::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_ascend( 0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_ascend( 1 ) );
/// itertools::assert_equal( vec![ (0,5), (1,7) ], (& matrix).view_minor_ascend( 2 ) );  
/// ```
pub struct VecOfVecSimpleViewMinorAscend< 'a, KeyMin, SnzVal > 
    where   KeyMin:     Clone,    
            SnzVal:     Clone,  
{
    vec_of_vec:             &'a VecOfVecSimple< KeyMin, SnzVal >,
    keymaj:                 usize,
    keymin:                 KeyMin,
    phantom_keymin:         PhantomData< KeyMin >,
}

// Iterator
impl < 'a, KeyMin, SnzVal > 

        Iterator for 

        VecOfVecSimpleViewMinorAscend< 'a, KeyMin, SnzVal >

    where   KeyMin:     Clone + PartialEq,    
            SnzVal:     Clone,          
{
    type Item = (usize, SnzVal);

    fn next( &mut self ) -> Option< Self::Item > {
        // extract the underlying vector of vectors
        let vecvec_data =  self.vec_of_vec.vec_of_vec_borrowed();                
        // (assuming the matrix is row-major) scan the rows of the matrix from bottom to top
        while self.keymaj < vecvec_data.len() {
            println!("for debuggin: keymaj: {:?}", self.keymaj);
            // get the row
            let view_major = & vecvec_data[ self.keymaj ];
            // scan the row to see if it contains an entry of form ( my_keymin, snzval )
            for ( keymin, snzval ) in view_major {
                // if it does, then return ( row_number, snzval )
                if keymin != & self.keymin { continue }
                else { 
                    // grow the row number by 1 (this is one of two branches where we do this)
                    self.keymaj += 1;
                    
                    // return the entry
                    return Some( ( self.keymaj - 1, snzval.clone() ) ) 
                }
            }
            // grow the row number by 1 (this is one of two branches where we do this)
            self.keymaj += 1;              
        }      

        // in this case the iterator is exhausted
        return None
    }
}            


//  OracleMinorDescend
impl < 'a, KeyMin, SnzVal > 

    OracleMinorDescend for 
    
    &'a VecOfVecSimple
        < KeyMin, SnzVal >

    where   KeyMin:     Clone + PartialEq,    
            SnzVal:     Clone,          
    
{
    type ViewMinorDescend = VecOfVecSimpleViewMinorDescend< 'a, KeyMin, SnzVal >;
    type ViewMinorDescendIntoIter = Self::ViewMinorDescend;
    type ViewMinorDescendEntry = (usize, SnzVal);

    fn view_minor_descend( &self, index: KeyMin ) -> VecOfVecSimpleViewMinorDescend< 'a, KeyMin, SnzVal > {
        VecOfVecSimpleViewMinorDescend{
            vec_of_vec:             self,
            keymaj:                 self.vec_of_vec.len(),
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }
}

//  OracleMinorAscend
impl < 'a, KeyMin, SnzVal > 

    OracleMinorAscend for 
    
    &'a VecOfVecSimple
        < KeyMin, SnzVal >

    where   KeyMin:     Clone + PartialEq,    
            SnzVal:     Clone,          
    
{
    type ViewMinorAscend = VecOfVecSimpleViewMinorAscend< 'a, KeyMin, SnzVal >;
    type ViewMinorAscendIntoIter = Self::ViewMinorAscend;
    type ViewMinorAscendEntry = (usize, SnzVal);

    fn view_minor_ascend( &self, index: KeyMin ) -> VecOfVecSimpleViewMinorAscend< 'a, KeyMin, SnzVal > {
        VecOfVecSimpleViewMinorAscend{
            vec_of_vec:             self,
            keymaj:                 0,
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }
}

//  OracleMinor
impl < 'a, KeyMin, SnzVal > 

    OracleMinor for 
    
    &'a VecOfVecSimple
        < KeyMin, SnzVal >

    where   KeyMin:     Clone + PartialEq,    
            SnzVal:     Clone,          
    
{
    type ViewMinor = VecOfVecSimpleViewMinorAscend< 'a, KeyMin, SnzVal >;
    type ViewMinorIntoIter = Self::ViewMinor;
    type ViewMinorEntry = (usize, SnzVal);

    fn view_minor( &self, index: KeyMin ) -> VecOfVecSimpleViewMinorAscend< 'a, KeyMin, SnzVal > {
        VecOfVecSimpleViewMinorAscend{
            vec_of_vec:             self,
            keymaj:                 0,
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }
}


//  ORACLE IMPLEMENTATIONS FOR &'b &'a VecOfVecSimple -- NOW DEPRECATED IN FAVOR OF THE TRAIT OracleRefInherit
//  -------------------------------------------------

// impl < 'a, 'b, KeyMin, SnzVal > 
    
//     OracleMajor
//     < usize, Cloned< Iter< 'b, (KeyMin, SnzVal) > > > 
    
//     for 
    
//     &'b &'a VecOfVecSimple 
//     < KeyMin, SnzVal > 

//     where   KeyMin:     Clone,    
//             SnzVal:     Clone,
// {
        
//     fn view_major( & self, index: usize ) -> Cloned< Iter< 'a, (KeyMin, SnzVal) > > {
//         return self.vec_of_vec[index].iter().cloned()
//     } 
// }

// impl < 'a, 'b, KeyMin, SnzVal > 
    
//     OracleMajorAscend
//     < usize, Cloned< Iter< 'b, (KeyMin, SnzVal) > >, >  for 
    
//     &'b &'a VecOfVecSimple 
//     < KeyMin, SnzVal > 

//     where   KeyMin:     Clone,    
//             SnzVal:     Clone,    
// {
//     /// Assumes that entries in each vector are sorted in ascending order.
//     fn view_major_ascend( & self, index: usize ) -> Cloned< Iter< 'a, (KeyMin, SnzVal) > > {
//         return self.view_major( index )
//     } 
// }

// impl < 'a, 'b, KeyMin, SnzVal > 
    
//     OracleMajorDescend
//     < usize, Cloned<Rev<std::slice::Iter<'a, (KeyMin, SnzVal)>>>, > 
    
//     for 
    
//     &'b &'a VecOfVecSimple < KeyMin, SnzVal > 

//     where   KeyMin:     Clone,    
//             SnzVal:     Clone,  
// {
//     /// Assumes that entries in each vector are sorted in ascending order.    
//     fn view_major_descend( & self, index: usize ) -> Cloned<Rev<std::slice::Iter<'a, (KeyMin, SnzVal)>>> {
//         return self.vec_of_vec[index].iter().rev().cloned()
//     } 
// }










//  VecOfVecSimpleFromBorrow: A simplified VecOfVec convenient for small examples / human readability
//  -------------------------------------------------------------------------------------------------------

/// Essentially the same as [`VecOfVecFromBorrow`], but with the added requriements that (i) index-coefficient pairs be stored
/// as tuples, and (ii) order on indices is the default implemented by the `PartialOrd` trait.  
/// These restrictions sometimes allow one to omit a large number of type parameters when working with [`VecOfVecSimpleFromBorrow`]
/// versus [`VecOfVecFromBorrow`].
/// 
/// Each of the internal vectors should have entries sorted in *strictly* asecneding order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::*;
/// use oat_rust::matrices::matrix_oracle_traits::*;
/// use std::marker::PhantomData;
/// 
/// // Create a new matrix.
/// let matrix  =   VecOfVecSimple::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                 );
/// 
/// 
/// 
/// ```
#[derive(Debug, Clone)]
pub struct  VecOfVecSimpleFromBorrow
            < 'a, KeyMin, SnzVal >

{
    vec_of_vec:         &'a Vec< Vec< ( KeyMin, SnzVal ) > >,
    order_comparator:   OrderComparatorAutoLtByKey < KeyMin, SnzVal, ( KeyMin, SnzVal ) >,
}

impl    < 'a, KeyMin, SnzVal >
        VecOfVecSimpleFromBorrow
        < 'a, KeyMin, SnzVal > 

{
    /// Make a new `VecOfVecSimple`.  Throws an error if one of the internal vectors is not sorted
    /// in strictly asecending order, by index.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
    /// 
    /// // This line should not panic.
    /// let matrix = VecOfVecSimple::new( vec![ vec![ (0,5), (1,6) ] ] ); // no panic
    /// 
    /// // This line should panic.
    /// // let matrix = VecOfVecSimple::new( vec![ vec![ (1,6), (0,5) ] ] );
    /// // And here is the test that confirms it.
    /// let result = std::panic::catch_unwind(|| {VecOfVecSimple::new( vec![ vec![ (1,6), (0,5) ] ] )} );
    /// assert!(result.is_err()); 
    /// ```
    pub fn new( vecvec: &'a Vec < Vec < ( KeyMin, SnzVal ) > > ) -> Self  
        where   KeyMin: PartialOrd,
                ( KeyMin, SnzVal ):     KeyValGet< KeyMin, SnzVal >, // if we comment this out then  we get an error sayint that `OrderComparatorAutoLtByKey` doesn't implement the `StrictlyLess` trait; probably this has something to do with the fact that `OrderComparatorAutoLtByKey` autoimplements `OrderComparatorAutoLtByKey` for structs that implement `KeyValGet`
    {
        for vec in vecvec.iter() {
            if ! is_sorted_strictly( & vec, &mut OrderComparatorAutoLtByKey::new() ) {
                panic!("Attempt to construct `VecOfVecSimple` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
            }
        }        

        VecOfVecSimpleFromBorrow{   
                    vec_of_vec:         vecvec,   
                    order_comparator:   OrderComparatorAutoLtByKey::new(),                 
                }
    }

    /// Returns a clone of the internally stored order comparator.
    pub fn clone_order_comparator( &self ) -> OrderComparatorAutoLtByKey < KeyMin, SnzVal, ( KeyMin, SnzVal ) > 
        where OrderComparatorAutoLtByKey < KeyMin, SnzVal, ( KeyMin, SnzVal ) > :   Clone
    { self.order_comparator.clone() }

    /// Returns a copy of the reference `& Vec< Vec< (KeyMin, SnzVal) > >` where the data of the matrix is stored.
    pub fn vec_of_vec_borrowed( &self ) -> &'a Vec< Vec< (KeyMin, SnzVal) > > 
    { self.vec_of_vec }    

}


// IndicesAndCoefficients
impl < 'a, KeyMin, SnzVal > 
    
    IndicesAndCoefficients  for 
    VecOfVecSimpleFromBorrow< 'a, KeyMin, SnzVal > 
{ type KeyMin = KeyMin; type KeyMaj = usize; type SnzVal = SnzVal; }    


// OracleMajor
impl < 'a, KeyMin, SnzVal > 
    
    OracleMajor  for 
    
    VecOfVecSimpleFromBorrow
        < 'a, KeyMin, SnzVal > 

    where   KeyMin:     Clone,    
            SnzVal:     Clone,
{
    type ViewMajor          =   Cloned< Iter< 'a, (KeyMin, SnzVal) > >;
    type ViewMajorIntoIter  =   Cloned< Iter< 'a, (KeyMin, SnzVal) > >;
    type ViewMajorEntry     =   (KeyMin, SnzVal);

    fn view_major( & self, index: usize ) -> Cloned< Iter< 'a, (KeyMin, SnzVal) > > {
        return self.vec_of_vec[index].iter().cloned()
    } 
}

impl < 'a, KeyMin, SnzVal > 
    
    OracleMajorAscend  for 
    
    VecOfVecSimpleFromBorrow 
        < 'a, KeyMin, SnzVal > 

    where   KeyMin:     Clone,    
            SnzVal:     Clone,    
{
    type ViewMajorAscend = Cloned< Iter< 'a, (KeyMin, SnzVal) > >;
    type ViewMajorAscendIntoIter = < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscendEntry = (KeyMin, SnzVal);

    /// Assumes that entries in each vector are sorted in ascending order.
    fn view_major_ascend( & self, index: usize ) -> Cloned< Iter< 'a, (KeyMin, SnzVal) > > {
        return self.view_major( index )
    } 
}

impl < 'a, KeyMin, SnzVal > 
    
    OracleMajorDescend  for 
    
    VecOfVecSimpleFromBorrow
        < 'a, KeyMin, SnzVal > 

    where   KeyMin:     Clone,    
            SnzVal:     Clone,  
{
    type ViewMajorDescend = Cloned<Rev<std::slice::Iter<'a, (KeyMin, SnzVal)>>>;
    type ViewMajorDescendIntoIter = < Self::ViewMajorDescend as IntoIterator >::IntoIter;
    type ViewMajorDescendEntry = (KeyMin, SnzVal);

    /// Assumes that entries in each vector are sorted in ascending order.    
    fn view_major_descend( & self, index: usize ) -> Cloned<Rev<std::slice::Iter<'a, (KeyMin, SnzVal)>>> {
        return self.vec_of_vec[index].iter().rev().cloned()
    } 
}






















#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_vec_of_vec_from_iterable() {
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        
        let iter = (0..2).map( |x| vec![(x,x)] );
        let vec_of_vec = VecOfVecSimple::from_iterable( iter );
        itertools::assert_equal( (& vec_of_vec).view_major_ascend(0), vec![(0,0)]);
        itertools::assert_equal( (& vec_of_vec).view_major_ascend(1), vec![(1,1)])   
    }

    #[test]
    fn test_vec_of_vec_construction() {

        let _matrix  =   VecOfVecUnsorted::new(
                                                vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ]
                                                    );                                      
    }

    #[test]    
    fn test_vec_of_vec_simple_descending_minor_view() {
        use crate::matrices::matrix_types::vec_of_vec::{VecOfVecSimple, VecOfVecSimpleViewMinorDescend};
        use crate::matrices::matrix_oracle_traits::OracleMinorDescend;
        
        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVecSimple::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
        itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_descend( 0 ) );
        itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_descend( 1 ) );
        itertools::assert_equal( vec![ (1,7), (0,5) ], (& matrix).view_minor_descend( 2 ) );        
    }

    #[test]     
    fn test_vec_of_vec_simple_ascending_minor_view() {
        use crate::matrices::matrix_types::vec_of_vec::{VecOfVecSimple, VecOfVecSimpleViewMinorAscend};
        use crate::matrices::matrix_oracle_traits::OracleMinorAscend;

        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVecSimple::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
        println!("waytpoint 1");
        itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_ascend( 0 ) );
        println!("waytpoint 2");        
        itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_ascend( 1 ) );
        println!("waytpoint 3");        
        itertools::assert_equal( vec![ (0,5), (1,7) ], (& matrix).view_minor_ascend( 2 ) ); 
    }

}
