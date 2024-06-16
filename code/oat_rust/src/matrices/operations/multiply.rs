//! Matrix-matrix and matrix-vector multiplication.
//! 
//! 
//! - to multiply a matrix with a matrix, use 
//!   - [ProductMatrixLazyMajorAscendSimplified::new](crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified)
//!   - [ProductMatrixLazyMajorAscendUnsimplified::new](crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendUnsimplified)
//! - to multiply a matrix with a vector, use 
//!   - [vector_matrix_multiply_major_ascend_simplified](    crate::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified)
//!   - [vector_matrix_multiply_major_ascend_unsimplified](  crate::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified)
//!   - [vector_matrix_multiply_minor_descend_simplified](   crate::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified)
//!   - [vector_matrix_multiply_minor_descend_unsimplified]( crate::matrices::operations::multiply::vector_matrix_multiply_minor_descend_unsimplified)


use std::{marker::PhantomData};
use crate::{matrices::{matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend, IndicesAndCoefficients}, matrix_types::oracle_ref::OracleRef}, rings::operator_traits::{Semiring}, utilities::{iterators::merge::heap_of_iterators::{HitMerge, hit_merge_by_predicate}, partial_order::{StrictlyLess, OrderComparatorReverse}}, vectors::{operations::{Scale, Transforms, Simplify}, linear_combinations::{LinearCombinationSimplified, LinearCombinationUnsimplified}}};
use crate::entries::{KeyValGet, KeyValSet};



//  MATRIX - VECTOR MULTIPLICATION
//  ===========================================================================


//  MAJOR ASCEND
//  ---------------------------------------------------------------------------

/// Returns the product of a matrix with a vector.
/// 
/// The resulting iterator is a linear combination of major views.  It returns entries in ascending (but **non-strictly** ascending) order of index, provided that
/// [OracleMajorAscend](crate::matrices::matrix_oracle_traits::OracleMajorAscend) returns entries in 
/// ascending order according to index.
/// 
/// **Note** the resulting iterator is not simplified -- it may return several entries with the same index.
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `vA`, a linear combination of **rows** of `A`.
/// If instead the representation is column-major, then this function returns `Av`.
/// 
/// ```
/// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::partial_order::OrderComparatorAutoLtByKey;
/// use oat_rust::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0usize, 1), (1usize, 1) ];
/// 
/// // the product
/// let product     =   vector_matrix_multiply_major_ascend_unsimplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderComparatorAutoLtByKey::new(), // this compares tuples according to their lexicographic order
///                     );
/// // the product should return the following sequence of entries: (0,1), (0,1), (1,1), (1,6)   
/// itertools::assert_equal( product, vec![(0isize,1usize), (0,1), (1,1), (1,6)] )
/// ```
pub fn vector_matrix_multiply_major_ascend_unsimplified 
            < Matrix, RingOperator, SparseVecIter, OrderComparator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_comparator:   OrderComparator,
            ) 
            ->
            LinearCombinationUnsimplified
                < Matrix::ViewMajorAscendIntoIter, Matrix::KeyMin, Matrix::SnzVal, RingOperator, OrderComparator >            
    where 
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        Matrix::SnzVal:                 Clone, // Clone is required for HitMerge
        OrderComparator:                Clone + StrictlyLess<  Matrix::ViewMajorAscendEntry >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < Matrix::KeyMaj, Matrix::SnzVal >,
{
    LinearCombinationUnsimplified{
        linear_combination_unsimplified:
            hit_merge_by_predicate(
                sparse_vec
                        .into_iter()
                        .map(   |x| 
                                matrix
                                    .view_major_ascend(x.key() )
                                    .into_iter()
                                    .scale( x.val(), ring_operator.clone() )
                        ),
                    order_comparator.clone()
            )        
    }
}


/// Returns the product of a matrix with a vector; the result is a linear combination of the 
/// matrix's major views.
/// 
/// The resulting iterator is a simplified linear combination of major views.  It returns entries in **strictly** ascending order of index, provided that
/// [OracleMajorAscend](crate::matrices::matrix_oracle_traits::OracleMajorAscend) returns entries in 
/// ascending order according to index.  In particular, no two consequtive indices have the same
/// index.
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `vA`, a linear combination of **rows** of `A`.
/// If instead the representation is column-major, then this function returns `Av`.
/// 
/// ```
/// use oat_rust::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;
/// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::partial_order::OrderComparatorAutoLtByKey;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0usize, 1), (1usize, 1) ];
/// 
/// // the product
/// let product     =   vector_matrix_multiply_major_ascend_simplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderComparatorAutoLtByKey::new(), // this compares tuples according to their lexicographic order
///                     );
/// // the product should equal vector [2 0]        
/// itertools::assert_equal( product, vec![(0isize, 2)] )
/// ```
pub fn vector_matrix_multiply_major_ascend_simplified 
            < Matrix, RingOperator, SparseVecIter, OrderComparator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_comparator:   OrderComparator,
            ) 
            ->
            LinearCombinationSimplified
                < Matrix::ViewMajorAscendIntoIter, Matrix::KeyMin, Matrix::SnzVal, RingOperator, OrderComparator >            
    where 
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone, // Clone is required for HitMerge
        Matrix::KeyMin:                 PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros        
        OrderComparator:                Clone + StrictlyLess<  Matrix::ViewMajorAscendEntry >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < Matrix::KeyMaj, Matrix::SnzVal >,
{
    let linear_combination_simplified = 
    vector_matrix_multiply_major_ascend_unsimplified( 
            sparse_vec,
            matrix,
            ring_operator.clone(),
            order_comparator,
        )
        .linear_combination_unsimplified
        .simplify(ring_operator);

        return LinearCombinationSimplified{ linear_combination_simplified: linear_combination_simplified }
}        



//  MINOR DESCEND
//  ---------------------------------------------------------------------------


/// Returns the product of a matrix with a vector; the resulting iterator is an unsimplified linear combination of the matrix's minor views.
/// 
/// The resulting iterator returns entries in descending (but **non-strictly** descending) order of index, provided that
/// [OracleMajorAscend](crate::matrices::matrix_oracle_traits::OracleMajorDescend) returns entries in 
/// descending order according to index.
/// 
/// **Note** the resulting iterator is not simplified -- it may return several entries with the same index.
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `A v`, a linear combination of **columns** of `A`.
/// If instead the representation is column-major, then this function returns `vA`.
/// 
/// ```
/// use oat_rust::matrices::operations::multiply::vector_matrix_multiply_minor_descend_unsimplified;
/// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::partial_order::OrderComparatorAutoAnyType;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0isize, 1), (1isize, 1) ];
/// 
/// // the product
/// let mut product     =   vector_matrix_multiply_minor_descend_unsimplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderComparatorAutoAnyType, // this compares tuples according to their lexicographic order
///                     );
/// // the product should return the following entries in sequence: (1,1), (1,1), (0,1), (0,6)       
/// itertools::assert_equal( product, vec![(1,1), (1,1), (0,6), (0,1)] )
/// ```
pub fn vector_matrix_multiply_minor_descend_unsimplified 
            < Matrix, KeyMaj, KeyMin, SnzVal, RingOperator, SparseVecIter, OrderComparator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_comparator:   OrderComparator,
            ) 
            ->
            LinearCombinationUnsimplified
                < Matrix::ViewMinorDescendIntoIter, KeyMaj, SnzVal, RingOperator, OrderComparatorReverse< OrderComparator > >            
    where 
        Matrix:                             OracleMinorDescend + IndicesAndCoefficients< KeyMin = KeyMin >,
        Matrix::ViewMinorDescendIntoIter:   Iterator,        
        < Matrix::ViewMinorDescendIntoIter as Iterator >::Item:     KeyValSet< KeyMaj, SnzVal >,
        RingOperator:                       Clone + Semiring< SnzVal >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        SnzVal:                             Clone, // Clone is required for HitMerge
        OrderComparator:                    Clone + StrictlyLess<  <Matrix::ViewMinorDescendIntoIter as Iterator>::Item >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                      IntoIterator,
        SparseVecIter::Item:                KeyValGet < KeyMin, SnzVal >,
{
    LinearCombinationUnsimplified{
        linear_combination_unsimplified:
            hit_merge_by_predicate(
                sparse_vec
                        .into_iter()
                        .map(   |x| 
                                matrix
                                    .view_minor_descend(x.key() )
                                    .into_iter()
                                    .scale( x.val(), ring_operator.clone() )
                        ),
                    OrderComparatorReverse::new( order_comparator.clone() )
            )        
    }
}


/// Returns the product of a matrix with a vector; the result is a linear combination of the matrix's minor views.
/// 
/// The resulting iterator is a simplified linear combinatino of minor views.  It returns entries in **strictly** descending order of index, provided that
/// [OracleMajorAscend](crate::matrices::matrix_oracle_traits::OracleMajorAscend) returns entries in 
/// ascending order according to index.  In particular, no two consequtive indices have the same
/// index.
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `A v`, a linear combination of **columns** of `A`.
/// If instead the representation is column-major, then this function returns `vA`.
/// 
/// ```
/// use oat_rust::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified;
/// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::partial_order::OrderComparatorAutoLtByKey;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0isize, 1), (1isize, 1) ];
/// 
/// // the product
/// let product     =   vector_matrix_multiply_minor_descend_simplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderComparatorAutoLtByKey::new(), // this compares order of entries according to their indices
///                     );
/// // the product should equal vector [0 1]        
/// itertools::assert_equal( product, vec![(1usize, 2)] )
/// ```
pub fn vector_matrix_multiply_minor_descend_simplified 
            < Matrix, KeyMaj, KeyMin, SnzVal, RingOperator, SparseVecIter, OrderComparator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_comparator:   OrderComparator,
            ) 
            ->
            LinearCombinationSimplified
                < Matrix::ViewMinorDescendIntoIter, KeyMaj, SnzVal, RingOperator, OrderComparatorReverse< OrderComparator > >

    where 
        Matrix:                         OracleMinorDescend< KeyMin = KeyMin, KeyMaj = KeyMaj >,
        Matrix::ViewMinorDescend:       IntoIterator,
        Matrix::ViewMinorDescendEntry:  KeyValGet < KeyMaj, SnzVal > + KeyValSet < KeyMaj, SnzVal >,
        RingOperator:                   Clone + Semiring< SnzVal >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        SnzVal:                         Clone, // Clone is required for HitMerge
        KeyMaj:                         std::cmp::PartialEq, // required by the struct that simplifies vectors (it has to compare the indices of different entries)
        OrderComparator:                Clone + StrictlyLess<  < Matrix::ViewMinorDescend as IntoIterator >::Item >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < KeyMin, SnzVal >
{
    let linear_combination_simplified = 
    vector_matrix_multiply_minor_descend_unsimplified( 
            sparse_vec,
            matrix,
            ring_operator.clone(),
            order_comparator,
        )
        .linear_combination_unsimplified
        .simplify(ring_operator);

        return LinearCombinationSimplified{ linear_combination_simplified: linear_combination_simplified }
}    


























//  MATRIX - MATRIX MULTIPLICATION
//  ===========================================================================

/// Lazy product of two matrix oracles, represented by another matrix oracle
/// 
/// This struct represnts the product of two matrices.  It contains only four essential pieces of data:
/// * a reference to a matrix, `A`,
/// * a reference to another matrix, `B`,
/// * an object that specifies the coefficient ring over which both matrices are defined
/// * an object that specifies a total order on indices (it should be compatible with the order in which `A` and `B` return entries)
/// 
/// When the user requests a *row*, the lazy object first (i) scales 
/// some rows of `B` by an appropriate factor, then (ii) merges those rows into a single iterator, `J`.  Iterator `J` returns
/// entries in sorted order (ascending), but does *not* simplify terms: in particular, it does not drop zero terms or combine terms that 
/// have the same coefficient, then (iii) wraps `J` in a [`Simplify`](`crate::vectors::transfomations::Simplify`) struct that
/// *does* combine terms and drop zeros.
/// 
/// When the user requests a *column*, the lazy object first (i) scales 
/// some columns of `B` by an appropriate factor, then (ii) merges those columns into a single iterator, `J`.  Iterator `J` returns
/// entries in sorted order (ascending), but does *not* simplify terms: in particular, it does not drop zero terms or combine terms that 
/// have the same coefficient, then (iii) wraps `J` in a [`Simplify`](`crate::vectors::transfomations::Simplify`) struct that
/// *does* combine terms and drop zeros.
/// 
/// **Major dimension matters**. The major dimension is the dimension you 
/// actually get when you call `A.view_major_ascend(..)`.  The reason for this is the following pair of facts,
/// (both facts follow from the prodecure by which rows/columns are constructed, as described above).
/// * If `A` and `B` are actually row-major, then the lazy product `matrix_multiply_major_ascend(A, B, ..)` represents `A * B`
/// * If `A` and `B` are actually col-major, then the lazy product `matrix_multiply_major_ascend(A, B, ..)` represents `B * A`
/// 
/// 
/// **Example**
/// Suppose that `X` represents the product of `A` and `B` (within `X`, we refer to `A` and `B` as `matrix_1` and `matrix_2`,
/// respectively).  If we call `X.view_major_ascend( i )`, then the following happens
/// 1) Let `(j1, a1)`, .., `(jn, an)` be the structural nonzero entries in the `i`th major view of `A`.
/// 2) Let `Jk` denote the iterator obtained by scaling row `jk` of matrix `B` by a factor of `ak`.
/// 3) Let `J` be the iterator obtained by merging `J1, .., Jn` via [hit_merge_by_predicate].
/// 4) Let `K` be the iterator obtained by combining all entries with the same index into a single entry, and dropping this entry if it is zero.
/// Then `K` is the output returned by `X.view_major_ascend( i )`.
/// 
/// # Examples
/// 
/// ```
/// // Import the necessary traits and structs
/// use oat_rust::rings::operator_structs::ring_native::DivisionRingNative;
/// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
/// use oat_rust::matrices::matrix_oracle_traits::{OracleMajorAscend};
/// use oat_rust::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;
/// use oat_rust::entries::KeyValGet;
/// use oat_rust::utilities::partial_order::OrderComparatorFnWrapper;
/// use std::iter::FromIterator;
/// // Define the operator for the coefficient ring
/// let ring_operator = DivisionRingNative::<f64>::new();     
/// // Define matrix A, a 2x2 matrix with 1's above the diagonal and 0's below
/// let matrix_1   =   VecOfVecSimple::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
/// // Define matrix B, a 2x2 matrix with 1's above the diagonal and 0's below
/// let matrix_2   =   VecOfVecSimple::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
/// // Define references to the two arrays (recall that OracleMajorAscend is only implemented on `& VecOfVecSimple`, not on `VecOfVecSimple`)
/// let matrix_1_ref   =   &matrix_1;
/// let matrix_2_ref    =   &matrix_2;        
/// // Define the lazy product of A and B
/// let product     =   ProductMatrixLazyMajorAscendSimplified::new( 
///                             matrix_1_ref,                 // matrix A
///                             matrix_2_ref,                 // matrix B
///                             ring_operator.clone(),      // ring operator
///                             OrderComparatorFnWrapper::new( |x:&(i32, f64), y:&(i32, f64)| x.key() < y.key() )   // defines a total order on indices
///                         );
///                     
/// // First, check the output iterators themselves (these are wrapped in `Simplify` structs).
/// // ---------------------------------------------------------------------------------------
///                     
/// // Use this object to create a vector-of-vectors that represents the product
/// let output =   Vec::from_iter(
///                                         (0..2)
///                                             .map(   |x| 
///                                                     Vec::from_iter(
///                                                             product
///                                                             .view_major_ascend(x)
///                                                         )
///                                                 )
///                                     );
///                                 
/// // Check that the answer is correct                                    
/// assert_eq!(     
///         output,
///         vec![ vec![(0,1.), (1,2.)], vec![(1,1.)]  ]
///     ); 
/// 
/// // Second, check the underlying, unsimplified iterators (these are are contained inside the `Simplify` structs).
/// // -------------------------------------------------------------------------------------------------------------
/// 
/// // Use this object to create a vector-of-vectors that represents the product
/// let output =   Vec::from_iter(
///                                         (0..2)
///                                             .map(   |x| 
///                                                     Vec::from_iter(
///                                                             product
///                                                             .view_major_ascend(x)
///                                                             .linear_combination_simplified
///                                                             .unsimplified // this unwraps the inner iterator
///                                                         )
///                                                 )
///                                     );
///                                 
/// // Check that the answer is correct                                    
/// assert_eq!(     
///         output,
///         vec![ vec![(0,1.), (1,1.), (1,1.)], vec![(1,1.)]  ]
///     ); 
///// ```
pub struct  ProductMatrixLazyMajorAscendSimplified < 
                    Matrix1, Matrix2,                
                    RingOperator,
                    OrderComparator,
                > 
            where 
                Matrix1:                            OracleMajorAscend + IndicesAndCoefficients,
                Matrix2:                            OracleMajorAscend + IndicesAndCoefficients< SnzVal = Matrix1::SnzVal, KeyMaj = Matrix1::KeyMin >, 
                Matrix1::ViewMajorAscend:           IntoIterator,
                Matrix2::ViewMajorAscend:           IntoIterator,
                Matrix1::ViewMajorAscendEntry:      KeyValGet < Matrix1::KeyMin, Matrix1::SnzVal >,
                Matrix2::ViewMajorAscendEntry:      KeyValGet < Matrix2::KeyMin, Matrix2::SnzVal >,  
                Matrix2::KeyMin:                    Clone,
                Matrix2::SnzVal:                    Clone,                          
                RingOperator:                       Clone + Semiring< Matrix1::SnzVal >,
                OrderComparator:                    Clone + StrictlyLess<  Matrix2::ViewMajorAscendEntry >,                                                           
{
    matrix_2:                   Matrix2,
    matrix_1:                   Matrix1, 
    ring_operator:              RingOperator,
    order_comparator:           OrderComparator,          
}

impl    < 
            Matrix1, 
            Matrix2,                
            RingOperator,
            OrderComparator,
        > 
    ProductMatrixLazyMajorAscendSimplified <
            Matrix1, 
            Matrix2,                
            RingOperator,
            OrderComparator,
        > 
    where 
        Matrix1:                            OracleMajorAscend + IndicesAndCoefficients,
        Matrix2:                            OracleMajorAscend + IndicesAndCoefficients< SnzVal = Matrix1::SnzVal, KeyMaj = Matrix1::KeyMin >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::ViewMajorAscendEntry:      KeyValGet < Matrix1::KeyMin, Matrix1::SnzVal >,
        Matrix2::ViewMajorAscendEntry:      KeyValGet < Matrix2::KeyMin, Matrix2::SnzVal >,  
        Matrix2::KeyMin:                    Clone,
        Matrix2::SnzVal:                    Clone,                          
        RingOperator:                       Clone + Semiring< Matrix1::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess<  Matrix2::ViewMajorAscendEntry >,    

{
    
    /// Generate a lazy lazy product of two matrix oracles
    /// 
    /// See documentation for [`ProductMatrixLazyMajorAscendSimplified`].  In that discussion, the function arguments `matrix_1` and `matrix_2`
    /// correspond to matrices `A` and `B`, respectively.
    pub fn new( 
                    matrix_1:           Matrix1, 
                    matrix_2:           Matrix2, 
                    ring_operator:      RingOperator, 
                    order_comparator:   OrderComparator 
                ) 
            -> Self {
        ProductMatrixLazyMajorAscendSimplified{
            matrix_1:                   matrix_1,
            matrix_2:                   matrix_2,
            ring_operator:              ring_operator,
            order_comparator:           order_comparator,         
        }
    }
  
}

// IndicesAndCoefficients
impl     < Matrix1, Matrix2, RingOperator, OrderComparator, > 

    IndicesAndCoefficients for

    ProductMatrixLazyMajorAscendSimplified< Matrix1, Matrix2, RingOperator, OrderComparator, > 

    where
        Matrix1:                            OracleMajorAscend + IndicesAndCoefficients,
        Matrix2:                            OracleMajorAscend + IndicesAndCoefficients< SnzVal = Matrix1::SnzVal, KeyMaj = Matrix1::KeyMin >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::ViewMajorAscendEntry:      KeyValGet < Matrix1::KeyMin, Matrix1::SnzVal >,
        Matrix2::ViewMajorAscendEntry:      KeyValGet < Matrix2::KeyMin, Matrix2::SnzVal >,  
        Matrix2::KeyMin:                    Clone,
        Matrix2::SnzVal:                    Clone,                          
        RingOperator:                       Clone + Semiring< Matrix1::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess<  Matrix2::ViewMajorAscendEntry >,         
{
    type KeyMin = Matrix2::KeyMin; type KeyMaj = Matrix1::KeyMaj; type SnzVal = Matrix1::SnzVal;
}   

// OracleMajorAscend
impl     < 
        Matrix1, 
        Matrix2, 
        RingOperator,
        OrderComparator,
        > 

    OracleMajorAscend for
    
    ProductMatrixLazyMajorAscendSimplified<
            Matrix1, 
            Matrix2, 
            RingOperator,
            OrderComparator,
        >     
    where
        Matrix1:                            OracleMajorAscend + IndicesAndCoefficients,
        Matrix2:                            OracleMajorAscend + IndicesAndCoefficients< SnzVal = Matrix1::SnzVal, KeyMaj = Matrix1::KeyMin >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::ViewMajorAscendEntry:      KeyValGet < Matrix1::KeyMin, Matrix1::SnzVal >,
        Matrix2::ViewMajorAscendEntry:      KeyValGet < Matrix2::KeyMin, Matrix2::SnzVal > + KeyValSet < Matrix2::KeyMin, Matrix2::SnzVal >,  
        Matrix2::KeyMin:                    Clone + PartialEq, // PartialEq is required by the struct that simplifies sparse vector iterators; it has to be able to compare the indices of different entries
        Matrix2::SnzVal:                    Clone,                          
        RingOperator:                       Clone + Semiring< Matrix1::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess<  Matrix2::ViewMajorAscendEntry >,                 

{   
    type ViewMajorAscend            =   LinearCombinationSimplified
                                            < Matrix2::ViewMajorAscendIntoIter, Matrix2::KeyMin, Matrix2::SnzVal, RingOperator, OrderComparator >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;
    type ViewMajorAscendEntry       =   Matrix2::ViewMajorAscendEntry;

    fn view_major_ascend( & self, index: Self::KeyMaj ) 
        -> 
        LinearCombinationSimplified< Matrix2::ViewMajorAscendIntoIter, Matrix2::KeyMin, Matrix2::SnzVal, RingOperator, OrderComparator >
    {

        vector_matrix_multiply_major_ascend_simplified( 
                self.matrix_1.view_major_ascend( index ),
                OracleRef::new( & self.matrix_2 ),
                self.ring_operator.clone(),
                self.order_comparator.clone(),
            )   

    }

}




































//  =========================================================================================================
//  TESTS
//  =========================================================================================================

//  ---------------------------------------------------------------------------
//  DOC-TEST DRAFTS
//  ---------------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_tests {
    use itertools::Itertools;

    use crate::utilities::partial_order::OrderComparatorFnWrapper;


    #[test]
    fn doc_test_vector_matrix_product_major_ascend_unsimplified() {
        use crate::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::partial_order::OrderComparatorAutoLtByKey;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0usize, 1), (1usize, 1) ];
        
        // the product
        let product     =   vector_matrix_multiply_major_ascend_unsimplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderComparatorAutoLtByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should return the following sequence of entries: (0,1), (0,1), (1,1), (1,6)   
        itertools::assert_equal( product, vec![(0isize,1usize), (0,1), (1,1), (1,6)] )
    }     

    #[test]
    fn doc_test_vector_matrix_product_major_ascend_simplified() {
        use crate::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::partial_order::OrderComparatorAutoLtByKey;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0usize, 1), (1usize, 1) ];
        
        // the product
        let product     =   vector_matrix_multiply_major_ascend_simplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderComparatorAutoLtByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should equal vector [2 0]        
        itertools::assert_equal( product, vec![(0isize, 2)] )
    }    

    #[test]
    fn doc_test_vector_matrix_product_minor_descend_unsimplified() {
        use crate::matrices::operations::multiply::vector_matrix_multiply_minor_descend_unsimplified;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::partial_order::OrderComparatorAutoAnyType;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0isize, 1), (1isize, 1) ];
        
        // the product
        let mut product     =   vector_matrix_multiply_minor_descend_unsimplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderComparatorAutoAnyType, // this compares tuples according to their lexicographic order
                            );
        // the product should return the following entries in sequence: (1,1), (1,1), (0,1), (0,6)       
        itertools::assert_equal( product, vec![(1,1), (1,1), (0,6), (0,1)] )
    }

    #[test]
    fn doc_test_vector_matrix_product_minor_descend_simplified() {
        use crate::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::partial_order::OrderComparatorAutoLtByKey;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVecSimple::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0isize, 1), (1isize, 1) ];
        
        // the product
        let product     =   vector_matrix_multiply_minor_descend_simplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderComparatorAutoLtByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should equal vector [0 2]        
        itertools::assert_equal( product, vec![(1usize, 2)] )
    }
    
    







    #[test]
    fn doc_test_matrix_product_lazy_major_ascend_unsimplified() {

        // Import the necessary traits and structs
        use crate::rings::operator_structs::ring_native::DivisionRingNative;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_oracle_traits::{OracleMajorAscend};
        use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;
        use crate::entries::KeyValGet;
        use std::iter::FromIterator;

        // Define the operator for the coefficient ring
        let ring_operator = DivisionRingNative::<f64>::new();     

        // Define matrix A, a 2x2 matrix with 1's above the diagonal and 0's below
        let matrix_1   =   VecOfVecSimple::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
        // Define matrix B, a 2x2 matrix with 1's above the diagonal and 0's below
        let matrix_2   =   VecOfVecSimple::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
        // Define references to the two arrays 
        let matrix_1_ref   =   &matrix_1;
        let matrix_2_ref    =   &matrix_2;        
        // Define the lazy product of A and B
        let product     =   ProductMatrixLazyMajorAscendSimplified::new( 
                                    matrix_1_ref,                 // matrix A
                                    matrix_2_ref,                 // matrix B
                                    ring_operator.clone(),      // ring operator
                                    OrderComparatorFnWrapper::new( |x:&(i32, f64), y:&(i32, f64)| x.key() < y.key() )   // total order on indices
                                );
                            
        // First, check the output iterators themselves (these are wrapped in `Simplify` structs).
        // ---------------------------------------------------------------------------------------
                            
        // Use this object to create a vector-of-vectors that represents the product
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .view_major_ascend(x)
                                                                )
                                                        )
                                            );
                                        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,2.)], vec![(1,1.)]  ]
            ); 
        
        // Second, check the underlying, unsimplified iterators (these are are contained inside the `Simplify` structs).
        // -------------------------------------------------------------------------------------------------------------
        
        // Use this object to create a vector-of-vectors that represents the product
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .view_major_ascend(x)
                                                                    .linear_combination_simplified // this unwraps the wrapper around the simplified linear combination
                                                                    .unsimplified // this unwraps the inner iterator
                                                                )
                                                        )
                                            );
                                        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,1.), (1,1.)], vec![(1,1.)]  ]
            ); 
    }

}



//  ---------------------------------------------------------------------------
//  TEST
//  ---------------------------------------------------------------------------


#[cfg(test)]
mod tests {
    use std::iter::FromIterator;

    use itertools::Itertools;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::matrices::matrix_types::vec_of_vec::{VecOfVecSimple};
    use crate::rings::operator_structs::ring_native::{SemiringNative, DivisionRingNative};
    use crate::utilities::partial_order::{OrderComparatorAutoLtByKey, OrderComparatorFnWrapper};


    // Test that lazy MATRIX-MATRIX multiplication works properly.
    #[test]
    pub fn matrix_by_matrix_multiply_major_ascend_test_1() {

        // Define the operator for the coefficient ring
        let ring_operator = DivisionRingNative::<f64>::new();     
        
        // Define the factor A, a 2x2 matrix with 1's above the diagonal and 0's below.
        let matrix_1   =   VecOfVecSimple::new(
                                                    // MajorDimension::Row,
                                             vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]
                                                );

        // Define the factor B, a 2x2 matrix with 1's above the diagonal and 0's below.
        let matrix_2   =   VecOfVecSimple::new(
                                                    // MajorDimension::Row,
                                             vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]
                                                );
        let matrix_1_ref    =   &matrix_1;
        let matrix_2_ref    =   &matrix_2;        

        // Define the lazy product of A and B
        let product     =   ProductMatrixLazyMajorAscendSimplified::new( 
                                                                        matrix_1_ref, // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                                        matrix_2_ref, 
                                                                        ring_operator.clone(), 
                                                                        OrderComparatorAutoLtByKey::new()
                                                                    );

        // CHECK THE MAJOR ITERATOR ITSELF

        // Use this object to create a vector-of-vectors that represents the product 
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .view_major_ascend(x)
                                                                )
                                                        )
                                            );
        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,2.)], vec![(1,1.)]  ]
            );     

        // CHECK THE UNDERLYING **UNSIMPLIFIED** ITERATOR   
        
        // Use this object to create a vector-of-vectors that represents the product 
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .view_major_ascend(x)
                                                                    .unwrap_simplified_lin_comb()
                                                                    .unsimplified
                                                                )
                                                        )
                                            );
        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,1.), (1,1.)], vec![(1,1.)]  ]
            );           
        
    }


    // Test that lazy MATRIX-MATRIX multiplication works properly - another way.
    #[test]
    pub fn matrix_by_matrix_multiply_major_ascend_test_2() {

        // Define the operator for the coefficient ring
        let ring_operator = SemiringNative::<i32>::new();     
        
        // Define factor A 
        let matrix_1 =   
            VecOfVecSimple::new(
                    vec![ 
                            vec![                  ],
                            vec![ (0,0)            ],                                                            
                            vec![ (0,4),   (1,-1)  ], 
                            vec![ (0,5),   (1,7)   ], 
                            vec![ (0,11),  (1,13)  ],                                                         
                    ]
                );
        // Define factor B 
        let matrix_2 =    
            VecOfVecSimple::new(
                    vec![                                                                                                                           
                            vec![ (0,1),   (1,2),   (2,3) ], 
                            vec![ (0,4),   (1,5),   (2,6) ],
                            vec![                         ],  
                            vec![ (0,0)                   ],                                                              
                        ]
                );
        let matrix_1_ref    =   &matrix_1;
        let matrix_2_ref    =   &matrix_2;  

        // This is the *actual* product A*B, if A and B are row-major.
        let product_true   =   
                    vec![ 
                            vec![                                   ], 
                            vec![                                   ],                                                             
                            vec![          (1,3),     (2,6)         ], 
                            vec![ (0,33),  (1,45),    (2,57)        ],
                            vec![ (0,63),  (1,87),    (2,111)       ],                                                             
                        ];                                              

        // Define the lazy product of A and B
        let product     =   ProductMatrixLazyMajorAscendSimplified::new( 
                                                                matrix_1_ref,  // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                                matrix_2_ref, 
                                                                ring_operator.clone(), 
                                                                OrderComparatorFnWrapper::new(|x:&(i32, i32), y:&(i32, i32)| x.key() < y.key() )  
                                                        );


        // Check that each row of the lazy MATRIX product agrees with each row of the true product.
        for row_index in 0 .. 5 {
            let vec_lazy   =   product
                                                .view_major_ascend( row_index )
                                                .collect_vec();
            let vec_true =  product_true[ row_index ].clone();
            assert_eq!( vec_lazy, vec_true )
        }   
        
    }  
    

   // Test that lazy MATRIX-VECTOR multiplication works properly.
   #[test]
   pub fn matrix_by_vector_multiply_major_ascend_test_2() {

       // Define the operator for the coefficient ring
       let ring_operator = SemiringNative::<i32>::new();     
       
       // Define factor A 
       let matrix_1 =   
           VecOfVecSimple::new(
                   vec![ 
                           vec![                  ],
                           vec![ (0,0)            ],                                                            
                           vec![ (0,4),   (1,-1)  ], 
                           vec![ (0,5),   (1,7)   ], 
                           vec![ (0,11),  (1,13)  ],                                                         
                   ]
               );
       // Define factor B 
       let matrix_2 =    
           VecOfVecSimple::new(
                   vec![                                                                                                                           
                           vec![ (0,1),   (1,2),   (2,3) ], 
                           vec![ (0,4),   (1,5),   (2,6) ],
                           vec![                         ],  
                           vec![ (0,0)                   ],                                                              
                       ]
                );
       // This is the *actual* product A*B, if A and B are row-major.
       let product_true   =   
            vec![ 
                    vec![                                   ], 
                    vec![                                   ],                                                             
                    vec![          (1,3),     (2,6)         ], 
                    vec![ (0,33),  (1,45),    (2,57)        ],
                    vec![ (0,63),  (1,87),    (2,111)       ],                                                             
                ];                                                                             


       // Check that each row of the lazy MATRIX product agrees with each row of the true product.
       for row_index in 0 .. 5 {
           let vec_lazy   =   vector_matrix_multiply_major_ascend_simplified(
                                                    (& matrix_1).view_major_ascend( row_index ).clone(),               
                                                    & matrix_2,  // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                    ring_operator.clone(),
                                                    OrderComparatorFnWrapper::new( |x:&(i32, i32), y:&(i32, i32)| x.key() < y.key()  )
                                                )
                                               .collect_vec();
           let vec_true =  product_true[ row_index ].clone();
           assert_eq!( vec_lazy, vec_true )
       }   
       
   } 

}