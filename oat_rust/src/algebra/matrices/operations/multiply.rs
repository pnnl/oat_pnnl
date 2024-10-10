//! Matrix-matrix and matrix-vector multiplication.
//! 
//! 
//! - to multiply a matrix with a matrix, use 
//!   - [ProductMatrix::new](crate::algebra::matrices::operations::multiply::ProductMatrix)
// //!   - [ProductMatrixLazyMajorAscendUnsimplified::new](crate::algebra::matrices::operations::multiply::ProductMatrixLazyMajorAscendUnsimplified)
//! - to multiply a matrix with a vector, use either the convenient [VectorOperations](crate::algebra::vectors::operations::VectorOperations) trait methods
//!   - [multiply_matrix](crate::algebra::vectors::operations::VectorOperations::multiply_matrix)
//!   - [multiply_matrix_major_ascend](crate::algebra::vectors::operations::VectorOperations::multiply_matrix_major_ascend)
//!   - [multiply_matrix_minor_descend](crate::algebra::vectors::operations::VectorOperations::multiply_matrix_minor_descend)
//!     
//!     or the following
//!   - [vector_matrix_multiply_major_ascend_simplified](    crate::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified)
//!   - [vector_matrix_multiply_major_ascend_unsimplified](  crate::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified)
//!   - [vector_matrix_multiply_minor_descend_simplified](   crate::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified)
//!   - [vector_matrix_multiply_minor_descend_unsimplified]( crate::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_unsimplified)



use crate::{algebra::matrices::{query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients}, }, algebra::rings::operator_traits::{Semiring}, utilities::{iterators::merge::hit::{hit_merge_by_predicate}, order::{JudgePartialOrder, ReverseOrder}}, algebra::vectors::{operations::{VectorOperations, LinearCombinationSimplified, LinearCombinationUnsimplified},}};
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet};



//  MATRIX - VECTOR MULTIPLICATION
//  ===========================================================================


//  MAJOR ASCEND
//  ---------------------------------------------------------------------------

/// Returns the product of a matrix with a vector.
/// 
/// The resulting iterator is a linear combination of major views.  It returns entries in ascending (but **non-strictly** ascending) order of index, provided that
/// [ViewRowAscend](crate::algebra::matrices::query::ViewRowAscend) returns entries in 
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
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// use oat_rust::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0usize, 1), (1usize, 1) ];
/// 
/// // the product
/// let product     =   vector_matrix_multiply_major_ascend_unsimplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
///                     );
/// // the product should return the following sequence of entries: (0,1), (0,1), (1,1), (1,6)   
/// itertools::assert_equal( product, vec![(0isize,1usize), (0,1), (1,1), (1,6)] )
/// ```
pub fn vector_matrix_multiply_major_ascend_unsimplified 
            < Matrix, RingOperator, SparseVecIter, OrderOperator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_operator:   OrderOperator,
            ) 
            ->
            LinearCombinationUnsimplified
                < Matrix::ViewMajorAscendIntoIter, Matrix::ColIndex, Matrix::Coefficient, RingOperator, OrderOperator >            
    where 
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:   KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        Matrix::Coefficient:                 Clone, // Clone is required for HitMerge
        OrderOperator:                Clone + JudgePartialOrder<  Matrix::EntryMajor >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < Matrix::RowIndex, Matrix::Coefficient >,
{
    // LinearCombinationUnsimplified{
    //     linear_combination_unsimplified:
            hit_merge_by_predicate(
                sparse_vec
                        .into_iter()
                        .map(   |x| 
                                matrix
                                    .view_major_ascend(x.key() )
                                    .into_iter()
                                    .scale( x.val(), ring_operator.clone() )
                        ),
                    order_operator
            )        
    // }
}


/// Returns the product of a matrix with a vector; the result is a linear combination of the 
/// matrix's major views.
/// 
/// The resulting iterator is a simplified linear combination of major views.  It returns entries in **strictly** ascending order of index, provided that
/// [ViewRowAscend](crate::algebra::matrices::query::ViewRowAscend) returns entries in 
/// ascending order according to index.  In particular, no two consequtive indices have the same
/// index.
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `vA`, a linear combination of **rows** of `A`.
/// If instead the representation is column-major, then this function returns `Av`.
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0usize, 1), (1usize, 1) ];
/// 
/// // the product
/// let product     =   vector_matrix_multiply_major_ascend_simplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
///                     );
/// // the product should equal vector [2 0]        
/// itertools::assert_equal( product, vec![(0isize, 2)] )
/// ```
pub fn vector_matrix_multiply_major_ascend_simplified 
            < Matrix, RingOperator, SparseVecIter, OrderOperator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_operator:   OrderOperator,
            ) 
            ->
            LinearCombinationSimplified
                < Matrix::ViewMajorAscendIntoIter, Matrix::ColIndex, Matrix::Coefficient, RingOperator, OrderOperator >            
    where 
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:   KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,
        Matrix::Coefficient:                 Clone, // Clone is required for HitMerge
        Matrix::ColIndex:                 PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros        
        OrderOperator:                Clone + JudgePartialOrder<  Matrix::EntryMajor >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < Matrix::RowIndex, Matrix::Coefficient >,
{
    vector_matrix_multiply_major_ascend_unsimplified( 
            sparse_vec,
            matrix,
            ring_operator.clone(),
            order_operator,
        )
        // .linear_combination_unsimplified
        .simplify(ring_operator)
}        



//  MINOR DESCEND
//  ---------------------------------------------------------------------------


/// Returns the product of a matrix with a vector; the resulting iterator is an unsimplified linear combination of the matrix's minor views.
/// 
/// The resulting iterator returns entries in descending (but **non-strictly** descending) order of index, provided that
/// [ViewRowAscend](crate::algebra::matrices::query::ViewRowDescend) returns entries in 
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
/// use oat_rust::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_unsimplified;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::order::OrderOperatorAuto;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0isize, 1), (1isize, 1) ];
/// 
/// // the product
/// let mut product     =   vector_matrix_multiply_minor_descend_unsimplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderOperatorAuto, // this compares tuples according to their lexicographic order
///                     );
/// // the product should return the following entries in sequence: (1,1), (1,1), (0,1), (0,6)       
/// itertools::assert_equal( product, vec![(1,1), (1,1), (0,6), (0,1)] )
/// ```
pub fn vector_matrix_multiply_minor_descend_unsimplified 
            < Matrix, RowIndex, ColIndex, Coefficient, RingOperator, SparseVecIter, OrderOperator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_operator:   OrderOperator,
            ) 
            ->
            LinearCombinationUnsimplified
                < Matrix::ViewMinorDescendIntoIter, RowIndex, Coefficient, RingOperator, ReverseOrder< OrderOperator > >            
    where 
        Matrix:                             ViewColDescend + IndicesAndCoefficients< ColIndex = ColIndex >,
        Matrix::ViewMinorDescendIntoIter:   Iterator,        
        < Matrix::ViewMinorDescendIntoIter as Iterator >::Item:     KeyValSet< RowIndex, Coefficient >,
        RingOperator:                       Clone + Semiring< Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        Coefficient:                             Clone, // Clone is required for HitMerge
        OrderOperator:                    Clone + JudgePartialOrder<  <Matrix::ViewMinorDescendIntoIter as Iterator>::Item >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                      IntoIterator,
        SparseVecIter::Item:                KeyValGet < ColIndex, Coefficient >,
{
    // LinearCombinationUnsimplified{
    //     linear_combination_unsimplified:
            hit_merge_by_predicate(
                sparse_vec
                        .into_iter()
                        .map(   |x| 
                                matrix
                                    .view_minor_descend(x.key() )
                                    .into_iter()
                                    .scale( x.val(), ring_operator.clone() )
                        ),
                    ReverseOrder::new( order_operator )
            )        
    // }
}


/// Returns the product of a matrix with a vector; the result is a linear combination of the matrix's minor views.
/// 
/// The resulting iterator is a simplified linear combinatino of minor views.  It returns entries in **strictly** descending order of index, provided that
/// [ViewRowAscend](crate::algebra::matrices::query::ViewRowAscend) returns entries in 
/// ascending order according to index.  In particular, no two consequtive indices have the same
/// index.
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `A v`, a linear combination of **columns** of `A`.
/// If instead the representation is column-major, then this function returns `vA`.
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// 
/// // a row-major matrix representation with usize major keys and isize minor keys
/// let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                                    vec![ (0isize, 1), (1isize, 1) ],     ]     );
/// // the vector [1, 1]
/// let vector      =   vec![ (0isize, 1), (1isize, 1) ];
/// 
/// // the product
/// let product     =   vector_matrix_multiply_minor_descend_simplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderFieldOperator::new(7), // the finite field of order 7
///                         OrderOperatorByKey::new(), // this compares order of entries according to their indices
///                     );
/// // the product should equal vector [0 1]        
/// itertools::assert_equal( product, vec![(1usize, 2)] )
/// ```
pub fn vector_matrix_multiply_minor_descend_simplified 
            < Matrix, RowIndex, ColIndex, Coefficient, RingOperator, SparseVecIter, OrderOperator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_operator:   OrderOperator,
            ) 
            ->
            LinearCombinationSimplified
                < Matrix::ViewMinorDescendIntoIter, RowIndex, Coefficient, RingOperator, ReverseOrder< OrderOperator > >

    where 
        Matrix:                         ViewColDescend< ColIndex = ColIndex, RowIndex = RowIndex >,
        Matrix::ViewMinorDescend:       IntoIterator,
        Matrix::EntryMinor:  KeyValGet < RowIndex, Coefficient > + KeyValSet < RowIndex, Coefficient >,
        RingOperator:                   Clone + Semiring< Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        Coefficient:                         Clone, // Clone is required for HitMerge
        RowIndex:                         std::cmp::PartialEq, // required by the struct that simplifies vectors (it has to compare the indices of different entries)
        OrderOperator:                Clone + JudgePartialOrder<  < Matrix::ViewMinorDescend as IntoIterator >::Item >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < ColIndex, Coefficient >
{

    vector_matrix_multiply_minor_descend_unsimplified( 
            sparse_vec,
            matrix,
            ring_operator.clone(),
            order_operator,
        )
        // .linear_combination_unsimplified
        .simplify(ring_operator)
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
/// have the same coefficient, then (iii) wraps `J` in a [`Simplify`](`crate::algebra::vectors::transfomations::Simplify`) struct that
/// *does* combine terms and drop zeros.
/// 
/// When the user requests a *column*, the lazy object first (i) scales 
/// some columns of `B` by an appropriate factor, then (ii) merges those columns into a single iterator, `J`.  Iterator `J` returns
/// entries in sorted order (ascending), but does *not* simplify terms: in particular, it does not drop zero terms or combine terms that 
/// have the same coefficient, then (iii) wraps `J` in a [`Simplify`](`crate::algebra::vectors::transfomations::Simplify`) struct that
/// *does* combine terms and drop zeros.
/// 
/// **Major dimension matters**. The major dimension is the dimension you 
/// actually get when you call `A.view_major_ascend(..)`.  The reason for this is the following pair of facts,
/// (both facts follow from the prodecure by which rows/columns are constructed, as described above).
/// * If `A` and `B` are row-major, then the lazy product `multiply_matrix_major_ascend(A, B, ..)` represents `A * B`
/// * If `A` and `B` are col-major, then the lazy product `multiply_matrix_major_ascend(A, B, ..)` represents `B * A`
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
/// // ---------------------------------------------------------------------------------------
/// 
/// use oat_rust::algebra::rings::operator_structs::ring_native::DivisionRingNative;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::{ViewRowAscend};
/// use oat_rust::algebra::matrices::operations::multiply::ProductMatrix;
/// use oat_rust::algebra::vectors::entries::KeyValGet;
/// use oat_rust::utilities::order::OrderOperatorByLessThan;
/// use std::iter::FromIterator;
/// 
/// // Initialize variables
/// // ---------------------------------------------------------------------------------------
/// 
/// // Define the operator for the coefficient ring
/// let ring_operator = DivisionRingNative::<f64>::new();     
/// // Define matrix A, a 2x2 matrix with 1's above the diagonal and 0's below
/// let matrix_1   =   VecOfVec::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
/// // Define matrix B, a 2x2 matrix with 1's above the diagonal and 0's below
/// let matrix_2   =   VecOfVec::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
/// // Define references to the two arrays (recall that ViewRowAscend is only implemented on `& VecOfVec`, not on `VecOfVec`)
/// let matrix_1_ref   =   &matrix_1;
/// let matrix_2_ref    =   &matrix_2;   
///      
/// // Define the lazy product of A and B
/// // ---------------------------------------------------------------------------------------
/// 
/// let product     =   ProductMatrix::new( 
///                             matrix_1_ref,                 // matrix A
///                             matrix_2_ref,                 // matrix B
///                             ring_operator.clone(),      // ring operator
///                             OrderOperatorByLessThan::new( |x:&(i32, f64), y:&(i32, f64)| x.key() < y.key() )   // defines a total order on indices
///                         );
///                     
/// // Check the output iterators themselves (these are wrapped in `Simplify` structs).
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
/// // Check the underlying, unsimplified iterators (these are are contained inside the `Simplify` structs).
/// // -------------------------------------------------------------------------------------------------------------
/// 
/// // Use this object to create a vector-of-vectors that represents the product
/// let output =   Vec::from_iter(
///                                         (0..2)
///                                             .map(   |x| 
///                                                     Vec::from_iter(
///                                                             product
///                                                             .view_major_ascend(x)
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
/// ```
/// 
/// # Limitations
/// 
/// This struct cannot return descending minor views of the product because it contains no order operator for the entries of minor views.
#[derive(Clone, Copy, Debug)]
pub struct ProductMatrix < 
                    Matrix1, 
                    Matrix2,                
                    RingOperator,
                    OrderOperator,
                > 
            where 
                Matrix1:                            ViewRowAscend + IndicesAndCoefficients,
                Matrix2:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Matrix1::Coefficient, RowIndex = Matrix1::ColIndex >, 
                Matrix1::ViewMajorAscend:           IntoIterator,
                Matrix2::ViewMajorAscend:           IntoIterator,
                Matrix1::EntryMajor:      KeyValGet < Matrix1::ColIndex, Matrix1::Coefficient >,
                Matrix2::EntryMajor:      KeyValGet < Matrix2::ColIndex, Matrix2::Coefficient >,  
                Matrix2::ColIndex:                    Clone,
                Matrix2::Coefficient:                    Clone,                          
                RingOperator:                       Clone + Semiring< Matrix1::Coefficient >,
                OrderOperator:                    Clone + JudgePartialOrder<  Matrix2::EntryMajor >,                                                           
{
    matrix_1:                   Matrix1,     
    matrix_2:                   Matrix2,
    ring_operator:              RingOperator,
    order_operator:           OrderOperator,          
}

impl    < 
            Matrix1, 
            Matrix2,                
            RingOperator,
            OrderOperator,
        > 
    ProductMatrix <
            Matrix1, 
            Matrix2,                
            RingOperator,
            OrderOperator,
        > 
    where 
        Matrix1:                            ViewRowAscend + IndicesAndCoefficients,
        Matrix2:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Matrix1::Coefficient, RowIndex = Matrix1::ColIndex >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::EntryMajor:                KeyValGet < Matrix1::ColIndex, Matrix1::Coefficient >,
        Matrix2::EntryMajor:                KeyValGet < Matrix2::ColIndex, Matrix2::Coefficient >,  
        Matrix2::ColIndex:                  Clone,
        Matrix2::Coefficient:              Clone,                          
        RingOperator:                       Clone + Semiring< Matrix1::Coefficient >,
        OrderOperator:                      Clone + JudgePartialOrder<  Matrix2::EntryMajor >,    

{
    
    /// Generate a lazy lazy product of two matrix oracles
    /// 
    /// See documentation for [`ProductMatrix`].  In that discussion, the function arguments `matrix_1` and `matrix_2`
    /// correspond to matrices `A` and `B`, respectively.
    pub fn new( 
                    matrix_1:           Matrix1, 
                    matrix_2:           Matrix2, 
                    ring_operator:      RingOperator, 
                    order_operator:     OrderOperator 
                ) 
            -> Self {
        ProductMatrix{
            matrix_1,
            matrix_2,
            ring_operator,
            order_operator,         
        }
    }
  
}

// IndicesAndCoefficients
impl     < Matrix1, Matrix2, RingOperator, OrderOperator, > 

    IndicesAndCoefficients for

    ProductMatrix< Matrix1, Matrix2, RingOperator, OrderOperator, > 

    where
        Matrix1:                            ViewRowAscend + IndicesAndCoefficients,
        Matrix2:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Matrix1::Coefficient, RowIndex = Matrix1::ColIndex >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::EntryMajor:      KeyValGet < Matrix1::ColIndex, Matrix1::Coefficient >,
        Matrix2::EntryMajor:      KeyValGet < Matrix2::ColIndex, Matrix2::Coefficient >,  
        Matrix2::ColIndex:                    Clone,
        Matrix2::Coefficient:                    Clone,                          
        RingOperator:                       Clone + Semiring< Matrix1::Coefficient >,
        OrderOperator:                    Clone + JudgePartialOrder<  Matrix2::EntryMajor >,         
{
    type EntryMajor = Matrix2::EntryMajor;
    type EntryMinor = Matrix1::EntryMinor;
    type ColIndex = Matrix2::ColIndex; 
    type RowIndex = Matrix1::RowIndex; 
    type Coefficient = Matrix1::Coefficient;
}   

// ViewRowAscend
impl     < 
        Matrix1, 
        Matrix2, 
        RingOperator,
        OrderOperator,
        > 

    ViewRowAscend for
    
    ProductMatrix<
            Matrix1, 
            Matrix2, 
            RingOperator,
            OrderOperator,
        >     
    where
        Matrix1:                            ViewRowAscend + IndicesAndCoefficients,
        Matrix2:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Matrix1::Coefficient, RowIndex = Matrix1::ColIndex >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::EntryMajor:      KeyValGet < Matrix1::ColIndex, Matrix1::Coefficient >,
        Matrix2::EntryMajor:      KeyValGet < Matrix2::ColIndex, Matrix2::Coefficient > + KeyValSet < Matrix2::ColIndex, Matrix2::Coefficient >,  
        Matrix2::ColIndex:                    Clone + PartialEq, // PartialEq is required by the struct that simplifies sparse vector iterators; it has to be able to compare the indices of different entries
        Matrix2::Coefficient:                    Clone,                          
        RingOperator:                       Clone + Semiring< Matrix1::Coefficient >,
        OrderOperator:                    Clone + JudgePartialOrder<  Matrix2::EntryMajor >,                 

{   
    type ViewMajorAscend            =   LinearCombinationSimplified
                                            < Matrix2::ViewMajorAscendIntoIter, Matrix2::ColIndex, Matrix2::Coefficient, RingOperator, OrderOperator >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( & self, index: Self::RowIndex ) 
        -> 
        LinearCombinationSimplified< Matrix2::ViewMajorAscendIntoIter, Matrix2::ColIndex, Matrix2::Coefficient, RingOperator, OrderOperator >
    {

        vector_matrix_multiply_major_ascend_simplified( 
                self.matrix_1.view_major_ascend( index ),
                & self.matrix_2,
                self.ring_operator.clone(),
                self.order_operator.clone(),
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
    

    use crate::utilities::order::OrderOperatorByLessThan;


    #[test]
    fn doc_test_vector_matrix_product_major_ascend_unsimplified() {
        use crate::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::order::OrderOperatorByKey;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0usize, 1), (1usize, 1) ];
        
        // the product
        let product     =   vector_matrix_multiply_major_ascend_unsimplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should return the following sequence of entries: (0,1), (0,1), (1,1), (1,6)   
        itertools::assert_equal( product, vec![(0isize,1usize), (0,1), (1,1), (1,6)] )
    }     

    #[test]
    fn doc_test_vector_matrix_product_major_ascend_simplified() {
        use crate::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::order::OrderOperatorByKey;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0usize, 1), (1usize, 1) ];
        
        // the product
        let product     =   vector_matrix_multiply_major_ascend_simplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should equal vector [2 0]        
        itertools::assert_equal( product, vec![(0isize, 2)] )
    }    

    #[test]
    fn doc_test_vector_matrix_product_minor_descend_unsimplified() {
        use crate::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_unsimplified;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::order::OrderOperatorAuto;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0isize, 1), (1isize, 1) ];
        
        // the product
        let product     =   vector_matrix_multiply_minor_descend_unsimplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderOperatorAuto, // this compares tuples according to their lexicographic order
                            );
        // the product should return the following entries in sequence: (1,1), (1,1), (0,1), (0,6)       
        itertools::assert_equal( product, vec![(1,1), (1,1), (0,6), (0,1)] )
    }

    #[test]
    fn doc_test_vector_matrix_product_minor_descend_simplified() {
        use crate::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::order::OrderOperatorByKey;

        // a row-major matrix representation with usize major keys and isize minor keys
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     );
        // the vector [1, 1]
        let vector      =   vec![ (0isize, 1), (1isize, 1) ];
        
        // the product
        let product     =   vector_matrix_multiply_minor_descend_simplified(
                                vector,
                                & matrix,
                                PrimeOrderFieldOperator::new(7), // the finite field of order 7
                                OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should equal vector [0 2]        
        itertools::assert_equal( product, vec![(1usize, 2)] )
    }
    
    







    #[test]
    fn doc_test_matrix_product_lazy_major_ascend_unsimplified() {

        // Import the necessary traits and structs
        use crate::algebra::rings::operator_structs::ring_native::DivisionRingNative;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::{ViewRowAscend};
        use crate::algebra::matrices::operations::multiply::ProductMatrix;
        use crate::algebra::vectors::entries::KeyValGet;
        use std::iter::FromIterator;

        // Define the operator for the coefficient ring
        let ring_operator = DivisionRingNative::<f64>::new();     

        // Define matrix A, a 2x2 matrix with 1's above the diagonal and 0's below
        let matrix_1   =   VecOfVec::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
        // Define matrix B, a 2x2 matrix with 1's above the diagonal and 0's below
        let matrix_2   =   VecOfVec::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    );
        // Define references to the two arrays 
        let matrix_1_ref   =   &matrix_1;
        let matrix_2_ref    =   &matrix_2;        
        // Define the lazy product of A and B
        let product     =   ProductMatrix::new( 
                                    matrix_1_ref,                 // matrix A
                                    matrix_2_ref,                 // matrix B
                                    ring_operator,      // ring operator
                                    OrderOperatorByLessThan::new( |x:&(i32, f64), y:&(i32, f64)| x.key() < y.key() )   // total order on indices
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
    use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    use crate::algebra::rings::operator_structs::ring_native::{SemiringNative, DivisionRingNative};
    use crate::utilities::order::{OrderOperatorByKey, OrderOperatorByLessThan};


    // Test that lazy MATRIX-MATRIX multiplication works properly.
    #[test]
    pub fn matrix_by_matrix_multiply_major_ascend_test_1() {

        // Define the operator for the coefficient ring
        let ring_operator = DivisionRingNative::<f64>::new();     
        
        // Define the factor A, a 2x2 matrix with 1's above the diagonal and 0's below.
        let matrix_1   =   VecOfVec::new(
                                                    // MajorDimension::Row,
                                             vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]
                                                );

        // Define the factor B, a 2x2 matrix with 1's above the diagonal and 0's below.
        let matrix_2   =   VecOfVec::new(
                                                    // MajorDimension::Row,
                                             vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]
                                                );
        let matrix_1_ref    =   &matrix_1;
        let matrix_2_ref    =   &matrix_2;        

        // Define the lazy product of A and B
        let product     =   ProductMatrix::new( 
                                                                        matrix_1_ref, // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                                        matrix_2_ref, 
                                                                        ring_operator, 
                                                                        OrderOperatorByKey::new()
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
            VecOfVec::new(
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
            VecOfVec::new(
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
        let product     =   ProductMatrix::new( 
                                                                matrix_1_ref,  // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                                matrix_2_ref, 
                                                                ring_operator, 
                                                                OrderOperatorByLessThan::new(|x:&(i32, i32), y:&(i32, i32)| x.key() < y.key() )  
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
           VecOfVec::new(
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
           VecOfVec::new(
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
                                                    ring_operator,
                                                    OrderOperatorByLessThan::new( |x:&(i32, i32), y:&(i32, i32)| x.key() < y.key()  )
                                                )
                                               .collect_vec();
           let vec_true =  product_true[ row_index ].clone();
           assert_eq!( vec_lazy, vec_true )
       }   
       
   } 

}