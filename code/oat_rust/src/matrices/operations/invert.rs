

// -----------------------------------------------------------------------------------


//! Inverse of a triangular matrix; see also [`solve`](crate::matrices::operations::solve).
//! 
//! # How to improve performance
//! If you only need to solve `Ax = b` for `x`, where `A` is a sparse matrix, then you will probably get
//! better results by using the [`solve`](crate::matrices::operations::solve) package instead.  See the blog post
//! [don't invert that matrix](https://www.johndcook.com/blog/2010/01/19/dont-invert-that-matrix/), by John Cook,
//! for general thoughts on inverting sparse matrices.
//! 
//! # Examples
//! 
//! ```
//! use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
//! use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
//! use oat_rust::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
//! use oat_rust::utilities::partial_order::OrderComparatorFnWrapper;
//! use oat_rust::matrices::display::print_indexed_major_views;
//! 
//! // Define the ring operator (prime field of order 1049)
//! let modulus = 1049;
//! let ring_operator  =   PrimeOrderFieldOperator::new(modulus);                  
//! 
//! // Define the order comparator
//! // This struct has a method that assigns a value of true/false to the statement
//! // that `entry (a,b) should precede entry (c,d)`.  We want (a,b) to preced (c,d)
//! // if a < c.
//! let order_comparator = 
//!     OrderComparatorFnWrapper::new( |x: &(usize, usize), y: &(usize, usize)| x.0 < y.0 ) ;
//! 
//! // Define a row-major upper triangular amtrix
//! let matrix  =   VecOfVecSimple::new(
//!                     vec![   
//!                         vec![ (0,1), (1,1), ], 
//!                         vec![        (1,1), ],
//!                     ],
//!                 );
//! 
//! // Define the inverse
//! let inverse =   InverseOfTriangularArrayLazyAscending::new( 
//!                     & matrix, 
//!                     ring_operator.clone(),
//!                     order_comparator,
//!                 );   
//! 
//! // Print the inverse
//! let row_indices = vec![ 0, 1 ];
//! print_indexed_major_views( &inverse, row_indices );
//! ```
//! 
//! This should print the following:
//! 
//! ```bash
//! $ major_view 0: [(0, 1), (1, 1048)]
//! $ major_view 1: [(1, 1)] 
//! ```
//! 
//! which is the correct solution, since the inverse of `matrix` is
//! ```bash 
//!       |  1   -1  |
//!       |  0    1  |
//! ```
//! and -1 = 1048, modulo 1049.


use crate::{entries::{KeyValGet, KeyValSet}, rings::operator_traits::{Semiring, Ring, DivisionRing}, matrices::{matrix_oracle_traits::{OracleMajorAscend, IndicesAndCoefficients}, matrix_types::oracle_ref::OracleRef}, utilities::{iterators::{merge::heap_of_iterators::{HitMerge, hit_bulk_insert, hit_merge_by_predicate}, general::PeekUnqualified}, partial_order::StrictlyLess}, vectors::operations::{Transforms, Scale, Simplify}};
use std::marker::PhantomData;



//  ---------------------------------------------------------------------------
//  INVERT A TRIANGULAR ARRAY


/// An iterator that represents a major ascending view of the inverse of a triangular array.
/// 
/// This struct will only iterate over the correct entries if the initial value of `next_entry_of_inv` and `entries_to_elim_simplified_heap`
/// are set correctly when the iterator is created. 
/// 
/// # Developer notes
/// 
/// This object is returned by `view_major_ascend( .. )` for [`InverseOfTriangularArrayLazyAscending`].
/// Alternatively, one could implement a different oracle for the same matrix (that is, for the inverse of a
/// triangular array) where `view_major_ascend( .. )` returns an object of type [`TriangularSolveAscend`](crate::matrices::operations::triangular_solve::TriangularSolveAscend).
/// In that implementation, `view_major_ascend( .. )` would simply return `TriangularSolveAscend::new( x, matrix, ring_operator, order_comparator)`,
/// where `x` is an entry iterator representing a standard unit vector.  However, there are some small differences between that approach and the one
/// implemented here, namely:
/// - the alternate approach wraps a (possibly large) number of iterators in wrappers of type [`TwoTypeIterator`](crate::utilities::iterators::general::IterTwoType)
/// - the approach which is actually applied here has a `head-tail` structure which allows peeking
/// 
/// It would be interesting to contrast the performance of these two approaches, to see if there is a meaningful difference.
pub struct ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                    Matrix,
                    RingOperator,
                    OrderComparator,
                > 
    where 
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,        
        Matrix::KeyMaj:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,                   

{
    ring_operator:                  RingOperator, // the operator for the coefficient ring_operator    
    matrix:                         Matrix, // the matrix to be inverted    
    next_entry_of_inv:              Option< Matrix::ViewMajorAscendEntry  >, // the next entry in this view of the inverse matrix     
    entries_to_elim_simplified_heap:  Simplify<
                                            HitMerge< 
                                                    // the HitMerge struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                    Scale< 
                                                            Matrix::ViewMajorAscendIntoIter,  // the iterators returned by the matrix
                                                            Matrix::KeyMaj,
                                                            RingOperator, // the ring_operator operator
                                                            Matrix::SnzVal,
                                                        >,
                                                    // the thing that declares whether one major key comes before of after another    
                                                    OrderComparator 
                                                >,
                                            Matrix::KeyMaj,
                                            RingOperator,
                                            Matrix::SnzVal,
                                        >, 
}                                   

impl    < 
        Matrix,
        RingOperator,
        OrderComparator,
        > 

        Iterator for    

        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderComparator,
            > 

    where 
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,        
        Matrix::KeyMaj:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               

{

    type Item = Matrix::ViewMajorAscendEntry;

    fn next( &mut self ) -> Option<Self::Item> { 

        match self.next_entry_of_inv.take() {
            // RETURN NONE IF `next_entry_of_inv` IS EMPTY
            None => None,  

            // OTHERWISE `next_entry_of_inv` IS NOT EMPTY, AND WE DO THE FOLLOWING
            Some( return_value ) => {

                // IF THE HEAD OF HEADTAIL IS NONEMPTY, IT BECOMES THE NEXT TARGET FOR ELIMINATION
                if let Some( mut entry_to_eliminate )       =   self.entries_to_elim_simplified_heap.next() {

                    // TO ELIMINATE THIS ENTRY, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX WE WANT TO INVERT; DENOTE THE UNSCLAED VECTOR BY v
                    let mut seed_of_eliminating_iterator    =   self.matrix.view_major_ascend( entry_to_eliminate.key() ).into_iter();
                    
                    // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
                    let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

                    // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
                    let scale_factor     =    self.ring_operator.negate(
                                                            self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                        );
                    
                    // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
                    // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_elim_simplified_heap
                    let eliminating_iterator    =   seed_of_eliminating_iterator.scale( scale_factor.clone(), self.ring_operator.clone() );

                    // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_elim_simplified_heap`
                    hit_bulk_insert( &mut self.entries_to_elim_simplified_heap.unsimplified, vec![eliminating_iterator] ); 

                    // NOW RE-USE `entry_to_eliminate` AS A CONTAINER FOR THE NEXT ENTRY IN THE INVERSE MATRIX
                    entry_to_eliminate.set_val( scale_factor );

                    // THE NEXT NONZERO ENTRY ON THE INVERSE MATRIX (AFTER THE ONE WE ARE ABOUT TO RETURN, NAMELY return_value) HAS INDEX EQUAL TO `entry_to_eliminate.key()` AND COEFFICIENT EQUAL TO `scale_factor`
                    self.next_entry_of_inv = Some( entry_to_eliminate );

                }

                // RETURN THE NEXT ENTRY OF THE INVERSE MATRIX
                // we could have done this much earlier, in principle, but we had to set up the next value of `self.next_entry_of_inv` before we returned `return_value`
                return Some( return_value )

            }
        }
    }
}         

impl    < 
        Matrix,
        RingOperator,
        OrderComparator,
        > 

        PeekUnqualified for    

        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderComparator,
            > 

    where 
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,        
        Matrix::KeyMaj:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               


{

    fn peek_unqualified(&mut self) -> std::option::Option<&<Self as Iterator>::Item> { 
        match &self.next_entry_of_inv {
            None => { None },
            Some(x) => { Some( x ) }
        }
    }
}

/// Given an invertible upper triangular row-major matrix `A`, returns a major view of the inverse of `A`.
pub fn view_major_of_inverse_of_triangular_array_lazy_ascend < 
            Matrix,
            RingOperator,
            OrderComparator,
        > 
        (   keymaj: Matrix::KeyMaj, 
            matrix: Matrix, 
            ring_operator: RingOperator,
            order_comparator: OrderComparator 
        ) 
        ->      
        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderComparator,
            > 
    where 
        Matrix:                 OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:  KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,        
        Matrix::KeyMaj:                 Clone + PartialEq,
        RingOperator:           Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:        Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >, 
{

    let mut unscaled_seed_of_entries_to_be_eliminated   =   matrix.view_major_ascend( keymaj ).into_iter();
    
    // CONSTRUCT THE FIRST ENTRY OF THE INVERSE ARRAY; NOTE THAT WE 
    //      (1) TAKE THE FIRST ENTRY OF THE ARRAY TO BE INVERTED
    //      (2) **INVERT THE SCALAR VALUE** OF THAT ENTRY
    let mut diagonal_entry_of_inverse                  =   unscaled_seed_of_entries_to_be_eliminated.next().unwrap();
    diagonal_entry_of_inverse.set_val(
            ring_operator.invert( diagonal_entry_of_inverse.val() )
        );

    // CONSTRUCT THE OBJECT THAT ITERATORS OVER ENTRIES TO ELIMINATE 
    // This is obtained by taking the row/column of the original matrix, removing its leading entry, then scaling by 1/(the value of the leading entry)
    let scaled_seed_of_tail_to_be_eliminated = unscaled_seed_of_entries_to_be_eliminated.scale( 
                    diagonal_entry_of_inverse.val(),
                    ring_operator.clone(),   
                );
    // let head_to_be_eliminated   =   scaled_seed_of_tail_to_be_eliminated.next();
    // let tail_to_be_eliminated   =   hit_merge_by_predicate(
    //                                         vec![ scaled_seed_of_tail_to_be_eliminated ],
    //                                         order_comparator,
    //                                     );

    let entries_to_elim_simplified_heap
            =   Simplify::new(
                        hit_merge_by_predicate(vec![ scaled_seed_of_tail_to_be_eliminated ], order_comparator),
                        ring_operator.clone(),                                                                            
                    );
    
    ViewMajorAscendOfInverseOfTriangularArrayLazy {
        ring_operator:                  ring_operator.clone(), // the operator for the coefficient ring_operator    
        matrix:                         matrix, // the matrix to be inverted    
        next_entry_of_inv:              Some( diagonal_entry_of_inverse ), // the next entry in this view of the inverse matrix    
        entries_to_elim_simplified_heap:  entries_to_elim_simplified_heap,
    }

 }


/// The inverse of a triangular matrix `M`, computed in a lazy fashion.
/// 
/// The matrix `M` must implement `OracleMajorAscend`, and the first entry of the sparse vector iterator 
/// `M.view_major_ascend( x )` must be a nonzero entry with index `x`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
/// use oat_rust::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
/// use oat_rust::utilities::partial_order::OrderComparatorFnWrapper;
/// use oat_rust::matrices::display::print_indexed_major_views;
/// 
/// // Define the ring operator
/// let modulus = 1049;
/// let ring_operator  =   PrimeOrderFieldOperator::new(modulus);                  
/// 
/// // Define a row-major upper triangular amtrix
/// let matrix  =   VecOfVecSimple::new(
///                     vec![   
///                         vec![ (0,1), (1,1), ], 
///                         vec![        (1,1), ],
///                     ],
///                 );
/// 
/// // Define the inverse
/// let inverse = InverseOfTriangularArrayLazyAscending::new( 
///     & matrix, 
///     ring_operator.clone(),
///     OrderComparatorFnWrapper::new( |x: &(usize, usize), y: &(usize, usize)| x.0 < y.0 ) 
/// );  
/// 
/// // Print the inverse
/// let row_indices = vec![ 0, 1 ];
/// print_indexed_major_views( &inverse, row_indices );
/// ```
/// 
/// This should print the following:
/// 
/// ```bash
/// $ major_view 0: [(0, 1), (1, 1048)]
/// $ major_view 1: [(1, 1)] 
/// ```
/// 
/// which is the correct solution, since the inverse of `matrix` is
/// ```bash 
///       |  1   -1  |
///       |  0    1  |
/// ```
/// and -1 = 1048, modulo 1049.
 #[derive(Clone, Debug)]
 pub struct InverseOfTriangularArrayLazyAscending
                < 
                Matrix,
                RingOperator,
                OrderComparator,
                > 
            where 
                Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
                Matrix::ViewMajorAscend:        IntoIterator,
                Matrix::ViewMajorAscendEntry:   KeyValSet< Matrix::KeyMaj, Matrix::SnzVal >,                 
                Matrix::KeyMaj:                 Clone + PartialEq,
                RingOperator:                   Clone + Semiring< Matrix::SnzVal >,
                Matrix::SnzVal:                 Clone,
                OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,                                              


{
    ring_operator:      RingOperator, // the operator for the coefficient ring_operator    
    matrix:             Matrix, // the matrix to be inverted    
    order_comparator:   OrderComparator,  
}



 impl < Matrix, RingOperator, OrderComparator, > 
 
    InverseOfTriangularArrayLazyAscending
        < Matrix, RingOperator, OrderComparator, > 
    where 
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValSet< Matrix::KeyMaj, Matrix::SnzVal >,                 
        Matrix::KeyMaj:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,                                              
{
    pub fn new( matrix: Matrix, ring_operator: RingOperator, order_comparator: OrderComparator ) -> Self {
        InverseOfTriangularArrayLazyAscending{
            ring_operator:          ring_operator, // the operator for the coefficient ring_operator    
            matrix:                 matrix, // the matrix to be inverted    
            order_comparator:       order_comparator,    
        }              
    }
}  




impl     < 
        Matrix,
        RingOperator,
        OrderComparator,
        > 
    
    IndicesAndCoefficients for

    InverseOfTriangularArrayLazyAscending< 
            Matrix,
            RingOperator,
            OrderComparator,
        > 

    where 
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValSet< Matrix::KeyMaj, Matrix::SnzVal >,                 
        Matrix::KeyMaj:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,  
{ type KeyMin = Matrix::KeyMin; type KeyMaj = Matrix::KeyMaj; type SnzVal = Matrix::SnzVal;}            


impl     < 
        Matrix,
        RingOperator,
        OrderComparator,
        > 

    OracleMajorAscend for

    InverseOfTriangularArrayLazyAscending< 
            Matrix,
            RingOperator,
            OrderComparator,
        >     

    where 
        Matrix:                         Copy + OracleMajorAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   KeyValSet< Matrix::KeyMaj, Matrix::SnzVal >,                 
        Matrix::KeyMaj:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,         

{   
    type ViewMajorAscend            =   ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                                                Matrix,
                                                RingOperator,
                                                OrderComparator,
                                            >;
    type ViewMajorAscendEntry       =   < Self::ViewMajorAscend as IntoIterator >::Item;
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;

    fn view_major_ascend( & self, index:Matrix::KeyMaj ) 
        -> 
        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderComparator,
            >            
    {        
        view_major_of_inverse_of_triangular_array_lazy_ascend
        (   index, 
            self.matrix, 
            self.ring_operator.clone(),
            self.order_comparator.clone(), 
        )         

    }

}


// pub fn inverse_of_triangular_array_lazy_ascend < 
//             Matrix,
//             RingOperator,
//             OrderComparator,
//         > 
//         (   
//             matrix: Matrix, 
//             ring_operator: RingOperator,
//             order_comparator: OrderComparator 
//         ) 
//         ->      
//         InverseOfTriangularArrayLazyAscending< 
//                 Matrix,
//                 RingOperator,
//                 OrderComparator,
//             > 
//     where 
//         Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
//         Matrix::ViewMajorAscend:        IntoIterator,
//         Matrix::ViewMajorAscendEntry:   KeyValSet< Matrix::KeyMaj, Matrix::SnzVal >,                 
//         Matrix::KeyMaj:                 Clone + PartialEq,
//         RingOperator:                   Clone + Semiring< Matrix::SnzVal >,
//         Matrix::SnzVal:                 Clone,
//         OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,                
// {
    
//     InverseOfTriangularArrayLazyAscending{
//         ring_operator:          ring_operator, // the operator for the coefficient ring_operator    
//         matrix:                 matrix, // the matrix to be inverted    
//         order_comparator:       order_comparator,    
//     }        
// }











#[cfg(test)]
mod doc_test_drafts {
    use crate::matrices::display;

    


    #[test]
    fn test_inverse_small() {
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
        use crate::utilities::partial_order::OrderComparatorFnWrapper;
        use crate::matrices::display::print_indexed_major_views;

        // Define the ring operator (prime field of order 1049)
        let modulus = 1049;
        let ring_operator  =   PrimeOrderFieldOperator::new(modulus);  
        
        // Define the order comparator
        // This struct has a method that assigns a value of true/false to the statement
        // that `entry (a,b) should precede entry (c,d)`.  We want (a,b) to preced (c,d)
        // if a < c.
        let order_comparator = 
            OrderComparatorFnWrapper::new( |x: &(usize, usize), y: &(usize, usize)| x.0 < y.0 ) ;

        // Define a row-major upper triangular amtrix
        let matrix  =   VecOfVecSimple::new(
                                                        vec![   
                                                            vec![ (0,1), (1,1), ], 
                                                            vec![        (1,1), ],
                                                        ],
        );

        // Define the inverse
        let inverse = InverseOfTriangularArrayLazyAscending::new( 
            & matrix, 
            ring_operator.clone(),
            order_comparator,
        );  

        // Print the inverse
        let row_indices = vec![ 0, 1 ];
        print_indexed_major_views( &inverse, row_indices );
        
        // This should print the following:
        //
        // $ major_view 0: [(0, 1), (1, 1048)]
        // $ major_view 1: [(1, 1)] 
        // 
        // This is the correct solution, since the inverse of `matrix` is
        //
        //      |  1   -1  |
        //      |  0    1  |
        //
        // and -1 = 1048, modulo 1049
    }

}





















#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use itertools::Itertools;
    use crate::{matrices::{matrix_types::vec_of_vec::VecOfVecSimple, random_constructors::random_upper_unitriangular}, matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified, rings::operator_structs::field_prime_order::PrimeOrderFieldOperator};

    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Computes the inverse of the given matrix, then checks that the product of the two matrices is identity.
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_inversion_of_input< SnzVal, RingOperator >( 
                // matrix:         ClonedEntries < 
                //                         'aa, 
                //                         VecOfVecSimple<'a, usize, Matrix::SnzVal >, 
                //                         & Vec< (usize, Matrix::SnzVal) >,
                //                     >,
                matrix:         & VecOfVecSimple< usize, SnzVal, >, 
                ring_operator:  RingOperator, 
                matrix_size:    usize,
            ) 
        where   SnzVal:         Clone + PartialEq + std::fmt::Debug,
                RingOperator:   Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal > + Clone,
    {
        use crate::utilities::partial_order::OrderComparatorFnWrapper;

        let inverse =   InverseOfTriangularArrayLazyAscending::new( 
                                                                    & matrix, 
                                                                    ring_operator.clone(),
                                                                    OrderComparatorFnWrapper::new( |x: &(usize, SnzVal), y: &(usize, SnzVal)| x.0 < y.0 ) 
                                                                );   
                                                                
        let product =   ProductMatrixLazyMajorAscendSimplified::new( 
                                                                                    matrix, 
                                                                                    & inverse, 
                                                                                    ring_operator.clone(), 
                                                                                    OrderComparatorFnWrapper::new( |x:&(usize, SnzVal), y:&(usize, SnzVal)| x.key() < y.key() )  
                                                                                );

        let row_vec = | p: usize | 
                                                    product.view_major_ascend( p )
                                                    .collect_vec();
                                     
        let one = <RingOperator as crate::rings::operator_traits::Semiring<SnzVal, >>::one();
                                                        
        for keymaj in 0 .. matrix_size {
            assert_eq!( vec![ ( keymaj, one.clone() ) ] , row_vec( keymaj ) );
        }                                                                
    }


    
    #[test]
    fn test_inversion_of_specific_matrices() {
        use num::rational::Ratio;        
        // use crate::matrices::matrix_oracle_traits::MajorDimension;
        use crate::rings::operator_structs::ring_native::DivisionRingNative;

        // Define the integer type for the rational numbers we'll use
        type IntegerType = isize;
        let modulus = 1049;

        // Define shorthand functions for creating new rational numbers
        let r = | n:IntegerType, d:IntegerType | Ratio::new( n, d );    
        let q = | n:IntegerType | Ratio::from_integer( n );   
        
        // Define the ring operators
        let ring_operator_q  =   DivisionRingNative::< Ratio<IntegerType> >::new();  
        let ring_operator_f  =   DivisionRingNative::< f64 >::new();          
        let ring_operator_p  =   PrimeOrderFieldOperator::new(modulus);                  

        let matrix  =   VecOfVecSimple::new(
                                                        vec![   
                                                            vec![ (0,1.), (1,1.), ], 
                                                            vec![         (1,1.), ],
                                                        ],
        );
        test_inversion_of_input( & matrix, ring_operator_f, 2 );        

        let matrix  =   VecOfVecSimple::new(
                            vec![   
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
        );     
        test_inversion_of_input( & matrix, ring_operator_q.clone(), 3 );        
        
        let matrix  =   VecOfVecSimple::new(
                            vec![   
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
        );        
        test_inversion_of_input( & matrix, ring_operator_q.clone(), 3 );           

        let matrix  =   VecOfVecSimple::new(
                                vec![   
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
        );     
        test_inversion_of_input( & matrix, ring_operator_q.clone(), 4 );         
        
        // MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        // use rand::Rng;        // we use this module to generate random elements
        // let matrix_size =   20;

        // let mut rng = rand::thread_rng(); // this generates random integers
        // let mut vec_of_vec = vec![];
        // for keymaj in 0 .. matrix_size {
        //     let coefficient_leading         =   rng.gen_range( 1 .. modulus );
        //     let mut new_vec     =   vec![ (keymaj, coefficient_leading) ]; // start with a nonzero diagonal element            
        //     for q in keymaj+1 .. matrix_size { // fill out the rest of this row of the matrix
        //         let coefficient   = rng.gen_range( 0 .. modulus );
        //         let flag = rng.gen_range(0usize .. 3usize);
        //         if      flag == 0 { new_vec.push( ( q, 0 )           ) }
        //         else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
        //         else              { continue }
        //     }
        //     vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        // }

        // let matrix  =   VecOfVecSimple::new( vec_of_vec ); // formally wrap the matrix in a VecOfVec struct

        let matrix_size = 20;
        let matrix = random_upper_unitriangular( matrix_size, modulus );
        test_inversion_of_input( & matrix, ring_operator_p.clone(), matrix_size );                 

    }
}

