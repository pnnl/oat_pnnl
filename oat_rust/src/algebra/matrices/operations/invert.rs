

// -----------------------------------------------------------------------------------


//! Invert a triangular matrix; see also [`solve`](crate::algebra::matrices::operations::solve).
//! 
//! # How to improve performance
//! If you only need to solve `Ax = b` for `x`, where `A` is a sparse matrix, then you will probably get
//! better results by using the [`solve`](crate::algebra::matrices::operations::solve) package instead.  See the blog post
//! [don't invert that matrix](https://www.johndcook.com/blog/2010/01/19/dont-invert-that-matrix/), by John Cook,
//! for general thoughts on inverting sparse matrices.
//! 
//! # Examples
//! 
//! ```
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
//! use oat_rust::utilities::order::OrderOperatorByLessThan;
//! use oat_rust::algebra::matrices::display::print_indexed_major_views;
//! 
//! // Define the ring operator (prime field of order 1049)
//! let modulus = 1049;
//! let ring_operator  =   PrimeOrderFieldOperator::new(modulus);                  
//! 
//! // Define the order comparator
//! // This struct has a method that assigns a value of true/false to the statement
//! // that `entry (a,b) should precede entry (c,d)`.  We want (a,b) to preced (c,d)
//! // if a < c.
//! let order_operator = 
//!     OrderOperatorByLessThan::new( |x: &(usize, usize), y: &(usize, usize)| x.0 < y.0 ) ;
//! 
//! // Define a row-major upper triangular amtrix
//! let matrix  =   VecOfVec::new(
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
//!                     order_operator,
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


use crate::{algebra::rings::operator_traits::{Semiring, Ring, DivisionRing}, algebra::matrices::{query::{ViewRowAscend, IndicesAndCoefficients}}, utilities::{iterators::{merge::hit::{HitMerge, hit_bulk_insert, hit_merge_by_predicate}, general::PeekUnqualified}, order::JudgePartialOrder}, algebra::vectors::operations::{VectorOperations, Scale, Simplify}};
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet};




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
/// triangular array) where `view_major_ascend( .. )` returns an object of type [`TriangularSolve`](crate::algebra::matrices::operations::triangular_solve::TriangularSolve).
/// In that implementation, `view_major_ascend( .. )` would simply return `TriangularSolve::new( x, matrix, ring_operator, order_operator)`,
/// where `x` is an entry iterator representing a standard unit vector.  However, there are some small differences between that approach and the one
/// implemented here, namely:
/// - the alternate approach wraps a (possibly large) number of iterators in wrappers of type [`TwoTypeIterator`](crate::utilities::iterators::general::IterTwoType)
/// - the approach which is actually applied here has a `head-tail` structure which allows peeking
/// 
/// It would be interesting to contrast the performance of these two approaches, to see if there is a meaningful difference.
pub struct ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                    Matrix,
                    RingOperator,
                    OrderOperator,
                > 
    where 
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:   KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,        
        Matrix::RowIndex:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient >,
        Matrix::Coefficient:                 Clone,
        OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,                   

{
    ring_operator:                  RingOperator, // the operator for the coefficient ring_operator    
    matrix:                         Matrix, // the matrix to be inverted    
    next_entry_of_inv:              Option< Matrix::EntryMajor  >, // the next entry in this view of the inverse matrix     
    entries_to_elim_simplified_heap:  Simplify<
                                            HitMerge< 
                                                    // the HitMerge struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                    Scale< 
                                                            Matrix::ViewMajorAscendIntoIter,  // the iterators returned by the matrix
                                                            Matrix::RowIndex,
                                                            RingOperator, // the ring_operator operator
                                                            Matrix::Coefficient,
                                                        >,
                                                    // the thing that declares whether one major key comes before of after another    
                                                    OrderOperator 
                                                >,
                                            Matrix::RowIndex,
                                            RingOperator,
                                            Matrix::Coefficient,
                                        >, 
}                                   

impl    < 
        Matrix,
        RingOperator,
        OrderOperator,
        > 

        Iterator for    

        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderOperator,
            > 

    where 
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:   KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,        
        Matrix::RowIndex:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                 Clone,
        OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,               

{

    type Item = Matrix::EntryMajor;

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
                Some( return_value )

            }
        }
    }
}         

impl    < 
        Matrix,
        RingOperator,
        OrderOperator,
        > 

        PeekUnqualified for    

        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderOperator,
            > 

    where 
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:   KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,        
        Matrix::RowIndex:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                 Clone,
        OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,               


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
            OrderOperator,
        > 
        (   keymaj: Matrix::RowIndex, 
            matrix: Matrix, 
            ring_operator: RingOperator,
            order_operator: OrderOperator 
        ) 
        ->      
        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderOperator,
            > 
    where 
        Matrix:                 ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:  KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,        
        Matrix::RowIndex:                 Clone + PartialEq,
        RingOperator:           Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                 Clone,
        OrderOperator:        Clone + JudgePartialOrder <  Matrix::EntryMajor >, 
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
    //                                         order_operator,
    //                                     );

    let entries_to_elim_simplified_heap
            =   Simplify::new(
                        hit_merge_by_predicate(vec![ scaled_seed_of_tail_to_be_eliminated ], order_operator),
                        ring_operator.clone(),                                                                            
                    );
    
    ViewMajorAscendOfInverseOfTriangularArrayLazy {
        ring_operator, // the operator for the coefficient ring_operator    
        matrix, // the matrix to be inverted    
        next_entry_of_inv:              Some( diagonal_entry_of_inverse ), // the next entry in this view of the inverse matrix    
        entries_to_elim_simplified_heap,
    }

 }


/// The inverse of a triangular matrix `M`, computed in a lazy fashion.
/// 
/// The matrix `M` must implement `ViewRowAscend`, and the first entry of the sparse vector iterator 
/// `M.view_major_ascend( x )` must be a nonzero entry with index `x`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
/// use oat_rust::utilities::order::OrderOperatorByLessThan;
/// use oat_rust::algebra::matrices::display::print_indexed_major_views;
/// 
/// // Define the ring operator
/// let modulus = 1049;
/// let ring_operator  =   PrimeOrderFieldOperator::new(modulus);                  
/// 
/// // Define a row-major upper triangular amtrix
/// let matrix  =   VecOfVec::new(
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
///     OrderOperatorByLessThan::new( |x: &(usize, usize), y: &(usize, usize)| x.0 < y.0 ) 
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
 #[derive(Clone, Copy, Debug)]
 pub struct InverseOfTriangularArrayLazyAscending
                < 
                Matrix,
                RingOperator,
                OrderOperator,
                > 
            where 
                Matrix:                         ViewRowAscend + IndicesAndCoefficients,
                Matrix::ViewMajorAscend:        IntoIterator,
                Matrix::EntryMajor:   KeyValSet< Matrix::RowIndex, Matrix::Coefficient >,                 
                Matrix::RowIndex:                 Clone + PartialEq,
                RingOperator:                   Clone + Semiring< Matrix::Coefficient >,
                Matrix::Coefficient:                 Clone,
                OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,                                              


{
    ring_operator:      RingOperator, // the operator for the coefficient ring_operator    
    matrix:             Matrix, // the matrix to be inverted    
    order_operator:   OrderOperator,  
}



 impl < Matrix, RingOperator, OrderOperator, > 
 
    InverseOfTriangularArrayLazyAscending
        < Matrix, RingOperator, OrderOperator, > 
    where 
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:   KeyValSet< Matrix::RowIndex, Matrix::Coefficient >,                 
        Matrix::RowIndex:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient >,
        Matrix::Coefficient:                 Clone,
        OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,                                              
{
    pub fn new( matrix: Matrix, ring_operator: RingOperator, order_operator: OrderOperator ) -> Self {
        InverseOfTriangularArrayLazyAscending{
            ring_operator, // the operator for the coefficient ring_operator    
            matrix, // the matrix to be inverted    
            order_operator,    
        }              
    }
}  




impl     < 
        Matrix,
        RingOperator,
        OrderOperator,
        > 
    
    IndicesAndCoefficients for

    InverseOfTriangularArrayLazyAscending< 
            Matrix,
            RingOperator,
            OrderOperator,
        > 

    where 
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:             KeyValSet< Matrix::RowIndex, Matrix::Coefficient >,                 
        Matrix::RowIndex:               Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient >,
        Matrix::Coefficient:           Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,  
{ 
    type EntryMajor = Matrix::EntryMajor;
    type EntryMinor = Matrix::EntryMinor;
    type ColIndex = Matrix::ColIndex; 
    type RowIndex = Matrix::RowIndex; 
    type Coefficient = Matrix::Coefficient;
}            


impl     < 
        Matrix,
        RingOperator,
        OrderOperator,
        > 

    ViewRowAscend for

    InverseOfTriangularArrayLazyAscending< 
            Matrix,
            RingOperator,
            OrderOperator,
        >     

    where 
        Matrix:                         Copy + ViewRowAscend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::EntryMajor:   KeyValSet< Matrix::RowIndex, Matrix::Coefficient >,                 
        Matrix::RowIndex:                 Clone + PartialEq,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                 Clone,
        OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,         

{   
    type ViewMajorAscend            =   ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                                                Matrix,
                                                RingOperator,
                                                OrderOperator,
                                            >;
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;

    fn view_major_ascend( & self, index:Matrix::RowIndex ) 
        -> 
        ViewMajorAscendOfInverseOfTriangularArrayLazy< 
                Matrix,
                RingOperator,
                OrderOperator,
            >            
    {        
        view_major_of_inverse_of_triangular_array_lazy_ascend
        (   index, 
            self.matrix, 
            self.ring_operator.clone(),
            self.order_operator.clone(), 
        )         

    }

}


// pub fn inverse_of_triangular_array_lazy_ascend < 
//             Matrix,
//             RingOperator,
//             OrderOperator,
//         > 
//         (   
//             matrix: Matrix, 
//             ring_operator: RingOperator,
//             order_operator: OrderOperator 
//         ) 
//         ->      
//         InverseOfTriangularArrayLazyAscending< 
//                 Matrix,
//                 RingOperator,
//                 OrderOperator,
//             > 
//     where 
//         Matrix:                         ViewRowAscend + IndicesAndCoefficients,
//         Matrix::ViewMajorAscend:        IntoIterator,
//         Matrix::EntryMajor:   KeyValSet< Matrix::RowIndex, Matrix::Coefficient >,                 
//         Matrix::RowIndex:                 Clone + PartialEq,
//         RingOperator:                   Clone + Semiring< Matrix::Coefficient >,
//         Matrix::Coefficient:                 Clone,
//         OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,                
// {
    
//     InverseOfTriangularArrayLazyAscending{
//         ring_operator:          ring_operator, // the operator for the coefficient ring_operator    
//         matrix:                 matrix, // the matrix to be inverted    
//         order_operator:       order_operator,    
//     }        
// }











#[cfg(test)]
mod doc_test_drafts {
    

    


    #[test]
    fn test_inverse_small() {
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
        use crate::utilities::order::OrderOperatorByLessThan;
        use crate::algebra::matrices::display::print_indexed_major_views;

        // Define the ring operator (prime field of order 1049)
        let modulus = 1049;
        let ring_operator  =   PrimeOrderFieldOperator::new(modulus);  
        
        // Define the order comparator
        // This struct has a method that assigns a value of true/false to the statement
        // that `entry (a,b) should precede entry (c,d)`.  We want (a,b) to preced (c,d)
        // if a < c.
        let order_operator = 
            OrderOperatorByLessThan::new( |x: &(usize, usize), y: &(usize, usize)| x.0 < y.0 ) ;

        // Define a row-major upper triangular amtrix
        let matrix  =   VecOfVec::new(
                                                        vec![   
                                                            vec![ (0,1), (1,1), ], 
                                                            vec![        (1,1), ],
                                                        ],
        );

        // Define the inverse
        let inverse = InverseOfTriangularArrayLazyAscending::new( 
            & matrix, 
            ring_operator,
            order_operator,
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
    use crate::{algebra::matrices::{types::vec_of_vec::sorted::VecOfVec, }, algebra::matrices::operations::multiply::ProductMatrix, algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator};

    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Computes the inverse of the given matrix, then checks that the product of the two matrices is identity.
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_inversion_of_input< Coefficient, RingOperator >( 
                // matrix:         ClonedRowEntries < 
                //                         'aa, 
                //                         VecOfVec<'a, usize, Matrix::Coefficient >, 
                //                         & Vec< (usize, Matrix::Coefficient) >,
                //                     >,
                matrix:         & VecOfVec< usize, Coefficient, >, 
                ring_operator:  RingOperator, 
                matrix_size:    usize,
            ) 
        where   Coefficient:         Clone + PartialEq + std::fmt::Debug,
                RingOperator:   Semiring< Coefficient > + Ring< Coefficient > + DivisionRing< Coefficient > + Clone,
    {
        use crate::utilities::order::OrderOperatorByLessThan;

        let inverse =   InverseOfTriangularArrayLazyAscending::new( 
                                                                    & matrix, 
                                                                    ring_operator.clone(),
                                                                    OrderOperatorByLessThan::new( |x: &(usize, Coefficient), y: &(usize, Coefficient)| x.0 < y.0 ) 
                                                                );   
                                                                
        let product =   ProductMatrix::new( 
                                                                                    matrix, 
                                                                                    & inverse, 
                                                                                    ring_operator, 
                                                                                    OrderOperatorByLessThan::new( |x:&(usize, Coefficient), y:&(usize, Coefficient)| x.key() < y.key() )  
                                                                                );

        let row_vec = | p: usize | 
                                                    product.view_major_ascend( p )
                                                    .collect_vec();
                                     
        let one = <RingOperator as crate::algebra::rings::operator_traits::Semiring<Coefficient, >>::one();
                                                        
        for keymaj in 0 .. matrix_size {
            assert_eq!( vec![ ( keymaj, one.clone() ) ] , row_vec( keymaj ) );
        }                                                                
    }


    
    #[test]
    fn test_inversion_of_specific_matrices() {
        use num::rational::Ratio;        
        // use crate::algebra::matrices::query::MajorDimension;
        use crate::algebra::rings::operator_structs::ring_native::DivisionRingNative;

        // Define the integer type for the rational numbers we'll use
        type IntegerType = isize;
        let modulus = 1049;

        // Define shorthand functions for creating new rational numbers
        let r = Ratio::new;    
        let q = Ratio::from_integer;   
        
        // Define the ring operators
        let ring_operator_q  =   DivisionRingNative::< Ratio<IntegerType> >::new();  
        let ring_operator_f  =   DivisionRingNative::< f64 >::new();          
        let ring_operator_p  =   PrimeOrderFieldOperator::new(modulus);                  

        let matrix  =   VecOfVec::new(
                                                        vec![   
                                                            vec![ (0,1.), (1,1.), ], 
                                                            vec![         (1,1.), ],
                                                        ],
        );
        test_inversion_of_input( & matrix, ring_operator_f, 2 );        

        let matrix  =   VecOfVec::new(
                            vec![   
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
        );     
        test_inversion_of_input( & matrix, ring_operator_q, 3 );        
        
        let matrix  =   VecOfVec::new(
                            vec![   
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
        );        
        test_inversion_of_input( & matrix, ring_operator_q, 3 );           

        let matrix  =   VecOfVec::new(
                                vec![   
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
        );     
        test_inversion_of_input( & matrix, ring_operator_q, 4 );         
        
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

        // let matrix  =   VecOfVec::new( vec_of_vec ); // formally wrap the matrix in a VecOfVec struct

        let matrix_size = 20;
        let matrix = VecOfVec::random_mod_p_upper_unitriangular( matrix_size, modulus );
        test_inversion_of_input( & matrix, ring_operator_p, matrix_size );                 

    }
}

