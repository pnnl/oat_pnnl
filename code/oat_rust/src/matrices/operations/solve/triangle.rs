//! Solve `Ax = b` for `x`, where `A` is a triangular matrix with nonzero diagonal entries.
//! 
//! Concretely, `A` is a matrix oracle; the major and minor keys of `A` must have the same type.
//! 
//! # Examples
//! 
//! ```
//! use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
//! use oat_rust::matrices::operations::solve::triangle::TriangularSolveAscend;
//! use oat_rust::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;        
//! use oat_rust::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::partial_order::OrderComparatorAutoLtByKey;
//! use itertools::Itertools;
//! 
//! // define the ring operator
//! let ring_operator = BooleanFieldOperator::new();
//! 
//! // build a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVecSimple::new(
//!                         vec![                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );    
//! 
//! // SOLVE xA = b FOR x
//! // --------------------------
//! 
//! // define a sparse vector b
//! let b = vec![ (0,true), (1,true) ];        
//! 
//! // define an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//! let order_comparator_major_ascend = OrderComparatorAutoLtByKey::new();          
//! 
//! // create a solver to solve xA = b for x
//! let solution =  TriangularSolveAscend::new(
//!                         b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
//!                         & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVecSimple, not on VecOfVecSimple itself
//!                         ring_operator.clone(), 
//!                         order_comparator_major_ascend.clone(),
//!                     );
//! 
//! // check the solution, i.e. check that xA = b
//! let product = vector_matrix_multiply_major_ascend_simplified(solution, &matrix, ring_operator, order_comparator_major_ascend);
//! assert_eq!( product.collect_vec(), b );    
//! ```

use std::marker::PhantomData;

use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
use crate::rings::operator_traits::{Semiring, Ring, DivisionRing, MinusOne};
use crate::utilities::iterators::merge::heap_of_iterators::{HitMerge, hit_bulk_insert, hit_merge_by_predicate};
use crate::utilities::partial_order::StrictlyLess;
use crate::entries::{KeyValSet, KeyValGet};
use crate::vectors::operations::{Simplify, Scale, Transforms};
use crate::utilities::iterators::general::IterTwoType;


// TRIANGULAR SOLUTION
// ---------------------------------------------------------------------------

/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, where `A` is invertible and triangular, and `A.view_major_ascend( i )` returns an iterator whose first entry index `i`.
/// 
/// - There are two types of triangular array `A` that satisfy this criterion: row-major upper triangular arrays, and column-major lower-triangular arrays.
/// - Iterator `A.view_major_ascend( i )` to be be *strictly* ascending; entries with repeat indices are alllowed.
pub struct TriangularSolveAscend< 
                    ProblemVector,
                    Matrix,
                    RingOperator,
                    OrderComparator,
                > 
    where 
        ProblemVector:                      IntoIterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             OracleMajorAscend,
        Matrix::KeyMaj:                     Clone + PartialEq,
        Matrix::SnzVal:                     Clone,   
        Matrix::ViewMajorAscend:            IntoIterator,        
        Matrix::ViewMajorAscendEntry:       KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                       Clone + Semiring< Matrix::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,                   

{
    ring_operator:                  RingOperator, // the operator for the coefficient ring_operator    
    matrix:                         Matrix, // the `A` in `Ax = b`
    entries_to_elim_simplified_heap:  Simplify<
                                            HitMerge< 
                                                    // the HitMerge struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                    Scale< 
                                                            IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                                                    < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                                                    ProblemVector::IntoIter,
                                                                    // Matrix::ViewMajorAscendEntry,
                                                                >,
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
            ProblemVector,
            Matrix,                  
            RingOperator,
            OrderComparator,
        > 

        Iterator for    

        TriangularSolveAscend< 
                ProblemVector,
                Matrix,                
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             OracleMajorAscend,
        Matrix::KeyMaj:                     Clone + PartialEq,
        Matrix::SnzVal:                     Clone,   
        Matrix::ViewMajorAscend:            IntoIterator,        
        Matrix::ViewMajorAscendEntry:       KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                       Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,  

{
    type Item = Matrix::ViewMajorAscendEntry;

    fn next( &mut self ) -> Option<Self::Item> { 

    // THIS ITERATOR PROCEED BY ELIMINATING ENTRIES IN `self.entries_to_elim_simplified_heap`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_elim_simplified_heap` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( mut entry_to_eliminate )       =   self.entries_to_elim_simplified_heap.next() {

            // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
            let mut seed_of_eliminating_iterator    =   self.matrix.view_major_ascend( entry_to_eliminate.key() ).into_iter();
            
            // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
            let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

            // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
            let scale_factor     =    self.ring_operator.negate(
                                                    self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                );
            
            // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
            // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_elim_simplified_heap
            let eliminating_iterator       
                    =   IterTwoType::Iter1( seed_of_eliminating_iterator ) 
                            .scale( scale_factor.clone(), self.ring_operator.clone() );

            // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_elim_simplified_heap`
            hit_bulk_insert( &mut self.entries_to_elim_simplified_heap.unsimplified, vec![eliminating_iterator] ); 

            // NOW RE-USE `entry_to_eliminate` AS A CONTAINER FOR THE NEXT ENTRY IN THE INVERSE MATRIX
            entry_to_eliminate.set_val( scale_factor );

            return Some( entry_to_eliminate ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function

        } else {
            return None
        }
    }
}  





impl    < 
            ProblemVector,
            Matrix,
            RingOperator,
            OrderComparator,
        >

        TriangularSolveAscend< 
                ProblemVector,
                Matrix,                 
                RingOperator,
                OrderComparator,
            >  
            
    where 
        ProblemVector:                      IntoIterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             OracleMajorAscend,
        Matrix::KeyMaj:                     Clone + PartialEq,
        Matrix::SnzVal:                     Clone,   
        Matrix::ViewMajorAscend:            IntoIterator,        
        Matrix::ViewMajorAscendEntry:       KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                       Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,  
{        
    /// Solve `Ax = b` for `x`, where `A` is triangular and `A.view_major_ascend( x )` returns an iterator whose first entry index `x`.
    /// 
    /// There are two types of triangular array `A` that satisfy this criterion: row-major upper triangular arrays, and column-major lower-triangular arrays.
    pub fn new (
                    b:              ProblemVector,                
                    a:              Matrix,
                    ring_operator:  RingOperator,
                    order_comparator:      OrderComparator,
                )
            ->
            TriangularSolveAscend< 
                    ProblemVector,
                    Matrix,                 
                    RingOperator,
                    OrderComparator,
                >  

    {

        let entries_to_elim_simplified_heap
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    IterTwoType::Iter2( b.into_iter() )   // wrap b in an IterTwoType enum so that it can be easily combined with  
                                        .scale( ring_operator.minus_one(), ring_operator.clone() ) 
                                ),
                            order_comparator,
                        )
                        .simplify( ring_operator.clone() );
        TriangularSolveAscend{
            ring_operator:                      ring_operator.clone() ,
            matrix:                             a,
            entries_to_elim_simplified_heap:    entries_to_elim_simplified_heap, 
        } 
    }     
}










//  DOCSTRING TESTS
//  ===========================================================================






#[cfg(test)]
mod doctring_tests {


    #[test]
    fn test_docstring_solve() {
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::operations::solve::triangle::TriangularSolveAscend; //echelon::{EchelonSolveMajorAscendWithMajorKeys, EchelonSolveMinorDescendWithMinorKeys};
        use crate::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;        
        use crate::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
        use crate::utilities::partial_order::OrderComparatorAutoLtByKey;
        use itertools::Itertools;

        // define the ring operator
        let ring_operator = BooleanFieldOperator::new();

        // build a matrix A with an invertible upper triangular submatrix
        let matrix  =   VecOfVecSimple::new(
                                vec![                                     
                                    vec![ (0,true),  (1,true), (2,true)       ], 
                                    vec![            (1,true), (2,true)       ],
                                    vec![                      (2,true)       ],                                       
                                ],
                            );    

        // SOLVE xA = b FOR x
        // --------------------------

        // define a sparse vector b
        let b = vec![ (0,true), (1,true) ];        

        // define an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
        let order_comparator_major_ascend = OrderComparatorAutoLtByKey::new();          

        // create a solver to solve xA = b for x
        let solution =  TriangularSolveAscend::new(
                                b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
                                & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVecSimple, not on VecOfVecSimple itself
                                ring_operator.clone(), 
                                order_comparator_major_ascend.clone(),
                            );
        
        // check the solution, i.e. check that xA = b
        let product = vector_matrix_multiply_major_ascend_simplified(solution, &matrix, ring_operator, order_comparator_major_ascend);
        assert_eq!( product.collect_vec(), b );     
        
    }

}





















//  TESTS
//  ===========================================================================


#[cfg(test)]
mod tests {

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use itertools::Itertools;
    use crate::{matrices::matrix_types::vec_of_vec::{VecOfVecSimple}, rings::operator_structs::field_prime_order::PrimeOrderFieldOperator};


    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Computes the inverse of the given matrix, then checks that the product of the two matrices is identity.
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_triangular_solve< SnzVal, RingOperator >( 
                matrix:         & VecOfVecSimple< usize, SnzVal >, 
                ring_operator:  RingOperator, 
                matrix_size:    usize 
            ) 
        where   SnzVal:         Clone + PartialEq + Ord + std::fmt::Debug,
                RingOperator:   Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal > + Clone,
    {
        use crate::{matrices::operations::multiply::{vector_matrix_multiply_major_ascend_unsimplified, vector_matrix_multiply_major_ascend_simplified}, utilities::partial_order::OrderComparatorAutoAnyType};

        
        // let order_comparator = |x:&(usize, SnzVal), y:&(usize, SnzVal)| x.key() < y.key();

        for num_nonzeros in [0, 1, 2, 3, matrix_size].iter().map(|x| std::cmp::min(x, &matrix_size) ) {
            for vector_support in (0 .. matrix_size).combinations( *num_nonzeros ) {
                let mut vec = vector_support
                                .iter()
                                .cloned()
                                .map(|x| ( x, RingOperator::one() ))
                                .collect_vec();
                vec.sort();
                let solution = TriangularSolveAscend::new(
                                        vec.iter().cloned(),                    
                                        & matrix,
                                        ring_operator.clone(),
                                        OrderComparatorAutoAnyType,                    
                                    )
                                    .peekable()
                                    .simplify( ring_operator.clone() )
                                    .collect_vec();
                assert_eq!(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ).collect_vec(),
                    vector_matrix_multiply_major_ascend_simplified( 
                            solution, matrix, 
                            ring_operator.clone(), 
                            OrderComparatorAutoAnyType
                        )
                        .collect_vec()
                );
            }
        }                                                               
    }


    
    #[test]
    fn test_triangular_solve_on_specific_matrices() {
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
        let ring_operator_p  =   PrimeOrderFieldOperator::new(modulus);                  

        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(1)), ], 
                                    vec![            (1,q(1)), ],
                                ],
        );
        test_triangular_solve( &matrix, ring_operator_q.clone(), 2 );

        // these lines result in an ERROR becuase the constructor `VecOfVecSimple::new` panics when it receives a vector of vectors where an internal vector is not sorted in *strictly* asecneding order
        // let matrix  =   VecOfVecSimple::new(
        //                         vec![ 
        //                             vec![ (0,q(1)), (0,q(1)),   (1,q(1)), (1,q(1)), ], 
        //                             vec![                       (1,q(1)), (1,q(1)), ],
        //                         ],
        // );
        // test_triangular_solve( &matrix, ring_operator_q.clone(), 2 );        

        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
        );     
        test_triangular_solve( &matrix, ring_operator_q.clone(), 3 );        
        
        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
        );        
        test_triangular_solve( &matrix, ring_operator_q.clone(), 3 );           

        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
        );     
        test_triangular_solve( &matrix, ring_operator_q.clone(), 4 );         
        
        //  MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        //  ----------------------------------------------------------------------

        use crate::matrices::random_constructors::random_upper_triangular_matrix_mod_p;
        let matrix_size =   20;
        let matrix = random_upper_triangular_matrix_mod_p( matrix_size, modulus );

        test_triangular_solve( &matrix, ring_operator_p.clone(), matrix_size );                 

    }
}