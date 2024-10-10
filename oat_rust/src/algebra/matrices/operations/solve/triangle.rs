//! Solve `Ax = b` for `x`, where `A` is a upper triangular with nonzero diagonal entries.
//! 
//! Concretely, `A` is a matrix oracle; the major and minor keys of `A` must have the same type.
//! 
//! # Alternate solutions
//! 
//! You can also solve the exact same problem, `Ax = b`, using the [echelon](crate::algebra::matrices::operations::solve::echelon) module.  That module doesn't assume that your matrix is triangular, so you have to use a function that "matches" the row and column indices of diagonal elements; for us, that function can just be [identity](crate::utilities::functions::evaluate::identity)



use std::marker::PhantomData;

use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients};
use crate::algebra::matrices::types::transpose::AntiTranspose;
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing, };
use crate::utilities::iterators::merge::hit::{HitMerge, hit_bulk_insert, hit_merge_by_predicate};
use crate::utilities::order::{JudgePartialOrder, ReverseOrder};
use crate::algebra::vectors::entries::{KeyValSet, KeyValGet};
use crate::algebra::vectors::operations::{Simplify, Scale, VectorOperations};
use crate::utilities::iterators::general::IterTwoType;


// TRIANGULAR SOLUTION -- ASCEND
// ---------------------------------------------------------------------------

/// Solve `xA = b` where `A` is  is upper-triangular and row-major
/// 
/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, **in strictly ascending order**.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::operations::solve::triangle::TriangularSolverMajorAscend;
/// use oat_rust::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;        
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// use itertools::Itertools;
/// 
/// // define the ring operator
/// let ring_operator = BooleanFieldOperator::new();
/// 
/// // build a matrix A with an invertible upper triangular submatrix
/// let matrix  =   VecOfVec::new(
///                         vec![                                     
///                             vec![ (0,true),  (1,true), (2,true)       ], 
///                             vec![            (1,true), (2,true)       ],
///                             vec![                      (2,true)       ],                                       
///                         ],
///                     );    
/// 
/// // SOLVE xA = b FOR x
/// // --------------------------
/// 
/// // define a sparse vector b
/// let b = vec![ (0,true), (1,true) ];        
/// 
/// // define an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
/// let order_operator_major_ascend = OrderOperatorByKey::new();          
/// 
/// // create a solver to solve xA = b for x
/// let solution =  TriangularSolverMajorAscend::solve(
///                         b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
///                         & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
///                         ring_operator.clone(), 
///                         order_operator_major_ascend.clone(),
///                     );
/// 
/// // check the solution, i.e. check that xA = b
/// let product = vector_matrix_multiply_major_ascend_simplified(solution, &matrix, ring_operator, order_operator_major_ascend);
/// assert_eq!( product.collect_vec(), b );    
/// ```
/// 
pub struct TriangularSolverMajorAscend< 
                    ProblemVector,
                    Matrix,
                    Key,
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                      IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             IndicesAndCoefficients<RowIndex=Key, ColIndex=Key> + ViewRowAscend,
        Matrix::RowIndex:                     Clone + PartialEq,
        Matrix::Coefficient:                     Clone,   
        Matrix::ViewMajorAscend:            IntoIterator,        
        Matrix::EntryMajor:       KeyValGet < Key, Matrix::Coefficient > + KeyValSet < Key, Matrix::Coefficient >,      
        RingOperator:                       Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                      Clone + JudgePartialOrder <  Matrix::EntryMajor >,                   

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
                                                                    // Matrix::EntryMajor,
                                                                >,
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
            ProblemVector,
            Matrix,
            Key,        
            RingOperator,
            OrderOperator,
        > 

    TriangularSolverMajorAscend< 
            ProblemVector,
            Matrix,       
            Key,                         
            RingOperator,
            OrderOperator,
        > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             IndicesAndCoefficients<RowIndex=Key, ColIndex=Key> + ViewRowAscend,
        Matrix::RowIndex:                     Clone + PartialEq,
        Matrix::Coefficient:                     Clone,   
        Matrix::ViewMajorAscend:            IntoIterator,        
        Matrix::EntryMajor:       KeyValGet < Key, Matrix::Coefficient > + KeyValSet < Key, Matrix::Coefficient >,      
        RingOperator:                       Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                      Clone + JudgePartialOrder <  Matrix::EntryMajor >,                   
 
{        
    /// Solve `Ax = b` for `x`
    /// 
    /// # Requirements
    /// 
    /// - `A.view_major_ascend( x )` returns an iterator whose first entry has index `x` and a nonzero coefficient.
    ///    There are two types of triangular array `A` that satisfy this criterion
    ///    - row-major upper triangular arrays,
    ///    - column-major lower-triangular arrays.
    /// - `b` must iterate over entries in strictly ascending order
    /// 
    /// 
    ///  See the [TriangularSolverMajorAscend] for an example.
    pub fn solve(
                    b:              ProblemVector,                
                    a:              Matrix,
                    ring_operator:  RingOperator,
                    order_operator:      OrderOperator,
                )
            ->
            TriangularSolverMajorAscend< 
                    ProblemVector,
                    Matrix,      
                    Key,           
                    RingOperator,
                    OrderOperator,
                >  

    {

        let entries_to_elim_simplified_heap
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    IterTwoType::Iter2( b.into_iter() )   // wrap b in an IterTwoType enum so that it can be easily combined with  
                                        .scale( ring_operator.minus_one(), ring_operator.clone() ) 
                                ),
                            order_operator,
                        )
                        .simplify( ring_operator.clone() );
        TriangularSolverMajorAscend{
            ring_operator ,
            matrix:                             a,
            entries_to_elim_simplified_heap, 
        } 
    }     
}






impl    < 
            ProblemVector,
            Matrix,
            Key,        
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        TriangularSolverMajorAscend< 
                ProblemVector,
                Matrix,       
                Key,                         
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             IndicesAndCoefficients<RowIndex=Key, ColIndex=Key> + ViewRowAscend,
        Matrix::RowIndex:                     Clone + PartialEq,
        Matrix::Coefficient:                     Clone,   
        Matrix::ViewMajorAscend:            IntoIterator,        
        Matrix::EntryMajor:       KeyValGet < Key, Matrix::Coefficient > + KeyValSet < Key, Matrix::Coefficient >,      
        RingOperator:                       Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                      Clone + JudgePartialOrder <  Matrix::EntryMajor >,                   

{
    type Item = Matrix::EntryMajor;

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

            Some( entry_to_eliminate ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function

        } else {
            None
        }
    }
}  








// TRIANGULAR SOLUTION -- DESCEND
// ---------------------------------------------------------------------------


/// Solve `Ax = b` where `A` is upper-triangular and row-major
/// 
/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, **in strictly descending order**.
/// 
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::operations::solve::triangle::TriangularSolverMinorDescend;
/// use oat_rust::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified;        
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// use itertools::Itertools;
/// 
/// // define the ring operator
/// let ring_operator = BooleanFieldOperator::new();
/// 
/// // build a matrix A with an invertible upper triangular submatrix
/// let matrix  =   VecOfVec::new(
///                         vec![                                     
///                             vec![ (0,true),  (1,true), (2,true)       ], 
///                             vec![            (1,true), (2,true)       ],
///                             vec![                      (2,true)       ],                                       
///                         ],
///                     );    
/// 
/// // SOLVE Ax = b FOR x
/// // --------------------------
/// 
/// // define a sparse vector b
/// let b = vec![ (2,true), (0,true) ];        
/// 
/// // define an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
/// let order_operator_major_ascend = OrderOperatorByKey::new();          
/// 
/// // create a solver to solve xA = b for x
/// let solution =  TriangularSolverMinorDescend::solve(
///                         b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
///                         & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
///                         ring_operator.clone(), 
///                         order_operator_major_ascend.clone(),
///                     );
/// 
/// // check the solution, i.e. check that xA = b
/// let product = vector_matrix_multiply_minor_descend_simplified(solution, &matrix, ring_operator, order_operator_major_ascend);
/// assert_eq!( product.collect_vec(), b );    
/// ```
pub struct TriangularSolverMinorDescend< 
                    ProblemVector,
                    Matrix,
                    Key,
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                      IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             IndicesAndCoefficients<RowIndex=Key,ColIndex=Key> + ViewColDescend,
        Matrix::ColIndex:                     Clone + PartialEq,
        Matrix::Coefficient:                     Clone,   
        Matrix::ViewMinorDescend:           IntoIterator,        
        Matrix::EntryMinor:      KeyValGet < Key, Matrix::Coefficient > + KeyValSet < Key, Matrix::Coefficient >,
        RingOperator:                       Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                      Clone + JudgePartialOrder <  Matrix::EntryMinor >,                   

{
    antitranspose_solution: TriangularSolverMajorAscend<
                            ProblemVector,
                            AntiTranspose< Matrix >,
                            Key,
                            RingOperator,
                            ReverseOrder< OrderOperator >
                        >,
    phantom_orderoperator:  PhantomData<OrderOperator>  
} 




impl    < 
            ProblemVector,
            Matrix,
            Key,        
            RingOperator,
            OrderOperator,
        > 

    TriangularSolverMinorDescend< 
            ProblemVector,
            Matrix,       
            Key,                         
            RingOperator,
            OrderOperator,
        > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             IndicesAndCoefficients<RowIndex=Key,ColIndex=Key> + ViewColDescend,
        Matrix::ColIndex:                     Clone + PartialEq,
        Matrix::Coefficient:                     Clone,   
        Matrix::ViewMinorDescend:           IntoIterator,        
        Matrix::EntryMinor:      KeyValGet < Key, Matrix::Coefficient > + KeyValSet < Key, Matrix::Coefficient >,
        RingOperator:                       Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                      Clone + JudgePartialOrder <  Matrix::EntryMinor >,                   
 
{        
    /// Solve `Ax = b` for `x`
    /// 
    /// # Requirements
    /// 
    /// - `A.view_minor_ascend( x )` returns an iterator whose first entry has index `x` and a nonzero coefficient.
    ///    There are two types of triangular array `A` that satisfy this criterion
    ///    - row-major upper triangular arrays,
    ///    - column-major lower-triangular arrays.
    /// - `b` must iterate over entries in strictly **descending** order
    /// 
    ///  See the [TriangularSolverMinorDescend] for an example.
    pub fn solve(
                    b:              ProblemVector,                
                    a:              Matrix,
                    ring_operator:  RingOperator,
                    order_operator:      OrderOperator,
                )
            ->
            TriangularSolverMinorDescend< 
                    ProblemVector,
                    Matrix,      
                    Key,           
                    RingOperator,
                    OrderOperator,
                >  

    {

        let antitranspose = AntiTranspose::new(a);
        let order_reverse = ReverseOrder::new( order_operator );
        let antitranspose_solution = TriangularSolverMajorAscend::solve(b, antitranspose, ring_operator, order_reverse );
        TriangularSolverMinorDescend{ antitranspose_solution, phantom_orderoperator: PhantomData }
    }     
}






impl    < 
            ProblemVector,
            Matrix,  
            Key,                
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        TriangularSolverMinorDescend< 
                ProblemVector,
                Matrix,         
                Key,       
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        Matrix:                             IndicesAndCoefficients<RowIndex=Key,ColIndex=Key> + ViewColDescend,
        Matrix::ColIndex:                     Clone + PartialEq,
        Matrix::Coefficient:                     Clone,   
        Matrix::ViewMinorDescend:           IntoIterator,        
        Matrix::EntryMinor:      KeyValGet < Key, Matrix::Coefficient > + KeyValSet < Key, Matrix::Coefficient >,
        RingOperator:                       Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                      Clone + JudgePartialOrder <  Matrix::EntryMinor >,                   

{
    type Item = Matrix::EntryMinor;

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitranspose_solution.next()
    }
}  









//  DOCSTRING TESTS
//  ===========================================================================






#[cfg(test)]
mod doctring_tests {


    #[test]
    fn test_docstring_solve() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::solve::triangle::TriangularSolverMajorAscend; //echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
        use crate::algebra::matrices::operations::multiply::vector_matrix_multiply_major_ascend_simplified;        
        use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        
        use crate::utilities::order::OrderOperatorByKey;
        use itertools::Itertools;

        // define the ring operator
        let ring_operator = BooleanFieldOperator::new();

        // build a matrix A with an invertible upper triangular submatrix
        let matrix  =   VecOfVec::new(
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
        let order_operator_major_ascend = OrderOperatorByKey::new();          

        // create a solver to solve xA = b for x
        let solution =  TriangularSolverMajorAscend::solve(
                                b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
                                & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                ring_operator, 
                                order_operator_major_ascend.clone(),
                            );
        
        // check the solution, i.e. check that xA = b
        let product = vector_matrix_multiply_major_ascend_simplified(solution, &matrix, ring_operator, order_operator_major_ascend);
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
    use crate::{algebra::matrices::types::vec_of_vec::sorted::VecOfVec, algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator};


    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Computes the inverse of the given matrix, then checks that the product of the two matrices is identity.
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_triangular_solve< Coefficient, RingOperator >( 
                matrix:         & VecOfVec< usize, Coefficient >, 
                ring_operator:  RingOperator, 
                matrix_size:    usize 
            ) 
        where   Coefficient:         Clone + PartialEq + Ord + std::fmt::Debug,
                RingOperator:   Semiring< Coefficient > + Ring< Coefficient > + DivisionRing< Coefficient > + Clone,
    {
        use crate::{algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified}, utilities::order::OrderOperatorAuto};

        
        // let order_operator = |x:&(usize, Coefficient), y:&(usize, Coefficient)| x.key() < y.key();

        for num_nonzeros in [0, 1, 2, 3, matrix_size].iter().map(|x| std::cmp::min(x, &matrix_size) ) {
            for vector_support in (0 .. matrix_size).combinations( *num_nonzeros ) {
                let mut vec = vector_support
                                .iter()
                                .cloned()
                                .map(|x| ( x, RingOperator::one() ))
                                .collect_vec();
                vec.sort();
                let solution = TriangularSolverMajorAscend::solve(
                                        vec.iter().cloned(),                    
                                        & matrix,
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .peekable()
                                    .simplify( ring_operator.clone() )
                                    .collect_vec();
                assert_eq!(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ).collect_vec(),
                    vector_matrix_multiply_major_ascend_simplified( 
                            solution, matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                        .collect_vec()
                );
            }
        }                                                               
    }


    
    #[test]
    fn test_triangular_solve_on_specific_matrices() {
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
        let ring_operator_p  =   PrimeOrderFieldOperator::new(modulus);                  

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(1)), ], 
                                    vec![            (1,q(1)), ],
                                ],
        );
        test_triangular_solve( &matrix, ring_operator_q, 2 );

        // these lines result in an ERROR becuase the constructor `VecOfVec::new` panics when it receives a vector of vectors where an internal vector is not sorted in *strictly* asecneding order
        // let matrix  =   VecOfVec::new(
        //                         vec![ 
        //                             vec![ (0,q(1)), (0,q(1)),   (1,q(1)), (1,q(1)), ], 
        //                             vec![                       (1,q(1)), (1,q(1)), ],
        //                         ],
        // );
        // test_triangular_solve( &matrix, ring_operator_q.clone(), 2 );        

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
        );     
        test_triangular_solve( &matrix, ring_operator_q, 3 );        
        
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
        );        
        test_triangular_solve( &matrix, ring_operator_q, 3 );           

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
        );     
        test_triangular_solve( &matrix, ring_operator_q, 4 );         
        
        //  MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        //  ----------------------------------------------------------------------

        let matrix_size =   20;
        let matrix = VecOfVec::random_mod_p_upper_triangular_invertible( matrix_size, modulus );

        test_triangular_solve( &matrix, ring_operator_p, matrix_size );                 

    }
}