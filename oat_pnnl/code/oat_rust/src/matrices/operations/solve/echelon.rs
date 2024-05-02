//! Solve `Ax = b` for `x`, where `A` is a matrix partially in echelon form.
//! 
//! Concretely, `A` is a matrix oracle and `b` is an iterator that runs over the
//! structural nonzero entries of a sparse vector.  
//! 
//! Here is a conceptual example of what is meant by "echelon solve".  Suppose we have
//! a matrix `A` and vector `b` with the sparse structure show below (`*` indicates a nonzero entry).
// //! It can be shown that `b` lies in the row space of `A` iff `b[1:3]` lies in the
// //! row space of the submatrix `A[2:4, 1:3]`.  
//! If there exists a 
//! vector `x` such that `xA = b`, then we can find such an `x` via "back substituion": first we add a scalar
//! multiple of `A[2,:]` to clear the first entry in `b`, resulting in a vector `c`.
//! Then we add a scalar multiple of `A[3,:]` to clear the first nonzero entry in `c` (which occurs in column 2).
//! This process continues until we have eliminated all nonzero entries.  The
//! solution `x` can then be constructed from the collection of scalar multipliers
//! used to perform the elimination.
//! 
//! 
//! ```ignore
//!              1  2  3  4  5  6  7
//!           -------------------------
//!         1 |                        |
//!         2 |  *  *  *  *  *  *  *   |
//!   A  =  3 |     *  *  *  *  *  *   |
//!         4 |        *  *  *  *  *   |
//!         5 |                        |
//!           -------------------------
//! 
//!              1  2  3  4  5  6  7
//!            -------------------------
//!   b  =    |  *  *  *  *  *  *  *   |
//!            -------------------------
//! ```
//! 
//! We can use an analogous method whenever we have the following data
//! 
//! - a matrix `A`
//! - a vector `b`
//! - a set of row indices `I` and a set of column indices `J` such that `A[I,J]` is upper triangular and invertible.
//! 
//! we require, in addition, that 
//! 
//! - if `[i,j]` is a diagonal entry of `A[I,J]`, then all entries to the left of `[i,j]` in `A` vanish.
//!   - this applies only in cases where we we add major views of `A` to `b`
//! - if `[i,j]` is a diagonal entry of `A[I,J]`, then all entries below `[i,j]`  in `A` vanish.
//!   - this applies only in cases where we we add minor views of `A` to `b`
//! 
//! In practice, when we use oat_rust to solve `xA = b`, we do not provide `I` and `J` explicitly;
//! rather we note that there are (mutually inverse) bijections `f: I -> J` and `g: J -> I` such that
//! `[i, f(i)]` and `[g(j), j]` lie on the diagonal of `A[I,J]`.  We pass these bijections to the solver, instead.
//! 
//! # Examples
//! 
//! ```
//! use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
//! use oat_rust::matrices::operations::solve::echelon::{EchelonSolveMajorAscendWithMajorKeys, EchelonSolveMinorDescendWithMinorKeys};
//! use oat_rust::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified, vector_matrix_multiply_minor_descend_simplified};        
//! use oat_rust::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::partial_order::{OrderComparatorAutoLtByKey, OrderComparatorAutoGtByKey};
//! use itertools::Itertools;
//! 
//! // define the ring operator
//! let ring_operator = BooleanFieldOperator::new();
//! 
//! // build a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVecSimple::new(
//!                         vec![ 
//!                             vec![                                     ], 
//!                             vec![                                     ],                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );    
//! 
//! // PART 1: SOLVE xA = b FOR x
//! // --------------------------
//! 
//! // define a sparse vector b
//! let b = vec![ (0,true), (1,true) ];        
//! 
//! // define the function that maps minor keys to matched major keys
//! // note: we have to place this function in a wrapper, for formatting reasons
//! let matching_from_keymin_to_keymaj = EvaluateFunctionFnMutWrapper::new(|x: usize| x+2 );        
//! 
//! // define an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//! let order_comparator_major_ascend = OrderComparatorAutoLtByKey::new();          
//! 
//! // create a solver to solve xA = b for x
//! let solution =  EchelonSolveMajorAscendWithMajorKeys::new(
//!                                                                 b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVecSimple, not on VecOfVecSimple itself
//!                                                                 matching_from_keymin_to_keymaj, 
//!                                                                 ring_operator.clone(), 
//!                                                                 order_comparator_major_ascend.clone(),
//!                                                             );
//! 
//! // check the solution, i.e. check that xA = b
//! let product = vector_matrix_multiply_major_ascend_simplified(solution, &matrix, ring_operator, order_comparator_major_ascend.clone());
//! assert_eq!( product.collect_vec(), b );
//! 
//! // PART 2: SOLVE Ay = c FOR y
//! // ---------------------------  
//! 
//! // define a sparse vector c
//! let c = vec![ (3,true), (2,true) ];        
//! 
//! // define the function that maps major keys to matched minor keys
//! // note: we have to place this function in a wrapper, for formatting reasons
//! let matching_from_keymaj_to_keymin = EvaluateFunctionFnMutWrapper::new(|x: usize| x-2 ); 
//! 
//! // define an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
//! let order_comparator_minor_descend = OrderComparatorAutoGtByKey::new();          
//! 
//! // create a solver to solve xA = b for x
//! let solution =  EchelonSolveMinorDescendWithMinorKeys::new(
//!                                                                 c.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVecSimple, not on VecOfVecSimple itself
//!                                                                 matching_from_keymaj_to_keymin, 
//!                                                                 ring_operator.clone(), 
//!                                                                 order_comparator_major_ascend.clone(), // note we do not invert the order comparator
//!                                                             );
//! 
//! // check the solution, i.e. check that Ay = c
//! let product = vector_matrix_multiply_minor_descend_simplified(solution, &matrix, ring_operator, order_comparator_minor_descend);
//! assert_eq!( product.collect_vec(), c );   
//! ```

use std::marker::PhantomData;

use crate::matrices::matrix_types::transpose::AntitransposeLazy;
use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend, IndicesAndCoefficients};
use crate::rings::operator_traits::{Semiring, Ring, DivisionRing, MinusOne};
use crate::utilities::iterators::merge::heap_of_iterators::{HitMerge, hit_bulk_insert, hit_merge_by_predicate};
use crate::utilities::partial_order::{StrictlyLess, OrderComparatorReverse};
use crate::utilities::functions::evaluate::EvaluateFunction;
use crate::entries::{KeyValSet, KeyValGet};
use crate::vectors::{operations::{Simplify, Scale, Transforms}, linear_combinations::LinearCombinationSimplified};
use crate::utilities::iterators::general::IterTwoType;
// use debugit::DebugIt as D;
use itertools::Itertools;








//  ======================================================================================
//  ======================================================================================
//  ======================================================================================
//  ASCENDING MAJOR SOLVE
//  ======================================================================================
//  ======================================================================================
//  ======================================================================================



//  ======================================================================================
//  ECHELON SOLVE ASCEND WITH MINOR KEYS
//  ======================================================================================




// ECHELON MAJOR ASCEND SOLVE W/ MINOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------

/// Iterates over the entries of the solution `x` to a matrix equation `xA = b`, where
/// (i) each entry of `x` representing a pair of form `(major_key, coefficient)` is replaced by an entry representing `(match(major_key), coefficient)`, where
/// `match(major_key)` is the matched minor key, and
/// (ii) entries appear in ascending order, according to minor index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the minor keys of `A`, and the entries of `b` appear in ascending order of index
/// - There is a partial bijection `match: keymin -> keymaj` from the minor keys of `A` to the major keys of `A`; concretely, this means a bijection from a subset of the minor keys to a subset of the major keys
/// - Whenever `keymaj = match( keymin )`, the leading entry of `A.view_major_ascend( keymaj )` has index `keymin`
/// - A solution to `xA = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymin, -alpha) )`.
pub struct EchelonSolveMajorAscendWithMinorKeys< 
                    ProblemVector,
                    MatchingFromKeyMinToKeyMaj,
                    Matrix,                   
                    RingOperator,
                    OrderComparator,
                > 
    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                         Clone + PartialEq,
        Matrix::SnzVal:                         Clone,   
        Matrix::ViewMajorAscend:                IntoIterator,        
        Matrix::ViewMajorAscendEntry:          KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,                   

{   
    ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
    matrix:                             Matrix, // the `A` in `Ax = b`
    entries_to_eliminate:               LinearCombinationSimplified<   
                                                IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                                        < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                                        ProblemVector,
                                                        // Matrix::ViewMajorAscendEntry,
                                                    >,
                                                Matrix::KeyMin, Matrix::SnzVal, RingOperator, OrderComparator 
                                            >, 
    matching_from_keymin_to_keymaj:     MatchingFromKeyMinToKeyMaj,
    // phantom_lifetime:                   PhantomData< Matrix::KeyMaj > 
} 




impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderComparator,
        >  

        EchelonSolveMajorAscendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,        
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                         Clone + PartialEq,
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                         Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               
        Matrix::ViewMajorAscendEntry:          KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,  

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolveMajorAscendWithMinorKeys`].
    pub fn new(
                b:                              ProblemVector,                
                A:                              Matrix,
                matching_from_keymin_to_keymaj: MatchingFromKeyMinToKeyMaj,
                ring_operator:                  RingOperator,
                order_comparator:               OrderComparator,
            )
            ->
            EchelonSolveMajorAscendWithMinorKeys<
                    ProblemVector,
                    MatchingFromKeyMinToKeyMaj,
                    Matrix,                    
                    RingOperator,
                    OrderComparator,
                >  
    {

        let entries_to_eliminate_unwrapped
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    IterTwoType::Iter2( b.into_iter() )   // wrap b in an IterTwoType enum so that it can be easily combined with  
                                        .scale( ring_operator.minus_one(), ring_operator.clone() ) 
                                ),
                            order_comparator,
                        )
                        .simplify( ring_operator.clone() );
        let entries_to_eliminate = LinearCombinationSimplified { linear_combination_simplified: entries_to_eliminate_unwrapped };
        EchelonSolveMajorAscendWithMinorKeys{
            ring_operator:                      ring_operator.clone(),
            matching_from_keymin_to_keymaj:     matching_from_keymin_to_keymaj,
            matrix:                             A,
            entries_to_eliminate:               entries_to_eliminate, 
            // phantom_lifetime:                   PhantomData,
        } 
    }




    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    pub fn remainder( self ) 
            -> 
            LinearCombinationSimplified<   
                    IterTwoType< // this enum allows us to treat two different iterator types as a single type
                            < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                            ProblemVector, // the vector we have tried to eliminate
                        >,
                    Matrix::KeyMin, Matrix::SnzVal, RingOperator, OrderComparator 
                >      
    {
        self.entries_to_eliminate
    }    
}





impl    <
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Iterator for    

        EchelonSolveMajorAscendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,        
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                         Clone + PartialEq,
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                         Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               
        Matrix::ViewMajorAscendEntry:          KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,  

{
    type Item = Matrix::ViewMajorAscendEntry;

    fn next( &mut self ) -> Option<Self::Item> { 

    // THIS ITERATOR PROCEED BY ELIMINATING ENTRIES IN `self.entries_to_eliminate`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_eliminate` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( mut entry_to_eliminate )       =   self.entries_to_eliminate.next() {

            println!("EchelonSolveMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: STARTING");

            // INDEX OF THE VECTOR THAT WILL EXECUTE THE ELIMINATION
            let keymaj_of_eliminating_viewmaj = self.matching_from_keymin_to_keymaj.evaluate_function( entry_to_eliminate.key() );

            println!("EchelonSolveMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: WAYPOINT 1");            

            // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
            let mut seed_of_eliminating_iterator    =   self.matrix.view_major_ascend( keymaj_of_eliminating_viewmaj ).into_iter();
            
            println!("EchelonSolveMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: WAYPOINT 2");            

            // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
            let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

            // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
            let scale_factor     =    self.ring_operator.negate(
                                                    self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                );
            
            // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
            // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate
            let eliminating_iterator       
                    =   IterTwoType::Iter1( seed_of_eliminating_iterator ) 
                            .scale( scale_factor.clone(), self.ring_operator.clone() );

            // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate`
            println!("EchelonSolveMajorAscendWithMinorKeys.next(): hit_bulk_insert: STARTING");
            hit_bulk_insert( &mut self.entries_to_eliminate.linear_combination_simplified.unsimplified, vec![eliminating_iterator] ); 
            println!("EchelonSolveMajorAscendWithMinorKeys.next(): hit_bulk_insert: ENDING");

            // NOW RE-USE `entry_to_eliminate` AS A CONTAINER FOR THE NEXT ENTRY IN THE INVERSE MATRIX
            entry_to_eliminate.set_val( scale_factor );

            return Some( entry_to_eliminate ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function

        } else {
            return None
        }
    }
}  







//  ======================================================================================
//  ECHELON MAJOR ASCEND SOLVE ASCEND WITH MAJOR KEYS
//  ======================================================================================




// ECHELON SOLVE W/ MAJOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------

/// Iterates over the entries of the solution `x` to a matrix equation `xA = b`, where entries appear in ascending order, according to (major) index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the minor keys of `A`, and the entries of `b` appear in ascending order of index
/// - There is a partial bijection `match: keymin -> keymaj` from the minor keys of `A` to the major keys of `A`; concretely, this means a bijection from a subset of the minor keys to a subset of the major keys
/// - Whenever `keymaj = match( keymin )`, the leading entry of `A.view_major_ascend( keymaj )` has index `keymin`
/// - A solution to `xA = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymaj, -alpha) )`.
pub struct EchelonSolveMajorAscendWithMajorKeys<
                    ProblemVector,
                    MatchingFromKeyMinToKeyMaj,
                    Matrix,                    
                    RingOperator,
                    OrderComparator,
                > 
    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                 Clone + PartialEq,
        Matrix::SnzVal:                 Clone,   
        Matrix::ViewMajorAscend:        IntoIterator,        
        Matrix::ViewMajorAscendEntry:   KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,                   

{   
    ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
    matrix:                             Matrix, // the `A` in `Ax = b`
    entries_to_eliminate:               LinearCombinationSimplified<   
                                                IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                                        < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                                        ProblemVector,
                                                        // Matrix::ViewMajorAscendEntry,
                                                    >,
                                                Matrix::KeyMin, Matrix::SnzVal, RingOperator, OrderComparator 
                                            >, 
    matching_from_keymin_to_keymaj:     MatchingFromKeyMinToKeyMaj,
} 


impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderComparator,
        >  

        EchelonSolveMajorAscendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,        
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                         Clone + PartialEq,
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                         Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               
        Matrix::ViewMajorAscendEntry:          KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,  

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolveMajorAscendWithMajorKeys`].
    pub fn new(
                b:                              ProblemVector,                
                A:                              Matrix,
                matching_from_keymin_to_keymaj: MatchingFromKeyMinToKeyMaj,
                ring_operator:                  RingOperator,
                order_comparator:               OrderComparator,
            )
            ->
            EchelonSolveMajorAscendWithMajorKeys<
                    ProblemVector,
                    MatchingFromKeyMinToKeyMaj,
                    Matrix,                    
                    RingOperator,
                    OrderComparator,
                >  
    {

        let entries_to_eliminate_unwrapped
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    IterTwoType::Iter2( b.into_iter() )   // wrap b in an IterTwoType enum so that it can be easily combined with  
                                        .scale( ring_operator.minus_one(), ring_operator.clone() ) 
                                ),
                            order_comparator,
                        )
                        .simplify( ring_operator.clone() );
        let entries_to_eliminate = LinearCombinationSimplified { linear_combination_simplified: entries_to_eliminate_unwrapped };
        EchelonSolveMajorAscendWithMajorKeys{
            ring_operator:                      ring_operator.clone(),
            matching_from_keymin_to_keymaj:     matching_from_keymin_to_keymaj,
            matrix:                             A,
            entries_to_eliminate:               entries_to_eliminate, 
        } 
    }


    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    pub fn remainder( self ) 
            -> 
            LinearCombinationSimplified<   
                    IterTwoType< // this enum allows us to treat two different iterator types as a single type
                            < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                            ProblemVector, // the vector we have tried to eliminate
                        >,
                    Matrix::KeyMin, Matrix::SnzVal, RingOperator, OrderComparator 
                >      
    {
        self.entries_to_eliminate
    }   
}





impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Iterator for    

        EchelonSolveMajorAscendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,        
        Matrix:                         OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                         Clone + PartialEq,
        Matrix::KeyMaj:                         Clone, // this is necessary to avoid a "move" error, since the following function uses a major key once to look up a major view (which consums the keymaj) and once to return the value of the keymaj to the user
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                         Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               
        Matrix::ViewMajorAscendEntry:          KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,  

{
    type Item = (Matrix::KeyMaj, Matrix::SnzVal);

    fn next( &mut self ) -> Option<Self::Item> { 

    // unsafe {
    //     println!("about to invoke EchelonSolveMajorAscendWithMajorKeys.nex(); problem vector =  {:?}", D(&self.entries_to_eliminate) );
    //     std::thread::sleep(std::time::Duration::from_secs(1));

    // }
    

    // THIS ITERATOR PROCEED BY ELIMINATING ENTRIES IN `self.entries_to_eliminate`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_eliminate` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( entry_to_eliminate )       =   self.entries_to_eliminate.next() {

            // unsafe {
            //     println!("entry to eliminate =  {:?}", D(&entry_to_eliminate) );
            // }
                            

            // INDEX OF THE VECTOR THAT WILL EXECUTE THE ELIMINATION
            let keymaj_of_eliminating_viewmaj = self.matching_from_keymin_to_keymaj.evaluate_function( entry_to_eliminate.key() );

            // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
            let mut seed_of_eliminating_iterator    =   self.matrix.view_major_ascend( keymaj_of_eliminating_viewmaj.clone() ).into_iter();
            
            // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
            let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

            // unsafe {
            //     println!("keymaj_of_eliminating_viewmaj =  {:?}", D(&keymaj_of_eliminating_viewmaj) );                
            //     println!("eliminating entry =  {:?}", D(&eliminating_entry) );
            //     let seed_of_eliminating_iterator_2    =   self.matrix.view_major_ascend( keymaj_of_eliminating_viewmaj.clone() ).into_iter().collect_vec();
            //     println!("eliminating iterator entries =  {:?}", D(&seed_of_eliminating_iterator_2) );
            // }            

            // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
            let scale_factor     =    self.ring_operator.negate(
                                                    self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                );
            
            // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
            // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate
            let eliminating_iterator       
                    =   IterTwoType::Iter1( seed_of_eliminating_iterator ) 
                            .scale( scale_factor.clone(), self.ring_operator.clone() );


            // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate`
            hit_bulk_insert( &mut self.entries_to_eliminate.linear_combination_simplified.unsimplified, vec![eliminating_iterator] ); 

            return Some( (keymaj_of_eliminating_viewmaj, scale_factor) ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function

        } else {
            return None
        }
    }
}  


















//  ======================================================================================
//  ======================================================================================
//  ======================================================================================
//  DESCENDING MINOR SOLVE
//  ======================================================================================
//  ======================================================================================
//  ======================================================================================







//  ======================================================================================
//  ECHELON MINOR DESCEND SOLVE WITH MAJOR KEYS
//  ======================================================================================




// ECHELON SOLVE W/ MAJOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------


/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, where
/// (i) each entry of `x`, which represents a pair of form `(minor_key, coefficient)`, is replaced by an entry representing `(match(minor_key), coefficient)`, where
/// `match(minor_key)` is the matched major key, and
/// (ii) entries appear in descending order, according to major index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the major keys of `A`, and the entries of `b` appear in descending order of index
/// - There is a partial bijection `match: keymaj -> keymin` from the major keys of `A` to the minor keys of `A`; concretely, this means a bijection from a subset of the major keys to a subset of the minor keys
/// - Whenever `keymin = match( keymaj )`, the leading entry of `A.view_minor_descend( keymin )` has index `keymaj`
/// - A solution to `Ax = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymaj, -alpha) )`.
/// 
/// Under the hood, `EchelonSolveMinorDescendWithMajorKeys` is simply a wrapper for a struct `EchelonSolveMajorAscendWithMinorKeys` which
/// contains a reference to the antitranspose of the matrix; see source code for full details.
pub struct EchelonSolveMinorDescendWithMajorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderComparator,
                > 
    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                         OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMaj:                 Clone + PartialEq,
        Matrix::SnzVal:                 Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::ViewMinorDescendEntry:         KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,                   

{   
    antitransposed_solver:      EchelonSolveMajorAscendWithMinorKeys< // the minor keys of the antitransposed matrix are the major keys of the un-antitransposed matrix
                                        ProblemVector,
                                        MatchingFromKeyMajToKeyMin,
                                        AntitransposeLazy< Matrix >,
                                        RingOperator,
                                        OrderComparatorReverse< OrderComparator >,
                                    >,
} 

impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderComparatorUnreversed,
        >  

        EchelonSolveMinorDescendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderComparatorUnreversed,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                         OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMaj:                 Clone + PartialEq,
        Matrix::SnzVal:                 Clone,   
        Matrix::ViewMinorDescend:       IntoIterator,        
        Matrix::ViewMinorDescendEntry: KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparatorUnreversed:      Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,   

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolveMajorAscendWithMinorKeys`].
    pub fn new(
                b:                              ProblemVector,                
                A:                              Matrix,
                matching_from_keymaj_to_keymin: MatchingFromKeyMajToKeyMin,
                ring_operator:                  RingOperator,
                order_comparator_unreversed:    OrderComparatorUnreversed,
            )
            ->
            EchelonSolveMinorDescendWithMajorKeys< 
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderComparatorUnreversed,
                >  
    {
        EchelonSolveMinorDescendWithMajorKeys{
            // Note that we use EchelonSolveMajorAscendWithMinorKeys, since the minor keys of the antitranspose are the major keys of the original matrix
            antitransposed_solver:   EchelonSolveMajorAscendWithMinorKeys::new(
                                            b,
                                            AntitransposeLazy::new( A ),
                                            matching_from_keymaj_to_keymin,
                                            ring_operator,
                                            OrderComparatorReverse::new( order_comparator_unreversed ),
                                        )
        }
    }




    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    pub fn remainder( self ) 
            -> 
            LinearCombinationSimplified<   
                    IterTwoType< // this enum allows us to treat two different iterator types as a single type
                            < Matrix::ViewMinorDescend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                            ProblemVector, // the vector we have tried to eliminate
                        >,
                    Matrix::KeyMaj, Matrix::SnzVal, RingOperator, OrderComparatorReverse< OrderComparatorUnreversed >
                >      
    {
        self.antitransposed_solver.remainder()
    }    
}





impl    <
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Iterator for    

        EchelonSolveMinorDescendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                         OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMaj:                         Clone + PartialEq,
        Matrix::SnzVal:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::ViewMinorDescendEntry:         KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,   

{
    type Item = Matrix::ViewMinorDescendEntry;

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitransposed_solver.next()
    }
}  










//  ======================================================================================
//  ECHELON MINOR DESCEND SOLVE WITH MINOR KEYS
//  ======================================================================================




// ECHELON SOLVE W/ MINOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------


/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, where entries appear in descending order, according to (minor) index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the major keys of `A`, and the entries of `b` appear in descending order of index
/// - There is a partial bijection `match: keymaj -> keymin` from the major keys of `A` to the minor keys of `A`; concretely, this means a bijection from a subset of the major keys to a subset of the minor keys
/// - Whenever `keymin = match( keymaj )`, the leading entry of `A.view_minor_descend( keymin )` has index `keymaj`
/// - A solution to `Ax = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymaj, -alpha) )`.
/// 
/// Under the hood, `EchelonSolveMinorDescendWithMajorKeys` is simply a wrapper for a struct `EchelonSolveMajorAscendWithMinorKeys` which
/// contains a reference to the antitranspose of the matrix; see source code for full details.
pub struct EchelonSolveMinorDescendWithMinorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderComparator,
                > 
    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                         OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMaj:                         Clone + PartialEq,
        Matrix::SnzVal:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::ViewMinorDescendEntry:         KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,                   

{   
    antitransposed_solver:      EchelonSolveMajorAscendWithMajorKeys< // the major keys of the antitransposed matrix are the minor keys of the un-antitransposed matrix
                                        ProblemVector,
                                        MatchingFromKeyMajToKeyMin,
                                        AntitransposeLazy< Matrix >,                 
                                        RingOperator,
                                        OrderComparatorReverse< OrderComparator >,
                                    >,
} 

impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderComparator,
        >  

        EchelonSolveMinorDescendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                         OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMaj:                         Clone + PartialEq,
        Matrix::SnzVal:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::ViewMinorDescendEntry:         KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,   

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolveMajorAscendWithMinorKeys`].
    pub fn new(
                b:                              ProblemVector,                
                A:                              Matrix,
                matching_from_keymaj_to_keymin: MatchingFromKeyMajToKeyMin,
                ring_operator:                  RingOperator,
                order_comparator:               OrderComparator,
            )
            ->
            EchelonSolveMinorDescendWithMinorKeys< 
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderComparator,
                >  
    {
        EchelonSolveMinorDescendWithMinorKeys{
            // Note that we use EchelonSolveMajorAscendWithMajorKeys, since the major keys of the antitranspose are the minor keys of the original matrix
            antitransposed_solver:   EchelonSolveMajorAscendWithMajorKeys::new(
                                            b,
                                            AntitransposeLazy::new( A ),
                                            matching_from_keymaj_to_keymin,
                                            ring_operator,
                                            OrderComparatorReverse::new( order_comparator ),
                                        )
        }
    }




    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    pub fn remainder( self ) 
            -> 
            LinearCombinationSimplified<   
                    IterTwoType< // this enum allows us to treat two different iterator types as a single type
                            < Matrix::ViewMinorDescend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                            ProblemVector, // the vector we have tried to eliminate
                        >,
                    Matrix::KeyMaj, Matrix::SnzVal, RingOperator, OrderComparatorReverse< OrderComparator >
                >      
    {
        self.antitransposed_solver.remainder()
    }    
}





impl    <
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Iterator for    

        EchelonSolveMinorDescendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                         OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMin:                         Clone,
        Matrix::KeyMaj:                         Clone + PartialEq,
        Matrix::SnzVal:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::ViewMinorDescendEntry:         KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,   

{
    type Item = (Matrix::KeyMin, Matrix::SnzVal);

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitransposed_solver.next()
    }
}  












//  =========================================================================
//  =========================================================================
//  =========================================================================
//  IMPLEMENTATIONS OF CLONE
//  =========================================================================
//  =========================================================================
//  =========================================================================


impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Clone for    

        EchelonSolveMajorAscendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Clone + Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     Clone + EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,        
        Matrix:                         Clone + OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                 Clone + PartialEq,
        Matrix::KeyMaj:                 Clone, // this is necessary to avoid a "move" error, since the following function uses a major key once to look up a major view (which consums the keymaj) and once to return the value of the keymaj to the user
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   Clone + KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,  
        Matrix::ViewMajorAscendIntoIter:Clone,        
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               
{
    fn clone(&self) -> Self {
        EchelonSolveMajorAscendWithMajorKeys{ 
                ring_operator: self.ring_operator.clone(),
                matrix: self.matrix.clone(),
                entries_to_eliminate: self.entries_to_eliminate.clone(),
                matching_from_keymin_to_keymaj: self.matching_from_keymin_to_keymaj.clone(),
            }
    }
}



impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Clone for    

        EchelonSolveMajorAscendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderComparator,
            > 

    where 
        ProblemVector:                  Clone + Iterator< Item = Matrix::ViewMajorAscendEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     Clone + EvaluateFunction< Matrix::KeyMin, Matrix::KeyMaj >,        
        Matrix:                         Clone + OracleMajorAscend + IndicesAndCoefficients,
        Matrix::KeyMin:                 Clone + PartialEq,
        Matrix::ViewMajorAscend:        IntoIterator,
        Matrix::ViewMajorAscendEntry:   Clone + KeyValGet < Matrix::KeyMin, Matrix::SnzVal > + KeyValSet < Matrix::KeyMin, Matrix::SnzVal >,  
        Matrix::ViewMajorAscendIntoIter:Clone,        
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        Matrix::SnzVal:                 Clone,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMajorAscendEntry >,               
{
    fn clone(&self) -> Self {
        EchelonSolveMajorAscendWithMinorKeys{ 
                ring_operator: self.ring_operator.clone(),
                matrix: self.matrix.clone(),
                entries_to_eliminate: self.entries_to_eliminate.clone(),
                matching_from_keymin_to_keymaj: self.matching_from_keymin_to_keymaj.clone(),
            }
    }
}




impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Clone for    

EchelonSolveMinorDescendWithMajorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderComparator,
                > 
    where 
        ProblemVector:                  Clone + Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     Clone + EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                         Clone + OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMaj:                 Clone + PartialEq,
        Matrix::SnzVal:                 Clone,   
        Matrix::ViewMinorDescend:       IntoIterator,        
        Matrix::ViewMinorDescendIntoIter:   Clone + IntoIterator,
        Matrix::ViewMinorDescendEntry:      Clone + KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                   Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,                    
{
    fn clone(&self) -> Self {
        EchelonSolveMinorDescendWithMajorKeys{ 
                antitransposed_solver: self.antitransposed_solver.clone()
            }
    }
}



impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderComparator,
        > 

        Clone for    

EchelonSolveMinorDescendWithMinorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderComparator,
                > 
    where 
        ProblemVector:                      Clone + Iterator< Item = Matrix::ViewMinorDescendEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:         Clone + EvaluateFunction< Matrix::KeyMaj, Matrix::KeyMin >,
        Matrix:                             Clone + OracleMinorDescend + IndicesAndCoefficients,
        Matrix::KeyMaj:                     Clone + PartialEq,
        Matrix::KeyMin:                     Clone, // this is necessary to avoid a "move" error, since the antitransposed solver uses a minor key once to look up a minor view (which consums the keymin) and once to return the value of the keymin to the user        
        Matrix::SnzVal:                     Clone,   
        Matrix::ViewMinorDescend:           IntoIterator,        
        Matrix::ViewMinorDescendIntoIter:   Clone,
        Matrix::ViewMinorDescendEntry:      Clone + KeyValGet < Matrix::KeyMaj, Matrix::SnzVal > + KeyValSet < Matrix::KeyMaj, Matrix::SnzVal >,      
        RingOperator:                       Clone + Semiring< Matrix::SnzVal > + Ring< Matrix::SnzVal > + DivisionRing< Matrix::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess <  Matrix::ViewMinorDescendEntry >,                    
{
    fn clone(&self) -> Self {
        EchelonSolveMinorDescendWithMinorKeys{ 
                antitransposed_solver: self.antitransposed_solver.clone()
            }
    }
}






























//  TESTS
//  ===========================================================================


#[cfg(test)]
mod doctring_tests {


    #[test]
    fn test_docstring_solve() {
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::operations::solve::echelon::{EchelonSolveMajorAscendWithMajorKeys, EchelonSolveMinorDescendWithMinorKeys};
        use crate::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified, vector_matrix_multiply_minor_descend_simplified};        
        use crate::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
        use crate::utilities::partial_order::{OrderComparatorAutoLtByKey, OrderComparatorAutoGtByKey};
        use itertools::Itertools;

        // define the ring operator
        let ring_operator = BooleanFieldOperator::new();

        // build a matrix A with an invertible upper triangular submatrix
        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![                                     ], 
                                    vec![                                     ],                                     
                                    vec![ (0,true),  (1,true), (2,true)       ], 
                                    vec![            (1,true), (2,true)       ],
                                    vec![                      (2,true)       ],                                       
                                ],
                            );    

        // PART 1: SOLVE xA = b FOR x
        // --------------------------

        // define a sparse vector b
        let b = vec![ (0,true), (1,true) ];        

        // define the function that maps minor keys to matched major keys
        // note: we have to place this function in a wrapper, for formatting reasons
        let matching_from_keymin_to_keymaj = EvaluateFunctionFnMutWrapper::new(|x: usize| x+2 );        

        // define an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
        let order_comparator_major_ascend = OrderComparatorAutoLtByKey::new();          

        // create a solver to solve xA = b for x
        let solution =  EchelonSolveMajorAscendWithMajorKeys::new(
                                                                        b.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVecSimple, not on VecOfVecSimple itself
                                                                        matching_from_keymin_to_keymaj, 
                                                                        ring_operator.clone(), 
                                                                        order_comparator_major_ascend.clone(),
                                                                    );
        
        // check the solution, i.e. check that xA = b
        let product = vector_matrix_multiply_major_ascend_simplified(solution, &matrix, ring_operator, order_comparator_major_ascend.clone());
        assert_eq!( product.collect_vec(), b );

        // PART 2: SOLVE Ay = c FOR y
        // ---------------------------  

        // define a sparse vector c
        let c = vec![ (3,true), (2,true) ];        

        // define the function that maps major keys to matched minor keys
        // note: we have to place this function in a wrapper, for formatting reasons
        let matching_from_keymaj_to_keymin = EvaluateFunctionFnMutWrapper::new(|x: usize| x-2 ); 
        
        // define an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
        let order_comparator_minor_descend = OrderComparatorAutoGtByKey::new();          

        // create a solver to solve xA = b for x
        let solution =  EchelonSolveMinorDescendWithMinorKeys::new(
                                                                        c.iter().cloned(), // the solver requires a bona-fide iterator, not an iterable
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVecSimple, not on VecOfVecSimple itself
                                                                        matching_from_keymaj_to_keymin, 
                                                                        ring_operator.clone(), 
                                                                        order_comparator_major_ascend.clone(), // note we do not invert the order comparator
                                                                    );
        
        // check the solution, i.e. check that Ay = c
        let product = vector_matrix_multiply_minor_descend_simplified(solution, &matrix, ring_operator, order_comparator_minor_descend);
        assert_eq!( product.collect_vec(), c );        
        
    }

}










//  TESTS
//  ===========================================================================


#[cfg(test)]
mod tests {

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use itertools::Itertools;
    use crate::{matrices::matrix_types::vec_of_vec::VecOfVecSimple, rings::operator_structs::field_prime_order::PrimeOrderFieldOperator, utilities::functions::evaluate::EvaluateFunctionFnMutWrapper};


    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Verifies that the output of a solver actually solves the problem Ax = b
    /// **IT IS IMPORTANT THAT THE MATRIX BE IN PARTIAL ECHECLON FORM**
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_echelon_solve_on_invertible_uppertriangular_matrix< 'a, MatchingFromKeyMinToKeyMaj, MatchingFromKeyMajToKeyMin, SnzVal, RingOperator >( 
                matrix:                             & VecOfVecSimple< usize, SnzVal >, 
                matching_from_keymin_to_keymaj:     MatchingFromKeyMinToKeyMaj,
                matching_from_keymaj_to_keymin:     MatchingFromKeyMajToKeyMin,
                ring_operator:                      RingOperator, 
                matrix_size:                        usize 
            ) 
        where   SnzVal:                         Clone + PartialEq + Ord + std::fmt::Debug,
                RingOperator:                   Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal > + Clone,
                MatchingFromKeyMinToKeyMaj:     Clone + EvaluateFunction< usize, usize >,
                MatchingFromKeyMajToKeyMin:     Clone + EvaluateFunction< usize, usize >,                
    {
        use crate::{matrices::operations::multiply::{vector_matrix_multiply_major_ascend_unsimplified, vector_matrix_multiply_major_ascend_simplified, vector_matrix_multiply_minor_descend_simplified}, utilities::partial_order::OrderComparatorAutoAnyType};

        


        // generate random vectors of 0's and 1's; try to get one vector each with 0, 1, 2, 3, and `matrix_size` nonzero entries in total
        for num_nonzeros in [0, 1, 2, 3, matrix_size].iter().map(|x| std::cmp::min(x, &matrix_size) ) {
            for vector_support in (0 .. matrix_size).combinations( *num_nonzeros ) {
                let mut vec = vector_support
                                .iter()
                                .cloned()
                                .map(|x| ( x, RingOperator::one() ))
                                .collect_vec();
                vec.sort();


                //  view-major solutions (corresponding to xA = b, where A is row-major)
                //  --------------------------------------------------------------------                

                // compute a view-major solution with minor keys
                let solution_major_with_minor_keys = EchelonSolveMajorAscendWithMinorKeys::new(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymin_to_keymaj.clone(),
                                        ring_operator.clone(),
                                        OrderComparatorAutoAnyType,                    
                                    )
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with minor keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_major_ascend_simplified( 
                            solution_major_with_minor_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderComparatorAutoAnyType
                        )
                );

                // compute a view-major solution with major keys
                let solution_major_with_major_keys = EchelonSolveMajorAscendWithMajorKeys::new(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymin_to_keymaj.clone(),
                                        ring_operator.clone(),
                                        OrderComparatorAutoAnyType,                    
                                    )
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with major keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_major_ascend_simplified( 
                            solution_major_with_major_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderComparatorAutoAnyType
                        )
                );    


                //  view-minor solutions (corresponding to Ax = b, where A is row-major)
                //  --------------------------------------------------------------------

                // REVERSE THE ORDER OF ENTRIES, SINCE WE ARE NOW DEALING WITH DESCENDING VIEWS
                vec.reverse(); 
                println!("REVERSED VECTOR FOR DESCENDING ECHELON SOLVERS: {:?}", & vec);

                // compute a view-minor solution with minor keys
                let solution_minor_with_minor_keys = EchelonSolveMinorDescendWithMinorKeys::new(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymaj_to_keymin.clone(),
                                        ring_operator.clone(),
                                        OrderComparatorAutoAnyType,                    
                                    )
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with minor keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_minor_descend_simplified( 
                            solution_minor_with_minor_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderComparatorAutoAnyType
                        )
                );

                println!("ECHELON SOLVE DESCEND: PASSED ESSENTIAL TEST 1");

                // compute a view-minor solution with major keys
                let solution_major_with_major_keys = EchelonSolveMinorDescendWithMajorKeys::new(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymaj_to_keymin.clone(),
                                        ring_operator.clone(),
                                        OrderComparatorAutoAnyType,                    
                                    )
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with major keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_minor_descend_simplified( 
                            solution_major_with_major_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderComparatorAutoAnyType
                        )
                );                    
                
                println!("ECHELON SOLVE DESCEND: PASSED ESSENTIAL TEST 2");                
                

            }
        }                                                               
    }


    
    #[test]
    fn test_echelon_solve_on_specific_matrices() {
        use num::rational::Ratio;        
        // use crate::matrices::matrix_oracle_traits::MajorDimension;
        use crate::rings::operator_structs::ring_native::DivisionRingNative;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;

        // Define the integer type for the rational numbers we'll use
        type IntegerType = isize;
        let modulus = 1049;

        // Define shorthand functions for creating new rational numbers
        let r = | n:IntegerType, d:IntegerType | Ratio::new( n, d );    
        let q = | n:IntegerType | Ratio::from_integer( n );   

        // Define matching function from Matrix::KeyMin (=usize) to Matrix::KeyMaj (=usize)
        let matching_keymin_to_keymaj_closure = | x: usize | -> usize { x };
        let matching_keymin_to_keymaj = EvaluateFunctionFnMutWrapper::new( matching_keymin_to_keymaj_closure );
        
        // Define the ring operators
        let ring_operator_q  =   DivisionRingNative::< Ratio<IntegerType> >::new();  
        let ring_operator_p  =   PrimeOrderFieldOperator::new(modulus);                  

        // Test individual matrices        
        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![ (0,q(1)) ]
                                ],
        );
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q.clone(), 1 );       

        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(1)), ], 
                                    vec![            (1,q(1)), ],
                                ],
        );
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q.clone(), 2 );       

        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
        );     
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q.clone(), 3 );        
        
        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
        );        
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q.clone(), 3 );           

        let matrix  =   VecOfVecSimple::new(
                                vec![ 
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
        );     
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q.clone(), 4 );         
        
        // MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        use rand::Rng;        // we use this module to generate random elements
        let matrix_size =   20;

        let mut rng = rand::thread_rng(); // this generates random integers
        let mut vec_of_vec = vec![];
        for majkey in 0 .. matrix_size {
            let coefficient_leading         =   rng.gen_range( 1 .. modulus );
            let mut new_vec     =   vec![ (majkey, coefficient_leading) ]; // start with a nonzero diagonal element            
            for q in majkey+1 .. matrix_size { // fill out the rest of this row of the matrix
                let coefficient   = rng.gen_range( 0 .. modulus );
                let flag = rng.gen_range(0usize .. 3usize);
                if      flag == 0 { new_vec.push( ( q, 0 )           ) }
                else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
                else              { continue }
            }
            vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        }

        let matrix  =   VecOfVecSimple::new(vec_of_vec); // formally wrap the matrix in a VecOfVec struct
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_p.clone(), matrix_size );                 

    }

    // #[test]
    // fn test_remainder() {
    //     !!! fill this in
    // }

}
