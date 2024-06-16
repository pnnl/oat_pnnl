//! U-match factorization for row-major oracles that implement `ViewMajorAscend` (does NOT give access to the columns of the matrices associated with the factorization).
//!
//! 
//! In some cases the user may find it onerous to implement the `OracleMinorDescend` trait for a new matrix class they
//! with to construct; however, they may still wish to obtain a U-match factorization (e.g. to solve a linear system of equations,
//! calculate rank, obtain a basis for the cokernel, etc.).
//! 
//! The purpose of this module is to provide such functionality.
//! 
//! 
//! This module is not maintained as regularly as [`umatch_rowmaor`](oat_rust::matrices::operations::umatch::umatch)
use itertools::Itertools;
use serde_json::ser::CompactFormatter;

use crate::entries::{KeyValSet, KeyValNew};
use crate::matrices::matrix_types::matching::{GeneralizedMatchingArrayWithMajorOrdinals, GeneralizedMatchingArray};
use crate::matrices::matrix_types::oracle_ref::OracleRef;
use crate::matrices::matrix_types::prepend_viewmaj::PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend;
use crate::matrices::matrix_types::scalar_matrices::ScalarMatrix;
use crate::matrices::matrix_types::vec_of_vec::{VecOfVec, VecOfVecSimple, VecOfVecSimpleFromBorrow};
use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMajor, OracleMinorDescend};
use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
use crate::utilities::iterators::general::{SkipUntil, IterTwoType, IterWrappedVec, OncePeekable};
use crate::utilities::iterators::merge::heap_of_iterators::{HitMerge, hit_bulk_insert, hit_merge_by_predicate};
use crate::entries::{KeyValGet};
use crate::rings::operator_traits::{Semiring, Ring, DivisionRing};
use crate::utilities::iterators::merge::two_iter_different_type::MergeTwoItersByOrderComparator;
use crate::utilities::partial_order::{StrictlyLess, OrderComparatorAutoLtByKey, OrderComparatorLtByKey};
use crate::utilities::sets::MapHasKeyOrSequenceHasElementRefWrapper;
use crate::vectors::linear_combinations::LinearCombinationSimplified;
use crate::vectors::operations::{Scale, Transforms, Simplify, OnlyIndicesInsideCollection, OnlyIndicesOutsideCollection, ChangeIndexSimple};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::{Cloned, Peekable, Once, Chain};
use std::marker::PhantomData;
use std::slice::{Iter};

use crate::matrices::operations::solve::echelon::{EchelonSolveMajorAscendWithMinorKeys};
use crate::matrices::operations::multiply::{ProductMatrixLazyMajorAscendSimplified, vector_matrix_multiply_major_ascend_simplified};
use crate::matrices::operations::transform_vector_wise::{OnlyKeyMinInsideCollection, OnlyKeyMinOutsideCollection};
use crate::matrices::operations::solve::triangle::TriangularSolveAscend;


// plan
// * remarks
//   * saving all the pivots (incl. obvious ones) need not be so bad in dim 1, since
//     there can't be any more pivots than edges
//   * why make SnzVal a parameter rather than an associated type: maybe it's hard to enforce the rule that Ring<>, Semiring<>, DivisionRing<> all have the same element type?
// * actions
//     * !!!! HAVEA **DEFAULT** UMATCH FACTORIZATION FUNCTION WHICH IS **SAFE** AND A SEPARATE FUNCTION WITH _unsafe THAT REMOVES SOME CHECKS
//     * !!!! IMPLEMENT THE FOLLOWING DESIGN DECISION: WE SHOULD ONLY ASK THE UMATCH STRUCT TO STORE ORDER COMPARATORS FOR MAJOR AND MINOR ENTRIES (NOT KEYS) BUT PROVIDE SAFE CONSTRUCTORS THAT ONLY TAKE ORDER COMPARATORS OF INDICES (NOT ENTRIES) AS INPUTS
//     * !!!! THERE ARE MORE EFFICIENT VERSIONS OF SOME FORMULAE TO RECONSTRUCT THE COMBS THAN REFLECTED IN THEOREM 6 OF THE UMATCH PAPER; for example, to get the matched rows of the of the codomain COMB, you can simply invert the matched block of in the inverse codomain COMB
//     * !!!! IMPLEMENT a version of MatrixOracleAscend on VecOfVec that returns references to keys and references to values; then use this together with a FilterChangeIndex that takes references for indices to eliminate unnecessary cloning (which happens right now) 
//     * !!!! the current implementation of CombDomain stores part of the output in a `IterWrappedVec`; there are lazier alternatives, e.g. creating a new struct that stores an iterator and a reference to the matching array, and scales coefficients as needed one at a time as they are returned
//     * !!!! rather than auto-implement OracleMajor for structs that implement OracleMajorAscend, create a wrapper object that does the auto-implementation; this avoids conflicting auto-implementations, e.g., if we ever wanted to auto-implement OracleMajor for somethng that implements OracleMajorDescend
//     * !!!! for C, first collect the elements of Row_c( A^{-1} ) into a Vec, and sort it; then use the elements of the vect to calculate the linera combination in the righthand block; then write a new vector for the lefthand block; this seems unavoidable due to the fact that the columns of A^{-1} and the columns of D have different key types
//     * !!!! create "OnlyPivot" and "OnlyNonpivot" objects that subselect pivot/nonpivot elements
//     * !!!! possibly definte wrapper objects for the major/minor views (just to cut down on the number of type parameters)
//     * !!!! ADD A CONVERSION FROM VECOFVEC TO VECOFVEC SIMPLE, AND BETWEEN THE REFERENCE VERSIONS
//     * !!!! Record this technical point: the HitMerge object stores "lower-in-order" iterators towards the front of the internally stored vector; therefore when placing a sequence of vectors in the heap (e.g. when takign a lienar combination of rows or R_{\rho \rho}^{-1}) it might be advantageous to do so in sorted order; this is relevant, for example, to the computation of the inverse of the codomain comb
//     * !!!! Since we must include an order comparator for major keys, it might be reasonable to include a safety check on the factorization algorithm to ensure that major keys appear in strictly descending order
//     * !!!! When computing a umatch decomposition (via the function codomain_comb_inv_off_diag_pivot_block_unsafecomparator), each matched row is sorted after its entries are computed; THERE MAY BE A WAY TO STORE THIS MATRIX IN UPPER TRIANGULAR FORM WITHOUT SORTYING -- FOR EXAMPLE, BY ORDERING ELEMENTS ACCORDING TO MINOR KEYS RATHER THAN MAJOR KEYS; THIS COULD SAVEA  LOT OF OPERATIONS, but I'm not sure if it's possible
//     * !!!! Discuss with Ryan -- when I write an implementation `impl OracleMajor for &'a T where T: OracleMajor` I get an error for &'a GeneralizedMatchingArrays but not for 
//     * search library for places where we make trait implementation requirements on a struct/trait implementation that could instead be levied on an individual function or functiosn within that implementation
//     * add section "Common challenges when implementing matrix oracles," where our approach to lifetimes vis-a-vis the vec-of-vec Rev lifetime problem
//     * resolve questions:
//     * remove lifetime parameter from oracle trait?
//       * are lifetime parameters really necessary for OracleMajorDescend?
//       * why do you need lifetimes for descend but not ascend?
//     * choose in a TRADEOFF:
//       * implement oracle trait for `&'a struct` -- this allows one to avoid putting the lifetime in the struct parameters
//       * be able to use auto_impl(&)
//     * (once trait specialization stabilizes) print the value that was unsuccessfully pushed to a matching array, in the resulting error message
//     * integer versus keymaj indexing of Ripiv
//       - integer indexing
//         - arguments in favor
//           - less memory when each entry takes a lot of memory (e.g. vecs)
//              - and entries only take a small amount of memory in very particular cases (namely low dimensions of simplicial complexes with relatively few vertices)
//              - also, tuples seem to take relatively little memory, c.f. https://stackoverflow.com/questions/65707981/why-is-the-size-of-a-tuple-or-struct-not-the-sum-of-the-members
//           - can use vec-of-vec data structure which is more efficient in time/memory than vec-of-vec 
//             - this makes for a more efficient triangular solve
//           - easier to transpose
//         - arguments against
//           - you have to choose between
//              - having inverted order of ordinals, or
//              - reversing order of ordinals; this entains reversing the order of 4 vectors and rewriting ordinals in a vec-of-vec as well as two hash-maps
//       - integer indexing with adaptor to reverse the indices of entries on the fly
//         - arguments in favor
//           - avoids the costly operation of reversing the entries in the whole matrix and reversing the vec-of-vec
//         - arguments against
//           - requires a custom data structure
//           - "more moving parts"
//       - keymaj indexing
//         - arguments in favor
//           - conceptual simplicity
//           - mathematically natural
//           - more efficient in the most common use case -- where key and val are stored in the same usize integer 
//           - same number of hashmaps stored, but no vectors (save 3 vectors worth of data)
//           - (weak) 
//         - arguments against
//           - **this necesitates storing a hashmap from KeyMaj to Vec.  However we also have to store a hashmap from KeyMaj to other things for the Matching Array**
//        
//   * the current implementation uses Peekable structs to pop out elts of iterators before they're added to hit_merge; this results in
//     essentially storing one extra entry per iterator; consider exploring ways to drop this extra baggage, e.g. by writing a new hit_merge
//     function that takes a sequence of `PutBack`s as input, instead of a sequence of iterators
//     * note also that since we put the peekable structs inside HeadTail structs for HitMerge, there really is some redundancy
//   * the current implementation clones the entries from a selected row of Ripiv each time it eliminates a leading entry; perhaps this can be avoided
//   * add unit tests for the matching arrays
//   * rename `OrderingPredicate` and/or `ordering_predicate`
//   * !!! rather than concatenating vectors to form parts of a row/col of R/Ri/C/Ci, use MERGE (can effect this by using an enum to effectively make a union of two different iterator types)
//   * add trait with auto-implementation for ring operators to generate 1 and 0 -- put this in the ring file
//   * add trait with auto-implementation for ring operators to generate minus 1 -- put this in the ring file
//   * consider imposing consistet rule re: order of arguments (e.g. high to low frequency)
//   * (in general) go back through existing code and remove SnzVal parametr (reducundant with the RingOperator parameter)
//   * (in general) make two versions of any given object: one private+safe, and one public+unsafe+modifiable
//   * (maybe) modify the inversion struct so that it doesn't constantly re-write the "head" entry
//   * (development goal) make it possible for people to write their own matrix oracle in python / export to exhact?
//   * (GOAL) make tutorials / demos a central, constant part of development goal
//   * enforce a rule that the structural nz entries in a VecOfVec must be nonzero?
// * structs
//   * Ri_piv -- probably want to make this a REFERENCE since some computations call for mulptiple copies
//   * VecOfVecUnitriangular
//   * ReductionMatrixConstructor
//      * iterates over "rows" of the reduction matrix in order
//      * each item is an iterator that iterates over a row
//      * after a row is generated, it's added to the internally stored VecOfVec, and the remaining 
//        coefficient gets saved to the internally-stored matching array
//      * potential differences with the ViewAscendingInverse struct
//        * modifies the matrix oracle from which it draws eliminating iterators
//        * writes its output to an internally stored attribute
//        * stores other, additional information at each step, too (e.g. the matching array)
//      
// DONE
//   * !!! consider switching order of KeyMin and KeyMaj



// -------
















//  =========================================================================================================
//  CALCULATION OF THE (COMPRESSED) MATCHING ARRAY AND CODOMAIN COMB
//  =========================================================================================================





//  HELPER FUNCTION 
//  ----------------------------------------------------------------------------------

/// Helper function for `try_writing_next_viewmaj_of_codomain_comb_to_buffer`.  That function
/// involves multiplying a row of a codomain COMB by a mapping array; sometimes entries
/// at the beginning of the resulting row vector can be deleted without penalty; the 
/// current function performs this deletion on a (scalar multiple) of an individual row of the 
/// mapping array.
fn codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array<
        Matrix,
        KeyMin,
        KeyMaj,                    
        SnzVal,
        ViewMajorAscend, 
        RingOperator,
        OrderComparatorMinor,
    > 
    (
        codomain_comb_inv_entry:            ( usize, SnzVal ),
        scale_factor:                       SnzVal,        
        truncation_limit:                   & ViewMajorAscend::Item,
        array_mapping:                      & Matrix,     
        array_matching:                     & GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,        
        ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
        mut order_comparator_minor:               OrderComparatorMinor,            
    ) 
    ->
    Peekable<
            Scale < 
                    ViewMajorAscend::IntoIter,  // a major view of the mapping array
                    KeyMin,
                    RingOperator, 
                    SnzVal,
                >, 
        >
    where 
        Matrix:                 OracleMajorAscend< KeyMaj, ViewMajorAscend >,
        KeyMin:                 Clone + Hash + PartialEq + std::cmp::Eq,
        KeyMaj:                 Clone + Hash + PartialEq + std::cmp::Eq,         
        SnzVal:                 Clone,               
        RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
        OrderComparatorMinor:        Clone + StrictlyLess <  ViewMajorAscend::Item >,  
        ViewMajorAscend:        IntoIterator,        
        ViewMajorAscend::Item:  KeyValSet < KeyMin, SnzVal >,  
{

    let ( codomain_comb_inv_entry_index, codomain_comb_inv_entry_coeff ) = codomain_comb_inv_entry;

    // we multiply the entry coefficient `codomain_comb_inv_entry.val()` by the `scale_factor`, 
    let scale_factor_aggregate = ring_operator.multiply( scale_factor, codomain_comb_inv_entry_coeff ); 
    
    // then we multiply `codomain_comb_inv_entry.val() * scale_factor_aggregate` by the view of the **mapping array**
    // (not matching array) that is indexed by `codomain_comb_inv_entry.key()`
    let iter = array_mapping // <-- note: mapping, **not matching**
                    .view_major_ascend( 
                            // this picks out the row index corresponding to the entry; this index might *not* be an integer, though `codomain_comb_inv_entry.key()` *is* an integer
                            array_matching.ordmaj_to_keymaj( codomain_comb_inv_entry_index )  
                        )
                    .into_iter()
                    .scale( scale_factor_aggregate, ring_operator.clone() )
                    // making the iterator peekable will help us prune some entries
                    .peekable() 
                    // we discard the entries of this iterator that have index <= `leading_entry_to_eliminate.key()` (because those entries would eventually be cancelled out, anyway -- this saves possibly a lot of heap insert + pop operations)
                    // in doing so, we must take care not to throw out too many entries by accident; that is why peekable is used                    
                    .skip_until( |x| order_comparator_minor.strictly_less( truncation_limit, x  )  );
    iter 
}


//  HELPER FUNCTION 
//  ----------------------------------------------------------------------------------


/// A major ascending view of an attempted solution to Ax = -b, calculated via cohomology
/// described in "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// 
/// The output of this function is best explained in terms of an example.  Let M be a matrix,
/// and suppose that we have partially executed the reduction algorithm described in 
/// "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek to 
/// obtain a U-match factorization.  This algorithm proceeds by reducing rows iteratively, 
/// from bottom to top.  To process a new row, the reduction procees by adding a sequence of
/// linear combinations of lower rows.  The addition of each linear combination eliminates a
/// new leading entry in the vector to be reduced.  The output of this reduction algorithm is a
/// sequence of index-coefficient pairs related to the sequence of linear combinations which we use to
/// reduce the vector.  The `ViewMajorAscendOfUmatchReducingVector` returns this sequence of 
/// entries in order, as an interator, assuming that 
/// - the `array_mapping` field stores matrix `M` in row-major form
/// - the portions of the matching array and inverse-codomain-COMB are stored in the fields
/// `array_matching` and `array_codomain_comb_inv_off_diag`.
/// - the field `order_comparator_minor` correctly determines which minor key precedes which 
/// 
/// # How to identify pivot indices
/// 
/// This struct modifies one of its private fields, the object `entries_to_elim_simplified_heap`.
/// The modified iterator represents the result of reducing `entries_to_elim_simplified_heap` via
/// clearing operations, as described in "U-match factorization, ..." by Hang, Ziegelmeier, Giusti,
/// and Henselman-Petrusek.
/// - if calling `entries_to_elim_simplified_heap.next()` returns `Some( minkey, coeff )` after
///  `ViewMajorAscendOfUmatchReducingVector` returns its last item,
///   then (assuming the struct is initialized correctly vis-a-vis the reduction algorithm
///   in  "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek) 
///   `minkey` is a pivot index, as is the corresponding `keymaj` (that is, the major key of
///   the mapping array, which indexes the vector that initialized `entries_to_elim_simplified_heap`).
/// - otherwise, `entries_to_elim_simplified_heap.next()` returns `None`, and we do not have a
///   new pivot index.
/// 


/// Needs new docs.
/// 
/// Note: entries are not inserted into the buffer vector in sorted order, and this function
/// does not sort the vector before returning it.
fn try_writing_next_viewmaj_of_codomain_comb_to_buffer< 
        Matrix,
        KeyMin,
        KeyMaj,                    
        SnzVal,
        ViewMajorAscend, 
        RingOperator,
        OrderComparatorMinor,
    > 
    (
        codomain_comb_inv_off_diag_view_buffer:      &mut Vec< ( usize, SnzVal ) >,          
        entries_to_elim_simplified_heap:    & mut Simplify <
                                                        HitMerge < 
                                                                Peekable<
                                                                        // the HitMerge struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                                        Scale < 
                                                                                ViewMajorAscend::IntoIter,  // a major view of the mapping array
                                                                                KeyMin,
                                                                                RingOperator, // the ring_operator operator
                                                                                SnzVal,
                                                                            >, 
                                                                    >,
                                                                // the thing that declares whether one major key comes before of after another    
                                                                OrderComparatorMinor                                                                
                                                            >,
                                                        KeyMin,
                                                        RingOperator,
                                                        SnzVal,
                                                    >,          
        array_mapping:                      & Matrix,    
        array_matching:                     & GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,
        array_codomain_comb_inv_off_diag:   & Vec< Vec< (usize, SnzVal) > >,  // the off-diagonal part of the (block indexed by pivot indices of) inverse of the codomain COMB
        // array_codomain_comb_inv_off_diag:   & VecOfVec<  // the off-diagonal part of the (block indexed by pivot indices of) inverse of the codomain COMB
        //                                                 ( usize, SnzVal ), 
        //                                                 OrderComparatorAutoLtByFirstTupleEntry, 
        //                                             >,          
        ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
        order_comparator_minor:             OrderComparatorMinor,            
    ) 
    ->

    Option< ViewMajorAscend::Item >

    where 
        Matrix:                 OracleMajorAscend< KeyMaj, ViewMajorAscend >,
        KeyMin:                 Clone + Hash + PartialEq + std::cmp::Eq + Debug, // !!! remove debug
        KeyMaj:                 Clone + Hash + PartialEq + std::cmp::Eq + Debug,      // !!! remove debug   
        SnzVal:                 Clone + Debug,                // remove debug
        RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
        OrderComparatorMinor:   Clone + StrictlyLess <  ViewMajorAscend::Item >,  
        ViewMajorAscend:        IntoIterator,        // !!!!! REMOVE THE DEBUG + CLONE
        ViewMajorAscend::Item:  KeyValSet < KeyMin, SnzVal > + Debug,       // !!!!! REMOVE THE DEBUG + CLONE REQUIREMENT WHEN DONE DEBUGGING!
        HitMerge<Peekable<Scale<<ViewMajorAscend as IntoIterator>::IntoIter, KeyMin, RingOperator, SnzVal>>, OrderComparatorMinor>: Clone, // !!! remove this later
{
    println!("initiating construction of next row");  // !!! delete later
    println!("array_matching {:?}", array_matching);  // !!! delete later    

    // we use this while loop because a for-loop results in a borrowing conflict (the conflict arrises because the for loop 
    // iteratos over `entries_to_elim_simplified_heap`, but we modify `entries_to_elim_simplified_heap` within the for-loop)
    while let Some( leading_entry_to_eliminate ) = entries_to_elim_simplified_heap.next() {

        println!("WHILE LOOP CHECKPOINT: leading_entry_to_eliminate: {:?}", &leading_entry_to_eliminate); // !!! REMOVE LATER
        println!("WHILE LOOP CHECKPOINT: entries_to_elim_simplified_heap: {:?}", entries_to_elim_simplified_heap.clone().collect_vec() ); // !!! REMOVE LATER        
    
        // IF THE MINOR INDEX OF THE ENTRY IS MATCHED, THEN WE CAN ELIMINATE
        if let Some( ordmaj_matched )       =   array_matching.keymin_to_ordmaj( & leading_entry_to_eliminate.key() ) {

            let scale_factor                =   ring_operator.negate(
                                                                ring_operator.divide(
                                                                        leading_entry_to_eliminate.val(),                                                                            
                                                                        array_matching.ordmaj_to_snzval( ordmaj_matched ),
                                                                    )
                                                            );
            

            // add the new (scaled and truncated) iterators to `entries_to_elim_simplified_heap`
            hit_bulk_insert( 
                    &mut entries_to_elim_simplified_heap.unsimplified,  
                    array_codomain_comb_inv_off_diag[ ordmaj_matched ]
                        .iter()
                        .cloned()
                        .chain(  // append a diagonal entry with coefficient 1 to the iterator
                                std::iter::once( 
                                        ( ordmaj_matched.clone(), RingOperator::one() ) 
                                    ) 
                            )
                        .map(   
                                | codomain_comb_inv_entry | 
                                {
                                println!(
                                        "codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array:      {:?}",
                                        codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array(
                                                codomain_comb_inv_entry.clone(),
                                                scale_factor.clone(),
                                                & leading_entry_to_eliminate, // truncation_limit:                   ViewMajorAscend::Item,
                                                array_mapping, //                     & Matrix,     
                                                array_matching, //                    & GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,        
                                                ring_operator.clone(), //                    RingOperator, // the operator for the coefficient ring_operator    
                                                order_comparator_minor.clone() //               OrderComparatorMinor,                                      
                                            )
                                            .collect_vec()
                                    );

                                codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array(
                                        codomain_comb_inv_entry,
                                        scale_factor.clone(),
                                        & leading_entry_to_eliminate, // truncation_limit:                   ViewMajorAscend::Item,
                                        array_mapping, //                     & Matrix,     
                                        array_matching, //                    & GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,        
                                        ring_operator.clone(), //                    RingOperator, // the operator for the coefficient ring_operator    
                                        order_comparator_minor.clone() //               OrderComparatorMinor,                                      
                                    )
                                }
                            )
                );                                     

            // push ( scale_factor * eliminating_row_with_diagonal_entry_re-added ) to the buffer; we do this in two steps
            // step 1: push the diagonal entry == (the scale factor and the index of the codomain COMB view that we used to perform the elimination) to the buffer
            codomain_comb_inv_off_diag_view_buffer.push(   ( ordmaj_matched, scale_factor.clone() )   );
            // step 1: push off-diagonal entries
            codomain_comb_inv_off_diag_view_buffer
                    .extend(  
                            array_codomain_comb_inv_off_diag[ ordmaj_matched ]
                                .iter()
                                .cloned() // note: we have to do this b/c `codomain_comb_inv_off_diag_view_buffer` is a `VecOfVec`, not a `VecOfVecSimple`
                                .scale( scale_factor, ring_operator.clone() )
                        );            

        } else {

            // REMARK: in this case we cannot clear `leading_entry_to_eliminate` via elementary operations; 
            // therefore `leading_entry_to_eliminate` is a matched (i.e. pivot) entry

            println!("terminating construction of next row"); // !!! delete later
            return Some( leading_entry_to_eliminate );
        }   
    }

    println!("terminating construction of next row"); // !!! delete later
    return None    
}





//  FUNCTION(S) TO EXTRACT Ripiv (the pivot portion of R^{-1}[pivot_indices, pivot_indices])
//  ------------------------------------------------------------------------------------------

/// Returns the block of the codomain COMB indexed by pivot indices, and a representation of the matching matrix.
/// 
/// For details on this factorization, see [this preprint](https://arxiv.org/pdf/2108.08831.pdf).
/// 
/// This U-match factorization is computed via the standard "cohomology algorithm."
/// 
/// **It is important that `order_comparator_minor` compares order of two entries based *solely* on
/// the value of the associated indices (not the associated coefficients).  This is easy
/// to get wrong, and it's hard to debug, so we are keeping the function private for now**
pub fn codomain_comb_inv_off_diag_pivot_block_unsafecomparator< 'a, Matrix, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, IterKeyMaj, OrderComparatorOfEntries >
    ( 
            array_mapping:                      &'a Matrix,  // use of a reference poses no problem for future calls to this function, since output has no references        
            iter_keymaj:                        IterKeyMaj,
            ring_operator:                      RingOperator,            
            order_comparator_of_entries:        OrderComparatorOfEntries,
            // mut order_comparator_of_keymaj:     OrderComparatorMajor,
    ) 
    -> 
    ( 
        VecOfVecSimple< usize, SnzVal >, 
        // VecOfVec< ( usize, SnzVal ), OrderComparatorAutoLtByFirstTupleEntry  >,         
        GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,
    )
    where   Matrix:                     OracleMajorAscend< KeyMaj, ViewMajorAscend >,
            IterKeyMaj:                 Iterator < Item = KeyMaj >,
            KeyMin:                     Clone + Hash + std::cmp::Eq + Debug, 
            KeyMaj:                     Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            SnzVal:                     Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:               Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
            ViewMajorAscend:            IntoIterator + Clone, // !!! remove clone
            ViewMajorAscend::Item:      KeyValSet< KeyMin, SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderComparatorOfEntries:   Clone + StrictlyLess <  ViewMajorAscend::Item >, // !!! remove clone
            // OrderComparatorMajor:       StrictlyLess<  KeyMaj >,
            HitMerge<Peekable<Scale<<ViewMajorAscend as IntoIterator>::IntoIter, KeyMin, RingOperator, SnzVal>>, OrderComparatorOfEntries>: Clone // !!!! remove this

{
    
    let mut entries_to_elim_simplified_heap    =   HitMerge::new( order_comparator_of_entries.clone() )
                                                                            .simplify( ring_operator.clone() ); // an empty (simplified) merger of iterators
    let mut array_matching          =   GeneralizedMatchingArrayWithMajorOrdinals::new(); // an all-zero generalized matching array
    let mut array_codomain_comb_inv_off_diag: Vec< Vec< ( usize, SnzVal ) > >  =   Vec::new();    
    // let mut array_codomain_comb_inv_off_diag: VecOfVec<(usize, SnzVal), OrderComparatorAutoLtByFirstTupleEntry >       =   VecOfVec::new(  vec![], OrderComparatorAutoLtByFirstTupleEntry::new()  );
    let mut codomain_comb_inv_off_diag_view_buffer   =   Vec::new();
    let mut codomain_comb_inv_off_diag_view_buffer_simplified   =   Vec::new();    

    // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
    // set up a variable to track the last major index produced by `iter_keymaj`; we will use this to ensure that the major keys returned by `iter_keymaj` appear in strictly ascending order
    // let mut prior_keymaj_opt = None;

    // build the (pivot block of the) codomain COMB row by row
    for keymaj in iter_keymaj {

        // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
        // check that this major key is strictly greater than the last
        // if let Some( ref prior_keyaj ) = prior_keymaj_opt {
        //     match order_comparator_of_keymaj.strictly_less( prior_keyaj, &keymaj ) {
        //         true => {  },
        //         false => { panic!("Umatch factorization cannot be completed because the sequence of major indices is not strictly ascending") }
        //     }               
        // }
        // prior_keymaj_opt.replace( keymaj.clone() );

        println!("keymaj: {:?}", &keymaj );

        // clear the collection of entries to eliminate
        entries_to_elim_simplified_heap.unsimplified.clear();        

        // insert the sequence of entries in row `keymaj`
        entries_to_elim_simplified_heap.unsimplified.insert_one_iter(
                array_mapping.view_major_ascend( keymaj.clone() )
                    .into_iter()
                    .scale( RingOperator::one(), ring_operator.clone() ) // !!! might be more efficient to use a two-type iter than to perform this multiplication, which only serves to make the iterator compatible with the HitMerge struct
                    .peekable()
        );

        codomain_comb_inv_off_diag_view_buffer.clear();
        let leading_entry_uneliminable_opt =
        try_writing_next_viewmaj_of_codomain_comb_to_buffer(
                & mut codomain_comb_inv_off_diag_view_buffer,
                & mut entries_to_elim_simplified_heap,        
                  array_mapping,
                & array_matching,
                & array_codomain_comb_inv_off_diag,
                ring_operator.clone(),
                order_comparator_of_entries.clone(),
            );

        match leading_entry_uneliminable_opt {
            // in this case we do not have a pivot
            Some( leading_entry_uneliminable )    => { 
                let pivot_keymin        =   leading_entry_uneliminable.key();
                let pivot_coeff         =   leading_entry_uneliminable.val();
                println!("leading_entry_uneliminable {:?}", &leading_entry_uneliminable);              
                // update the matching array
                array_matching.push( pivot_keymin, keymaj, pivot_coeff ); 
                // sort the buffer vector
                codomain_comb_inv_off_diag_view_buffer.sort_by( |a, b| b.0.cmp( &a.0 ) ); // note this yields DESCENDING ORDER -- which is what we want, because the order of indices is currently inverted (the greatest keymaj gets the lowest ordinal, and the least keymaj gets the highest ordinal; we will correct this later)
                println!("THE BUFFER:      {:?}", &codomain_comb_inv_off_diag_view_buffer);
                // simplify the sequence of entries in the buffer (note that the buffer itself could have entries with duplicate indices)
                // !!! hypothetically there's a opportunity for a small optimization here, where we simultaneously simplify and sort at the same time; this would avoid a few sorting operations
                let mut iter_to_push     =   codomain_comb_inv_off_diag_view_buffer
                                                            .iter()
                                                            .cloned()
                                                            .peekable()
                                                            .simplify( ring_operator.clone() );
                
                // write the simplified sequence to a buffer vector (this avoids having to reallocate, later, since we can't know the length of the simplified sequence 'till we iterate over it)
                codomain_comb_inv_off_diag_view_buffer_simplified.clear();
                codomain_comb_inv_off_diag_view_buffer_simplified.extend( iter_to_push );
                
                // write the simplified sequence to a new vector of exactly the right length
                let num_entries = codomain_comb_inv_off_diag_view_buffer_simplified.len();
                let mut vec_to_push     =   Vec::with_capacity( num_entries );
                vec_to_push.extend( codomain_comb_inv_off_diag_view_buffer_simplified.drain( 0..num_entries ) );
                
                // update the vector that stores the (off diagonal entries of the) "pivot" part of the inverse of the codomain comb
                array_codomain_comb_inv_off_diag.push( vec_to_push );
            }, 

            // in this case we don't have a pivot entry, so there is nothing to do
            None   => {}
        }
    }
    
    // if there are no pairs, return empty matrices (otherwise we will need to perform some transformations)
    if array_codomain_comb_inv_off_diag.len() == 0 { return ( VecOfVecSimple::new(array_codomain_comb_inv_off_diag), array_matching ) }

    // remove excess capacity
    array_codomain_comb_inv_off_diag.shrink_to_fit();
    

    // reverse the order of rows (because they are currently inverted, since we started with the bottom row and worked up)
    array_codomain_comb_inv_off_diag.reverse();

    // invert the ordinal used to index each entry
    let num_matched_pairs_minus_one = array_codomain_comb_inv_off_diag.len() - 1;
    for row_vec in array_codomain_comb_inv_off_diag.iter_mut() {  
        for entry in row_vec.iter_mut() {
            entry.set_key( num_matched_pairs_minus_one - entry.key() )
        }
    }

    // invert the ordinals of entries in the matching array
    array_matching.reverse();        

    // return the (off-diagonal entries of the pivot block of the) codomain COMB
    return ( VecOfVecSimple::new(array_codomain_comb_inv_off_diag), array_matching )

}


/// Calculates two quantities: (1) the square submatrix of the inverse of the codomain comb 
/// indexed by row pivot indices, and (2) the matching array.
/// 
/// # Design notes
/// 
/// We use `VecOfVecSimple< usize, SnzVal >` to store the pivot block of the inverse of 
/// the codomain comb rather than a `VecOfVec< ... >` for the following reasons: (1) we typically
/// need 
/// 
pub fn get_codomain_comb_inv_off_diag_pivot_block< 'a, Matrix, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, IterKeyMaj, OrderComparatorMinor, OrderComparatorMajor >
    ( 
            array_mapping:                      &'a Matrix,             
            iter_keymaj:                        IterKeyMaj,
            ring_operator:                      RingOperator,            
            order_comparator_for_keymin:        OrderComparatorMinor,
            order_comparator_for_keymaj:        OrderComparatorMajor,            
    ) 
    -> 
    ( 
        // VecOfVec< ( usize, SnzVal ), OrderComparatorAutoLtByFirstTupleEntry  >, 
        VecOfVecSimple< usize, SnzVal >,
        GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,
    )
    where   Matrix:                 OracleMajorAscend< KeyMaj, ViewMajorAscend >,
            IterKeyMaj:             Iterator < Item = KeyMaj >,
            KeyMin:                 Clone + Hash + std::cmp::Eq + Debug, 
            KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            SnzVal:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
            ViewMajorAscend:        IntoIterator + Clone, // !!! remove clone
            ViewMajorAscend::Item:  KeyValSet< KeyMin, SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderComparatorMinor:   Clone + StrictlyLess <  KeyMin >, // !!! remove clone
            OrderComparatorMajor:   StrictlyLess <  KeyMaj >, // !!! remove clone            
            HitMerge<Peekable<Scale<<ViewMajorAscend as IntoIterator>::IntoIter, KeyMin, RingOperator, SnzVal>>, OrderComparatorLtByKey< KeyMin, SnzVal, ViewMajorAscend::Item, OrderComparatorMinor>>: Clone // !!!! remove this

{
    // let order_comparator_for_entries_with_minor_keys = 
    //         |x: &ViewMajorAscend::Item, y: &ViewMajorAscend::Item | 
    //             order_comparator_for_keymin.strictly_less( &x.key(), &y.key() );

    let order_comparator_for_entries_with_minor_keys = OrderComparatorLtByKey::new( order_comparator_for_keymin );
    // let order_comparator_for_entries_with_major_keys = OrderComparatorLtByKey::new( order_comparator_for_keymaj );
    
    let ( array_codomain_comb_inv_off_diag_pivot_block, array_matching ) = 
    codomain_comb_inv_off_diag_pivot_block_unsafecomparator(
            array_mapping,
            iter_keymaj,
            ring_operator,
            order_comparator_for_entries_with_minor_keys,
            // order_comparator_for_keymaj,            
    );

    ( array_codomain_comb_inv_off_diag_pivot_block, array_matching )
}





















//  =========================================================================================================
//  U-MATCH OBJECT
//  =========================================================================================================





/// A (compressed representation of) a U-match factorization; all matrices concerned are represented in row-major format.
/// 
/// Internally, this object stores a copy of the matrix `R_{\rho \rho}` defined in Hang et al., 
/// "Umatch factorization: ..." 2021 (paper link)[https://arxiv.org/abs/2108.08831] **with its diagonal elements deleted**.  
/// More precisely, it contains a `VecOfVecSimple< usize, SnzVal >`, whose `k`th row contains the off-diagonal elements of
/// `R_{\rho \rho}`; where we replace each index `\rho_i` with the corresponding integer `i`, so that elements of the 
/// `VecOfVecSimple` are tuples of form `(usize, SnzVal)`.
/// 
/// 
/// # Design notes
/// 
/// **Why keep `ViewMajorAscend` as a generic type parameter?**  
/// 
/// **Why store the "off diagonal seed" of the codomain COMB as a `VecOfVecSimple` instead of a `VecOfVec`?**  Because we want to wrap this
/// struct in a `PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend` struct; this wrapper suppies the missing diagonal entries of,
/// the seed, thus forming the complete seed of the codomain COMB, which is a primitive and fundamental object for this 
/// computational library.  If we wanted to create an oracle for the seed that exposed entries as *references* to tuples,
/// then we would have to store the diagonal entries in memory.  However,  [`numerical experiments`](https://arxiv.org/pdf/2108.08831.pdf) 
/// have shown that the number of diagonal entries often eclipses the number of off-diagonal elements.  Thus storing 
/// the diagonal entries might incur a nontrivial memory cost.
/// 
/// **Why store the "off diagonal seed" of the codomain COMB as a `VecOfVecSimple` instead of a `Vec< Vec< (usize, SnzVal) > >`?**
/// Because `VecOfVecSimple` comes with certain guarantees about order of entries, which allows us to safely implement 
/// the `OracleMajorAscend` and `OracleMajorDescend` traits.
/// 
/// **Remark** One can always obtain a `VecOfVecFromBorrow` from a `VecOfVecSimple` via the
/// [`from_vecofvecsimple`](crate::matrices::matrix_types::vec_of_vec::VecOfVecFromBorrow::from_vecofvecsimple) method.
pub struct UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where   
        KeyMin:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        // OrderComparatorAutoLtByKey<usize, SnzVal, (usize, SnzVal)>: StrictlyLess< (usize, SnzVal)>, // this seems extraneous but the compiler seems to want it
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    array_mapping:                      ArrayMapping,
    array_matching:                     GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,
    array_comb_codomain_inv_matched_block_off_diagonal:   VecOfVecSimple< usize, SnzVal >,
    // array_comb_codomain_inv_matched_block_off_diagonal: VecOfVec< ( usize, SnzVal ), OrderComparatorAutoLtByFirstTupleEntry >,
    // array_comb_codomain_inv_matched_block_off_diagonal: Vec< Vec< ( usize, SnzVal ) > >,
    ring_operator:                      RingOperator,
    order_comparator_minor:             OrderComparatorMinor,
    order_comparator_major:             OrderComparatorMajor,    
    phantom_viewmajorascend:            PhantomData< ViewMajorAscend >, // required b/c otherwise the compiler complains that the type parameter `ViewMajorAscend` is unused
}

     


impl < ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  

    UmatchRowMajor 
    < ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  
    
    where   
        KeyMin:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >,
        ViewMajorAscend:        IntoIterator,
        ViewMajorAscend::Item:  KeyValGet< KeyMin, SnzVal >,
        // OrderComparatorAutoLtByKey<usize, SnzVal, (usize, SnzVal)>: StrictlyLess< (usize, SnzVal)>
{

    /// Generate a new U-match factorization.
    pub fn new
            < IterKeyMaj > ( 
                array_mapping:              ArrayMapping, 
                iter_keymaj:                IterKeyMaj,
                ring_operator:              RingOperator,
                order_comparator_minor:     OrderComparatorMinor,
                order_comparator_major:     OrderComparatorMajor,                
            ) 
        -> 
        UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 

    where   ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >,
            IterKeyMaj:             Iterator < Item = KeyMaj >,
            KeyMin:                 Clone + Hash + std::cmp::Eq + Debug, 
            KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            SnzVal:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
            ViewMajorAscend:        IntoIterator + Clone, // !!! remove clone
            ViewMajorAscend::Item:  KeyValSet< KeyMin, SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderComparatorMinor:   Clone + StrictlyLess <  KeyMin >, // !!! remove clone
            OrderComparatorMajor:   Clone + StrictlyLess<  KeyMaj >,
            HitMerge<Peekable<Scale<<ViewMajorAscend as IntoIterator>::IntoIter, KeyMin, RingOperator, SnzVal>>, OrderComparatorLtByKey< KeyMin, SnzVal, ViewMajorAscend::Item, OrderComparatorMinor>>: Clone // !!!! remove this        
            
    {
        
        let ( array_codomain_comb_inv_off_diag_pivot_block, array_matching ) : ( VecOfVecSimple<usize, SnzVal>, GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal > )
            = get_codomain_comb_inv_off_diag_pivot_block( 
                    & array_mapping, 
                    iter_keymaj, 
                    ring_operator.clone(),
                    order_comparator_minor.clone(),
                    order_comparator_major.clone(),
                );
        
        UmatchRowMajor{ 
                array_mapping, 
                array_matching, 
                array_comb_codomain_inv_matched_block_off_diagonal:     array_codomain_comb_inv_off_diag_pivot_block,   
                ring_operator:                                          ring_operator,
                order_comparator_minor:                                 order_comparator_minor,
                order_comparator_major:                                 order_comparator_major,                
                phantom_viewmajorascend:                                PhantomData,
            }
        
    }
  
// }








// //  =========================================================================================================
// //  U-MATCH REF OBJECT
// //  =========================================================================================================



// #[derive(Debug)]
// pub struct UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
//     where   
//         KeyMin:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ViewMajorAscend:            IntoIterator,
//         ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
// {
//     array_mapping:                                          &'a ArrayMapping,
//     array_matching:                                         &'b GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,
//     array_comb_codomain_inv_matched_block_off_diagonal:     &'b VecOfVecSimple< usize, SnzVal >,
//     ring_operator:                                          RingOperator,
//     order_comparator_minor:                                 OrderComparatorMinor,
//     order_comparator_major:                                 OrderComparatorMajor,    
//     phantom_viewmajorascend:                                PhantomData< ViewMajorAscend >, // required b/c otherwise the compiler complains that the type parameter `ViewMajorAscend` is unused
// }

// //  Implement
// //  ---------------------------------------------------------------------------------------------------------

// impl < 'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  

//     UmatchRowMajorWithRefs
//     < 'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  
    
//     where   
//         KeyMin:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
//         ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >,
//         ViewMajorAscend:        IntoIterator,
//         ViewMajorAscend::Item:  KeyValGet< KeyMin, SnzVal >,


// {
//     pub fn new
//         ( umatch: &'b UmatchRowMajor < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  ) 
//         ->
//         UmatchRowMajorWithRefs
//              <'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
//         where
//             RingOperator:       Clone,
//             OrderComparatorMinor: Clone,
//             OrderComparatorMajor: Clone,
//     {
//         UmatchRowMajorWithRefs{ 
//                 array_mapping:                                            umatch.array_mapping, 
//                 array_matching:                                         & umatch.array_matching,
//                 array_comb_codomain_inv_matched_block_off_diagonal:     & umatch.array_comb_codomain_inv_matched_block_off_diagonal,
//                 ring_operator:                                            umatch.ring_operator.clone(),
//                 order_comparator_minor:                                   umatch.order_comparator_minor.clone(),
//                 order_comparator_major:                                   umatch.order_comparator_major.clone(),                
//                 phantom_viewmajorascend:                                  PhantomData,
//             }
//     }

    /// Returns a copy of the ring operator
    pub fn ring_operator( &self ) -> RingOperator 
        where
            RingOperator:   Clone,
    { self.ring_operator.clone() }

    /// Returns a copy of the order comparator for index-value pairs whose index is a `KeyMin`.
    pub fn order_comparator_minor( &self ) -> OrderComparatorMinor 
        where
            OrderComparatorMinor:   Clone,
    { self.order_comparator_minor.clone() }    

    /// Returns a copy of the order comparator for index-value pairs whose index is a `KeyMin`.
    pub fn order_comparator_major( &self ) -> OrderComparatorMajor 
        where
            OrderComparatorMajor:   Clone,
    { self.order_comparator_major.clone() }        

    /// Returns the (row-major) codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn array_comb_codomain< 'a >( &'a self ) -> UmatchRowMajorCombCodomain< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  {
        UmatchRowMajorCombCodomain{ umatch: self }
    }

    /// Returns the (row-major) inverse of the codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn array_comb_codomain_inv< 'a >( &'a self ) -> UmatchRowMajorCombCodomainInv< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  {
        UmatchRowMajorCombCodomainInv{ umatch: self }
    }  
    
    /// Returns the (row-major) domain COMB, indexed by `KeyMin`
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn array_comb_domain< 'a >( &'a self ) -> UmatchRowMajorCombDomain< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  {
        UmatchRowMajorCombDomain{ umatch: self }
    }

    /// Returns the (row-major) inverse of the domain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn array_comb_domain_inv< 'a >( &'a self ) -> UmatchRowMajorCombDomainInv< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  {
        UmatchRowMajorCombDomainInv{ umatch: self }
    }      

    /// Returns a reference to the matching array of the internally stored  U-match factorization.
    pub fn array_matching_ref( &self ) -> & GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal > { & self.array_matching }

    /// Returns a reference to the mapping array of the internally stored U-match factorization.
    pub fn array_mapping_ref( &self ) -> & ArrayMapping { & self.array_mapping }    

    /// The column submatrix of the mapping array indexed by matched column indices.
    /// 
    /// Returns a wrapper that implements `OracleMajorAscend`.
    /// 
    /// # Design notes
    /// 
    /// It seems necessary to use the lifetime parameter 'b because (i) the matching array in the U-match decomposition
    /// is stored by value within each `UmatchRowMajor` object, (ii) this is the object that contains the hashmap
    /// which we use to store the set of matched indices, however (iii) the object `OnlyKeyMinOutsideCollection`
    /// requires a *reference* to a hashmap, and (iv) said reference has to live long enough to do interesting things.  
    /// The most evident way to get that "long-enough" lifetime parameter, at present, is by passing the `UmatchRowMajor`
    /// object by reference.
    pub fn array_mapping_matched_cols_only< 'a >( &'a self ) -> OnlyKeyMinInsideCollection< OracleRef< 'a, ArrayMapping >, ViewMajorAscend, &'a HashMap< KeyMin, usize >, KeyMin, KeyMaj, SnzVal >
        // where 'a: 'b,
    {
        OnlyKeyMinInsideCollection::<_,_,_,_,_,_>::new( OracleRef::new( &self.array_mapping ), self.array_matching.bimap_min_ref().hashmap_val_to_ord() )
    }    
    // fn array_mapping_matched_cols_only<'c>( &'c self ) -> OnlyKeyMinInsideCollection<'c,'c, ArrayMapping, ViewMajorAscend, HashMap< KeyMin, usize >, KeyMin, KeyMaj, SnzVal >
    //     where 'a: 'c, 'b: 'c,
    // {
    //     OnlyKeyMinInsideCollection::<'c, 'c,_,_,_,_,_,_>::new( self.array_mapping, self.array_matching.bimap_min().hashmap_val_to_ord() )
    // }    

    /// The column submatrix of the mapping array indexed by unmatched column indices.
    /// 
    /// Returns a wrapper that implements `OracleMajorAscend`.
    /// 
    /// # Design notes
    /// 
    /// It seems necessary to use the lifetime parameter 'b because (i) the matching array in the U-match decomposition
    /// is stored by value within each `UmatchRowMajor` object, (ii) this is the object that contains the hashmap
    /// which we use to store the set of matched indices, however (iii) the object `OnlyKeyMinOutsideCollection`
    /// requires a *reference* to a hashmap, and (iv) said reference has to live long enough to do interesting things.  
    /// The most evident way to get that "long-enough" lifetime parameter, at present, is by passing the `UmatchRowMajor`
    /// object by reference.
    pub fn array_mapping_matchless_cols_only< 'a >( &'a self ) -> OnlyKeyMinOutsideCollection< OracleRef< 'a, ArrayMapping >, ViewMajorAscend, &'a HashMap< KeyMin, usize >, KeyMin, KeyMaj, SnzVal >
    {
        OnlyKeyMinOutsideCollection::<_,_,_,_,_,_>::new( OracleRef::new( &self.array_mapping ), & self.array_matching.bimap_min_ref().hashmap_val_to_ord() )
    }    

    /// Returns a reference to the internally stored compressed representation of the inverse of the codomain COMB;
    /// this representation consists of a `VecOfVecSimple` which encodes the off-diagonal entries of the square submatrix of the 
    /// inverse of the codomain comb which is indexed by matched (i.e. pivot) indices.
    /// 
    /// Concretely, the represented matrix is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn array_comb_codomain_inv_matched_block_off_diagonal_ref< 'a >( &'a self ) 
        -> 
        &'a VecOfVecSimple< usize, SnzVal >
        // & VecOfVec< (usize, SnzVal), OrderComparatorAutoLtByFirstTupleEntry > 
        { & self.array_comb_codomain_inv_matched_block_off_diagonal }   

        
    /// Returns a nested double reference to the internally stored compressed representation of the inverse of the codomain COMB;
    /// this representation consists of a `VecOfVecSimple` which encodes the off-diagonal entries of the square submatrix of the 
    /// inverse of the codomain comb which is indexed by matched (i.e. pivot) indices.
    /// 
    /// Concretely, the represented matrix is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).    
    /// 
    /// # Design note
    /// 
    /// This function exists because many parts of the `oat_rust` library use *references* to objects that implement
    /// matrix oracle traits.  A `VecOfVec` simple does not implement oracle traits in general, but a reference
    /// `& VecOfVecSimple` does.  Therefore we often need to work with objects of form `&'b &'b VecOfVecSimple`.  In
    /// practice, we find that Rust is prone to inferring the wrong lifetime if we simply write 
    /// `& self.array_comb_codomain_inv_matched_block_off_diagonal_ref()` (for example, one finds errors alluding to
    /// dropped temprorary values).  This function has succeeded in sidestepping such errors in the past; please 
    /// let us know if it fails to do so successfully in future examples.
    // pub fn array_comb_codomain_inv_matched_block_off_diagonal_ref_ref( &'b self ) 
    //     -> 
    //     &'b &'b VecOfVecSimple< usize, SnzVal >
    //     // & VecOfVec< (usize, SnzVal), OrderComparatorAutoLtByFirstTupleEntry > 
    //     { & self.array_comb_codomain_inv_matched_block_off_diagonal }           
    


    /// The square, block submatrix of the inverse of the codomain COMB which is indexed by matched row indices.
    /// 
    /// Concretely, this object is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn array_comb_codomain_inv_matched_block< 'a >( &'a self ) 
        -> 
        PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend<
                &'a VecOfVecSimple< usize, SnzVal >,
                // &'a VecOfVec< (usize, SnzVal), OrderComparatorAutoLtByFirstTupleEntry >,                
                Cloned< Iter< 'a, (usize, SnzVal) > >,
                // &'a [(usize, SnzVal)],
                usize,
                SnzVal,
            >
        where
            RingOperator:   Semiring< SnzVal >,
        {   

            let prepended : PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< &'a VecOfVecSimple<usize, SnzVal>, Cloned< Iter< 'a, (usize, SnzVal) >>, usize, SnzVal >
            = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( 
                        self.array_comb_codomain_inv_matched_block_off_diagonal_ref(), 
                        RingOperator::one() 
                    );  
            
            return prepended
        }      


    pub fn array_comb_codomain_inv_times_mapping_array_matched_block< 'a >( &'a self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlock< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
        where 
            RingOperator: Clone,
            OrderComparatorMinor: Clone,
            OrderComparatorMajor: Clone,           
    {
        CombCodomainInvTimesMappingMatchedBlock{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }


    pub fn array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_row_ordinals< 'a >( &'a self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
        where 
            RingOperator: Clone,
            OrderComparatorMinor: Clone,
            OrderComparatorMajor: Clone,
    {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }    
}

//  Implement clone
//  ---------------------------------------------------------------------------------------------------------

// impl < 'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >   

//     Clone for

//     UmatchRowMajorWithRefs
//         < 'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
//     where   
//         KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ViewMajorAscend:            IntoIterator,
//         ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
//         RingOperator:               Clone,
//         OrderComparatorMinor:       Clone,
//         OrderComparatorMajor:       Clone,        
// {
//     fn clone( &self ) -> Self {
//         UmatchRowMajorWithRefs{
//             array_mapping:                                          self.array_mapping,
//             array_matching:                                         self.array_matching,
//             array_comb_codomain_inv_matched_block_off_diagonal:     self.array_comb_codomain_inv_matched_block_off_diagonal,
//             ring_operator:                                          self.ring_operator.clone(),
//             order_comparator_minor:                                 self.order_comparator_minor.clone(),
//             order_comparator_major:                                 self.order_comparator_major.clone(),            
//             phantom_viewmajorascend:                                PhantomData, // required b/c otherwise the compiler complains that the type parameter `ViewMajorAscend` is unused
//         }
//     }
// }



    







//  =========================================================================================================
//  COMPRESSED CODOMAIN COMB
//  =========================================================================================================


// ///
// /// 
// /// # Design note
// /// 
// /// It makes sense to talk about a row-major matrix in this setting because rows vs. columns of upper
// /// triangular matrices can be differentiated by their sparcity patterns (in particular, whether the
// /// diagonal element is the first or the last nonzero entry)
// struct RowMajorCombCodCompressed< 
//                 'a, KeyMin, KeyMaj, SnzVal, ArrayMapping, ViewMajorAscend, RingOperator, OrderComparatorMinor,
//             > 
//     where   
//         KeyMin:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct            
// {
//     comb_cod_compressed_ordinal:    &'a VecOfVecSimple< usize, SnzVal >,
//     matching:                       &'a GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,
//     mapping:                        &'a ArrayMapping, // the mapping array
//     ring_opeartor:                  RingOperator,
//     order_comparator_minor:               OrderComparatorMinor,
//     phandom_viewmajorascend:        PhantomData< ViewMajorAscend >,
// }

// impl < 'a, KeyMin, KeyMaj, SnzVal, ArrayMapping, ViewMajorAscend, RingOperator, OrderComparatorMinor, OrderComparatorMajor >

//         OracleMajorAscend< 
//                 KeyMaj, Cloned< Iter< '_, (KeyMin, SnzVal) > > 
//             > for
    
//         RowMajorCombCodCompressed< 
//                 'a, KeyMin, KeyMaj, SnzVal, ArrayMapping, ViewMajorAscend, RingOperator, OrderComparatorMinor,
//             >     
//     where
//         KeyMin:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct            
//         SnzVal:                 Clone,
//         ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >, 
//         RingOperator:           Semiring< SnzVal >,
//         OrderComparatorMinor:        Clone + StrictlyLess<  KeyMaj >,
//         ViewMajorAscend:        IntoIterator,
//         ViewMajorAscend::Item:  KeyValGet< KeyMin, SnzVal >,

// {

//     fn view_major_ascend(&self, keymaj: KeyMaj) -> Cloned< Iter< '_, (KeyMin, SnzVal) > >  { 

//         let iter_diag   =   std::iter::once(  ( keymaj.clone(), )  );
        
//         let vec_masked  =   self
//                                             .mapping
//                                             .view_major_ascend( keymaj )
//                                             .into_iter()
//                                             .filter_map(
//                                                     |x| 
//                                                     self
//                                                         .matching
//                                                         .keymin_to_ordmaj(keymin)
//                                                 )

//         let iter_over_diag  =   hit_merge_by(


//                 );
        
//     }
// }
        














//  =========================================================================================================
//  COMB'S (COLUMNAR ORDERED MATCHING BASES -- UNCOMPRESSED)
//  =========================================================================================================



/// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
/// to construct rows of the codomain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombCodomain< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where   
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > ,
}

/// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
/// to construct rows of the inverse of the codomain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombCodomainInv< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where   
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > ,
}

/// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
/// to construct rows of the domain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombDomain< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where   
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > ,
}

/// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
/// to construct rows of the inverse of the domain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombDomainInv< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where   
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > ,
}




//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: DOMAIN
//  ---------------------------------------------------------------------------------------------------------

//  FOR UNIT TESTS, RUN A TEXT SEARCH FOR <#CombDomainImplementOracleMajorAscend>

impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >

    OracleMajorAscend<  
            KeyMin,  // !! NOTE THAT THE MINOR KEYS OF THE MAPPING ARRAY ARE ALSO MAJOR KEYS OF THE DOMAIN COMB !!
            IterTwoType< 
                Once< ViewMajorAscend::Item >,
                MergeTwoItersByOrderComparator<
                        Peekable< IterWrappedVec< ViewMajorAscend::Item > >, 
                        Peekable< 
                                LinearCombinationSimplified<
                                        OnlyIndicesOutsideCollection< ViewMajorAscend::IntoIter, &'a HashMap<KeyMin, usize>, KeyMin, SnzVal>, 
                                        KeyMin, SnzVal, RingOperator, OrderComparatorMinor
                                    >, 
                            >,
                        OrderComparatorMinor
                    >                                     
            >  
        > for 

    UmatchRowMajorCombDomain< 
            'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor,
        >

    where 
        KeyMin:                 Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        KeyMaj:                 Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        SnzVal:                 Clone,
        RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
        OrderComparatorMinor:   Clone + StrictlyLess<  ViewMajorAscend::Item >,
        OrderComparatorMajor:   Clone + StrictlyLess<  ( KeyMaj, SnzVal ) >,        
        ViewMajorAscend:        IntoIterator,        
        ViewMajorAscend::Item:  Clone + KeyValSet< KeyMin, SnzVal > + KeyValNew< KeyMin, SnzVal >, // KeyValNew is used (i) to create the diagonal elements of the array, and (ii) to create the 'problem vector' b in a linear problem xA = b whose solution is used in our constructor
        ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >
{
    fn view_major_ascend( &self, keymin: KeyMin ) 
        ->  
        IterTwoType< 
                Once< ViewMajorAscend::Item >,
                MergeTwoItersByOrderComparator<
                        Peekable< IterWrappedVec< ViewMajorAscend::Item > >, 
                        Peekable< 
                                LinearCombinationSimplified<
                                        OnlyIndicesOutsideCollection< ViewMajorAscend::IntoIter, &'a HashMap<KeyMin, usize>, KeyMin, SnzVal>, 
                                        KeyMin, SnzVal, RingOperator, OrderComparatorMinor
                                    >, 
                            >,
                        OrderComparatorMinor
                    >                                     
            >        

    {
        match self.umatch.array_matching.contains_keymin( &keymin ) { 
            false => { 

                println!("MADE IT THROUGH END OF UmatchRowMajorCombDomain.view_major_ascend (UNMATCHED INDEX)");                   
                IterTwoType::Iter1(
                        std::iter::once(  ViewMajorAscend::Item::new( keymin, RingOperator::one() )  ) 
                    )
            }
            true => { 
                
                println!("STARTING CONSTRUCTION FOR UmatchRowMajorCombDomain.view_major_ascend (UNMATCHED INDEX)");                   

                // The matrix A from Hang et al., "Umatch factorization ...", with rows indexed by major key ordinals
                // Note: this struct only contains a reference
                let seed = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_row_ordinals();

                // Struct that encodes the matching from minor keys to major key ordinals
                let array_matching_ref = self.umatch.array_matching_ref();
                let keymin_to_ordmaj = | keymin: KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
                let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );

                println!("WAYPOINT 1");  

                // obtain a row of A^{-1}
                let seed_inv_row = 
                    EchelonSolveMajorAscendWithMinorKeys::new(
                            std::iter::once( ViewMajorAscend::Item::new( keymin.clone(), RingOperator::one() ) ), // the standard unit vector supported on `keymin`
                            seed, // matrix A
                            keymin_to_ordmaj_wrapped,
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_comparator_minor.clone(),
                        );

                println!("WAYPOINT 2");  

                let mut seed_inv_row_vec = seed_inv_row.collect_vec();

                println!("WAYPOINT 3"); 

                seed_inv_row_vec.shrink_to_fit();               

                // obtain a copy of R_{\rho \rho}
                let comb_codomain_inv_matched_block = self.umatch.array_comb_codomain_inv_matched_block();

                // obtain a copy of the column submatrix of the mapping array indexed by unmatched column indices
                let array_mapping_npcols = self.umatch.array_mapping_matchless_cols_only();

                // y = - seed_inv_row_vec * R_{\rho \rho}^{-1} * D_{ \rho \bar \kappa}
                let x 
                    =   vector_matrix_multiply_major_ascend_simplified(
                                seed_inv_row_vec
                                    .iter()
                                    .map(   |x| 
                                            (   
                                                self.umatch.array_matching_ref().keymin_to_ordmaj(& x.key() ).unwrap(),  // replace integer indices (major ordinals) with major key indices
                                                self.umatch.ring_operator.negate ( x.val() ) // multiply by -1 (recally that we are computing *MINUS* seed_inv_row_vec * R_{\rho \rho}^{-1} * D_{ \rho \bar \kappa})
                                            ) 
                                        ),
                                comb_codomain_inv_matched_block,
                                self.umatch.ring_operator.clone(),
                                OrderComparatorAutoLtByKey::new(),
                            );

                let y
                    =   vector_matrix_multiply_major_ascend_simplified(
                                x.map( 
                                        | ( ordmaj, snzval ) |
                                        ( 
                                                array_matching_ref.ordmaj_to_keymaj( ordmaj ),
                                                snzval
                                            )
                                    ),
                                array_mapping_npcols,
                                self.umatch.ring_operator.clone(),
                                self.umatch.order_comparator_minor.clone(),                                
                            );

                // rescale entries of seed_inv_row_vec, transforming it into row `keymin` of `A^{-1} M_{\rho \kappa}`
                // NB: this has to occur AFTER we've used `seed_inv_row_vec` to construct `y`
                let mut seed_inv_row_vec_times_matching = seed_inv_row_vec;
                for entry in seed_inv_row_vec_times_matching.iter_mut() {
                    entry.set_val( 
                            self.umatch.ring_operator.multiply( 
                                    entry.val(),  
                                    self.umatch.array_matching_ref().keymin_to_snzval( &entry.key() )
                                )                         
                        )
                }

                // merge seed_inv_row_vec
                let z = MergeTwoItersByOrderComparator::new( 
                                                            IterWrappedVec::new( seed_inv_row_vec_times_matching ).peekable(), 
                                                            y.peekable(), 
                                                            self.umatch.order_comparator_minor.clone()
                                                        );


                println!("MADE IT THROUGH END OF UmatchRowMajorCombDomain.view_major_ascend (MATCHED INDEX)");                
                return  IterTwoType::Iter2( z );                

            }
        }
    }
}



//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: DOMAIN INVERSE
//  ---------------------------------------------------------------------------------------------------------

//  FOR UNIT TESTS, RUN A TEXT SEARCH FOR <#CombDomainInvImplementOracleMajorAscend>

impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >

    OracleMajorAscend<  
            KeyMin,  // !! NOTE THAT THE MINOR KEYS OF THE MAPPING ARRAY ARE ALSO MAJOR KEYS OF THE DOMAIN COMB !!
            IterTwoType< 
                    Once
                        < ViewMajorAscend::Item >,
                    LinearCombinationSimplified
                        < ViewMajorAscend::IntoIter, KeyMin, SnzVal, RingOperator, OrderComparatorMinor >,
                    // ViewMajorAscend::Item,                        
                >,
        > for 

    UmatchRowMajorCombDomainInv< 
            'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor,
        >

    where 
        KeyMin:                 Clone + Hash + std::cmp::Eq,
        KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        SnzVal:                 Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
        OrderComparatorMinor:        Clone + StrictlyLess<  ViewMajorAscend::Item >,
        ViewMajorAscend:        IntoIterator,        
        ViewMajorAscend::Item:  KeyValSet< KeyMin, SnzVal > + KeyValNew< KeyMin, SnzVal >,
        ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >,
{
    fn view_major_ascend( &self, keymin: KeyMin )  // this function signature looks strange b/c the major keys of the domain comb have type KeyMin!  
        ->  
        IterTwoType< 
            Once
                < ViewMajorAscend::Item >,
            LinearCombinationSimplified 
                < ViewMajorAscend::IntoIter, KeyMin, SnzVal, RingOperator, OrderComparatorMinor >,           
            // ViewMajorAscend::Item,                        
        >        

    {
        match self.umatch.array_matching.keymin_to_ordmaj( &keymin ) { // this function call looks strange b/c the major keys of the domain comb have type KeyMin!  
            None => { 
                IterTwoType::Iter1(
                        std::iter::once(  ViewMajorAscend::Item::new( keymin, RingOperator::one() )  ) 
                    )
            }
            Some( ordmaj ) => {  
                let scalar  =   self.umatch.array_matching.ordmaj_to_snzval( ordmaj );
                let scalar_inv      =   self.umatch.ring_operator.invert( scalar );
                // let key_to_codomain_comb    =   self.umatch.array_matching.ordmaj_to_keymaj(ordmaj);

                // this equals row `keymaj` of M_{rk}^{-1} * R^{-1}_{\rho \rho}
                let lefthand_factor_vec     =   
                        (& self.umatch.array_comb_codomain_inv_matched_block() )                                                    
                            .view_major_ascend( ordmaj )
                            .scale( scalar_inv, self.umatch.ring_operator.clone() )
                            .map(   |(x,y)|   // re-index the entries, so that indices align with the row indices of the mapping array
                                    (self.umatch.array_matching.ordmaj_to_keymaj(x), y) 
                                );

                println!("LEFTHAND FACTOR VEC: {:?}", lefthand_factor_vec.clone().collect_vec() );

                let iter = vector_matrix_multiply_major_ascend_simplified ( 
                            lefthand_factor_vec,
                            OracleRef::new( &self.umatch.array_mapping ),
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_comparator_minor.clone(),
                        );

                IterTwoType::Iter2(
                    iter
                )
            }
        }
    }
}



//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: CODOMAIN
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >

    OracleMajorAscend<  
            KeyMaj,  
            IterTwoType< 
                    Chain<
                            Once
                                < (KeyMaj, SnzVal) >,
                            IterWrappedVec
                                < (KeyMaj, SnzVal) >,
                        >,
                    IterWrappedVec
                        < (KeyMaj, SnzVal) >,
                >,
        > for 

    UmatchRowMajorCombCodomain< 
            'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor,
        >

    where 
        KeyMin:                 Clone + Hash + std::cmp::Eq,
        KeyMaj:                 Clone + Hash + std::cmp::Eq,
        SnzVal:                 Clone,
        RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
        OrderComparatorMinor:   Clone + StrictlyLess<  ViewMajorAscend::Item >,
        OrderComparatorMajor:   Clone + StrictlyLess<  (KeyMaj, SnzVal) >,
        ViewMajorAscend:        IntoIterator,        
        ViewMajorAscend::Item:  KeyValSet< KeyMin, SnzVal > + KeyValNew< KeyMin, SnzVal >,
        ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >
{
    /// The developers are unsure how to design a lazy ascending sparse vector iterator for a solution to `xA = b`, where `A = R_{\rho \rho}^{-1} * D_{\rho \kappa}` is the 
    /// matrix defined vis-a-vis inner identities in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf).  In particular, while there are
    /// several direct means to compute a solution `x`, it is unclear how to generate the entries of `x` in ascending order of index, in a lazy 
    /// fashion.
    /// The reason for this difficulty can be traced to the matrix `A = R_{\rho \rho}^{-1} * D_{\rho \kappa}`:
    /// If we arrange the columns of `A` to match the order of `\rho` (concretely, this is the same as taking the matrix `R_{\rho \rho}^{-1} * D_{\rho \kappa^*}`)
    /// then the resulting matrix is lower triangular.  On the other hand, if we arrange entries accoring the the ordering of
    /// `\kappa` (that is, if we take the matrix `R_{\rho^* \rho} * D_{\rho \kappa}`) then the resulting matrix is upper triangular, its
    /// entries follow the order of `\kappa` rather than `\rho`.  For the time being, 
    /// we simply generate every entry of `x` (in ascending order according to `\kappa`), store these entries in a `Vec`, reindex according to 
    /// `\rho`, sort the resulting sequence by index, and return the sorted vector, wrapped in struct that allows one to iterate. 
    fn view_major_ascend
        ( 
            &self, 
            keymaj: KeyMaj 
        ) 
        -> 
        IterTwoType< 
                Chain<
                        Once
                            < (KeyMaj, SnzVal) >,
                        IterWrappedVec
                            < (KeyMaj, SnzVal) >,
                    >,
                IterWrappedVec
                    < (KeyMaj, SnzVal) >,
            >
    {

        // define the matrix A that will fit into an equation xA = b
        let seed_with_integer_indexed_rows = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_row_ordinals();

        // define the problem vector b that will fit into an equation xA = b
        let array_mapping_matched_cols_only = self.umatch.array_mapping_matched_cols_only();
        let problem_vector = array_mapping_matched_cols_only.view_major_ascend( keymaj.clone() );

        // Struct that encodes the matching from minor keys of A to major key ordinals of A
        let array_matching_ref = self.umatch.array_matching_ref();
        let keymin_to_ordmaj = | keymin: KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
        let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

        // Solve xA = b for x
        let mut solution_vec = 
        EchelonSolveMajorAscendWithMinorKeys::new( // solution to xA = b
                problem_vector, // mabrix b
                seed_with_integer_indexed_rows, // matrix A
                keymin_to_ordmaj_wrapped,
                self.umatch.ring_operator.clone(),
                self.umatch.order_comparator_minor.clone(),
            )
            .collect_vec();

        // Sort solution vector x
        println!("REDO THIS SECTION OF CODE AFTER REFACTORING <CompareOrder>");
        let mut order_comparator_minor_clone = self.umatch.order_comparator_minor.clone(); // we have to clone the order comparator in order to compare order, since the `strictly_less` method takes `& mut self` as an argument -- but to get a mutable reference in the present context we would have to borrow the umatch struct as mutable
        let compare_order_for_vector_sort 
            = | x: &ViewMajorAscend::Item, y: &ViewMajorAscend::Item |-> Ordering { // we have "repackage" order_comparator_major so that it returns an std::cmp::Ordering rather than a bool
            match order_comparator_minor_clone.strictly_less( x, y ) {
                true => { return Ordering::Less },
                false => { 
                    match order_comparator_minor_clone.strictly_less( y, x ) {
                        true => { return Ordering::Greater },
                        false => { return Ordering:: Equal }
                    }
                }
            }
        };
        solution_vec.sort_by( compare_order_for_vector_sort );

        // Reindex solution vector x
        let mut solution_vec_reindexed = 
            solution_vec.into_iter()
                .map(   | item |
                        ( 
                            self.umatch.array_matching.keymin_to_keymaj( &item.key() ).unwrap(),
                            item.val(),
                        )
                    )
                .collect_vec();

        // a function that sends (a,b) to `Less` if a < b and to `Greater` otherwise
        let mut order_comparator_major = self.umatch.order_comparator_major.clone(); // we have to make a clone because, as currently written, `self` is behind a mutable reference and `order_comparator_major` may mutate itself when it compares two objects
        solution_vec_reindexed.sort_by( 
                |a,b| {  
                    if order_comparator_major.strictly_less( a, b ) { Ordering::Less }
                    else { Ordering::Greater }
                    }
            ); 
        let solution_vec_reindexed_and_sorted = IterWrappedVec::new( solution_vec_reindexed );                            

        // Possibly append an entry with coefficient 1
        match self.umatch.array_matching_ref().contains_keymaj( & keymaj ){
            false =>    { // if keymaj is UNMATCHED then merge in a copy of ( keymaj, 1 )
                return  IterTwoType::Iter1(
                                std::iter::once( (keymaj, RingOperator::one() ) )
                                    .chain( solution_vec_reindexed_and_sorted ) 
                            )

            }
            true =>     { // otherwise just return the reindexed, sorted solution vector
                return IterTwoType::Iter2(
                        solution_vec_reindexed_and_sorted
                    )
            }
        }
    }   
}




//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------



//  IMPLEMENT ORACLE MAJOR ASCEND FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------

impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >

    OracleMajorAscend<  
            KeyMaj,  
            IterTwoType< 
                    MergeTwoItersByOrderComparator<
                            Peekable<
                                    ChangeIndexSimple<
                                            LinearCombinationSimplified
                                                <Chain<Once<(usize, SnzVal)>, Cloned<Iter<'a, (usize, SnzVal)>>>, usize, SnzVal, RingOperator, OrderComparatorAutoLtByKey<usize, SnzVal, (usize, SnzVal)>>, 
                                            &'a Vec<KeyMaj>, 
                                            usize, 
                                            KeyMaj, 
                                            SnzVal
                                        >
                                >, 
                            OncePeekable<(KeyMaj, SnzVal)>, 
                            OrderComparatorMajor
                        >,                        
                    ChangeIndexSimple<
                            Chain<
                                    Once<(usize, SnzVal)>, 
                                    Cloned<Iter<'a, (usize, SnzVal)>>
                                >, 
                            &'a Vec<KeyMaj>, usize, KeyMaj, SnzVal
                        >
                >
        > for 

    UmatchRowMajorCombCodomainInv< 
            'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor,
        >

    where 
        KeyMin:                 Clone + Hash + std::cmp::Eq,
        KeyMaj:                 Clone + Hash + std::cmp::Eq,
        SnzVal:                 Clone,
        RingOperator:           Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
        OrderComparatorMinor:   Clone + StrictlyLess<  ViewMajorAscend::Item >,
        OrderComparatorMajor:   Clone + StrictlyLess<  ( KeyMaj, SnzVal ) >,
        ViewMajorAscend:        IntoIterator,        
        ViewMajorAscend::Item:  KeyValSet< KeyMin, SnzVal > + KeyValNew< KeyMin, SnzVal >,
        ArrayMapping:           OracleMajorAscend< KeyMaj, ViewMajorAscend >
{
    /// The developers are unsure how to design a lazy ascending sparse vector iterator for a solution to `xA = b`, where `A = R_{\rho \rho}^{-1} * D_{\rho \kappa}` is the 
    /// matrix defined vis-a-vis inner identities in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf).  In particular, while there are
    /// several direct means to compute a solution `x`, it is unclear how to generate the entries of `x` in ascending order of index, in a lazy 
    /// fashion.
    /// The reason for this difficulty can be traced to the matrix `A = R_{\rho \rho}^{-1} * D_{\rho \kappa}`:
    /// If we arrange the columns of `A` to match the order of `\rho` (concretely, this is the same as taking the matrix `R_{\rho \rho}^{-1} * D_{\rho \kappa^*}`)
    /// then the resulting matrix is lower triangular.  On the other hand, if we arrange entries accoring the the ordering of
    /// `\kappa` (that is, if we take the matrix `R_{\rho^* \rho} * D_{\rho \kappa}`) then the resulting matrix is upper triangular, its
    /// entries follow the order of `\kappa` rather than `\rho`.  For the time being, 
    /// we simply generate every entry of `x` (in ascending order according to `\kappa`), store these entries in a `Vec`, reindex according to 
    /// `\rho`, sort the resulting sequence by index, and return the sorted vector, wrapped in struct that allows one to iterate. 
    fn view_major_ascend
        ( 
            &self, 
            keymaj: KeyMaj 
        ) 
        -> 
            IterTwoType< 
                    MergeTwoItersByOrderComparator<
                            Peekable<
                                    ChangeIndexSimple<
                                            LinearCombinationSimplified
                                                <Chain<Once<(usize, SnzVal)>, Cloned<Iter<'a, (usize, SnzVal)>>>, usize, SnzVal, RingOperator, OrderComparatorAutoLtByKey<usize, SnzVal, (usize, SnzVal)>>, 
                                            &'a Vec<KeyMaj>, 
                                            usize, 
                                            KeyMaj, 
                                            SnzVal
                                        >
                                >, 
                            OncePeekable<(KeyMaj, SnzVal)>, 
                            OrderComparatorMajor
                        >,                        
                    ChangeIndexSimple<
                            Chain<
                                    Once<(usize, SnzVal)>, 
                                    Cloned<Iter<'a, (usize, SnzVal)>>
                                >, 
                            &'a Vec<KeyMaj>, usize, KeyMaj, SnzVal
                        >,
                >
    {

        match self.umatch.array_matching_ref().keymaj_to_ordmaj( & keymaj ) {
            // unmatched row index
            None => {
                // define the matrix A that will fit into an equation xA = b
                let seed_with_integer_indexed_rows = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_row_ordinals();
                // let seed_with_integer_indexed_rows_ref = & seed_with_integer_indexed_rows;

                // define the problem vector b that will fit into an equation xA = b
                let array_mapping_matched_cols_only = self.umatch.array_mapping_matched_cols_only();
                let problem_vector 
                    = array_mapping_matched_cols_only
                        .view_major_ascend( keymaj.clone() )
                        .negate( self.umatch.ring_operator().clone() );  // recall that there is a MINUS SIGN in the formula in Theorem 6 of Hang et al., "Umatch factorization: ..."

                // Struct that encodes the matching from minor keys of A to major key ordinals of A
                let array_matching_ref = self.umatch.array_matching_ref();
                let keymin_to_ordmaj = | keymin: KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
                let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

                // Solve xA = b for x
                let mut solution_vec = 
                EchelonSolveMajorAscendWithMinorKeys::new( // solution to xA = b
                        problem_vector, // mabrix b
                        seed_with_integer_indexed_rows, // matrix A
                        keymin_to_ordmaj_wrapped,
                        self.umatch.ring_operator.clone(),
                        self.umatch.order_comparator_minor.clone(),
                    );
                let solution_vec_integer_indexed 
                    = ChangeIndexSimple::new( solution_vec, self.umatch.array_matching_ref().bimap_min_ref().hashmap_val_to_ord() );
                    
                // Multiply the solution x by matrix R_{\rho \rho}^{-1}
                let comb_codomain_inv_matched_block = self.umatch.array_comb_codomain_inv_matched_block();

                // Compute the portion of desired matrix view indexed by matched major keys -- however, this vector will actually have usize indices which we must change next
                let vec_with_matched_indices 
                    = vector_matrix_multiply_major_ascend_simplified( 
                            solution_vec_integer_indexed, 
                            comb_codomain_inv_matched_block, 
                            self.umatch.ring_operator(),
                            OrderComparatorAutoLtByKey::new(),
                        );

                // reindex the preceding vector so that it is indeed indexed by major keys
                let vec_with_matched_indices_reindexed 
                    = ChangeIndexSimple::new( 
                            vec_with_matched_indices, 
                            self.umatch.array_matching_ref().bimap_maj_ref().vec_ord_to_val() 
                        );
                
                // Compute the portion of desired matrix view indexed by unmatched major keys -- this portion of the vector has exactly one nonzero entry
                let vec_with_unmatched_indices = OncePeekable::new(  (keymaj, RingOperator::one() )  );

                // Merge the two portions of the vector together, to form a whole
                let merged  
                    =   MergeTwoItersByOrderComparator::new(
                                vec_with_matched_indices_reindexed.peekable(),
                                vec_with_unmatched_indices,
                                self.umatch.order_comparator_major.clone(),
                            );
                
                return IterTwoType::Iter1( merged )
            }
            // matched row index
            Some( ordmaj ) => {
                let reindexed_iter = 
                    ChangeIndexSimple::new(  
                            self.umatch.array_comb_codomain_inv_matched_block().view_major_ascend( ordmaj ),
                            self.umatch.array_matching_ref().bimap_maj_ref().vec_ord_to_val(),
                        );
                return  IterTwoType::Iter2( reindexed_iter )
            }
        }

        // ON 2022/04/29 THIS CODE SEEMS TO BE UNCESSSARY/DEPRECATED; CONSIDER DELETING
        // // define the matrix A that will fit into an equation xA = b
        // let seed_with_integer_indexed_rows = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_row_ordinals();

        // // define the problem vector b that will fit into an equation xA = b
        // let array_mapping_matched_cols_only = self.umatch.array_mapping_matched_cols_only();
        // let problem_vector = array_mapping_matched_cols_only.view_major_ascend( keymaj.clone() );

        // // Struct that encodes the matching from minor keys of A to major key ordinals of A
        // let array_matching_ref = self.umatch.array_matching_ref();
        // let keymin_to_ordmaj = | keymin: KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
        // let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

        // // Solve xA = b for x
        // let mut solution_vec = 
        // EchelonSolveMajorAscendWithMinorKeys::new( // solution to xA = b
        //         problem_vector, // mabrix b
        //         && seed_with_integer_indexed_rows, // matrix A
        //         keymin_to_ordmaj_wrapped,
        //         self.umatch.ring_operator.clone(),
        //         self.umatch.order_comparator_minor.clone(),
        //     )
        //     .collect_vec();

        // // Sort solution vector x
        // let mut order_comparator_minor_clone = self.umatch.order_comparator_minor.clone(); // we have to clone the order comparator in order to compare order, since the `strictly_less` method takes `& mut self` as an argument -- but to get a mutable reference in the present context we would have to borrow the umatch struct as mutable
        // let compare_order_for_vector_sort = | x: &ViewMajorAscend::Item, y: &ViewMajorAscend::Item |-> Ordering {
        //     match order_comparator_minor_clone.strictly_less( x, y ) {
        //         true => { return Ordering::Less },
        //         false => { return Ordering::Greater }
        //     }
        // };
        // solution_vec.sort_by( compare_order_for_vector_sort );

        // // Reindex solution vector x
        // let mut solution_vec_reindexed = 
        //     solution_vec.into_iter()
        //         .map(   | item |
        //                 ( 
        //                     self.umatch.array_matching.keymin_to_keymaj( &item.key() ).unwrap(),
        //                     item.val(),
        //                 )
        //             )
        //         .collect_vec();
        // let solution_vec_reindexed_and_sorted = IterWrappedVec::new( solution_vec_reindexed );                            

        // // Possibly append an entry with coefficient 1
        // match self.umatch.array_matching_ref().contains_keymaj( & keymaj ){
        //     false =>    { // if keymaj is UNMATCHED then merge in a copy of ( keymaj, 1 )
        //         return  IterTwoType::Iter1(
        //                         std::iter::once( (keymaj, RingOperator::one() ) )
        //                             .chain( solution_vec_reindexed_and_sorted ) 
        //                     )

        //     }
        //     true =>     { // otherwise just return the reindexed, sorted solution vector
        //         return IterTwoType::Iter2(
        //                 solution_vec_reindexed_and_sorted
        //             )
        //     }
        // }
    }   
}



















//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY)
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain comb) * (the pivot block of the matching array).
#[derive(Copy, Clone)]
pub struct CombCodomainInvTimesMappingMatchedBlock< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where     
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    umatch_ref:     & 'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > ,  
}

impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    
    CombCodomainInvTimesMappingMatchedBlock
        < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where   
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
        RingOperator:               Clone,
        OrderComparatorMinor:       Clone,
        OrderComparatorMajor:       Clone,        
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlock`].
    pub fn new( umatch_ref: &'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlock{ umatch_ref }
    }
}


//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT OracleMajorAscend FOR < KEYMAJ, KEYMIN >
//  --------------------------------------------------------------------------------------------------------------

impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 

    OracleMajorAscend< 
            KeyMaj, 
            // LinearCombinationSimplified< ViewMajorAscend::IntoIter, KeyMin, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
            LinearCombinationSimplified< 
                    OnlyIndicesInsideCollection< ViewMajorAscend::IntoIter, &'a HashMap<KeyMin, usize>, KeyMin, SnzVal>, 
                    KeyMin, SnzVal, RingOperator, OrderComparatorMinor 
                >  
        > for

    &'a CombCodomainInvTimesMappingMatchedBlock
        < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 

    where   
        ArrayMapping:               OracleMajorAscend< KeyMaj, ViewMajorAscend >,
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        SnzVal:                     Clone,
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal > + KeyValSet< KeyMin, SnzVal >, // KeyValSet is required in order to construct a `Simplify` struct
        RingOperator:               Clone + Semiring< SnzVal >,
        OrderComparatorMinor:            Clone + StrictlyLess<  ViewMajorAscend::Item >,
        // &'a VecOfVecSimple<usize, SnzVal>: OracleMajorAscend<usize, Cloned<Iter<'a, (usize, SnzVal)>> >,
        // &'a VecOfVecSimple<usize, SnzVal>: OracleMajorAscend<usize, Cloned<std::slice::Iter<'a, (KeyMin, SnzVal)>>>,
        // OrderComparatorMinor: 'a, RingOperator: 'a, SnzVal: 'a, KeyMaj: 'a, KeyMin: 'a, ViewMajorAscend: 'a,            

{
    fn view_major_ascend( &self, keymaj: KeyMaj ) 
        -> 
        LinearCombinationSimplified< 
                OnlyIndicesInsideCollection< ViewMajorAscend::IntoIter, &'a HashMap<KeyMin, usize>, KeyMin, SnzVal>, 
                KeyMin, 
                SnzVal, 
                RingOperator, 
                OrderComparatorMinor 
            >  
        {
        
        // define a row vector
        // NB: we have to translate the index `keymaj` which has type `KeyMaj` into the index `ordmaj` which has type `usize`, because `self.umatch.array_comb_codomain_inv_matched_block_off_diagonal` is indexed by unsigned integers 
        let ordmaj  =   self.umatch_ref.array_matching.keymaj_to_ordmaj( &keymaj ).unwrap();

        let combining_coefficients 
            = self.umatch_ref.array_comb_codomain_inv_matched_block().view_major_ascend( ordmaj );     

        // the matched columns of the mapping array
        let matched_cols_of_mapping_array : OnlyKeyMinInsideCollection< OracleRef< 'a, ArrayMapping >, ViewMajorAscend, &'a HashMap<KeyMin, usize>, KeyMin, KeyMaj, SnzVal>
            =   self.umatch_ref.array_mapping_matched_cols_only();
        // let matched_cols_of_mapping_array // : OnlyKeyMinInsideCollection<'a, 'b, ArrayMapping, ViewMajorAscend, HashMap<KeyMin, usize>, KeyMin, KeyMaj, SnzVal>
        //     =   self.umatch.array_mapping_ref();            

        // the terms whose sum equals the product of the vector with the matrix
        let iter_over_scaled_views = 
            combining_coefficients
                    // .iter()
                    .map(   
                            |( ordmaj, snzval)|  
                            matched_cols_of_mapping_array.view_major_ascend( 
                                    self.umatch_ref.array_matching.ordmaj_to_keymaj( ordmaj ) 
                                )
                                .into_iter()
                                .scale( snzval, self.umatch_ref.ring_operator.clone() )                                
                        );                 
                              

        println!("CombCodomainInvTimesMappingMatchedBlock: WAYPOINT 1");                                

        // sum the terms
        let product_vector  
            =   hit_merge_by_predicate( iter_over_scaled_views, self.umatch_ref.order_comparator_minor.clone() )
                    .simplify( self.umatch_ref.ring_operator.clone() );

        println!("CombCodomainInvTimesMappingMatchedBlock: WAYPOINT 2");                                

        // wrap in a LinearCombinationSimplified struct
        return LinearCombinationSimplified{ linear_combination_simplified: product_vector }

    }
}





//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- ROWS INDEXED BY ORDMAJ
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain comb) * (the pivot block of the matching array).
/// 
/// This struct is almost identical to [`CombCodomainInvTimesMappingMatchedBlock`].  The only difference is 
/// that the corresponding matrix oracle has rows indexed by the integer ordinals of the matched major keys, rather than 
/// by the matched major keys themselves.
#[derive(Copy, Clone)]
pub struct CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where     
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    umatch_ref:    &'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > ,  
}

impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    
    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    where   
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
        RingOperator:               Clone,
        OrderComparatorMinor:       Clone,
        OrderComparatorMajor:       Clone,        
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlock`].
    fn new( umatch_ref: &'a UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj{ umatch_ref }
    }
}


//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT ORACLE FOR < usize, KEYMIN >
//  --------------------------------------------------------------------------------------------------------------

impl < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 

    OracleMajorAscend< 
            usize, 
            // LinearCombinationSimplified< ViewMajorAscend::IntoIter, KeyMin, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
            LinearCombinationSimplified< 
                    OnlyIndicesInsideCollection< ViewMajorAscend::IntoIter, &'a HashMap<KeyMin, usize>, KeyMin, SnzVal>, 
                    KeyMin, SnzVal, RingOperator, OrderComparatorMinor 
                >  
        > for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 

    where   
        ArrayMapping:               OracleMajorAscend< KeyMaj, ViewMajorAscend >,
        KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        SnzVal:                     Clone,
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal > + KeyValSet< KeyMin, SnzVal >, // KeyValSet is required in order to construct a `Simplify` struct
        RingOperator:               Clone + Semiring< SnzVal >,
        OrderComparatorMinor:            Clone + StrictlyLess<  ViewMajorAscend::Item >,
        // &'a VecOfVecSimple<usize, SnzVal>: OracleMajorAscend<usize, Cloned<Iter<'a, (usize, SnzVal)>> >,
        // &'a VecOfVecSimple<usize, SnzVal>: OracleMajorAscend<usize, Cloned<std::slice::Iter<'a, (KeyMin, SnzVal)>>>,
        // OrderComparatorMinor: 'a, RingOperator: 'a, SnzVal: 'a, KeyMaj: 'a, KeyMin: 'a, ViewMajorAscend: 'a,            

{
    fn view_major_ascend( &self, ordmaj: usize ) 
        -> 
        LinearCombinationSimplified< 
                OnlyIndicesInsideCollection< ViewMajorAscend::IntoIter, &'a HashMap<KeyMin, usize>, KeyMin, SnzVal>, 
                KeyMin, 
                SnzVal, 
                RingOperator, 
                OrderComparatorMinor 
            >  
        {
        
        // define a row vector
        let combining_coefficients 
            = self.umatch_ref.array_comb_codomain_inv_matched_block().view_major_ascend( ordmaj );     

        // the matched columns of the mapping array
        let matched_cols_of_mapping_array : OnlyKeyMinInsideCollection< OracleRef< 'a, ArrayMapping >, ViewMajorAscend, &'a HashMap<KeyMin, usize>, KeyMin, KeyMaj, SnzVal>
            =   self.umatch_ref.array_mapping_matched_cols_only();
        // let matched_cols_of_mapping_array // : OnlyKeyMinInsideCollection<'a, 'b, ArrayMapping, ViewMajorAscend, HashMap<KeyMin, usize>, KeyMin, KeyMaj, SnzVal>
        //     =   self.umatch.array_mapping_ref();            

        // the terms whose sum equals the product of the vector with the matrix
        let iter_over_scaled_views = 
            combining_coefficients
                    // .iter()
                    .map(   
                            |( ordmaj, snzval)|  
                            matched_cols_of_mapping_array.view_major_ascend( 
                                    self.umatch_ref.array_matching.ordmaj_to_keymaj( ordmaj ) 
                                )
                                .into_iter()
                                .scale( snzval, self.umatch_ref.ring_operator.clone() )                                
                        );

        // sum the terms
        let product_vector  
            =   hit_merge_by_predicate( iter_over_scaled_views, self.umatch_ref.order_comparator_minor.clone() )
                    .simplify( self.umatch_ref.ring_operator.clone() );

        // wrap in a LinearCombinationSimplified struct
        return LinearCombinationSimplified{ linear_combination_simplified: product_vector }

    }
}




//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- INDEXED BY MINOR KEY ORDINAL INTEGERS
//  =========================================================================================================


//  !!! UNDER CONSTRUCTION
//  !!! AN IMPORTANT CHALLENGE IS THAT MATCHING MATRICES ARE NOT CURRENTLY DESIGNED TO STORE MINOR KEY ORDINALS



/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain comb) * (the pivot block of the matching array).
/// 
/// This marix is indexed by **integers** (its major and minor keys are `usize`).  Specifically, row `i`
/// corresponds to `\rho_i` and column `j` corresponds to `\kappa_j`, in the notation of "U-match factorization, ...".
// pub struct CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal< 'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
//     where     
//         KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ViewMajorAscend:            IntoIterator,
//         ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
// {
//     umatch_ref:     UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > ,  
// }

// impl < 'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
    
//     CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal
//         < 'a, 'b, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
//     where   
//         KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ViewMajorAscend:            IntoIterator,
//         ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
//         RingOperator:               Clone,
//         OrderComparatorMinor:            Clone,
// {
//     // Make a new [`CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal`].
//     fn new( umatch_ref: & UmatchRowMajor< ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >  ) -> Self {
//         CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal{ umatch_ref: (*umatch_ref).clone() }
//     }
// }

    


//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- INDEXED BY MAJOR KEY ORDINAL INTEGERS
//  =========================================================================================================


//  !!! UNDER CONSTRUCTION
//
//  !!! ISSUE: THE OBJECT array_mapping DOES NOT RETURN ENTRIES IN ASCENDING ORDER OF MAJOR KEY ORDINAL; WOULD HAVE TO GENERATE ALL ENTRIES AND RETURN AS A VECTOR






















//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
}    


//  ---------------------------------------------------------------------
//  Unit tests
//  ---------------------------------------------------------------------

#[cfg(test)]
mod unit_tests {
    use std::{iter::Cloned, array};

    use itertools::{Product, assert_equal, Itertools};

    use crate::matrices::{operations::{multiply::vector_matrix_multiply_major_ascend_simplified, display::print_indexed_major_views, umatch::row_major_only::{UmatchRowMajor}}, matrix_types::oracle_ref::OracleRef};
    use crate::matrices::random_constructors::random_vec_of_vec_simple;
    use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
    use crate::utilities::partial_order::OrderComparatorAutoLtByKey; 
    use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    use crate::rings::operator_traits::MinusOne;
    use crate::matrices::random_constructors::random_m_by_n_matrix;



    //  ===================================================================================
    //  CONSTRUCTION OF PIVOT BLOCK OF INVERSE OF THE CODOMAIN COMB
    //  ===================================================================================    

    

    /// This test targets the initial computation of the pivot block of the inverse of the codomain COMB.
    #[cfg(test)]
    fn test_initial_decomposition() {
        use crate::matrices::operations::umatch::row_major_only::get_codomain_comb_inv_off_diag_pivot_block;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::partial_order::OrderComparatorAutoAnyType;

        let array_mapping       =   VecOfVecSimple::new(
                                            vec![
                                                vec![ (0, 1), (1, 1) ],
                                                vec![ (0, 1), (1, 1) ],
                                                vec![ (0, 1), (1, 1) ],                                                
                                            ]
                                        );
        let iter_keymaj         = (0..2).rev();
        let ring_operator       =   PrimeOrderFieldOperator::new(5);
        
        let codomain_comb_inv_pivot_block = get_codomain_comb_inv_off_diag_pivot_block(
                                            & (& array_mapping),
                                            iter_keymaj,
                                            ring_operator,
                                            OrderComparatorAutoAnyType,
                                            OrderComparatorAutoAnyType,
                                        );       
                                        
        println!("{:?}", & codomain_comb_inv_pivot_block.0);
        println!("{:?}", & codomain_comb_inv_pivot_block.1);        
    }


    /// This test targets the initial computation of the pivot block of the inverse of the codomain COMB.    
    #[test]
    fn test_initial_decomposition_another_example() {
        use itertools::Itertools;

        use crate::matrices::operations::umatch::row_major_only::get_codomain_comb_inv_off_diag_pivot_block;
        use crate::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified;        
        use crate::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
        use crate::matrices::operations::transform_vector_wise::VecWiseTransformed;
        use crate::matrices::random_constructors::random_upper_unitriangular;
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::partial_order::{OrderComparatorAutoLtByKey, OrderComparatorAutoAnyType};
        use crate::vectors::operations::Transforms;

        use std::slice::Iter;

        let matrix_size = 5;
        let modulus = 7;
        let array_mapping = VecOfVecSimple::new(
                                    vec![
                                        vec![(0, 1), (1, 2),                 (4, 0)],
                                        vec![        (1, 1), (2, 0),         (4, 1)],
                                        vec![                (2, 1), (3, 0), (4, 0)],
                                        vec![                        (3, 1), (4, 0)],
                                        vec![                                (4, 1)],
                                    ]
                                );
        let array_array_mapping_ref = & array_mapping;
        let iter_keymaj         = (0 .. matrix_size).rev();
        let ring_operator       =   PrimeOrderFieldOperator::new( modulus );

        // compute the inverse
        let inverse =   InverseOfTriangularArrayLazyAscending::new( 
            & array_array_mapping_ref,
            ring_operator.clone(),
            OrderComparatorAutoLtByKey::new(),
        );    
        
        // let `array_mapping_transformed` be the matrix obtained by reversing the order of columns of `array_mapping`
        let vector_transformer = |x: Cloned< Iter< '_, (usize,usize) > >| 
                                                                    x
                                                                        .rev()
                                                                        .map(   |(a,b)| 
                                                                                ( matrix_size - 1 - a, b) 
                                                                            )
                                                                        .collect_vec() ;

        let array_mapping_transformed: VecWiseTransformed< _, Cloned< Iter< '_, (usize,usize) > >, _>      =   VecWiseTransformed::new( 
                                            array_array_mapping_ref,
                                            vector_transformer,
                                        );

        // compute the codomain COMB of `array_mapping_transformed`
        //
        // NOTE: the codomain COMB of `array_mapping_transformed` is the inverse of `array_mapping`, where -- by contrast -- 
        //       the codomain COMB of `array_mapping` is the identity matrix (which is less interesting)
        let (array_comb, array_match) = get_codomain_comb_inv_off_diag_pivot_block(
            & array_mapping_transformed, // the matrix array_mapping_transformed ought to implement Copy
            iter_keymaj,
            ring_operator,
            OrderComparatorAutoAnyType,
            OrderComparatorAutoAnyType,
        );          


        for keymaj in 0 .. matrix_size {

            // obtain a row of the inverse
            let inv_row     =    inverse.view_major_ascend(keymaj).collect_vec();

            // obtain a row of the codomain COMB
            let ordmaj    =   array_match.keymaj_to_ordmaj( &keymaj ).unwrap();

            println!("KEYMAJ: {:?}", keymaj);
            println!("ORDMAJ: {:?}", ordmaj);       
            println!("view_major_ascend: {:?}", (& array_comb).view_major_ascend(ordmaj) );     

            let comb_off_diag_view =    (& array_comb)
                                        .view_major_ascend( ordmaj )
                                        // .iter()
                                        .map(   | (x,y)| 
                                                ( array_match.ordmaj_to_keymaj( x ), y ) // reindex the row from ordinals to major keys
                                            );
            let mut comb_view   =   comb_off_diag_view.collect_vec();
            comb_view.push( (keymaj, 1) );
            comb_view.sort();       
            
            // compute products
            let product_inv=   vector_matrix_multiply_major_ascend_simplified(
                                                                        inv_row.clone(),
                                                                        OracleRef::new( & array_mapping_transformed ),
                                                                        ring_operator.clone(),
                                                                        OrderComparatorAutoLtByKey::new()
                                                                    );
            let product_umatch=   vector_matrix_multiply_major_ascend_simplified(
                                                                        comb_view.clone(),
                                                                        OracleRef::new( & array_mapping_transformed ),
                                                                        ring_operator.clone(),
                                                                        OrderComparatorAutoLtByKey::new()
                                                                    );                                                                    
                                                                

            // compare
            println!("MATRIX:");
            for k in 0 .. matrix_size { println!("{:?}", array_array_mapping_ref.view_major_ascend(k).collect_vec()  ) }
            println!("MATRIX TRANSFORMED:");
            for k in 0 .. matrix_size { println!("{:?}", array_mapping_transformed.view_major_ascend(k)  ) }            
            println!("product_inv: {:?}", product_inv.clone().collect_vec() );
            println!("product_umatch: {:?}", product_umatch.clone().collect_vec());            
            println!("comb row len: {:?}, inv row len{:?}", comb_view.len(), inv_row.len());
            assert_eq!(comb_view, inv_row);
                            
        } 
        
        
    }    


    /// This test targets the initial computation of the pivot block of the inverse of the codomain COMB.
    ///
    /// The key idea of this test is the following fact: let M be a square upper-unitriangular matrix,
    /// and let N be the matrix obtained by reversing the order of columns of M.  Then the standard
    /// "cohomology algorithm," applied to N, produces a codomain COMB equal to M^{-1}.
    /// 
    /// This test applies the standard cohomology algorithm to compute a codomain COMB of N.  We 
    /// check to ensure that this codomain COMB equals M^{-1}.
    #[test]
    fn test_initial_decomposition_larger() {
        use itertools::Itertools;

        use crate::matrices::operations::umatch::row_major_only::get_codomain_comb_inv_off_diag_pivot_block;
        use crate::matrices::operations::multiply::vector_matrix_multiply_major_ascend_unsimplified;        
        use crate::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
        use crate::matrices::operations::transform_vector_wise::VecWiseTransformed;
        use crate::matrices::random_constructors::random_upper_unitriangular;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::partial_order::{OrderComparatorAutoLtByKey, OrderComparatorAutoAnyType};
        use crate::vectors::operations::Transforms;

        use std::slice::Iter;

        let matrix_size     =   10;
        let modulus             =   7;
        
        let array_mapping = random_upper_unitriangular( matrix_size, modulus );
        let array_array_mapping_ref = & array_mapping;
        let iter_keymaj         = (0 .. matrix_size).rev();
        let ring_operator       =   PrimeOrderFieldOperator::new( modulus );

        // compute the inverse
        let inverse =   InverseOfTriangularArrayLazyAscending::new( 
            & array_array_mapping_ref,
            ring_operator.clone(),
            OrderComparatorAutoLtByKey::new(),
        );    
        
        // let `array_mapping_transformed` be the matrix obtained by reversing the order of columns of `array_mapping`        
        let vector_transformer = |x: Cloned< Iter< '_, (usize,usize) > >| 
                                                                    x
                                                                        .rev()
                                                                        .map(   |(a,b)| 
                                                                                ( matrix_size - 1 - a, b) 
                                                                            )
                                                                        .collect_vec() ;

        let array_mapping_transformed: VecWiseTransformed< _, Cloned< Iter< '_, (usize,usize) > >, _>      =   VecWiseTransformed::new( 
                                            array_array_mapping_ref,
                                            vector_transformer,
                                        );

        // compute the codomain COMB of `array_mapping_transformed`
        //
        // NOTE: the codomain COMB of `array_mapping_transformed` is the inverse of `array_mapping`, where -- by contrast -- 
        //       the codomain COMB of `array_mapping` is the identity matrix (which is less interesting)        
        let (array_comb, array_match) = get_codomain_comb_inv_off_diag_pivot_block(
            & array_mapping_transformed,
            iter_keymaj,
            ring_operator,
            OrderComparatorAutoAnyType,
            OrderComparatorAutoAnyType,
        );          


        for keymaj in 0 .. matrix_size {

            // obtain a row of the inverse
            let inv_row     =    inverse.view_major_ascend(keymaj).collect_vec();

            // obtain a row of the codomain COMB
            let ordmaj    =   array_match.keymaj_to_ordmaj( &keymaj ).unwrap();

            let comb_off_diag_view =    (& array_comb)
                                        .view_major_ascend( ordmaj )
                                        // .iter()
                                        .map(   | (x,y)| 
                                                ( array_match.ordmaj_to_keymaj( x ), y ) // reindex the row from ordinals to major keys
                                            );
            let mut comb_view   =   comb_off_diag_view.collect_vec();
            comb_view.push( (keymaj, 1) );
            comb_view.sort();       
            
            // compute products
            let product_inv=   vector_matrix_multiply_major_ascend_simplified(
                                                                        inv_row.clone(),
                                                                        OracleRef::new( & array_mapping_transformed ),
                                                                        ring_operator.clone(),
                                                                        OrderComparatorAutoLtByKey::new()
                                                                    );
            let product_umatch=   vector_matrix_multiply_major_ascend_simplified(
                                                                        comb_view.clone(),
                                                                        OracleRef::new( & array_mapping_transformed ),
                                                                        ring_operator.clone(),
                                                                        OrderComparatorAutoLtByKey::new()
                                                                    );                                                                    
                                                                

            // compare
            println!("MATRIX:");
            for k in 0 .. matrix_size { println!("{:?}", array_array_mapping_ref.view_major_ascend(k).collect_vec()  ) }
            println!("MATRIX TRANSFORMED:");
            for k in 0 .. matrix_size { println!("{:?}", array_mapping_transformed.view_major_ascend(k)  ) }            
            println!("product_inv: {:?}", product_inv.clone().collect_vec() );
            println!("product_umatch: {:?}", product_umatch.clone().collect_vec());            
            println!("comb row len: {:?}, inv row len{:?}", comb_view.len(), inv_row.len());
            assert_eq!(comb_view, inv_row);
                            
        } 
        
        
    }    



    //  ===================================================================================
    //  RECOVERY OF COMBS -- COMPREHENSIVE
    //  ===================================================================================    


    /// Checks that Umatch decomposition is correct (using a small example matrix, D) in the following sense:
    /// R^{-1} * R = I
    /// C^{-1} * C = I
    /// R^{-1} * D * C = M
    /// And the rows of R, R^{-1}, C, and C^{-1} appear in strictly ascending order
    #[test]
    fn test_umatchrowmajor_comprehensive_small() {

        use crate::matrices::random_constructors::random_vec_of_vec_simple;
        use crate::utilities::partial_order::OrderComparatorAutoAnyType;
        use crate::matrices::debug::verify_that_product_is_identity;
        use crate::matrices::operations::umatch::row_major_only::ProductMatrixLazyMajorAscendSimplified;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::utilities::iterators::is_sorted::IsSortedBy;

        let num_indices_major           =   2;
        let num_indices_minor           =   3;
        // let approximate_density         =   0.3;
        let modulus                     =   3;
        // let allow_nonstructural_zero         =   true;

        let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
        let array_mapping_data       =   VecOfVecSimple::new( vec![ vec![(0,1), (1,2), (2,0)], vec![(2,0)]  ] );
        let array_mapping                               =   & array_mapping_data;

        let umatch 
            =   UmatchRowMajor::new( 
                    array_mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator.clone(), 
                    OrderComparatorAutoAnyType, 
                    OrderComparatorAutoAnyType,
                );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
        let array_matching = umatch.array_matching_ref();

        let array_comb_codomain = umatch.array_comb_codomain();
        let array_comb_codomain_inv = umatch.array_comb_codomain_inv();        
        let array_comb_domain = umatch.array_comb_domain();        
        let array_comb_domain_inv = umatch.array_comb_domain_inv(); 

        let array_comb_codomain_ref         =   OracleRef::new( & array_comb_codomain );
        let array_comb_codomain_inv_ref         =   OracleRef::new( & array_comb_codomain_inv );
        let array_comb_domain_ref         =   OracleRef::new( & array_comb_domain );
        let array_comb_domain_inv_ref         =   OracleRef::new( & array_comb_domain_inv );                        
        
        
        let product_domain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_domain_ref, array_comb_domain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );
        let product_codomain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_ref, array_comb_codomain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );        
        let product_codomain_comb_inv_times_mapping = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_inv_ref, array_mapping, ring_operator.clone(), OrderComparatorAutoAnyType );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrixLazyMajorAscendSimplified::new( product_codomain_comb_inv_times_mapping, array_comb_domain_ref, ring_operator.clone(), OrderComparatorAutoAnyType );                


        println!("array_mapping:");
        print_indexed_major_views( & array_mapping, 0 .. num_indices_major );
        println!("array_matching:");
        print_indexed_major_views( & array_matching, 0 .. num_indices_major );        
        println!("comb_domain:");
        print_indexed_major_views( & array_comb_domain, 0 .. num_indices_minor );        
        println!("comb_domain_inv:");
        print_indexed_major_views( & array_comb_domain_inv, 0 .. num_indices_minor );     
        println!("comb_codomain:");
        print_indexed_major_views( & array_comb_codomain, 0 .. num_indices_major );        
        println!("comb_codomain_inv:");
        print_indexed_major_views( & array_comb_codomain_inv, 0 .. num_indices_major );    
        println!("comb_codomain_inv * mapping * comb_domain:");
        print_indexed_major_views( & product_codomain_comb_inv_times_mapping_times_domain_comb, 0 .. num_indices_major );                                
        for column_index in 0 .. num_indices_minor {
            println!("{:?}", product_domain.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( product_domain.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }


        // check that the product of the domain comb with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, 1) ]
            ) 
        }

        // check that the product of the codomain comb with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, 1) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                array_matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }     
        
        // check that rows are sorted in strictly ascending order
        for keymaj in 0 .. num_indices_major { 
            assert!(    array_mapping.view_major_ascend( keymaj.clone()             ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_codomain.view_major_ascend( keymaj.clone()       ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_codomain_inv.view_major_ascend( keymaj.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_domain.view_major_ascend( keymaj.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_domain_inv.view_major_ascend( keymaj.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }                
    }



    /// Checks that Umatch decomposition is correct (using a random example matrix, D) in the following sense:
    /// R^{-1} * R = I
    /// C^{-1} * C = I
    /// R^{-1} * D * C = M   
    /// And the rows of R, R^{-1}, C, and C^{-1} appear in strictly ascending order 
    #[test]
    fn test_umatchrowmajor_comprehensive() {

        use crate::matrices::random_constructors::random_vec_of_vec_simple;
        use crate::utilities::partial_order::OrderComparatorAutoAnyType;
        use crate::matrices::debug::verify_that_product_is_identity;
        use crate::matrices::operations::umatch::row_major_only::ProductMatrixLazyMajorAscendSimplified;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::utilities::iterators::is_sorted::IsSortedBy;        

        let num_indices_major           =   50;
        let num_indices_minor           =   100;
        let approximate_density         =   0.05;
        let modulus                     =   3;
        let allow_nonstructural_zero         =   true;

        let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
        let array_mapping_data       =   random_vec_of_vec_simple( num_indices_major, num_indices_minor, approximate_density, modulus, allow_nonstructural_zero );
        let array_mapping                               =   & array_mapping_data;

        let umatch 
            =   UmatchRowMajor::new( 
                    array_mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator.clone(), 
                    OrderComparatorAutoAnyType, 
                    OrderComparatorAutoAnyType,
                );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
        let array_matching = umatch.array_matching_ref();

        let array_comb_codomain = umatch.array_comb_codomain();
        let array_comb_codomain_inv = umatch.array_comb_codomain_inv();        
        let array_comb_domain = umatch.array_comb_domain();        
        let array_comb_domain_inv = umatch.array_comb_domain_inv();    

        let array_comb_codomain_ref         =   OracleRef::new( & array_comb_codomain );
        let array_comb_codomain_inv_ref         =   OracleRef::new( & array_comb_codomain_inv );
        let array_comb_domain_ref         =   OracleRef::new( & array_comb_domain );
        let array_comb_domain_inv_ref         =   OracleRef::new( & array_comb_domain_inv );            
        

        let product_domain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_domain_ref, array_comb_domain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );
        let product_codomain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_ref, array_comb_codomain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );        
        let product_codomain_comb_inv_times_mapping = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_inv_ref, array_mapping, ring_operator.clone(), OrderComparatorAutoAnyType );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrixLazyMajorAscendSimplified::new( product_codomain_comb_inv_times_mapping, array_comb_domain_ref, ring_operator.clone(), OrderComparatorAutoAnyType );                        


        println!("array_mapping:");
        print_indexed_major_views( & array_mapping, 0 .. num_indices_major );
        println!("array_matching:");
        print_indexed_major_views( & array_matching, 0 .. num_indices_major );        
        println!("comb_domain:");
        print_indexed_major_views( & array_comb_domain, 0 .. num_indices_minor );        
        println!("comb_domain_inv:");
        print_indexed_major_views( & array_comb_domain_inv, 0 .. num_indices_minor );     
        println!("comb_codomain:");
        print_indexed_major_views( & array_comb_codomain, 0 .. num_indices_major );        
        println!("comb_codomain_inv:");
        print_indexed_major_views( & array_comb_codomain_inv, 0 .. num_indices_major );                        
        println!("comb_codomain_inv * mapping * comb_domain:");
        print_indexed_major_views( & product_codomain_comb_inv_times_mapping_times_domain_comb, 0 .. num_indices_major );                                        
        for column_index in 0 .. num_indices_minor {
            println!("{:?}", product_domain.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( product_domain.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }



        // check that the product of the domain comb with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, 1) ]
            ) 
        }

        // check that the product of the codomain comb with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, 1) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                array_matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }       
        
        // check that rows are sorted in strictly ascending order
        for keymaj in 0 .. num_indices_major { 
            assert!(    array_mapping.view_major_ascend( keymaj.clone()             ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_codomain.view_major_ascend( keymaj.clone()       ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_codomain_inv.view_major_ascend( keymaj.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_domain.view_major_ascend( keymaj.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    array_comb_domain_inv.view_major_ascend( keymaj.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }           

    }



    //  ===================================================================================
    //  RECOVERY OF COMBS -- TARGETTED AT SPECIFIC POINTS IN THE PROCESS
    //  ===================================================================================

    
    //  COMB DOMAIN INV (SMALL + LARGE)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #CombDomainInvImplementOracleMajorAscend    

    #[test]
    fn test_retreival() {
        use itertools::Itertools;

        use crate::matrices::operations::umatch::row_major_only::{get_codomain_comb_inv_off_diag_pivot_block, UmatchRowMajor, CombCodomainInvTimesMappingMatchedBlock};
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::rings::operator_structs::ring_native::RingNative;
        use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorAutoLtByKey};

        let array_mapping: VecOfVecSimple< usize, usize >   =   
        VecOfVecSimple::new(  
                vec![
                    vec![ (0, 1), (1, 1), (2, 2) ],
                    vec![ (0, 1),         (2, 1) ],
                    vec![         (1, 1), (2, 1) ],
                ]
            );
        let array_mapping_ref   =   & array_mapping; // the oracle trait is only implemented for *references* to a `VecOfVecSimple` object
        
        let order_comparator_minor                =   OrderComparatorAutoAnyType;     
        let order_comparator_major                =   OrderComparatorAutoAnyType;                
        let ring_operator       =   PrimeOrderFieldOperator::new( 13 );

        let array_mapping_ref   =   & array_mapping;
        let umatch  
            =   UmatchRowMajor::new( 
                        array_mapping_ref, 
                        (0..3).rev(), 
                        ring_operator.clone(),
                        order_comparator_minor.clone(),
                        order_comparator_major.clone(),                        
                    );
        let umatch_ref = & umatch;
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( umatch_ref);
        // let umatch_with_refs_ref = & umatch_with_refs;
        
        //  the "seed" matrix A (equal to the pivot block of the inverse of the codomain COMB times the pivot block of the matching array)
        let A   =   CombCodomainInvTimesMappingMatchedBlock::new( & umatch );
        let A_ref = &A;
        // println!("{:?}", umatch_with_refs.array_matching_ref());

        for keymaj in 1..3 { println!( "{:?}", A_ref.view_major_ascend(keymaj).into_iter().collect_vec() ) }

        //  the domain comb
        let comb_domain_inv = umatch.array_comb_domain_inv();
        let comb_domain_inv_ground_truth
                =   VecOfVecSimple::new(
                            vec![
                                vec![ (0, 1),         (2, 1) ],
                                vec![         (1, 1), (2, 1) ],
                                vec![                 (2, 1 )                        ],
                            ]
                        );
        let comb_domain_inv_ground_truth_ref = & comb_domain_inv_ground_truth;
        for keymaj in 0 .. 3 {
            // println!("GROUND TRUTH  : {:?}", comb_domain_inv_ground_truth_ref.view_major_ascend( keymaj ).collect_vec() );
            // println!("UNPACKED      : {:?}", comb_domain_inv.view_major_ascend( keymaj ).into_iter().collect_vec() );   
            // println!("SCALE FACTORS : {:?}", umatch_with_refs.array_matching.vec_snzval_ref() );    
            // println!("keymaj        : {:?}", keymaj );                                
            itertools::assert_equal(
                    comb_domain_inv_ground_truth_ref.view_major_ascend( keymaj ),
                    comb_domain_inv.view_major_ascend( keymaj ).into_iter(),
                )
        }
    }





    //  COMB DOMAIN (SMALL)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #CombDomainImplementOracleMajorAscend

    #[test]
    fn test_umatchrowmajor_comb_domain_small_example() {

        use crate::matrices::operations::umatch::row_major_only::{get_codomain_comb_inv_off_diag_pivot_block, UmatchRowMajor, CombCodomainInvTimesMappingMatchedBlock};
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::rings::operator_structs::ring_native::RingNative;
        use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorAutoLtByKey};   
        use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;     

        // let num_rows = 2; let num_cols = 2; let modulus = 7;
        // let array_mapping = VecOfVecSimple::new( vec![ vec![ (0usize,5usize), (1,5)], vec![ (1,6)]] );
        let num_rows = 1; let num_cols = 4; let modulus = 7;
        let array_mapping = VecOfVecSimple::new( vec![ vec![ (2usize, 6usize), (3,1)], ]  );        
        //  NOTE: array_mapping can be regarded as     [  0  0  6  1  ]
        let array_mapping_ref = & array_mapping;
        let ring_operator = PrimeOrderFieldOperator::new( modulus );
        let order_comparator_minor = OrderComparatorAutoAnyType;
        let order_comparator_major = OrderComparatorAutoAnyType;        

        let umatch_root = 
                UmatchRowMajor::new( 
                        array_mapping_ref, 
                        (0 .. num_rows).rev(), 
                        ring_operator.clone(),
                        order_comparator_minor.clone(),
                        order_comparator_major.clone(),
                    );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch_root );
        // let umatch_ref = & umatch_with_refs;

        let comb_domain = umatch_root.array_comb_domain();
        let comb_domain_inv = umatch_root.array_comb_domain_inv();     
        let array_mapping_matched_cols_only = umatch_root.array_mapping_matched_cols_only();                
                
        // check that C * C^{-1} = identity
        let c_times_c_inv = 
            ProductMatrixLazyMajorAscendSimplified::new( 
                    OracleRef::new( &comb_domain ), 
                    OracleRef::new( &comb_domain_inv ), 
                    ring_operator.clone(), 
                    OrderComparatorAutoLtByKey::new() 
                );

        println!("array_mapping:");
        print_indexed_major_views( & array_mapping_ref, 0 .. num_rows );
        println!("array_mapping_matched_cols_only:");
        print_indexed_major_views( & array_mapping_matched_cols_only, 0 .. num_rows );        
        println!("array_matching:");
        print_indexed_major_views( & umatch_root.array_matching_ref(), 0 .. num_rows );    
        println!("array_comb_codomain_inv_times_mapping_array_matched_block (recall that num_rows = {:?}):", num_rows);        
        print_indexed_major_views( && umatch_root.array_comb_codomain_inv_times_mapping_array_matched_block(), 0 .. num_rows );
        println!("comb_domain (recall that num_cols = {:?}) (THIS FUNCTION CALL SEEMS TO BREAK DOWN INTERNALLY WHERE array_comb_codomain_inv_times_mapping_array_matched_block IS CALLED):", num_cols);
        print_indexed_major_views( & comb_domain, 0 .. num_cols );        
        println!("comb_domain_inv:");
        print_indexed_major_views( & comb_domain_inv, 0 .. num_cols );                
        for column_index in 0 .. num_cols {
            println!("{:?}", c_times_c_inv.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( c_times_c_inv.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }

        // check D * C is right-reduced

    }


    //  COMB DOMAIN (LARGER + RANDOM)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #CombDomainImplementOracleMajorAscend)    


    #[test]
    fn test_umatchrowmajor_comb_domain() {

        use crate::matrices::operations::umatch::row_major_only::{get_codomain_comb_inv_off_diag_pivot_block, UmatchRowMajor, CombCodomainInvTimesMappingMatchedBlock};
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::rings::operator_structs::ring_native::RingNative;
        use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorAutoLtByKey};   
        use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;     

        let num_rows = 10; let num_cols = 10; let modulus = 7;
        let array_mapping = random_m_by_n_matrix(num_rows, num_cols, modulus);
        let array_mapping_ref = & array_mapping;
        let ring_operator = PrimeOrderFieldOperator::new( modulus );
        let order_comparator_minor = OrderComparatorAutoAnyType;
        let order_comparator_major = OrderComparatorAutoAnyType;        

        let umatch_root = 
                UmatchRowMajor::new( 
                        array_mapping_ref, 
                        (0 .. num_rows).rev(), 
                        ring_operator.clone(),
                        order_comparator_minor.clone(),
                        order_comparator_major.clone(),                        
                    );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch_root );
        // let umatch_ref = & umatch_with_refs;

        let comb_domain = umatch_root.array_comb_domain();
        let comb_domain_inv = umatch_root.array_comb_domain_inv();     
                
        // check that C * C^{-1} = identity
        let c_times_c_inv = 
            ProductMatrixLazyMajorAscendSimplified::new( 
                    OracleRef::new( & comb_domain ), 
                    OracleRef::new( & comb_domain_inv ), 
                    ring_operator.clone(), 
                    OrderComparatorAutoLtByKey::new() 
                );

        println!("array_mapping:");
        print_indexed_major_views( & array_mapping_ref, 0 .. num_rows );
        println!("array_matching:");
        print_indexed_major_views( & umatch_root.array_matching_ref(), 0 .. num_rows );        
        println!("comb_domain:");
        print_indexed_major_views( & comb_domain, 0 .. num_cols );        
        println!("comb_domain_inv:");
        print_indexed_major_views( & comb_domain_inv, 0 .. num_cols );                
        for column_index in 0 .. num_cols {
            println!("{:?}", c_times_c_inv.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( c_times_c_inv.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }

        // check D * C is right-reduced

    }





}    


























//  TESTS REGARDING WHAT REFERENCES CAN / CANNOT BE RETURNED BY A FUNCTION
//  ====================================================================================================




fn test_oracles() {
    let matrix_unprepended = VecOfVecSimple::new( vec![ vec![ (1,1)] ]);
    let matrix_unprepended_ref: & VecOfVecSimple< i32, i32 > = & matrix_unprepended;

    let a = matrix_unprepended_ref.view_major_ascend(0);

    let prepended : PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< & VecOfVecSimple<i32, i32>, Cloned< Iter< '_, (i32, i32) >>, i32, i32 >
    = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( matrix_unprepended_ref, 1 );

    // let a = prepended.view_major_ascend(0);

    // < usize, Cloned< Iter< 'a, (KeyMin, SnzVal) > >, > 

    //-------------------------------------------------------------------------------------------------------------------------

    let matrix_scalar = ScalarMatrix::new( 1 );

    let prepended_scalar 
    = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( & matrix_scalar, 1 );    

    let b = prepended_scalar.view_major_ascend(0);

    //-------------------------------------------------------------------------------------------------------------------------

    let data = vec![vec![(1usize,1usize)]];
    let matrix_vecvecfromborrow = VecOfVecSimpleFromBorrow::new( &data );

    let prepended_vecvecfromborrow : PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< VecOfVecSimpleFromBorrow<'_, usize, usize>, Cloned< Iter< '_, (usize, usize) >>, usize, usize > 
    = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( matrix_vecvecfromborrow, 1 );    

    let b = prepended_vecvecfromborrow.view_major_ascend(0usize);    


}

fn test_oracle_with_lifetime< 'a >( matrix_unprepended_ref_ref: &'a &'a VecOfVecSimple< usize, i32 > ) {

    let prepended : PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< &'a VecOfVecSimple<usize, i32>, Cloned< Iter< 'a, (usize, i32) >>, usize, i32 >
    = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( matrix_unprepended_ref_ref, 1 ); 

}

// fn test_reference< 'a >( variable_ref: &'a usize ) -> &'a usize { variable_ref } 


// CHECK WHERE REFERENCES CAN LIE
// ----------------------------------------------------------------------------------------------------

pub struct UmatchRowMajorVAR< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor > 
    where   
        KeyMin:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        KeyMaj:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        // OrderComparatorAutoLtByKey<usize, SnzVal, (usize, SnzVal)>: StrictlyLess< (usize, SnzVal)>, // this seems extraneous but the compiler seems to want it
        ViewMajorAscend:            IntoIterator,
        ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,
{
    array_mapping:                      &'a ArrayMapping,
    array_matching:                     GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >,
    array_comb_codomain_inv_matched_block_off_diagonal: &'a VecOfVecSimple< usize, SnzVal >,
    // array_comb_codomain_inv_matched_block_off_diagonal: VecOfVec< ( usize, SnzVal ), OrderComparatorAutoLtByFirstTupleEntry >,
    // array_comb_codomain_inv_matched_block_off_diagonal: Vec< Vec< ( usize, SnzVal ) > >,
    ring_operator:                      RingOperator,
    order_comparator_minor:                   OrderComparatorMinor,
    phantom_viewmajorascend:            PhantomData< ViewMajorAscend >, // required b/c otherwise the compiler complains that the type parameter `ViewMajorAscend` is unused
}

// fn test_where_refs_can_live
//         < 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor > 
//         (  umatch_rowmajor_var: UmatchRowMajorVAR< 'a, ArrayMapping, ViewMajorAscend, KeyMin, KeyMaj, SnzVal, RingOperator, OrderComparatorMinor, OrderComparatorMajor >, scalar: SnzVal  ) 
//         ->
//         PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< 'a, &'a VecOfVecSimple<usize, SnzVal>, Cloned< Iter< 'a, (usize, SnzVal) >>, usize, SnzVal >
//     where
//         KeyMin:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         KeyMaj:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ViewMajorAscend:            IntoIterator,
//         ViewMajorAscend::Item:      KeyValGet< KeyMin, SnzVal >,  
//         SnzVal:     Clone,              
// {
//     let prepended : PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< 'a, &'a VecOfVecSimple<usize, SnzVal>, Cloned< Iter< 'a, (usize, SnzVal) >>, usize, SnzVal >
//     = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( & umatch_rowmajor_var.array_comb_codomain_inv_matched_block_off_diagonal, scalar );  

//     let _a = <  
//         PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend<
//                 'a,
//                 &'a VecOfVecSimple< usize, SnzVal >,
//                 // &'a VecOfVec< (usize, SnzVal), OrderComparatorAutoLtByFirstTupleEntry >,                
//                 Cloned< Iter< 'a, (usize, SnzVal) > >,
//                 // &'a [(usize, SnzVal)],
//                 usize,
//                 SnzVal,
//             >
//         as OracleMajorAscend < usize, std::iter::Chain < std::iter::Once < (usize, SnzVal) >, Cloned< Iter< 'a, (usize, SnzVal) > >> >>::
//         view_major_ascend( &prepended, 0 ) ;

//     let _b = prepended.view_major_ascend(0);

    
//     return prepended
// }

