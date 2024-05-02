//! U-match factorization for row-major oracles.
//!
//! # To develop
//! 
//! - command to get barcode (or birth/death pairs) in dimension d
//!   - take two iterators as input
//!   - return two iterators as output, excluding pairs with equal birth/death time
//! - command to get element of jordan basis
//!   - if index is paired row index, multiply the paired col of C by D
//!   - otherwise return the column of C
//! - command to get element of dual basis
//!   - ?
//! - make new constructor function for boundary matrices that skips matched entries
//! - place COMB's in their own submodule
//! - place technical components in their own submodule (so people can use things that are tested)
//! - (possibly) add minor keyval pair type parameter to umatch, and create 3 contstructor functions
//!   - note that if you incorporate lazy pairs then you need to have an order comparator for the key-val pairs with major keys
//! 
//! # To-do (primary)
//! 
//! - modify Umatch struct to include entry type for minor views (this alleviates need to implement oracle minor descend in order to get an initial decomposition)
//! - make the element type of a ring operator an ASSOCIATED type
//! - make `one` and `zero` into methods on ring operators
//! - get Ryan's help to write a small demo library crate and small binary crate 
//! - try removing all the Debug requirements
//! 
//! # To-do
//! 
//! - run experiments to see what's faster: multiplicatoin with inverse or triangle solve
//! - determinant
//! - solve operations: Ax = b, where
//!   - A = D
//!     - umatch.solve_Dx_equal_b
//!     - umatch.solve_xD_equal_b
//!   - A = basis for kernel
//!     - umatch.solve_Kx_equal_b
//!     - umatch.solve_xK_equal_b
//!   - A = basis for image
//!     - umatch.solve_Ix_equal_b
//!     - umatch.solve_xI_equal_b
//! - compute basis for kernel, image
//!   - specifically, iterator that runs over all kernel / image vectors
//!     - umatch.kernel_basis_iter( col_index_iter )
//!     - umatch.image_basis_iter( col_index_iter )
//! - (separate, requires chain complex) jordan basis
//! 
//! # Example
//! 
//! ```
//! use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
//! use oat_rust::matrices::matrix_types::{vec_of_vec::VecOfVecSimple, oracle_ref::OracleRef};
//! use oat_rust::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorLtByKey, OrderComparatorAutoLt};
//! use oat_rust::matrices::debug::verify_that_product_is_identity;
//! use oat_rust::matrices::operations::umatch::row_major::{UmatchRowMajor, new_umatchrowmajor};
//! use oat_rust::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;
//! use oat_rust::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend};
//! use oat_rust::utilities::iterators::is_sorted::IsSortedBy;
//! use itertools::Itertools;
//! 
//! // define the coefficient ring
//! let modulus                     =   5;
//! let ring_operator                   =   PrimeOrderFieldOperator::new( modulus );        
//! 
//! // define the matrix we wish to factor
//! let num_indices_major           =   2;
//! let num_indices_minor           =   3;
//! let array_mapping_data          =   VecOfVecSimple::new( 
//!                                         vec![   
//!                                                     vec![(0,1), (1,2), (2,3)], 
//!                                                     vec![              (2,1)]  
//!                                         ] 
//!                                     );
//! let array_mapping               =   & array_mapping_data;
//! 
//! // compute the U-match factorization
//! let umatch
//!     =   new_umatchrowmajor(
//!             array_mapping,  // the matrix we wish to factor
//!             (0..num_indices_major).rev(), // an iterator that runs over all row indices, from bottom to top
//!             ring_operator.clone(), // the operator for the coefficient ring
//!             OrderComparatorAutoLt::<usize>::new(), // order comparator for the entries in each row
//!             OrderComparatorAutoLt::<usize>::new(), // order comparator for the entries in each column
//!         );
//! 
//! // extract R, R^{-1}, C, C^{-1}, and M
//! let array_matching              =   umatch.array_matching_ref();
//! let array_comb_codomain         =   umatch.array_comb_codomain();
//! let array_comb_codomain_inv     =   umatch.array_comb_codomain_inv();        
//! let array_comb_domain           =   umatch.array_comb_domain();        
//! let array_comb_domain_inv       =   umatch.array_comb_domain_inv(); 
//! 
//! // get references to R, R^{-1}, C, C^{-1}, and M        
//! let array_comb_codomain_ref     =   OracleRef::new( & array_comb_codomain );
//! let array_comb_codomain_inv_ref =   OracleRef::new( & array_comb_codomain_inv );
//! let array_comb_domain_ref       =   OracleRef::new( & array_comb_domain );
//! let array_comb_domain_inv_ref   =   OracleRef::new( & array_comb_domain_inv );   
//! 
//! // compute some products
//! let product_domain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_domain_ref, array_comb_domain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );
//! let product_codomain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_ref, array_comb_codomain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );        
//! let product_codomain_comb_inv_times_mapping = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_inv_ref, array_mapping, ring_operator.clone(), OrderComparatorAutoAnyType );      
//! let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrixLazyMajorAscendSimplified::new( product_codomain_comb_inv_times_mapping, array_comb_domain_ref, ring_operator.clone(), OrderComparatorAutoAnyType );                
//! 
//! 
//! // check that the product of the domain comb with its inverse is identity: C * C^{-1} = I
//! for keymin in 0 .. num_indices_minor { 
//!     assert_eq!(
//!         product_domain.view_major_ascend( keymin ).collect_vec(),
//!         vec![ (keymin, 1) ]
//!     ) 
//! }
//! 
//! // check that the product of the codomain comb with its inverse is identity: R * R^{-1} = I
//! for keymaj in 0 .. num_indices_major { 
//!     assert_eq!(
//!         product_codomain.view_major_ascend( keymaj ).collect_vec(),
//!         vec![ (keymaj, 1) ]
//!     ) 
//! }    
//! 
//! // check the factorization: R^{-1} * D * C = M
//! for keymaj in 0 .. num_indices_major { 
//!     assert_eq!(
//!         product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
//!         array_matching.view_major_ascend( keymaj ).collect_vec()
//!     ) 
//! }    
//! ```
//! 
//! 


use itertools::Itertools;
use serde_json::ser::CompactFormatter;

use crate::entries::{KeyValSet, KeyValNew};
use crate::matrices::operations::transform_entry_wise::ReindexSquareMatrix;
use crate::matrices::matrix_types::matching::{GeneralizedMatchingArrayWithMajorOrdinals, GeneralizedMatchingArray};
use crate::matrices::matrix_types::oracle_ref::OracleRef;
use crate::matrices::matrix_types::prepend_viewmaj::PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend;
use crate::matrices::matrix_types::scalar_matrices::ScalarMatrix;
use crate::matrices::matrix_types::vec_of_vec::{VecOfVec, VecOfVecSimple, VecOfVecSimpleFromBorrow, VecOfVecSimpleViewMinorDescend};
use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMajor, OracleMinorDescend, IndicesAndCoefficients};
use crate::utilities::functions::evaluate::{ IdentityFunction, EvaluateFunctionFnMutWrapper };
use crate::utilities::functions::compose::ComposeFunctions;
use crate::utilities::iterators::general::{SkipUntil, IterTwoType, IterWrappedVec, OncePeekable, MapByTransform};
use crate::utilities::iterators::merge::heap_of_iterators::{HitMerge, hit_bulk_insert, hit_merge_by_predicate};
use crate::entries::{KeyValGet, ReindexEntry};
use crate::rings::operator_traits::{Semiring, Ring, DivisionRing};
use crate::utilities::iterators::merge::two_iter_different_type::MergeTwoItersByOrderComparator;
use crate::utilities::partial_order::{StrictlyLess, OrderComparatorAutoLtByKey, OrderComparatorLtByKey, OrderComparatorReverse, InferTotalOrderFromStrictlyLess, DecideOrder};
use crate::utilities::sets::MapHasKeyOrSequenceHasElementRefWrapper;
use crate::vectors::linear_combinations::LinearCombinationSimplified;
use crate::vectors::operations::{Scale, Transforms, Simplify, OnlyIndicesInsideCollection, OnlyIndicesOutsideCollection, ChangeIndexSimple, Negate, ChangeEntryType};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::{Cloned, Peekable, Once, Chain};
use std::marker::PhantomData;
use std::slice::{Iter};

use crate::matrices::operations::solve::echelon::{EchelonSolveMajorAscendWithMinorKeys, EchelonSolveMinorDescendWithMinorKeys, EchelonSolveMinorDescendWithMajorKeys};
use crate::matrices::operations::multiply::{ProductMatrixLazyMajorAscendSimplified, vector_matrix_multiply_major_ascend_simplified, vector_matrix_multiply_minor_descend_simplified};
use crate::matrices::operations::transform_vector_wise::{OnlyKeyMinInsideCollection, OnlyKeyMinOutsideCollection, OnlyKeyMajInsideCollection, OnlyKeyMajOutsideCollection};
use crate::matrices::operations::solve::triangle::TriangularSolveAscend;






// -------
















//  =========================================================================================================
//  CALCULATION OF THE (COMPRESSED) MATCHING ARRAY AND CODOMAIN COMB
//  =========================================================================================================




pub trait ParetoShortCircuit<T> {
    fn pareto_short_circuit(& self) -> Option< T >;
}



//  HELPER FUNCTION 
//  ----------------------------------------------------------------------------------

/// Helper function for `try_writing_next_viewmaj_of_codomain_comb_to_buffer`.  That function
/// involves multiplying a row of a codomain COMB by a mapping array; sometimes entries
/// at the beginning of the resulting row vector can be deleted without penalty; the 
/// current function performs this deletion on a (scalar multiple) of an individual row of the 
/// mapping array.
fn codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array<
        ArrayMapping,
        RingOperator,
        OrderComparatorViewMajorAscendEntry,
    > 
    (
        codomain_comb_inv_entry:            ( usize, ArrayMapping::SnzVal ),
        scale_factor:                       ArrayMapping::SnzVal,        
        truncation_limit:                   & ArrayMapping::ViewMajorAscendEntry,
        array_mapping:                      & ArrayMapping,     
        array_matching:                     & GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,        
        ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
        order_comparator_viewmajorascendentry:               OrderComparatorViewMajorAscendEntry,            
    ) 
    ->
    Peekable<
            Scale < 
                    ArrayMapping::ViewMajorAscendIntoIter,  // a major view of the mapping array
                    ArrayMapping::KeyMin,
                    RingOperator, 
                    ArrayMapping::SnzVal,
                >, 
        >
    where 
        ArrayMapping:                 OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                 Clone + Hash + PartialEq + std::cmp::Eq,
        ArrayMapping::KeyMaj:                 Clone + Hash + PartialEq + std::cmp::Eq,         
        ArrayMapping::SnzVal:                 Clone,               
        RingOperator:           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess <  ArrayMapping::ViewMajorAscendEntry >,  
        ArrayMapping::ViewMajorAscend:        IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:  KeyValSet < ArrayMapping::KeyMin, ArrayMapping::SnzVal >,  
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
                    .skip_until( |x| order_comparator_viewmajorascendentry.strictly_less( truncation_limit, x  )  );
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
/// - the field `order_comparator_viewmajorascendentry` correctly determines which minor key precedes which 
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
        ArrayMapping,
        RingOperator,
        OrderComparatorViewMajorAscendEntry,
    > 
    (
        codomain_comb_inv_off_diag_view_buffer:      &mut Vec< ( usize, ArrayMapping::SnzVal ) >,          
        entries_to_elim_simplified_heap:    & mut Simplify <
                                                        HitMerge < 
                                                                Peekable<
                                                                        // the HitMerge struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                                        Scale < 
                                                                                ArrayMapping::ViewMajorAscendIntoIter,  // a major view of the mapping array
                                                                                ArrayMapping::KeyMin,
                                                                                RingOperator, // the ring_operator operator
                                                                                ArrayMapping::SnzVal,
                                                                            >, 
                                                                    >,
                                                                // the thing that declares whether one major key comes before of after another    
                                                                OrderComparatorViewMajorAscendEntry                                                                
                                                            >,
                                                        ArrayMapping::KeyMin,
                                                        RingOperator,
                                                        ArrayMapping::SnzVal,
                                                    >,          
        array_mapping:                      & ArrayMapping,    
        array_matching:                     & GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
        array_codomain_comb_inv_off_diag:   & Vec< Vec< (usize, ArrayMapping::SnzVal) > >,  // the off-diagonal part of the (block indexed by pivot indices of) inverse of the codomain COMB
        // array_codomain_comb_inv_off_diag:   & VecOfVec<  // the off-diagonal part of the (block indexed by pivot indices of) inverse of the codomain COMB
        //                                                 ( usize, ArrayMapping::SnzVal ), 
        //                                                 OrderComparatorAutoLtByFirstTupleEntry, 
        //                                             >,          
        ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
        order_comparator_viewmajorascendentry:             OrderComparatorViewMajorAscendEntry,            
    ) 
    ->

    Option< ArrayMapping::ViewMajorAscendEntry >

    where 
        ArrayMapping:                 OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                 Clone + Hash + PartialEq + std::cmp::Eq + Debug, // !!! remove debug
        ArrayMapping::KeyMaj:                 Clone + Hash + PartialEq + std::cmp::Eq + Debug,      // !!! remove debug   
        ArrayMapping::SnzVal:                 Clone + Debug,                // remove debug
        RingOperator:           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:   Clone + StrictlyLess <  ArrayMapping::ViewMajorAscendEntry >,  
        ArrayMapping::ViewMajorAscend:        IntoIterator,        // !!!!! REMOVE THE DEBUG + CLONE
        ArrayMapping::ViewMajorAscendEntry:  KeyValSet < ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug,       // !!!!! REMOVE THE DEBUG + CLONE REQUIREMENT WHEN DONE DEBUGGING!
        HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorViewMajorAscendEntry>: Clone, // !!! remove this later
{
    // println!("initiating construction of next row");  // !!! delete later
    // println!("array_matching {:?}", array_matching);  // !!! delete later    

    // we use this while loop because a for-loop results in a borrowing conflict (the conflict arrises because the for loop 
    // iteratos over `entries_to_elim_simplified_heap`, but we modify `entries_to_elim_simplified_heap` within the for-loop)
    while let Some( leading_entry_to_eliminate ) = entries_to_elim_simplified_heap.next() {

        // if entries_to_elim_simplified_heap.unsimplified.len() > 10 { println!("styx[[{:?}]]", entries_to_elim_simplified_heap.unsimplified.len()) }

        // println!("WHILE LOOP CHECKPOINT: leading_entry_to_eliminate: {:?}", &leading_entry_to_eliminate); // !!! REMOVE LATER
        // println!("WHILE LOOP CHECKPOINT: entries_to_elim_simplified_heap: {:?}", entries_to_elim_simplified_heap.clone().collect_vec() ); // !!! REMOVE LATER        
    
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
                                // println!(
                                //         "codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array:      {:?}",
                                //         codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array(
                                //                 codomain_comb_inv_entry.clone(),
                                //                 scale_factor.clone(),
                                //                 & leading_entry_to_eliminate, // truncation_limit:                   ArrayMapping::ViewMajorAscendEntry,
                                //                 array_mapping, //                     & ArrayMapping,     
                                //                 array_matching, //                    & GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,        
                                //                 ring_operator.clone(), //                    RingOperator, // the operator for the coefficient ring_operator    
                                //                 order_comparator_viewmajorascendentry.clone() //               OrderComparatorViewMajorAscendEntry,                                      
                                //             )
                                //             .collect_vec()
                                //     );

                                codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array(
                                        codomain_comb_inv_entry,
                                        scale_factor.clone(),
                                        & leading_entry_to_eliminate, // truncation_limit:                   ArrayMapping::ViewMajorAscendEntry,
                                        array_mapping, //                     & ArrayMapping,     
                                        array_matching, //                    & GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,        
                                        ring_operator.clone(), //                    RingOperator, // the operator for the coefficient ring_operator    
                                        order_comparator_viewmajorascendentry.clone() //               OrderComparatorViewMajorAscendEntry,                                      
                                    )
                                }
                            )
                );            
            // if entries_to_elim_simplified_heap.unsimplified.len() > 10 { println!("styx2[[{:?}]]", entries_to_elim_simplified_heap.unsimplified.len()) }                                         

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
            // if entries_to_elim_simplified_heap.unsimplified.len() > 10 { println!("styxdone[[{:?}]]", entries_to_elim_simplified_heap.unsimplified.len()) }                                         
        } else {

            // REMARK: in this case we cannot clear `leading_entry_to_eliminate` via elementary operations; 
            // therefore `leading_entry_to_eliminate` is a matched (i.e. pivot) entry

            return Some( leading_entry_to_eliminate );
        }   
    }

    // println!("terminating construction of next row"); // !!! delete later
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
/// **It is important that `order_comparator_viewmajorascendentry` compares order of two entries based *solely* on
/// the value of the associated indices (not the associated coefficients).  This is easy
/// to get wrong, and it's hard to debug, so we are keeping the function private for now**
pub fn codomain_comb_inv_off_diag_pivot_block_unsafecomparator< 'a, ArrayMapping, RingOperator, IterKeyMaj, OrderComparatorOfEntries >
    ( 
            array_mapping:                      &'a ArrayMapping,  // use of a reference poses no problem for future calls to this function, since output has no references        
            iter_keymaj:                        IterKeyMaj,
            ring_operator:                      RingOperator,            
            order_comparator_of_entries:        OrderComparatorOfEntries,
            // mut order_comparator_of_keymaj:     OrderComparatorViewMinorDescendEntry,
    ) 
    -> 
    ( 
        VecOfVecSimple< usize, ArrayMapping::SnzVal >, 
        // VecOfVec< ( usize, ArrayMapping::SnzVal ), OrderComparatorAutoLtByFirstTupleEntry  >,         
        GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
    )
    where   ArrayMapping:                     OracleMajorAscend + IndicesAndCoefficients,
            IterKeyMaj:                 Iterator < Item = ArrayMapping::KeyMaj >,
            ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq + Debug, 
            ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            ArrayMapping::SnzVal:                     Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
            ArrayMapping::ViewMajorAscend:            IntoIterator + Clone, // !!! remove clone
            ArrayMapping::ViewMajorAscendEntry:      KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderComparatorOfEntries:   Clone + StrictlyLess <  ArrayMapping::ViewMajorAscendEntry >, // !!! remove clone
            // OrderComparatorViewMinorDescendEntry:       StrictlyLess<  ArrayMapping::KeyMaj >,
            HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorOfEntries>: Clone // !!!! remove this

{
    
    // Initialize some objects
    let mut entries_to_elim_simplified_heap    =   HitMerge::new( order_comparator_of_entries.clone() )
                                                                            .simplify( ring_operator.clone() ); // an empty (simplified) merger of iterators
    let mut array_matching          =   GeneralizedMatchingArrayWithMajorOrdinals::new(); // an all-zero generalized matching array
    let mut array_codomain_comb_inv_off_diag: Vec< Vec< ( usize, ArrayMapping::SnzVal ) > >  =   Vec::new();    
    // let mut array_codomain_comb_inv_off_diag: VecOfVec<(usize, ArrayMapping::SnzVal), OrderComparatorAutoLtByFirstTupleEntry >       =   VecOfVec::new(  vec![], OrderComparatorAutoLtByFirstTupleEntry::new()  );
    let mut codomain_comb_inv_off_diag_view_buffer   =   Vec::new();
    let mut codomain_comb_inv_off_diag_view_buffer_simplified   =   Vec::new();    

    // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
    // set up a variable to track the last major index produced by `iter_keymaj`; we will use this to ensure that the major keys returned by `iter_keymaj` appear in strictly ascending order
    // let mut prior_keymaj_opt = None;

    // let mut counter = 0;
    // build the (pivot block of the) codomain COMB row by row
    for keymaj in iter_keymaj {

        // counter +=1;
        // print!("row({:?})", counter);

        // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
        // check that this major key is strictly greater than the last
        // if let Some( ref prior_keyaj ) = prior_keymaj_opt {
        //     match order_comparator_of_keymaj.strictly_less( prior_keyaj, &keymaj ) {
        //         true => {  },
        //         false => { panic!("Umatch factorization cannot be completed because the sequence of major indices is not strictly ascending") }
        //     }               
        // }
        // prior_keymaj_opt.replace( keymaj.clone() );

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
                // println!("leading_entry_uneliminable {:?}", &leading_entry_uneliminable);              
                // update the matching array
                array_matching.push( pivot_keymin, keymaj, pivot_coeff ); 
                // sort the buffer vector
                codomain_comb_inv_off_diag_view_buffer.sort_by( |a, b| b.0.cmp( &a.0 ) ); // note this yields DESCENDING ORDER -- which is what we want, because the order of indices is currently inverted (the greatest keymaj gets the lowest ordinal, and the least keymaj gets the highest ordinal; we will correct this later)
                // println!("THE BUFFER:      {:?}", &codomain_comb_inv_off_diag_view_buffer);
                // simplify the sequence of entries in the buffer (note that the buffer itself could have entries with duplicate indices)
                // !!! hypothetically there's a opportunity for a small optimization here, where we simultaneously simplify and sort at the same time; this would avoid a few sorting operations
                let iter_to_push     =   codomain_comb_inv_off_diag_view_buffer
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










/// Returns (i) the block of the codomain COMB indexed by pivot indices, and (ii)) the matching matrix; output is correct ONLY if pivot column indices and pivot row indices are disjoint.
/// 
/// For details on this factorization, see [this preprint](https://arxiv.org/pdf/2108.08831.pdf).
/// 
/// This U-match factorization is computed via the standard "cohomology algorithm,"
/// using the clear/compress/twist optimization.  This optimization skips over rows of the matrix indexed by keys
/// that already index pivot columns.
/// 
/// 
/// **It is important that `order_comparator_viewmajorascendentry` compares order of two entries based *solely* on
/// the value of the associated indices (not the associated coefficients).  This is easy
/// to get wrong, and it's hard to debug, so we are keeping the function private for now**
pub fn codomain_comb_inv_off_diag_pivot_block_unsafecomparator_skipmatched< 'a, ArrayMapping, RingOperator, IterKeyMaj, KeyBoth, OrderComparatorOfEntries >
    ( 
            array_mapping:                      &'a ArrayMapping,  // use of a reference poses no problem for future calls to this function, since output has no references        
            iter_keymaj:                        IterKeyMaj,
            ring_operator:                      RingOperator,            
            order_comparator_of_entries:        OrderComparatorOfEntries,
            // mut order_comparator_of_keymaj:     OrderComparatorViewMinorDescendEntry,
    ) 
    -> 
    ( 
        VecOfVecSimple< usize, ArrayMapping::SnzVal >, 
        // VecOfVec< ( usize, ArrayMapping::SnzVal ), OrderComparatorAutoLtByFirstTupleEntry  >,         
        GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
    )
    where   ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients< KeyMaj = KeyBoth, KeyMin = KeyBoth >,
            IterKeyMaj:                             Iterator < Item = ArrayMapping::KeyMaj >,
            ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq + Debug, 
            ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            ArrayMapping::SnzVal:                   Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
            ArrayMapping::ViewMajorAscend:          IntoIterator + Clone + ParetoShortCircuit<ArrayMapping::ViewMajorAscendEntry>, // !!! remove clone
            ArrayMapping::ViewMajorAscendEntry:     KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderComparatorOfEntries:               Clone + StrictlyLess <  ArrayMapping::ViewMajorAscendEntry >, // !!! remove clone
            // OrderComparatorViewMinorDescendEntry:       StrictlyLess<  ArrayMapping::KeyMaj >,
            HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorOfEntries>: Clone // !!!! remove this

{
    
    // Initialize some objects
    let mut entries_to_elim_simplified_heap    =   HitMerge::new( order_comparator_of_entries.clone() )
                                                                            .simplify( ring_operator.clone() ); // an empty (simplified) merger of iterators
    let mut array_matching          =   GeneralizedMatchingArrayWithMajorOrdinals::new(); // an all-zero generalized matching array
    let mut array_codomain_comb_inv_off_diag: Vec< Vec< ( usize, ArrayMapping::SnzVal ) > >  =   Vec::new();    
    // let mut array_codomain_comb_inv_off_diag: VecOfVec<(usize, ArrayMapping::SnzVal), OrderComparatorAutoLtByFirstTupleEntry >       =   VecOfVec::new(  vec![], OrderComparatorAutoLtByFirstTupleEntry::new()  );
    let mut codomain_comb_inv_off_diag_view_buffer   =   Vec::new();
    let mut codomain_comb_inv_off_diag_view_buffer_simplified   =   Vec::new();    

    // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
    // set up a variable to track the last major index produced by `iter_keymaj`; we will use this to ensure that the major keys returned by `iter_keymaj` appear in strictly ascending order
    // let mut prior_keymaj_opt = None;

    // let mut counter = 0;
    // let pb = indicatif::ProgressBar::new(100000);
    // let mut sc_counter = 0;

    // build the (pivot block of the) codomain COMB row by row
    for keymaj in iter_keymaj {


        // print!("row={}",counter);
        // pb.println(format!("row = {}", counter));
        // counter +=1;        
        // print!("out");        
        // print!("simplex={:?}",keymaj.clone());
        // use std::{thread, time};
        // thread::sleep(time::Duration::from_millis(10));
    
        if array_matching.contains_keymin( & keymaj ) { 
            continue }      
          

        // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
        // check that this major key is strictly greater than the last
        // if let Some( ref prior_keyaj ) = prior_keymaj_opt {
        //     match order_comparator_of_keymaj.strictly_less( prior_keyaj, &keymaj ) {
        //         true => {  },
        //         false => { panic!("Umatch factorization cannot be completed because the sequence of major indices is not strictly ascending") }
        //     }               
        // }
        // prior_keymaj_opt.replace( keymaj.clone() );

        // clear the collection of entries to eliminate
        entries_to_elim_simplified_heap.unsimplified.clear();        

        // if possible, short circuit
        let next_iter = array_mapping.view_major_ascend( keymaj.clone() );
        if let Some( pivot_entry ) = next_iter.pareto_short_circuit() {
            let pivot_keymin = pivot_entry.key();   let pivot_coeff = pivot_entry.val();
            array_matching.push( pivot_keymin, keymaj, pivot_coeff );
            array_codomain_comb_inv_off_diag.push( vec![] );
            // print!("sc: {}", sc_counter );
            // sc_counter += 1;
            continue
        }

        // insert the sequence of entries in row `keymaj`
        entries_to_elim_simplified_heap.unsimplified.insert_one_iter(
                next_iter
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
                // println!("leading_entry_uneliminable {:?}", &leading_entry_uneliminable);              
                // update the matching array
                array_matching.push( pivot_keymin, keymaj, pivot_coeff ); 
                // sort the buffer vector
                codomain_comb_inv_off_diag_view_buffer.sort_by( |a, b| b.0.cmp( &a.0 ) ); // note this yields DESCENDING ORDER -- which is what we want, because the order of indices is currently inverted (the greatest keymaj gets the lowest ordinal, and the least keymaj gets the highest ordinal; we will correct this later)
                // println!("THE BUFFER:      {:?}", &codomain_comb_inv_off_diag_view_buffer);
                // simplify the sequence of entries in the buffer (note that the buffer itself could have entries with duplicate indices)
                // !!! hypothetically there's a opportunity for a small optimization here, where we simultaneously simplify and sort at the same time; this would avoid a few sorting operations
                let iter_to_push     =   codomain_comb_inv_off_diag_view_buffer
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
/// We use `VecOfVecSimple< usize, ArrayMapping::SnzVal >` to store the pivot block of the inverse of 
/// the codomain comb rather than a `VecOfVec< ... >` for the following reasons: (1) we typically
/// need 
/// 
pub fn get_codomain_comb_inv_off_diag_pivot_block< 'a, ArrayMapping, RingOperator, IterKeyMaj, OrderComparatorViewMajorAscendEntry, > // OrderComparatorViewMinorDescendEntry >
    ( 
            array_mapping:                      &'a ArrayMapping,             
            iter_keymaj:                        IterKeyMaj,
            ring_operator:                      RingOperator,            
            order_comparator_for_keymin:        OrderComparatorViewMajorAscendEntry,
            // order_comparator_for_keymaj:        OrderComparatorViewMinorDescendEntry,            
    ) 
    -> 
    ( 
        // VecOfVec< ( usize, ArrayMapping::SnzVal ), OrderComparatorAutoLtByFirstTupleEntry  >, 
        VecOfVecSimple< usize, ArrayMapping::SnzVal >,
        GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
    )
    where   ArrayMapping:                 OracleMajorAscend + IndicesAndCoefficients,
            IterKeyMaj:             Iterator < Item = ArrayMapping::KeyMaj >,
            ArrayMapping::KeyMin:                 Clone + Hash + std::cmp::Eq + Debug, 
            ArrayMapping::KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            ArrayMapping::SnzVal:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
            ArrayMapping::ViewMajorAscend:        IntoIterator + Clone, // !!! remove clone
            ArrayMapping::ViewMajorAscendEntry:  KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderComparatorViewMajorAscendEntry:   Clone + StrictlyLess <  ArrayMapping::KeyMin >, // !!! remove clone
            // OrderComparatorViewMinorDescendEntry:   StrictlyLess <  ArrayMapping::KeyMaj >, // !!! remove clone            
            HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorViewMajorAscendEntry>>: Clone // !!!! remove this

{
    // let order_comparator_for_entries_with_minor_keys = 
    //         |x: &ArrayMapping::ViewMajorAscendEntry, y: &ArrayMapping::ViewMajorAscendEntry | 
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



/// Sames as [get_codomain_comb_inv_off_diag_pivot_block], but applies the clearing optimization.
pub fn get_codomain_comb_inv_off_diag_pivot_block_with_clearing< 'a, ArrayMapping, RingOperator, IterKeyMaj, KeyBoth, OrderComparatorViewMajorAscendEntry, > // OrderComparatorViewMinorDescendEntry >
    ( 
            array_mapping:                      &'a ArrayMapping,             
            iter_keymaj:                        IterKeyMaj,
            ring_operator:                      RingOperator,            
            order_comparator_for_keymin:        OrderComparatorViewMajorAscendEntry,
            // order_comparator_for_keymaj:        OrderComparatorViewMinorDescendEntry,            
    ) 
    -> 
    ( 
        // VecOfVec< ( usize, ArrayMapping::SnzVal ), OrderComparatorAutoLtByFirstTupleEntry  >, 
        VecOfVecSimple< usize, ArrayMapping::SnzVal >,
        GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
    )
    where   ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients<KeyMin = KeyBoth, KeyMaj=KeyBoth >,
            IterKeyMaj:                             Iterator < Item = ArrayMapping::KeyMaj >,
            ArrayMapping::KeyMin:                 Clone + Hash + std::cmp::Eq + Debug, 
            ArrayMapping::KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            ArrayMapping::SnzVal:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
            ArrayMapping::ViewMajorAscend:        IntoIterator + Clone + ParetoShortCircuit<ArrayMapping::ViewMajorAscendEntry>, // !!! remove clone
            ArrayMapping::ViewMajorAscendEntry:  KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderComparatorViewMajorAscendEntry:   Clone + StrictlyLess <  ArrayMapping::KeyMin >, // !!! remove clone
            // OrderComparatorViewMinorDescendEntry:   StrictlyLess <  ArrayMapping::KeyMaj >, // !!! remove clone            
            HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorViewMajorAscendEntry>>: Clone // !!!! remove this

{
    // let order_comparator_for_entries_with_minor_keys = 
    //         |x: &ArrayMapping::ViewMajorAscendEntry, y: &ArrayMapping::ViewMajorAscendEntry | 
    //             order_comparator_for_keymin.strictly_less( &x.key(), &y.key() );

    let order_comparator_for_entries_with_minor_keys = OrderComparatorLtByKey::new( order_comparator_for_keymin );
    // let order_comparator_for_entries_with_major_keys = OrderComparatorLtByKey::new( order_comparator_for_keymaj );
    
    let ( array_codomain_comb_inv_off_diag_pivot_block, array_matching ) = 
    codomain_comb_inv_off_diag_pivot_block_unsafecomparator_skipmatched(
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
/// More precisely, it contains a `VecOfVecSimple< usize, ArrayMapping::SnzVal >`, whose `k`th row contains the off-diagonal elements of
/// `R_{\rho \rho}`; where we replace each index `\rho_i` with the corresponding integer `i`, so that elements of the 
/// `VecOfVecSimple` are tuples of form `(usize, ArrayMapping::SnzVal)`.
/// 
/// 
/// # Design notes
/// 
/// **Why keep `ArrayMapping::ViewMajorAscend` as a generic type parameter?**  
/// 
/// **Why store the "off diagonal seed" of the codomain COMB as a `VecOfVecSimple` instead of a `VecOfVec`?**  Because we want to wrap this
/// struct in a `PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend` struct; this wrapper suppies the missing diagonal entries of,
/// the seed, thus forming the complete seed of the codomain COMB, which is a primitive and fundamental object for this 
/// computational library.  If we wanted to create an oracle for the seed that exposed entries as *references* to tuples,
/// then we would have to store the diagonal entries in memory.  However,  [`numerical experiments`](https://arxiv.org/pdf/2108.08831.pdf) 
/// have shown that the number of diagonal entries often eclipses the number of off-diagonal elements.  Thus storing 
/// the diagonal entries might incur a nontrivial memory cost.
/// 
/// **Why store the "off diagonal seed" of the codomain COMB as a `VecOfVecSimple` instead of a `Vec< Vec< (usize, ArrayMapping::SnzVal) > >`?**
/// Because `VecOfVecSimple` comes with certain guarantees about order of entries, which allows us to safely implement 
/// the `OracleMajorAscend` and `OracleMajorDescend` traits.
/// 
/// **Remark** One can always obtain a `VecOfVecFromBorrow` from a `VecOfVecSimple` via the
/// [`from_vecofvecsimple`](crate::matrices::matrix_types::vec_of_vec::VecOfVecFromBorrow::from_vecofvecsimple) method.
pub struct UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        // OrderComparatorAutoLtByKey<usize, ArrayMapping::SnzVal, (usize, ArrayMapping::SnzVal)>: StrictlyLess< (usize, ArrayMapping::SnzVal)>, // this seems extraneous but the compiler seems to want it
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    array_mapping:                      ArrayMapping,
    array_matching:                     GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
    array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:   VecOfVecSimple< usize, ArrayMapping::SnzVal >,
    // array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal: VecOfVec< ( usize, ArrayMapping::SnzVal ), OrderComparatorAutoLtByFirstTupleEntry >,
    // array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal: Vec< Vec< ( usize, ArrayMapping::SnzVal ) > >,
    ring_operator:                      RingOperator,
    order_comparator_viewmajorascendentry:      OrderComparatorViewMajorAscendEntry,
    order_comparator_viewminordescendentry:     OrderComparatorViewMinorDescendEntry,    
    // phantom_viewmajorascend:            PhantomData< ArrayMapping::ViewMajorAscend >, // required b/c otherwise the compiler complains that the type parameter `ArrayMapping::ViewMajorAscend` is unused
    // phantom_viewminordescend:           PhantomData< ViewMinorDescend >, // required b/c otherwise the compiler complains that the type parameter `ViewMinorDescend` is unused    
}


//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- CONSTRUCTOR
//  ---------------------------------------------------------------------------------------------------------


/// Generate a new U-match factorization.
pub fn new_umatchrowmajor
        < ArrayMapping, IterKeyMaj, RingOperator, OrderComparatorKeyMin, OrderComparatorKeyMaj > ( 
            array_mapping:              ArrayMapping, 
            iter_keymaj:                IterKeyMaj,
            ring_operator:              RingOperator,
            order_comparator_keymin:    OrderComparatorKeyMin,
            order_comparator_keymaj:    OrderComparatorKeyMaj,                
        ) 
    -> 
    UmatchRowMajor< 
            ArrayMapping, 
            RingOperator, 
            OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin >, 
            OrderComparatorLtByKey< ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ArrayMapping::ViewMinorDescendEntry, OrderComparatorKeyMaj >, 
        > 

where   ArrayMapping:           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients, // we require OracleMinorDescend for the function `new(..)`, but not for the `UmatchRowMajor` struct itself.
        IterKeyMaj:             Iterator < Item = ArrayMapping::KeyMaj >,
        ArrayMapping::KeyMin:                 Clone + Hash + std::cmp::Eq + Debug, 
        ArrayMapping::KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        ArrayMapping::SnzVal:                 Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        ArrayMapping::ViewMajorAscend:        Clone, // !!! remove clone
        ArrayMapping::ViewMajorAscendEntry:  KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
        OrderComparatorKeyMin:   Clone + StrictlyLess <  ArrayMapping::KeyMin >, // !!! remove clone
        OrderComparatorKeyMaj:   Clone + StrictlyLess <  ArrayMapping::KeyMaj >,
        HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin>>: Clone // !!!! remove this        
        
{
    
    let ( array_codomain_comb_inv_off_diag_pivot_block, array_matching ) : ( VecOfVecSimple<usize, ArrayMapping::SnzVal>, GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal > )
        = get_codomain_comb_inv_off_diag_pivot_block( 
                & array_mapping, 
                iter_keymaj, 
                ring_operator.clone(),
                order_comparator_keymin.clone(),
                // order_comparator_keymaj.clone(),
            );
    
    UmatchRowMajor{ 
            array_mapping, 
            array_matching, 
            array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     array_codomain_comb_inv_off_diag_pivot_block,   
            ring_operator:                              ring_operator,
            order_comparator_viewmajorascendentry:      OrderComparatorLtByKey
                                                            ::<ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin, >
                                                            ::new( order_comparator_keymin ),
            order_comparator_viewminordescendentry:     OrderComparatorLtByKey
                                                            ::<ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ArrayMapping::ViewMinorDescendEntry, OrderComparatorKeyMaj, >
                                                            ::new( order_comparator_keymaj ),                                                                
            // order_comparator_viewminordescendentry:     OrderComparatorLtByKey::new( order_comparator_keymaj ),                
            // phantom_viewmajorascend:        PhantomData,
        }
    
}


/// Same as [new_umatchrowmajor], but applies the clearning optimization.
pub fn new_umatchrowmajor_with_clearing
        < ArrayMapping, IterKeyMaj, KeyBoth, RingOperator, OrderComparatorKeyMin, OrderComparatorKeyMaj > ( 
            array_mapping:              ArrayMapping, 
            iter_keymaj:                IterKeyMaj,
            ring_operator:              RingOperator,
            order_comparator_keymin:    OrderComparatorKeyMin,
            order_comparator_keymaj:    OrderComparatorKeyMaj,                
        ) 
    -> 
    UmatchRowMajor< 
            ArrayMapping, 
            RingOperator, 
            OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin >, 
            OrderComparatorLtByKey< ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ArrayMapping::ViewMinorDescendEntry, OrderComparatorKeyMaj >, 
        > 

where   ArrayMapping:           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients< KeyMin=KeyBoth, KeyMaj=KeyBoth >, // we require OracleMinorDescend for the function `new(..)`, but not for the `UmatchRowMajor` struct itself.
        IterKeyMaj:             Iterator < Item = ArrayMapping::KeyMaj >,
        ArrayMapping::KeyMin:                 Clone + Hash + std::cmp::Eq + Debug, 
        ArrayMapping::KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        ArrayMapping::SnzVal:                 Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        ArrayMapping::ViewMajorAscend:        Clone  + ParetoShortCircuit<ArrayMapping::ViewMajorAscendEntry>, // !!! remove clone
        ArrayMapping::ViewMajorAscendEntry:  KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
        OrderComparatorKeyMin:   Clone + StrictlyLess <  ArrayMapping::KeyMin >, // !!! remove clone
        OrderComparatorKeyMaj:   Clone + StrictlyLess <  ArrayMapping::KeyMaj >,
        HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin>>: Clone // !!!! remove this        
        
{
    
    let ( array_codomain_comb_inv_off_diag_pivot_block, array_matching ) : ( VecOfVecSimple<usize, ArrayMapping::SnzVal>, GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal > )
        = get_codomain_comb_inv_off_diag_pivot_block_with_clearing( 
                & array_mapping, 
                iter_keymaj, 
                ring_operator.clone(),
                order_comparator_keymin.clone(),
                // order_comparator_keymaj.clone(),
            );
    
    UmatchRowMajor{ 
            array_mapping, 
            array_matching, 
            array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     array_codomain_comb_inv_off_diag_pivot_block,   
            ring_operator:                              ring_operator,
            order_comparator_viewmajorascendentry:      OrderComparatorLtByKey
                                                            ::<ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin, >
                                                            ::new( order_comparator_keymin ),
            order_comparator_viewminordescendentry:     OrderComparatorLtByKey
                                                            ::<ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ArrayMapping::ViewMinorDescendEntry, OrderComparatorKeyMaj, >
                                                            ::new( order_comparator_keymaj ),                                                                
            // order_comparator_viewminordescendentry:     OrderComparatorLtByKey::new( order_comparator_keymaj ),                
            // phantom_viewmajorascend:        PhantomData,
        }
    
}


//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- GENERAL IMPLEMENTATIONS (THERE ARE SPECIFIC IMPLEMENTATIONS FOR KEYMIN=KEMAJ BELOW)
//  ---------------------------------------------------------------------------------------------------------

     
impl < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >  

    UmatchRowMajor 
    < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >  
    
    where   
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        // OrderComparatorAutoLtByKey<usize, ArrayMapping::SnzVal, (usize, ArrayMapping::SnzVal)>: StrictlyLess< (usize, ArrayMapping::SnzVal)>
{
  
// }


    /// Generate a new U-match factorization (this constructor usually throws an error due to a technial issue re: type inference, so it is usually more convenient to use the function `new_umatchrowmajor`).
    pub fn new
            < IterKeyMaj, OrderComparatorKeyMin, OrderComparatorKeyMaj, ViewMinorDescendEntry > ( 
                array_mapping:              ArrayMapping, 
                iter_keymaj:                IterKeyMaj,
                ring_operator:              RingOperator,
                order_comparator_keymin:    OrderComparatorKeyMin,
                order_comparator_keymaj:    OrderComparatorKeyMaj,                
            ) 
        -> 
        UmatchRowMajor< 
                ArrayMapping, 
                RingOperator, 
                OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin >, 
                OrderComparatorLtByKey< ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ViewMinorDescendEntry, OrderComparatorKeyMaj >, 
            > 

    where   ArrayMapping:           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients, // we require OracleMinorDescend for the function `new(..)`, but not for the `UmatchRowMajor` struct itself.
            IterKeyMaj:             Iterator < Item = ArrayMapping::KeyMaj >,
            ArrayMapping::KeyMin:                 Clone + Hash + std::cmp::Eq + Debug, 
            ArrayMapping::KeyMaj:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            ArrayMapping::SnzVal:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
            ArrayMapping::ViewMajorAscend:        Clone, // !!! remove clone
            ArrayMapping::ViewMajorAscendEntry:  KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
            OrderComparatorKeyMin:   Clone + StrictlyLess <  ArrayMapping::KeyMin >, // !!! remove clone
            OrderComparatorKeyMaj:   Clone + StrictlyLess <  ArrayMapping::KeyMaj >,
            HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin>>: Clone // !!!! remove this        
            
    {
        
        let ( array_codomain_comb_inv_off_diag_pivot_block, array_matching ) : ( VecOfVecSimple<usize, ArrayMapping::SnzVal>, GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal > )
            = get_codomain_comb_inv_off_diag_pivot_block( 
                    & array_mapping, 
                    iter_keymaj, 
                    ring_operator.clone(),
                    order_comparator_keymin.clone(),
                    // order_comparator_keymaj.clone(),
                );
        
        UmatchRowMajor{ 
                array_mapping, 
                array_matching, 
                array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     array_codomain_comb_inv_off_diag_pivot_block,   
                ring_operator:                  ring_operator,
                order_comparator_viewmajorascendentry:      OrderComparatorLtByKey
                                                                ::<ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin, >
                                                                ::new( order_comparator_keymin ),
                order_comparator_viewminordescendentry:     OrderComparatorLtByKey
                                                                ::<ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ViewMinorDescendEntry, OrderComparatorKeyMaj, >
                                                                ::new( order_comparator_keymaj ),                                                                
                // order_comparator_viewminordescendentry:     OrderComparatorLtByKey::new( order_comparator_keymaj ),                
                // phantom_viewmajorascend:        PhantomData,
            }
        
    }


// //  =========================================================================================================
// //  U-MATCH REF OBJECT
// //  =========================================================================================================



// #[derive(Debug)]
// pub struct UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
//     where   
//         ArrayMapping::KeyMin:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::KeyMaj:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::ViewMajorAscend:            IntoIterator,
//         ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
// {
//     array_mapping:                                          &'a ArrayMapping,
//     array_matching:                                         &'b GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
//     array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     &'b VecOfVecSimple< usize, ArrayMapping::SnzVal >,
//     ring_operator:                                          RingOperator,
//     order_comparator_viewmajorascendentry:                                 OrderComparatorViewMajorAscendEntry,
//     order_comparator_viewminordescendentry:                                 OrderComparatorViewMinorDescendEntry,    
//     phantom_viewmajorascend:                                PhantomData< ArrayMapping::ViewMajorAscend >, // required b/c otherwise the compiler complains that the type parameter `ArrayMapping::ViewMajorAscend` is unused
// }

// //  Implement
// //  ---------------------------------------------------------------------------------------------------------

// impl < 'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  

//     UmatchRowMajorWithRefs
//     < 'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  
    
//     where   
//         ArrayMapping::KeyMin:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::KeyMaj:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
//         ArrayMapping:           OracleMajorAscend + IndicesAndCoefficients,
//         ArrayMapping::ViewMajorAscend:        IntoIterator,
//         ArrayMapping::ViewMajorAscendEntry:  KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,


// {
//     pub fn new
//         ( umatch: &'b UmatchRowMajor < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  ) 
//         ->
//         UmatchRowMajorWithRefs
//              <'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
//         where
//             RingOperator:       Clone,
//             OrderComparatorViewMajorAscendEntry: Clone,
//             OrderComparatorViewMinorDescendEntry: Clone,
//     {
//         UmatchRowMajorWithRefs{ 
//                 array_mapping:                                            umatch.array_mapping, 
//                 array_matching:                                         & umatch.array_matching,
//                 array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     & umatch.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal,
//                 ring_operator:                                            umatch.ring_operator.clone(),
//                 order_comparator_viewmajorascendentry:                                   umatch.order_comparator_viewmajorascendentry.clone(),
//                 order_comparator_viewminordescendentry:                                   umatch.order_comparator_viewminordescendentry.clone(),                
//                 phantom_viewmajorascend:                                  PhantomData,
//             }
//     }

    /// Returns a copy of the ring operator
    pub fn ring_operator( &self ) -> RingOperator 
        where
            RingOperator:   Clone,
    { self.ring_operator.clone() }

    /// Returns a copy of the order comparator for index-value pairs whose index is a `ArrayMapping::KeyMin`.
    pub fn order_comparator_viewmajorascendentry( &self ) -> OrderComparatorViewMajorAscendEntry 
        where
            OrderComparatorViewMajorAscendEntry:   Clone,
    { (self.order_comparator_viewmajorascendentry).clone() } 
    
    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `ArrayMapping::KeyMin`.
    pub fn order_comparator_viewmajorascendentry_reverse( &self ) -> OrderComparatorReverse< OrderComparatorViewMajorAscendEntry >
        where
            OrderComparatorViewMajorAscendEntry:   Clone,
    { OrderComparatorReverse::new( (self.order_comparator_viewmajorascendentry).clone() ) }    
    
    
    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `ArrayMapping::KeyMin`, wrapped in a struct that implements additional ordering traits.
    pub fn order_comparator_viewmajorascendentry_reverse_extended_functionality( &self ) -> InferTotalOrderFromStrictlyLess< OrderComparatorReverse< OrderComparatorViewMajorAscendEntry > >
        where
            OrderComparatorViewMajorAscendEntry:   Clone,
    { 
        InferTotalOrderFromStrictlyLess::new( // we have to construct an order comparator that spits out values of type Ord rather than type bool
            self.order_comparator_viewmajorascendentry_reverse() 
        )         
    }        

    /// Returns a copy of the order comparator for index-value pairs whose index is a `ArrayMapping::KeyMaj`.
    pub fn order_comparator_viewminordescendentry( &self ) -> OrderComparatorViewMinorDescendEntry 
        where
            OrderComparatorViewMinorDescendEntry:   Clone,
    { self.order_comparator_viewminordescendentry.clone() }  
    
    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `ArrayMapping::KeyMaj`.
    pub fn order_comparator_viewminordescendentry_reverse( &self ) -> OrderComparatorReverse< OrderComparatorViewMinorDescendEntry >
        where
            OrderComparatorViewMinorDescendEntry:   Clone,
    { OrderComparatorReverse::new( self.order_comparator_viewminordescendentry.clone() ) }   


    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `ArrayMapping::KeyMaj`, wrapped in a struct that implements additional ordering traits.
    pub fn order_comparator_viewminordescendentry_reverse_extended_functionality( &self ) -> InferTotalOrderFromStrictlyLess< OrderComparatorReverse< OrderComparatorViewMinorDescendEntry > >
        where
            OrderComparatorViewMinorDescendEntry:   Clone,
    { 
        InferTotalOrderFromStrictlyLess::new( // we have to construct an order comparator that spits out values of type Ord rather than type bool
            self.order_comparator_viewminordescendentry_reverse() 
        )         
    }       
    
    
   

    /// Returns the (row-major) codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn array_comb_codomain< 'a >( &'a self ) -> UmatchRowMajorCombCodomain< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  {
        UmatchRowMajorCombCodomain{ umatch: self }
    }

    /// Returns the (row-major) inverse of the codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn array_comb_codomain_inv< 'a >( &'a self ) -> UmatchRowMajorCombCodomainInv< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  {
        UmatchRowMajorCombCodomainInv{ umatch: self }
    }  
    
    /// Returns the (row-major) domain COMB, indexed by `ArrayMapping::KeyMin`
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn array_comb_domain< 'a >( &'a self ) -> UmatchRowMajorCombDomain< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  {
        UmatchRowMajorCombDomain{ umatch: self }
    }

    /// Returns the (row-major) inverse of the domain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn array_comb_domain_inv< 'a >( &'a self ) -> UmatchRowMajorCombDomainInv< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  {
        UmatchRowMajorCombDomainInv{ umatch: self }
    }      

    /// Returns a reference to the matching array of the internally stored  U-match factorization.
    pub fn array_matching_ref( &self ) -> & GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal > { & self.array_matching }

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
    pub fn array_mapping_matched_cols_only< 'a >( &'a self ) -> OnlyKeyMinInsideCollection< OracleRef< 'a, ArrayMapping >, &'a HashMap< ArrayMapping::KeyMin, usize >, >
        // where 'a: 'b,
    {
        OnlyKeyMinInsideCollection::new( OracleRef::new( &self.array_mapping ), self.array_matching.bimap_min_ref().hashmap_val_to_ord() )
    }    
    // fn array_mapping_matched_cols_only<'c>( &'c self ) -> OnlyKeyMinInsideCollection<'c,'c, ArrayMapping, ArrayMapping::ViewMajorAscend, HashMap< ArrayMapping::KeyMin, usize >, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >
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
    pub fn array_mapping_matchless_cols_only< 'a >( &'a self ) -> OnlyKeyMinOutsideCollection< OracleRef< 'a, ArrayMapping >, &'a HashMap< ArrayMapping::KeyMin, usize >, >
    {
        OnlyKeyMinOutsideCollection::new( OracleRef::new( &self.array_mapping ), & self.array_matching.bimap_min_ref().hashmap_val_to_ord() )
    }    

    /// The row submatrix of the mapping array indexed by matched row indices.
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
    pub fn array_mapping_matched_rows_only< 'a >( &'a self ) -> OnlyKeyMajInsideCollection< OracleRef< 'a, ArrayMapping >, &'a HashMap< ArrayMapping::KeyMaj, usize >, >
        // where 'a: 'b,
    {
        OnlyKeyMajInsideCollection::new( OracleRef::new( &self.array_mapping ), self.array_matching.bimap_maj_ref().hashmap_val_to_ord() )
    }    
    // fn array_mapping_matched_cols_only<'c>( &'c self ) -> OnlyKeyMinInsideCollection<'c,'c, ArrayMapping, ArrayMapping::ViewMajorAscend, HashMap< ArrayMapping::KeyMin, usize >, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >
    //     where 'a: 'c, 'b: 'c,
    // {
    //     OnlyKeyMinInsideCollection::<'c, 'c,_,_,_,_,_,_>::new( self.array_mapping, self.array_matching.bimap_min().hashmap_val_to_ord() )
    // }    

    /// The row submatrix of the mapping array indexed by unmatched row indices.
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
    pub fn array_mapping_matchless_rows_only< 'a >( &'a self ) -> OnlyKeyMajOutsideCollection< OracleRef< 'a, ArrayMapping >, &'a HashMap< ArrayMapping::KeyMaj, usize >, >
    {
        OnlyKeyMajOutsideCollection::new( OracleRef::new( &self.array_mapping ), & self.array_matching.bimap_maj_ref().hashmap_val_to_ord() )
    }    

    /// Returns a reference to the internally stored compressed representation of the inverse of the codomain COMB;
    /// this representation consists of a `VecOfVecSimple` which encodes the off-diagonal entries of the square submatrix of the 
    /// inverse of the codomain comb which is indexed by matched (i.e. pivot) indices.
    /// 
    /// Concretely, the represented matrix is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal_ref< 'a >( &'a self ) 
        -> 
        &'a VecOfVecSimple< usize, ArrayMapping::SnzVal >
        // & VecOfVec< (usize, ArrayMapping::SnzVal), OrderComparatorAutoLtByFirstTupleEntry > 
        { & self.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal }   

        
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
    /// `& self.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal_ref()` (for example, one finds errors alluding to
    /// dropped temprorary values).  This function has succeeded in sidestepping such errors in the past; please 
    /// let us know if it fails to do so successfully in future examples.
    // pub fn array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal_ref_ref( &'b self ) 
    //     -> 
    //     &'b &'b VecOfVecSimple< usize, ArrayMapping::SnzVal >
    //     // & VecOfVec< (usize, ArrayMapping::SnzVal), OrderComparatorAutoLtByFirstTupleEntry > 
    //     { & self.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal }           
    


    /// The square, block submatrix of the inverse of the codomain COMB which is indexed by matched row indices.
    /// 
    /// Concretely, this object is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal< 'a >( &'a self ) 
        -> 
        PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend<
                &'a VecOfVecSimple< usize, ArrayMapping::SnzVal >,
            >
        where
            RingOperator:   Semiring< ArrayMapping::SnzVal >,
        {   

            let prepended : PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< &'a VecOfVecSimple<usize, ArrayMapping::SnzVal> >
            = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( 
                        self.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal_ref(), 
                        RingOperator::one() 
                    );  
            
            return prepended
        }  
        
    /// The square, block submatrix of the inverse of the codomain COMB which is indexed by matched row indices.
    /// 
    /// Concretely, this object is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn array_comb_codomain_inv_matched_block_indexed_by_keymaj< 'a >( &'a self ) 
        -> 
        ReindexSquareMatrix< 
                PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend<
                        &'a VecOfVecSimple< usize, ArrayMapping::SnzVal >,
                    >,
                &'a Vec< ArrayMapping::KeyMaj >,
                &'a HashMap< ArrayMapping::KeyMaj, usize >,
                usize,
                ArrayMapping::KeyMaj,
                ArrayMapping::ViewMinorDescendEntry,
            >
        where
            RingOperator:   Semiring< ArrayMapping::SnzVal >,
            ArrayMapping:   OracleMinorDescend, // this constraint exists only to guarantee that there is a well-defined type ArrayMapping::ViewMinorDescendEntry
        {   

            let matrix_integer_indexed = self.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal();
            ReindexSquareMatrix::new(
                    matrix_integer_indexed,
                    self.array_matching.bimap_maj_ref().vec_ord_to_val(),
                    self.array_matching.bimap_maj_ref().hashmap_val_to_ord(),                    
                )
        }          


    /// Returns a matrix with rows indexed by `KeyMaj` and columns indexed by `KeyMin`
    pub fn array_comb_codomain_inv_times_mapping_array_matched_block< 'a >( &'a self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlock< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
        where 
            RingOperator: Clone,
            OrderComparatorViewMajorAscendEntry: Clone,
            OrderComparatorViewMinorDescendEntry: Clone,           
    {
        CombCodomainInvTimesMappingMatchedBlock{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }


    pub fn array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_ordmaj< 'a >( &'a self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
        where 
            RingOperator: Clone,
            OrderComparatorViewMajorAscendEntry: Clone,
            // OrderComparatorViewMinorDescendEntry: Clone,
    {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }    

    pub fn array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_keymin< 'a >( &'a self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
        where 
            RingOperator: Clone,
            OrderComparatorViewMajorAscendEntry: Clone,
            // OrderComparatorViewMinorDescendEntry: Clone,
    {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }


    

}


//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- IMPLEMENTATIONS FOR KEYMIN=KEMAJ
//  ---------------------------------------------------------------------------------------------------------


pub struct Barcode< Entry >{
    intervals:          Vec< Vec< (f64, f64)> >,
    representatives:    Option< Vec< Vec< Vec< Entry > > > >,
    bounding_chains:    Option< Vec< Vec< Vec< Entry > > > >,
}

impl < Entry >
    Barcode< Entry > {

    pub fn intervals(&self) -> & Vec< Vec< (f64, f64) > > { &self.intervals }
    pub fn representatives(&self) -> & Option< Vec< Vec< Vec< Entry > > > > { &self.representatives }  
    pub fn bounding_chains(&self) -> & Option< Vec< Vec< Vec< Entry > > > > { &self.bounding_chains }       
    /// Consumes the barcode, returning its component parts
    pub fn unwrap(self) -> 
        (   Vec< Vec< (f64, f64)> >, 
            Option< Vec< Vec< Vec< Entry > > > >, 
            Option< Vec< Vec< Vec< Entry > > > >    ) 
        { ( self.intervals, self.representatives, self.bounding_chains ) }
}

impl < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, KeyBoth, EntryBoth >  

    UmatchRowMajor 
    < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMajorAscendEntry, >  
    
    where   
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        ArrayMapping:                           OracleMajorAscend<ViewMajorAscendEntry = EntryBoth> + 
                                                OracleMinorDescend<ViewMinorDescendEntry = EntryBoth> + 
                                                IndicesAndCoefficients< KeyMin=KeyBoth, KeyMaj=KeyBoth >,
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        // OrderComparatorAutoLtByKey<usize, ArrayMapping::SnzVal, (usize, ArrayMapping::SnzVal)>: StrictlyLess< (usize, ArrayMapping::SnzVal)>
{
    
    /// The betti numbers of a chain complex.
    /// 
    /// The betti numbers of a chain complex with boundary matrix `D` can be computed as 
    /// follows, c.f. [U-match Factorization](https://arxiv.org/abs/2108.08831) (recall
    /// that `D` has rows and columns for chains in every dimension): (i) obtain
    /// a U-match factorization of `D` with matching matrix `M`, (ii) the `k`th betti number
    /// is the number of indices `i` such that `M[i,:]=0` and `M[:,i]=0`, and index `i`
    /// corrsponds to a chain of dimension `k`.
    /// 
    /// This function computes betti numbers according to this formula, assuming that `Self`
    /// is the matching matrix of some boundary matrix `D`.  Argument `I` is an iterator
    /// that runs over all the row (equivalently column) indices of `D`.  If you only need
    /// homology up through dimension `d`, you only need to include indices for chains
    /// of dimension `d` and below.  Argument `dim_fn`
    /// is a function that returns the dimension of the chain associated with each index.
    /// 
    /// *Remark* under the hood this method is identical to [unmatched_histo].  It is 
    /// included as a separate function primarily for the purpose of readability.
    pub fn betti_numbers< I, F >( &self, iter_keyboth: I, dim_fn: F ) 
            -> 
            Vec< usize > 
        where
            I: Iterator<Item=KeyBoth>, 
            F: FnMut(KeyBoth)->usize  
        {
        return self.array_matching_ref().unmatched_histo(iter_keyboth, dim_fn)
    }        


    /// The betti numbers of a chain complex.
    /// 
    /// The betti numbers of a chain complex with boundary matrix `D` can be computed as 
    /// follows, c.f. [U-match Factorization](https://arxiv.org/abs/2108.08831) (recall
    /// that `D` has rows and columns for chains in every dimension): (i) obtain
    /// a U-match factorization of `D` with matching matrix `M`, (ii) the `k`th betti number
    /// is the number of indices `i` such that `M[i,:]=0` and `M[:,i]=0`, and index `i`
    /// corrsponds to a chain of dimension `k`.
    /// 
    /// This function computes betti numbers according to this formula, assuming that `Self`
    /// is the matching matrix of some boundary matrix `D`.  Argument `I` is an iterator
    /// that runs over all the row (equivalently column) indices of `D`.  If you only need
    /// homology up through dimension `d`, you only need to include indices for chains
    /// of dimension `d` and below.  Argument `dim_fn`
    /// is a function that returns the dimension of the chain associated with each index.
    /// 
    /// *Remark* under the hood this method is identical to [unmatched_histo].  It is 
    /// included as a separate function primarily for the purpose of readability.
    pub fn barcode< I, DimFn, FilFn >( 
                &self, iter_keymaj: I, 
                mut dim_fn: DimFn, 
                mut fil_fn: FilFn, 
                return_representatives: bool,
                return_bounding_chains: bool,
            ) 
            -> 
            Barcode< ArrayMapping::ViewMinorDescendEntry > 
        where
            I: Iterator<Item=KeyBoth>, 
            DimFn: FnMut( & KeyBoth)->usize,
            FilFn: FnMut( & KeyBoth)->f64,      
            ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,     
            KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
            ArrayMapping::SnzVal:                       Clone + Debug,
            ArrayMapping::ViewMajorAscend:              IntoIterator,        
            ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, 
            ArrayMapping::ViewMinorDescend:             IntoIterator,        
            ArrayMapping::ViewMinorDescendEntry:        Clone + Debug + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,  // !!!! try to delete debug
            OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry  > + StrictlyLess< ArrayMapping::ViewMinorDescendEntry>,
            RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,        
            EntryBoth:                                  std::cmp::PartialEq
        {
        println!("we will check that bounding chains indeed bound");
        println!("see that we are checking?");

        let mut bd: Vec< Vec<(f64, f64)>>  =   Vec::new();
        let mut reps = Vec::new(); // holder for representatives
        let mut bchs = Vec::new();    // holder for bounding chains
        let matching = self.array_matching_ref();
        let mut dim; let mut birth; let mut death; let mut death_index; let mut rep; let mut bch;

        for keymaj in iter_keymaj {
            if matching.contains_keymin( & keymaj ) { continue } // in this case we ignore the key

            death_index = matching.keymaj_to_keymin( & keymaj );
            death = death_index.map_or( f64::INFINITY, |x| fil_fn( & x) );
            birth = fil_fn( & keymaj );
            
            if birth == death { continue }

            dim = dim_fn( & keymaj );            
            while bd.len() < dim + 1 { bd.push( vec![] ) } // extend bd if necessary
            if return_representatives { while reps.len() < dim + 1 { reps.push( vec![] ) } }  // reps if necessary  
            if return_bounding_chains { while bchs.len() < dim + 1 { bchs.push( vec![] ) } }  // reps if necessary                  
            bd[dim].push( (birth, death) );            

            if return_representatives {
                rep = {
                    if let Some( keymin ) = matching.keymaj_to_keymin( & keymaj ) {
                        let mut boundary_vec = self.array_comb_codomain().view_minor_descend( keymaj.clone() ).into_iter().collect_vec();
                        let scalar = self.array_matching_ref().keymaj_to_snzval(&keymaj);
                        let boundary_vec_scaled =  boundary_vec.iter().cloned().scale( scalar, self.ring_operator() ).collect_vec();
                        boundary_vec_scaled
                    } else {
                        self.array_comb_domain().view_minor_descend( keymaj.clone() ).into_iter().collect_vec()
                    }
                };
                reps[ dim ].push( rep );
            };

            if return_bounding_chains {
                bch = {
                    if let Some( keymin ) = matching.keymaj_to_keymin( & keymaj ) {
                        let bounding_chain = self.array_comb_domain().view_minor_descend( keymin ).into_iter().collect_vec();
                        let boundary_vec =  vector_matrix_multiply_minor_descend_simplified(
                                bounding_chain, 
                                OracleRef::new(self.array_mapping_ref()), 
                                self.ring_operator(), 
                                self.order_comparator_viewmajorascendentry(),
                            ).collect_vec();
                        assert_eq!( reps[dim].last(), Some(&boundary_vec) );
                        boundary_vec
                    } else {
                        Vec::with_capacity(0)
                    }
                };                
                bchs[ dim ].push( bch );
            }
        }

        // tighten up the rep list, if relevant
        reps.shrink_to_fit();
        for rep_list in reps.iter_mut() {
            rep_list.shrink_to_fit();
            for rep in rep_list.iter_mut() { rep.shrink_to_fit(); }
        }

        let representatives = match return_representatives {
            true    =>  Some( reps ),
            false   =>  None,
        };

        let bounding_chains = match return_bounding_chains {
            true    =>  Some( bchs ),
            false   =>  None,
        };        

        return Barcode{ intervals: bd, representatives, bounding_chains }
    }    

}






//  Implement clone
//  ---------------------------------------------------------------------------------------------------------

// impl < 'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >   

//     Clone for

//     UmatchRowMajorWithRefs
//         < 'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
//     where   
//         ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::ViewMajorAscend:            IntoIterator,
//         ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
//         RingOperator:               Clone,
//         OrderComparatorViewMajorAscendEntry:       Clone,
//         OrderComparatorViewMinorDescendEntry:       Clone,        
// {
//     fn clone( &self ) -> Self {
//         UmatchRowMajorWithRefs{
//             array_mapping:                                          self.array_mapping,
//             array_matching:                                         self.array_matching,
//             array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     self.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal,
//             ring_operator:                                          self.ring_operator.clone(),
//             order_comparator_viewmajorascendentry:                                 self.order_comparator_viewmajorascendentry.clone(),
//             order_comparator_viewminordescendentry:                                 self.order_comparator_viewminordescendentry.clone(),            
//             phantom_viewmajorascend:                                PhantomData, // required b/c otherwise the compiler complains that the type parameter `ArrayMapping::ViewMajorAscend` is unused
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
//                 'a, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ArrayMapping, ArrayMapping::ViewMajorAscend, RingOperator, OrderComparatorViewMajorAscendEntry,
//             > 
//     where   
//         ArrayMapping::KeyMin:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::KeyMaj:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct            
// {
//     comb_cod_compressed_ordinal:    &'a VecOfVecSimple< usize, ArrayMapping::SnzVal >,
//     matching:                       &'a GeneralizedMatchingArrayWithMajorOrdinals< ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
//     mapping:                        &'a ArrayMapping, // the mapping array
//     ring_opeartor:                  RingOperator,
//     order_comparator_viewmajorascendentry:               OrderComparatorViewMajorAscendEntry,
//     phandom_viewmajorascend:        PhantomData< ArrayMapping::ViewMajorAscend >,
// }

// impl < 'a, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ArrayMapping, ArrayMapping::ViewMajorAscend, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

//         OracleMajorAscend< 
//                 ArrayMapping::KeyMaj, Cloned< Iter< '_, (ArrayMapping::KeyMin, ArrayMapping::SnzVal) > > 
//             > for
    
//         RowMajorCombCodCompressed< 
//                 'a, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, ArrayMapping, ArrayMapping::ViewMajorAscend, RingOperator, OrderComparatorViewMajorAscendEntry,
//             >     
//     where
//         ArrayMapping::KeyMin:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::KeyMaj:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct            
//         ArrayMapping::SnzVal:                 Clone,
//         ArrayMapping:           OracleMajorAscend + IndicesAndCoefficients, 
//         RingOperator:           Semiring< ArrayMapping::SnzVal >,
//         OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::KeyMaj >,
//         ArrayMapping::ViewMajorAscend:        IntoIterator,
//         ArrayMapping::ViewMajorAscendEntry:  KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,

// {

//     fn view_major_ascend(&self, keymaj: ArrayMapping::KeyMaj) -> Cloned< Iter< '_, (ArrayMapping::KeyMin, ArrayMapping::SnzVal) > >  { 

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



// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the codomain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombCodomain< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:            IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the inverse of the codomain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombCodomainInv< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:            IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the domain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombDomain< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the inverse of the domain COMB in a lazy fashion.
#[derive(Copy, Clone)]
pub struct UmatchRowMajorCombDomainInv< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients, 
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    umatch:     &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,
}




//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: DOMAIN
//  ---------------------------------------------------------------------------------------------------------



//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    IndicesAndCoefficients for 

    UmatchRowMajorCombDomain
            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,            
            
{   type KeyMaj = ArrayMapping::KeyMin; type KeyMin = ArrayMapping::KeyMin; type SnzVal = ArrayMapping::SnzVal;  }            



//  ORACLE MAJOR ASCEND
//  ---------------------------------------------------------------------------------------------------------


//  FOR UNIT TESTS, RUN A TEXT SEARCH FOR <#CombDomainImplementOracleMajorAscend>

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    OracleMajorAscend for 

    UmatchRowMajorCombDomain< 
            'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry,
        >

    where 
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::SnzVal:                       Clone,
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        OrderComparatorViewMinorDescendEntry:       Clone + StrictlyLess<  ( ArrayMapping::KeyMaj, ArrayMapping::SnzVal ) >,        
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValNew is used (i) to create the diagonal elements of the array, and (ii) to create the 'problem vector' b in a linear problem xA = b whose solution is used in our constructor
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients
{
    type ViewMajorAscend            =   UmatchRowMajorCombDomainViewMajorAscend< 
                                                'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, 
                                            >;
    type ViewMajorAscendEntry       =   ArrayMapping::ViewMajorAscendEntry;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymin: Self::KeyMaj ) -> Self::ViewMajorAscend       

    {
        match self.umatch.array_matching.contains_keymin( &keymin ) { 
            false => { 

                return  UmatchRowMajorCombDomainViewMajorAscend{
                                iter_unwrapped:         IterTwoType::Iter1(
                                                                std::iter::once(  ArrayMapping::ViewMajorAscendEntry::new( keymin, RingOperator::one() )  ) 
                                                            ),
                                phantom_arraymapping:   PhantomData,
                            };  
            }
            true => { 

                // The matrix A from Hang et al., "Umatch factorization ...", with rows indexed by major key ordinals
                // Note: this struct only contains a reference
                let seed = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_ordmaj();

                // Struct that encodes the matching from minor keys to major key ordinals
                let array_matching_ref = self.umatch.array_matching_ref();
                let keymin_to_ordmaj = | keymin: ArrayMapping::KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
                let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );

                // obtain a row of A^{-1}
                let seed_inv_row = 
                    EchelonSolveMajorAscendWithMinorKeys::new(
                            std::iter::once( ArrayMapping::ViewMajorAscendEntry::new( keymin.clone(), RingOperator::one() ) ), // the standard unit vector supported on `keymin`
                            seed, // matrix A
                            keymin_to_ordmaj_wrapped,
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_comparator_viewmajorascendentry(),
                        );

                let mut seed_inv_row_vec = seed_inv_row.collect_vec();

                seed_inv_row_vec.shrink_to_fit();               

                // obtain a copy of R_{\rho \rho}
                let comb_codomain_inv_matched_block = self.umatch.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal();

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
                                self.umatch.order_comparator_viewmajorascendentry.clone(),                                
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
                                                            self.umatch.order_comparator_viewmajorascendentry.clone()
                                                        );


                // println!("MADE IT THROUGH END OF UmatchRowMajorCombDomain.view_major_ascend (MATCHED INDEX)");                
                return  UmatchRowMajorCombDomainViewMajorAscend{
                                iter_unwrapped:         IterTwoType::Iter2( z ),
                                phantom_arraymapping:   PhantomData,
                            };                

            }
        }
    }
}


pub struct UmatchRowMajorCombDomainViewMajorAscend< 
                    'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, 
                >

    where 
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::SnzVal:                       Clone,
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValNew is used (i) to create the diagonal elements of the array, and (ii) to create the 'problem vector' b in a linear problem xA = b whose solution is used in our constructor
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients                

{
    iter_unwrapped:     IterTwoType< 
                                Once< ArrayMapping::ViewMajorAscendEntry >,
                                MergeTwoItersByOrderComparator<
                                        Peekable< IterWrappedVec< ArrayMapping::ViewMajorAscendEntry > >, 
                                        Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection< ArrayMapping::ViewMajorAscendIntoIter, &'a HashMap<ArrayMapping::KeyMin, usize>, ArrayMapping::KeyMin, ArrayMapping::SnzVal>, 
                                                        ArrayMapping::KeyMin, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry
                                                    >, 
                                            >,
                                        OrderComparatorViewMajorAscendEntry
                                    >,
                            >,
    phantom_arraymapping:   PhantomData< ArrayMapping >
}  


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, >

    Iterator for

    UmatchRowMajorCombDomainViewMajorAscend< 
            'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, 
        >   

    where 
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::SnzVal:                       Clone,
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValNew is used (i) to create the diagonal elements of the array, and (ii) to create the 'problem vector' b in a linear problem xA = b whose solution is used in our constructor
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients                        

{
    type Item = ArrayMapping::ViewMajorAscendEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}         



//  ORACLE MINOR DESCEND
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    OracleMinorDescend for 

    UmatchRowMajorCombDomain< 
            'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry,
        >
    where
        ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,        
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::SnzVal:                       Clone,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, 
        ArrayMapping::ViewMinorDescend:             IntoIterator,        
        ArrayMapping::ViewMinorDescendEntry:        Clone + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, 
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry  >,
        OrderComparatorViewMinorDescendEntry:       Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,                
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,        
{
    type ViewMinorDescend           =   UmatchRowMajorCombDomainViewMinorDescend< 
                                                'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry,
                                            >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;
    type ViewMinorDescendEntry      =   ArrayMapping::ViewMajorAscendEntry;
    
    fn  view_minor_descend( &self, index: ArrayMapping::KeyMin ) -> Self::ViewMinorDescend {

        match self.umatch.array_matching.bimap_min_ref().ord( & index ) {

            Some( ordmaj ) => {                     
                let matched_snzval  =   self.umatch.array_matching.ordmaj_to_snzval(ordmaj);
                // let matched_keymaj = self.umatch.array_matching.ordmaj_to_keymaj(ordmaj);
                let unit_vec = OncePeekable::new( ArrayMapping::ViewMajorAscendEntry::new( index, matched_snzval ) );
                let col = 
                    EchelonSolveMinorDescendWithMinorKeys::new(
                            unit_vec,
                            self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_keymin(), // matrix A
                            IdentityFunction{},
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_comparator_viewmajorascendentry.clone(), // the constructor function `EchelonSolveMajorAscendWithMinorKeys::new(` takes care of reversing the order comparator, so we don't have to
                        );   
                let col = ChangeEntryType::new( col ); // change entries from tuples to ArrayMappingViewMajorAscendEntry
                return UmatchRowMajorCombDomainViewMinorDescend{ 
                                iter_unwrapped:     IterTwoType::Iter1( col ) 
                            }
            }
            None => {
                // a column c of D[ matched_rows, : ] 
                let col_matched_rows_only 
                    =   self.umatch.array_mapping_matched_rows_only().view_minor_descend( index.clone() )
                            .negate( self.umatch.ring_operator() ); // this multiplies the column by -1
                            // .map( | entry | (self.umatch.array_matching.bimap_maj.ord( entry.key() ), entry.val()) ); // reindex so that entries are indexed by ordinals of matched major keys
                // get an object that represents the REVERSED total order on minor keys
                let order_comparator_viewmajorascendentry_reverse_total = self.umatch.order_comparator_viewmajorascendentry_reverse_extended_functionality();
                // a linear combination v of columns of R_{\rho \rho}^{-1}; we reindex by minor key, then sort by minor key
                let mut lin_comb_r 
                    =   vector_matrix_multiply_minor_descend_simplified(
                                col_matched_rows_only, 
                                self.umatch.array_comb_codomain_inv_matched_block_indexed_by_keymaj(), 
                                self.umatch.ring_operator(), 
                                self.umatch.order_comparator_viewminordescendentry()
                            )
                            .map(   |x| //|(key,val)| 
                                    ArrayMapping::ViewMajorAscendEntry::new(
                                            self.umatch.array_matching.keymaj_to_keymin(&x.key() ).unwrap(),
                                            x.val(),
                                        )
                                )
                            .collect_vec(); // collect into a vector
                lin_comb_r.sort_by( |a,b| order_comparator_viewmajorascendentry_reverse_total.decide_cmp( a, b ) );  // sort the vector according to minor key
                // let lin_comb_r = lin_comb_r.into_iter();
                // A^{-1} v  <-- this is equal to the matched part of the column we want to construct
                let matched_part = 
                    EchelonSolveMinorDescendWithMinorKeys::new(
                            IterWrappedVec::new( lin_comb_r ),
                            self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_keymin(), // matrix A
                            IdentityFunction{},
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_comparator_viewmajorascendentry.clone(), // the constructor function `EchelonSolveMajorAscendWithMinorKeys::new(` takes care of reversing the order comparator, so we don't have to
                        );
                // change the entry type of this iterator from (KeyMaj, SnzVal) to ArrayMapping::ViewMinorDescendEntry
                let matched_part =   ChangeEntryType::new( matched_part, ); // rust infers the new entry type that we want, automatically
                let col =   MergeTwoItersByOrderComparator::new(
                                    matched_part.peekable(),
                                    OncePeekable::new( ArrayMapping::ViewMajorAscendEntry::new( index, RingOperator::one() ) ),
                                    self.umatch.order_comparator_viewmajorascendentry_reverse(), // recall that we want entries returned in *descending* order
                                );                             
                return      UmatchRowMajorCombDomainViewMinorDescend{ 
                                    iter_unwrapped:     IterTwoType::Iter2( col )
                                }                                
            }
        }
    }
}




pub struct UmatchRowMajorCombDomainViewMinorDescend< 
                    'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry,
                >
    where
        ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,        
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::SnzVal:                       Clone,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, 
        ArrayMapping::ViewMinorDescend:             IntoIterator,        
        ArrayMapping::ViewMinorDescendEntry:        Clone + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, 
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry  >,
        OrderComparatorViewMinorDescendEntry:       Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,                
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,                        
{
    iter_unwrapped:    IterTwoType<     
                                ChangeEntryType<
                                        EchelonSolveMinorDescendWithMinorKeys<
                                                OncePeekable< ArrayMapping::ViewMajorAscendEntry >, 
                                                IdentityFunction, 
                                                CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
                                                    <'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry>, 
                                                RingOperator, 
                                                OrderComparatorViewMajorAscendEntry
                                            >,
                                        ArrayMapping::ViewMajorAscendEntry,
                                        ArrayMapping::KeyMin,
                                        ArrayMapping::SnzVal,                              
                                    >,                    
                                MergeTwoItersByOrderComparator<
                                        Peekable< 
                                                ChangeEntryType<
                                                        EchelonSolveMinorDescendWithMinorKeys<
                                                                IterWrappedVec< ArrayMapping::ViewMajorAscendEntry >, 
                                                                IdentityFunction, 
                                                                CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
                                                                    <'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry>, 
                                                                RingOperator, 
                                                                OrderComparatorViewMajorAscendEntry
                                                            >, 
                                                        ArrayMapping::ViewMajorAscendEntry,
                                                        ArrayMapping::KeyMin,
                                                        ArrayMapping::SnzVal,                                    
                                                    >,
                                            >,
                                        OncePeekable<ArrayMapping::ViewMajorAscendEntry>, 
                                        OrderComparatorReverse< OrderComparatorViewMajorAscendEntry >
                                    >,        
                            >,
}

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    Iterator for 

    UmatchRowMajorCombDomainViewMinorDescend< 
            'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry,
        >

    where
        ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,        
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        ArrayMapping::SnzVal:                       Clone,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, 
        ArrayMapping::ViewMinorDescend:             IntoIterator,        
        ArrayMapping::ViewMinorDescendEntry:        Clone + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, 
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry  >,
        OrderComparatorViewMinorDescendEntry:       Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,                
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,                        
{
    type Item = ArrayMapping::ViewMajorAscendEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}            










//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: DOMAIN INVERSE
//  ---------------------------------------------------------------------------------------------------------

//  FOR UNIT TESTS, RUN A TEXT SEARCH FOR <#CombDomainInvImplementOracleMajorAscend>

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    IndicesAndCoefficients for 

    UmatchRowMajorCombDomainInv
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >   
        
    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients, 
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        
{   type KeyMin = ArrayMapping::KeyMin; type KeyMaj = ArrayMapping::KeyMin; type SnzVal = ArrayMapping::SnzVal;}        


//  ORACLE MAJOR ASCEND
//  ---------------------------------------------------------------------------------------------------------

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    OracleMajorAscend for 

    UmatchRowMajorCombDomainInv
            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    where 
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        ArrayMapping::SnzVal:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients,
{
    type ViewMajorAscend            =   UmatchRowMajorCombDomainInvViewMajorAscend
                                            < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, >;
    type ViewMajorAscendEntry       =   ArrayMapping::ViewMajorAscendEntry;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymin: ArrayMapping::KeyMin )  // this function signature looks strange b/c the major keys of the domain comb have type ArrayMapping::KeyMin!  
        -> Self::ViewMajorAscend      

    {
        match self.umatch.array_matching.keymin_to_ordmaj( &keymin ) { // this function call looks strange b/c the major keys of the domain comb have type ArrayMapping::KeyMin!  
            None => { 
                return  UmatchRowMajorCombDomainInvViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter1(
                                                            std::iter::once(  ArrayMapping::ViewMajorAscendEntry::new( keymin, RingOperator::one() )  ) 
                                                        )
                    }
            }
            Some( ordmaj ) => {  
                let scalar  =   self.umatch.array_matching.ordmaj_to_snzval( ordmaj );
                let scalar_inv      =   self.umatch.ring_operator.invert( scalar );
                // let key_to_codomain_comb    =   self.umatch.array_matching.ordmaj_to_keymaj(ordmaj);

                // this equals row `keymaj` of M_{rk}^{-1} * R^{-1}_{\rho \rho}
                let lefthand_factor_vec     =   
                        (& self.umatch.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal() )                                                    
                            .view_major_ascend( ordmaj )
                            .scale( scalar_inv, self.umatch.ring_operator.clone() )
                            .map(   |(x,y)|   // re-index the entries, so that indices align with the row indices of the mapping array
                                    (self.umatch.array_matching.ordmaj_to_keymaj(x), y) 
                                );

                let iter = vector_matrix_multiply_major_ascend_simplified ( 
                            lefthand_factor_vec,
                            OracleRef::new( &self.umatch.array_mapping ),
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_comparator_viewmajorascendentry.clone(),
                        );

                return  UmatchRowMajorCombDomainInvViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter2(  iter  )
                            }                        
            }
        }
    }
}


pub struct  UmatchRowMajorCombDomainInvViewMajorAscend
                < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, >
    where 
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        ArrayMapping::SnzVal:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients,                
{
    iter_unwrapped:     IterTwoType< 
                                Once
                                    < ArrayMapping::ViewMajorAscendEntry >,
                                LinearCombinationSimplified
                                    < ArrayMapping::ViewMajorAscendIntoIter, ArrayMapping::KeyMin, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry >,
                            >,
}

impl < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, >

    Iterator for

    UmatchRowMajorCombDomainInvViewMajorAscend
        < ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, >
        
    where 
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        ArrayMapping::SnzVal:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients,                
        
{
    type Item = ArrayMapping::ViewMajorAscendEntry;
    
    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}



//  ORACLE MINOR DESCEND
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    OracleMinorDescend for 

    UmatchRowMajorCombDomainInv
            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    where 
        ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        ArrayMapping::SnzVal:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        ArrayMapping::ViewMajorAscend:              IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:         KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping::ViewMinorDescend:             IntoIterator,        
        ArrayMapping::ViewMinorDescendEntry:        KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,        
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,        
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,        
        OrderComparatorViewMinorDescendEntry:       Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,                
{
    type ViewMinorDescend           =   Vec< ArrayMapping::ViewMajorAscendEntry >;
    type ViewMinorDescendIntoIter   =   < Self::ViewMinorDescend as IntoIterator >::IntoIter;    
    type ViewMinorDescendEntry      =   ArrayMapping::ViewMajorAscendEntry;

    fn   view_minor_descend( &self, index: Self::KeyMin ) -> Self::ViewMinorDescend {
        let col     =   self.umatch.array_mapping_matched_rows_only().view_minor_descend( index.clone() );
        let mut solution    =   
            vector_matrix_multiply_minor_descend_simplified(
                    col, 
                    self.umatch.array_comb_codomain_inv_matched_block_indexed_by_keymaj(), 
                       self.umatch.ring_operator(), 
                    self.umatch.order_comparator_viewminordescendentry(),
                )
                .map(
                    |x|
                    {
                        let ordmaj = self.umatch.array_matching.keymaj_to_ordmaj( &x.key() ).unwrap();
                        let matched_keymin = self.umatch.array_matching.ordmaj_to_keymin( ordmaj );
                        let matched_snzval = self.umatch.array_matching.ordmaj_to_snzval( ordmaj );
                        let coefficient = self.umatch.ring_operator.divide( x.val(), matched_snzval );
                        ArrayMapping::ViewMajorAscendEntry::new( matched_keymin, coefficient )
                    }
                )
                .collect_vec();
        // if `index` is not a matched minor key, then append an entry of form (index, 1)
        if ! self.umatch.array_matching.contains_keymin( & index ) {
            solution.push( ArrayMapping::ViewMajorAscendEntry::new( index, RingOperator::one() ) );
        }

        // sort the entries in the solution in descending order
        let order_comparator_viewmajorascendentry_reverse_extended = self.umatch.order_comparator_viewmajorascendentry_reverse_extended_functionality();
        solution.sort_by( |a,b| order_comparator_viewmajorascendentry_reverse_extended.decide_cmp(a, b) ); // sort the entries of the vector
        
        return solution

    }
}


//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: CODOMAIN
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    IndicesAndCoefficients for

    UmatchRowMajorCombCodomain
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    where
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,            

{   type KeyMin = ArrayMapping::KeyMaj; type KeyMaj = ArrayMapping::KeyMaj; type SnzVal = ArrayMapping::SnzVal; }        



//  IMPLEMENT ORACLE MAJOR ASCEND FOR COMB: CODOMAIN
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    OracleMajorAscend for 

    UmatchRowMajorCombCodomain
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    where 
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::SnzVal:                   Clone,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,
        ArrayMapping::ViewMajorAscend:          IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:     KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping::ViewMinorDescendEntry:    Clone + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,        
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients
{
    type ViewMajorAscend            =   UmatchRowMajorCombCodomainViewMajorAscend< ArrayMapping::ViewMinorDescendEntry >;
    type ViewMajorAscendEntry       =   ArrayMapping::ViewMinorDescendEntry;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

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
            keymaj: ArrayMapping::KeyMaj 
        ) 
        -> Self::ViewMajorAscend
    {

        // define the matrix A that will fit into an equation xA = b
        let seed_with_integer_indexed_rows = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_ordmaj();

        // define the problem vector b that will fit into an equation xA = b
        let array_mapping_matched_cols_only = self.umatch.array_mapping_matched_cols_only();
        let problem_vector = array_mapping_matched_cols_only.view_major_ascend( keymaj.clone() );

        // Struct that encodes the matching from minor keys of A to major key ordinals of A
        let array_matching_ref = self.umatch.array_matching_ref();
        let keymin_to_ordmaj = | keymin: ArrayMapping::KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
        let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

        // Solve xA = b for x
        let mut solution_vec = 
        EchelonSolveMajorAscendWithMinorKeys::new( // solution to xA = b
                problem_vector, // mabrix b
                seed_with_integer_indexed_rows, // matrix A
                keymin_to_ordmaj_wrapped,
                self.umatch.ring_operator.clone(),
                self.umatch.order_comparator_viewmajorascendentry.clone(),
            )
            .collect_vec();

        // Sort solution vector x
        // println!("REDO THIS SECTION OF CODE AFTER REFACTORING <CompareOrder>");
        let mut order_comparator_viewmajorascendentry_clone = self.umatch.order_comparator_viewmajorascendentry.clone(); // we have to clone the order comparator in order to compare order, since the `strictly_less` method takes `& mut self` as an argument -- but to get a mutable reference in the present context we would have to borrow the umatch struct as mutable
        let compare_order_for_vector_sort 
            = | x: &ArrayMapping::ViewMajorAscendEntry, y: &ArrayMapping::ViewMajorAscendEntry |-> Ordering { // we have "repackage" order_comparator_viewminordescendentry so that it returns an std::cmp::Ordering rather than a bool
            match order_comparator_viewmajorascendentry_clone.strictly_less( x, y ) {
                true => { return Ordering::Less },
                false => { 
                    match order_comparator_viewmajorascendentry_clone.strictly_less( y, x ) {
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
                        ArrayMapping::ViewMinorDescendEntry::new( 
                            self.umatch.array_matching.keymin_to_keymaj( &item.key() ).unwrap(),
                            item.val(),
                        )
                    )
                .collect_vec();

        // a function that sends (a,b) to `Less` if a < b and to `Greater` otherwise
        let mut order_comparator_viewminordescendentry = self.umatch.order_comparator_viewminordescendentry.clone(); // we have to make a clone because, as currently written, `self` is behind a mutable reference and `order_comparator_viewminordescendentry` may mutate itself when it compares two objects
        solution_vec_reindexed.sort_by( 
                |a,b| {  
                    if order_comparator_viewminordescendentry.strictly_less( a, b ) { Ordering::Less }
                    else { Ordering::Greater }
                    }
            ); 
        let solution_vec_reindexed_and_sorted = IterWrappedVec::new( solution_vec_reindexed );                            

        // Possibly append an entry with coefficient 1
        match self.umatch.array_matching_ref().contains_keymaj( & keymaj ){
            false =>    { // if keymaj is UNMATCHED then merge in a copy of ( keymaj, 1 )
                return  UmatchRowMajorCombCodomainViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter1(
                                                            std::iter::once( ArrayMapping::ViewMinorDescendEntry::new(keymaj, RingOperator::one() ) )
                                                                .chain( solution_vec_reindexed_and_sorted ) 
                                                        )
                            }

            }
            true =>     { // otherwise just return the reindexed, sorted solution vector
                return  UmatchRowMajorCombCodomainViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter2(
                                                            solution_vec_reindexed_and_sorted
                                                        )
                            }
            }
        }
    }   
}


pub struct UmatchRowMajorCombCodomainViewMajorAscend< T: Clone >
{
    iter_unwrapped:     IterTwoType< 
                                Chain<
                                        Once
                                            < T >,
                                        IterWrappedVec
                                            < T >,
                                    >,
                                IterWrappedVec
                                    < T >,
                            >
}

impl < T: Clone >

    Iterator for

    UmatchRowMajorCombCodomainViewMajorAscend< T >
{
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}






//  IMPLEMENT ORACLE MINOR DESCEND FOR COMB: CODOMAIN
//  ---------------------------------------------------------------------------------------------------------

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    OracleMinorDescend for

    UmatchRowMajorCombCodomain
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector      
        ArrayMapping::ViewMinorDescendEntry:    Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector     
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,
        OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        ArrayMapping::KeyMin:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        ArrayMapping::SnzVal:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING        
{
    type ViewMinorDescend           =   UmatchRowMajorCombCodomainViewMinorDescend
                                            < ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >;
    type ViewMinorDescendEntry      =   ArrayMapping::ViewMinorDescendEntry;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;


    fn   view_minor_descend( &self, keymaj: ArrayMapping::KeyMaj ) -> Self::ViewMinorDescend {


        match self.umatch.array_matching_ref().keymaj_to_keymin( & keymaj ) {

            // follow this branch if keymaj is matched to a minor key, keymin
            Some( keymin ) => {

                // 
                let unit_vector_keymin = std::iter::once( ArrayMapping::ViewMajorAscendEntry::new( keymin.clone(), RingOperator::one() ) );   // the standard unit vector supported on `keymaj`                                

                // the matched block of R^{-1} D, indexed by minor keys along rows and columns
                let seed = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_keymin();

                // DEPRECATED (OK TO DELETE)
                // the matching relation from rows of A to columns of A
                // let keymaj_to_keymin = self.umatch.array_matching_ref().bijection_keymaj_to_keymin();
                
                // a column c of A^{-1}
                let column_of_A_inv = EchelonSolveMinorDescendWithMinorKeys::new(
                        unit_vector_keymin,
                        seed, // matrix A
                        EvaluateFunctionFnMutWrapper::new( |x| x ),
                        self.umatch.ring_operator.clone(),
                        self.umatch.order_comparator_viewmajorascendentry.clone(), // the constructor function `EchelonSolveMajorAscendWithMinorKeys::new(` takes care of reversing the order comparator, so we don't have to
                    );
                
// ------------------ diagnostic: begin

                // // MAJOR KEYS

                // let unit_vector_keymin_copy = std::iter::once( ArrayMapping::ViewMajorAscendEntry::new( keymin.clone(), RingOperator::one() ) );   // the standard unit vector supported on `keymaj`                                                
                // let seed_copy = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_keymin();
                // let mut column_of_A_inv_copy = EchelonSolveMinorDescendWithMinorKeys::new(
                //     unit_vector_keymin_copy,
                //     seed_copy, // matrix A
                //     EvaluateFunctionFnMutWrapper::new( |x| x ),
                //     self.umatch.ring_operator.clone(),
                //     self.umatch.order_comparator_viewmajorascendentry.clone(), // the constructor function `EchelonSolveMajorAscendWithMinorKeys::new(` takes care of reversing the order comparator, so we don't have to
                // );

                // for p in 0 .. 3 { println!("column_of_A_inv_copy_majorkeys.next() = {:?}", column_of_A_inv_copy.next()) }            
                // println!("LENGTH( column_of_A_inv_copy_majorkeys ) = (next line)");
                // println!("{:?}", column_of_A_inv_copy.count() );


                // // MINOR KEYS
                // // let keymaj_to_keymin_copy = self.umatch.array_matching_ref().bijection_keymaj_to_keymin();
                // // let unit_vector_copy = std::iter::once( ArrayMapping::ViewMinorDescendEntry::new( keymaj.clone(), RingOperator::one() ) );   // the standard unit vector supported on `keymaj`
                // // let mut column_of_A_inv_copy = EchelonSolveMinorDescendWithMinorKeys::new(
                // //         unit_vector_copy,
                // //         & seed, // matrix A
                // //         keymaj_to_keymin_copy,
                // //         self.umatch.ring_operator.clone(),
                // //         self.umatch.order_comparator_viewminordescendentry.clone(), // the constructor function `EchelonSolveMajorAscendWithMinorKeys::new(` takes care of reversing the order comparator, so we don't have to
                // //     );    
                // // for p in 0 .. 3 { println!("column_of_A_inv_copy.next() = {:?}", column_of_A_inv_copy.next()) }            
                // // println!("LENGTH( column_of_A_inv ) = (next line)");
                // // println!("{:?}", column_of_A_inv_copy.count() );

// ------------------ diagnostic: end


                // the product vector D_{m k} * c  (refer to the matching identities in the Umatch facotrization paper for explanation)
                let column_of_comb  =   vector_matrix_multiply_minor_descend_simplified(
                                    column_of_A_inv,
                                                OracleRef::new( self.umatch.array_mapping_ref() ),
                                                self.umatch.ring_operator.clone(),
                                                self.umatch.order_comparator_viewminordescendentry.clone(), // the constructor handles reversing the order for us, so we don't have to
                                            );                                          
                
                // wrap the product vector in an enum
                return UmatchRowMajorCombCodomainViewMinorDescend{
                                iter_unwrapped:     IterTwoType::Iter1( column_of_comb )
                            }

            }

            // follow this branch if keymaj is not matched to a minor key
            None => {

                let unit_vector_keymaj = std::iter::once( ArrayMapping::ViewMinorDescendEntry::new( keymaj.clone(), RingOperator::one() ) );   // the standard unit vector supported on `keymaj`                
                return UmatchRowMajorCombCodomainViewMinorDescend{
                                iter_unwrapped:     IterTwoType::Iter2( unit_vector_keymaj )
                            }                                
            }
        }


        
        
    }
        
}




pub struct UmatchRowMajorCombCodomainViewMinorDescend
                < ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector      
        ArrayMapping::ViewMinorDescendEntry:    Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector     
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,
        ArrayMapping::KeyMin:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        ArrayMapping::SnzVal:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING                        
{
    iter_unwrapped:     IterTwoType<
                                LinearCombinationSimplified
                                    < ArrayMapping::ViewMinorDescendIntoIter, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorReverse< OrderComparatorViewMinorDescendEntry > >,
                                std::iter::Once< ArrayMapping::ViewMinorDescendEntry >,
                            >
}             

impl    < ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

        Iterator for 

        UmatchRowMajorCombCodomainViewMinorDescend
                < ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector      
        ArrayMapping::ViewMinorDescendEntry:    Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector     
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,
        ArrayMapping::KeyMin:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        ArrayMapping::SnzVal:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING                        
{
    type Item = ArrayMapping::ViewMinorDescendEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}





//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    IndicesAndCoefficients for

    UmatchRowMajorCombCodomainInv
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >    

    where
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,            

{   type KeyMin = ArrayMapping::KeyMaj; type KeyMaj = ArrayMapping::KeyMaj; type SnzVal = ArrayMapping::SnzVal; }  



//  IMPLEMENT ORACLE MAJOR ASCEND FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------



impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    OracleMajorAscend for 

    UmatchRowMajorCombCodomainInv 
            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry, >

    where 
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::SnzVal:                   Clone,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ( ArrayMapping::KeyMaj, ArrayMapping::SnzVal ) >,
        ArrayMapping::ViewMajorAscend:          IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:     KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
{

    type ViewMajorAscend            =   UmatchRowMajorCombCodomainInvViewMajorAscend
                                            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >;
                                        // ChangeEntryType<   
                                        //         IterTwoType< 
                                        //                 MergeTwoItersByOrderComparator<
                                        //                         Peekable<
                                        //                                 ChangeIndexSimple<
                                        //                                         LinearCombinationSimplified
                                        //                                             <Chain<Once<(usize, ArrayMapping::SnzVal)>, Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>>>, usize, ArrayMapping::SnzVal, RingOperator, OrderComparatorAutoLtByKey<usize, ArrayMapping::SnzVal, (usize, ArrayMapping::SnzVal)>>, 
                                        //                                         &'a Vec<ArrayMapping::KeyMaj>, 
                                        //                                         usize, 
                                        //                                         ArrayMapping::KeyMaj, 
                                        //                                         ArrayMapping::SnzVal
                                        //                                     >
                                        //                             >, 
                                        //                         OncePeekable<(ArrayMapping::KeyMaj, ArrayMapping::SnzVal)>, 
                                        //                         OrderComparatorViewMinorDescendEntry
                                        //                     >,                        
                                        //                 ChangeIndexSimple<
                                        //                         Chain<
                                        //                                 Once<(usize, ArrayMapping::SnzVal)>, 
                                        //                                 Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>>
                                        //                             >, 
                                        //                         &'a Vec<ArrayMapping::KeyMaj>, usize, ArrayMapping::KeyMaj, ArrayMapping::SnzVal
                                        //                     >
                                        //             >, 
                                        //         ArrayMapping::ViewMinorDescendEntry,
                                        //         ArrayMapping::KeyMaj, 
                                        //         ArrayMapping::SnzVal,    
                                        //     >;
    type ViewMajorAscendEntry       =   ArrayMapping::ViewMinorDescendEntry;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;


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
            keymaj: ArrayMapping::KeyMaj 
        ) 
        -> 
        Self::ViewMajorAscend
    {

        match self.umatch.array_matching_ref().keymaj_to_ordmaj( & keymaj ) {
            // unmatched row index
            None => {
                // define the matrix A that will fit into an equation xA = b
                let seed_with_integer_indexed_rows = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_ordmaj();
                // let seed_with_integer_indexed_rows_ref = & seed_with_integer_indexed_rows;

                // define the problem vector b that will fit into an equation xA = b
                let array_mapping_matched_cols_only = self.umatch.array_mapping_matched_cols_only();
                let problem_vector 
                    = array_mapping_matched_cols_only
                        .view_major_ascend( keymaj.clone() )
                        .negate( self.umatch.ring_operator().clone() );  // recall that there is a MINUS SIGN in the formula in Theorem 6 of Hang et al., "Umatch factorization: ..."

                // Struct that encodes the matching from minor keys of A to major key ordinals of A
                let array_matching_ref = self.umatch.array_matching_ref();
                let keymin_to_ordmaj = | keymin: ArrayMapping::KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
                let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

                // Solve xA = b for x
                let mut solution_vec = 
                EchelonSolveMajorAscendWithMinorKeys::new( // solution to xA = b
                        problem_vector, // mabrix b
                        seed_with_integer_indexed_rows, // matrix A
                        keymin_to_ordmaj_wrapped,
                        self.umatch.ring_operator.clone(),
                        self.umatch.order_comparator_viewmajorascendentry.clone(),
                    );
                let solution_vec_integer_indexed 
                    = ChangeIndexSimple::new( solution_vec, self.umatch.array_matching_ref().bimap_min_ref().hashmap_val_to_ord() );
                    
                // Multiply the solution x by matrix R_{\rho \rho}^{-1}
                let comb_codomain_inv_matched_block = self.umatch.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal();

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
                                self.umatch.order_comparator_viewminordescendentry.clone(),
                            );
                
                return  UmatchRowMajorCombCodomainInvViewMajorAscend{ iter_unwrapped: ChangeEntryType::new( IterTwoType::Iter1( merged ) ) }
            }
            // matched row index
            Some( ordmaj ) => {
                let reindexed_iter = 
                    ChangeIndexSimple::new(  
                            self.umatch.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal().view_major_ascend( ordmaj ),
                            self.umatch.array_matching_ref().bimap_maj_ref().vec_ord_to_val(),
                        );
                return  UmatchRowMajorCombCodomainInvViewMajorAscend{ iter_unwrapped: ChangeEntryType::new( IterTwoType::Iter2( reindexed_iter ) ) }
            }
        }

        // ON 2022/04/29 THIS CODE SEEMS TO BE UNCESSSARY/DEPRECATED; CONSIDER DELETING
        // // define the matrix A that will fit into an equation xA = b
        // let seed_with_integer_indexed_rows = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_ordmaj();

        // // define the problem vector b that will fit into an equation xA = b
        // let array_mapping_matched_cols_only = self.umatch.array_mapping_matched_cols_only();
        // let problem_vector = array_mapping_matched_cols_only.view_major_ascend( keymaj.clone() );

        // // Struct that encodes the matching from minor keys of A to major key ordinals of A
        // let array_matching_ref = self.umatch.array_matching_ref();
        // let keymin_to_ordmaj = | keymin: ArrayMapping::KeyMin | -> usize { array_matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
        // let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

        // // Solve xA = b for x
        // let mut solution_vec = 
        // EchelonSolveMajorAscendWithMinorKeys::new( // solution to xA = b
        //         problem_vector, // mabrix b
        //         && seed_with_integer_indexed_rows, // matrix A
        //         keymin_to_ordmaj_wrapped,
        //         self.umatch.ring_operator.clone(),
        //         self.umatch.order_comparator_viewmajorascendentry.clone(),
        //     )
        //     .collect_vec();

        // // Sort solution vector x
        // let mut order_comparator_viewmajorascendentry_clone = self.umatch.order_comparator_viewmajorascendentry.clone(); // we have to clone the order comparator in order to compare order, since the `strictly_less` method takes `& mut self` as an argument -- but to get a mutable reference in the present context we would have to borrow the umatch struct as mutable
        // let compare_order_for_vector_sort = | x: &ArrayMapping::ViewMajorAscendEntry, y: &ArrayMapping::ViewMajorAscendEntry |-> Ordering {
        //     match order_comparator_viewmajorascendentry_clone.strictly_less( x, y ) {
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





pub struct UmatchRowMajorCombCodomainInvViewMajorAscend
            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    where 
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::SnzVal:                   Clone,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        // OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ( ArrayMapping::KeyMaj, ArrayMapping::SnzVal ) >,
        ArrayMapping::ViewMajorAscend:          IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:     KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,


{
    iter_unwrapped:     ChangeEntryType<   
                                IterTwoType< 
                                        MergeTwoItersByOrderComparator<
                                                Peekable<
                                                        ChangeIndexSimple<
                                                                LinearCombinationSimplified
                                                                    <Chain<Once<(usize, ArrayMapping::SnzVal)>, Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>>>, usize, ArrayMapping::SnzVal, RingOperator, OrderComparatorAutoLtByKey<usize, ArrayMapping::SnzVal, (usize, ArrayMapping::SnzVal)>>, 
                                                                &'a Vec<ArrayMapping::KeyMaj>, 
                                                                usize, 
                                                                ArrayMapping::KeyMaj, 
                                                                ArrayMapping::SnzVal
                                                            >
                                                    >, 
                                                OncePeekable<(ArrayMapping::KeyMaj, ArrayMapping::SnzVal)>, 
                                                OrderComparatorViewMinorDescendEntry
                                            >,                        
                                        ChangeIndexSimple<
                                                Chain<
                                                        Once<(usize, ArrayMapping::SnzVal)>, 
                                                        Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>>
                                                    >, 
                                                &'a Vec<ArrayMapping::KeyMaj>, usize, ArrayMapping::KeyMaj, ArrayMapping::SnzVal
                                            >
                                    >, 
                                ArrayMapping::ViewMinorDescendEntry,
                                ArrayMapping::KeyMaj, 
                                ArrayMapping::SnzVal,    
                            >,
}

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    UmatchRowMajorCombCodomainInvViewMajorAscend
            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    where 
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::SnzVal:                   Clone,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        // OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ( ArrayMapping::KeyMaj, ArrayMapping::SnzVal ) >,
        ArrayMapping::ViewMajorAscend:          IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:     KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
{
    /// Returns the wrapped iterator, and consumes the wrapper.
    pub fn unwrap_iter( self ) ->     ChangeEntryType<   
                                    IterTwoType< 
                                            MergeTwoItersByOrderComparator<
                                                    Peekable<
                                                            ChangeIndexSimple<
                                                                    LinearCombinationSimplified
                                                                        <Chain<Once<(usize, ArrayMapping::SnzVal)>, Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>>>, usize, ArrayMapping::SnzVal, RingOperator, OrderComparatorAutoLtByKey<usize, ArrayMapping::SnzVal, (usize, ArrayMapping::SnzVal)>>, 
                                                                    &'a Vec<ArrayMapping::KeyMaj>, 
                                                                    usize, 
                                                                    ArrayMapping::KeyMaj, 
                                                                    ArrayMapping::SnzVal
                                                                >
                                                        >, 
                                                    OncePeekable<(ArrayMapping::KeyMaj, ArrayMapping::SnzVal)>, 
                                                    OrderComparatorViewMinorDescendEntry
                                                >,                        
                                            ChangeIndexSimple<
                                                    Chain<
                                                            Once<(usize, ArrayMapping::SnzVal)>, 
                                                            Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>>
                                                        >, 
                                                    &'a Vec<ArrayMapping::KeyMaj>, usize, ArrayMapping::KeyMaj, ArrayMapping::SnzVal
                                                >
                                        >, 
                                    ArrayMapping::ViewMinorDescendEntry,
                                    ArrayMapping::KeyMaj, 
                                    ArrayMapping::SnzVal,    
                                > 
    {
        self.iter_unwrapped
    }
}


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    Iterator for

    UmatchRowMajorCombCodomainInvViewMajorAscend
            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    where 
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq,
        ArrayMapping::SnzVal:                   Clone,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        // OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess<  ( ArrayMapping::KeyMaj, ArrayMapping::SnzVal ) >,
        ArrayMapping::ViewMajorAscend:          IntoIterator,        
        ArrayMapping::ViewMajorAscendEntry:     KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >,
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,            
{
    type Item = ArrayMapping::ViewMinorDescendEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}            



//  IMPLEMENT ORACLE MINOR DESCEND FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------



impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    OracleMinorDescend for

    UmatchRowMajorCombCodomainInv
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    
    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendIntoIter: Iterator,        
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal>,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess< ArrayMapping::ViewMinorDescendEntry >,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess< ArrayMapping::ViewMajorAscendEntry >,
        
{
    type ViewMinorDescend                   =   UmatchRowMajorCombCodomainInvViewMinorDescend
                                                    < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry, >;
    type ViewMinorDescendEntry              =   ArrayMapping::ViewMinorDescendEntry;
    type ViewMinorDescendIntoIter           =   Self::ViewMinorDescend;

    fn view_minor_descend
    ( 
        &self, 
        keymaj: ArrayMapping::KeyMaj
    ) 
    -> 
    Self::ViewMinorDescend
    {
        match self.umatch.array_matching.contains_keymaj( &keymaj ) {
            true => {
                // get the ordinal of the matched major key
                let ordmaj = self.umatch.array_matching.keymaj_to_ordmaj( &keymaj ).unwrap();
                // a column, c, of (R_{\rho \rho}^{-1})  [THIS NOTATION COMES FROM INNER IDENTITIES, C.F. UMATCH FACTORIZATION]
                // we have to reindex this column with minor keys, then sort it according to the order on minor keys, because our next step will be to run a triangular solve, and the triangular matrix is only triangular with respect to the order on minor keys                                
                let mut comb_codomain_inv_matched_block_col     = 
                    self.umatch.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal()
                        .view_minor_descend( ordmaj )
                        .map(   |(key,val)|   // here we convert indices to minor keys
                                ArrayMapping::ViewMajorAscendEntry::new(
                                        self.umatch.array_matching.ordmaj_to_keymin(key), 
                                        val
                                    )  
                            )
                        .collect_vec(); // collect entries into a vector to simplify sorting
                // sort the entries of the column
                let cmp_style_comparator = InferTotalOrderFromStrictlyLess::new( // we have to construct an order comparator that spits out values of type Ord rather than type bool
                        self.umatch.order_comparator_viewmajorascendentry() 
                    );    
                comb_codomain_inv_matched_block_col.sort_by( |x,y| cmp_style_comparator.decide_cmp( x, y ));
                let A = self.umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_keymin();
                // A^{-1} * c
                let echelon_solution    =   EchelonSolveMinorDescendWithMinorKeys::new(
                                                    comb_codomain_inv_matched_block_col.into_iter(),
                                                    A,
                                                    EvaluateFunctionFnMutWrapper::new(|x|  x),
                                                    self.umatch.ring_operator(),
                                                    self.umatch.order_comparator_viewmajorascendentry(),
                                                );
                let echelon_solution_minus = echelon_solution.negate( self.umatch.ring_operator() );
                let unmatched_part_of_solution = vector_matrix_multiply_minor_descend_simplified(
                        echelon_solution_minus,
                        self.umatch.array_mapping_matchless_rows_only(),
                        self.umatch.ring_operator(),
                        self.umatch.order_comparator_viewminordescendentry(),
                    );

                // this equals the variable comb_codomain_inv_matched_block_col defined above; we need a second copy, and calling the constructor a second time avoids the need to impose "Clone" on the struct
                let matched_part_of_solution     = 
                    self.umatch.array_comb_codomain_inv_matched_block_indexed_by_keymaj()
                        .view_minor_descend( keymaj );
                      
                // let   matched_part_of_solution = IterTwoType::Iter1( matched_part_of_solution );
                // let unmatched_part_of_solution = IterTwoType::Iter2( unmatched_part_of_solution );

                let solution = MergeTwoItersByOrderComparator::new(
                            matched_part_of_solution.peekable(),
                            unmatched_part_of_solution.peekable(),
                            self.umatch.order_comparator_viewminordescendentry_reverse(),
                );
                let solution = ChangeEntryType::new( IterTwoType::Iter1( solution ) );

                // wrap the output in a struct that has a simpler type signature
                return UmatchRowMajorCombCodomainInvViewMinorDescend{ iter_unwrapped: solution }
                
            }
            false => {
                let solution = IterTwoType::Iter2( 
                        std::iter::once( 
                                Self::ViewMinorDescendEntry::new( keymaj, RingOperator::one()  ) 
                            ),
                    );
                let solution = ChangeEntryType::new( solution );
                // wrap the output in a struct that has a simpler type signature                    
                return UmatchRowMajorCombCodomainInvViewMinorDescend{ iter_unwrapped: solution }
            }
        }
    }
}


pub struct UmatchRowMajorCombCodomainInvViewMinorDescend
                < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry, >
    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendIntoIter: Iterator,        
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal>,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess< ArrayMapping::ViewMinorDescendEntry >,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        // OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess< ArrayMapping::ViewMajorAscendEntry >,                
{
    iter_unwrapped: 
                    // <   EntryIter, EntryNew, Index, RingElement, >
                    ChangeEntryType<
                            IterTwoType<
                                    MergeTwoItersByOrderComparator<
                                            Peekable<
                                                    MapByTransform<
                                                            Chain<
                                                                    Once<(usize, ArrayMapping::SnzVal)>, 
                                                                    VecOfVecSimpleViewMinorDescend<'a, usize, ArrayMapping::SnzVal>
                                                                >, 
                                                            ArrayMapping::ViewMinorDescendEntry, 
                                                            ReindexEntry<(usize, ArrayMapping::SnzVal), ArrayMapping::ViewMinorDescendEntry, usize, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, &'a Vec<ArrayMapping::KeyMaj>>
                                                        >,                                                                
                                                >,
                                            Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection<
                                                                ArrayMapping::ViewMinorDescendIntoIter, 
                                                                &'a HashMap< ArrayMapping::KeyMaj, usize>, 
                                                                ArrayMapping::KeyMaj, 
                                                                ArrayMapping::SnzVal, 
                                                            >,
                                                            ArrayMapping::KeyMaj, 
                                                            ArrayMapping::SnzVal, 
                                                            RingOperator, 
                                                            OrderComparatorReverse<OrderComparatorViewMinorDescendEntry>
                                                    >,                                                                                                                             
                                                >,                                                        
                                            OrderComparatorReverse< 
                                                    OrderComparatorViewMinorDescendEntry
                                                >,                                                                
                                        >,
                                    std::iter::Once< ArrayMapping::ViewMinorDescendEntry >,
                                >,
                            ArrayMapping::ViewMinorDescendEntry,
                            ArrayMapping::KeyMaj,
                            ArrayMapping::SnzVal,
                        >,
}

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry, >

    UmatchRowMajorCombCodomainInvViewMinorDescend
                < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry, >
    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendIntoIter: Iterator,        
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal>,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess< ArrayMapping::ViewMinorDescendEntry >,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        // OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess< ArrayMapping::ViewMajorAscendEntry >,                
{
    pub fn unwrap_iter( self ) -> 
                    ChangeEntryType<
                            IterTwoType<
                                    MergeTwoItersByOrderComparator<
                                            Peekable<
                                                    MapByTransform<
                                                            Chain<
                                                                    Once<(usize, ArrayMapping::SnzVal)>, 
                                                                    VecOfVecSimpleViewMinorDescend<'a, usize, ArrayMapping::SnzVal>
                                                                >, 
                                                            ArrayMapping::ViewMinorDescendEntry, 
                                                            ReindexEntry<(usize, ArrayMapping::SnzVal), ArrayMapping::ViewMinorDescendEntry, usize, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, &'a Vec<ArrayMapping::KeyMaj>>
                                                        >,                                                                
                                                >,
                                            Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection<
                                                                ArrayMapping::ViewMinorDescendIntoIter, 
                                                                &'a HashMap< ArrayMapping::KeyMaj, usize>, 
                                                                ArrayMapping::KeyMaj, 
                                                                ArrayMapping::SnzVal, 
                                                            >,
                                                            ArrayMapping::KeyMaj, 
                                                            ArrayMapping::SnzVal, 
                                                            RingOperator, 
                                                            OrderComparatorReverse<OrderComparatorViewMinorDescendEntry>
                                                    >,                                                                                                                             
                                                >,                                                        
                                            OrderComparatorReverse< 
                                                    OrderComparatorViewMinorDescendEntry
                                                >,                                                                
                                        >,
                                    std::iter::Once< ArrayMapping::ViewMinorDescendEntry >,
                                >,
                            ArrayMapping::ViewMinorDescendEntry,
                            ArrayMapping::KeyMaj,
                            ArrayMapping::SnzVal,
                        >
    {
        self.iter_unwrapped
    }    
}



impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry, >

    Iterator for

    UmatchRowMajorCombCodomainInvViewMinorDescend
                < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry, >
    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendIntoIter: Iterator,        
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal> + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal>,
        OrderComparatorViewMinorDescendEntry:   Clone + StrictlyLess< ArrayMapping::ViewMinorDescendEntry >,
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        // OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess< ArrayMapping::ViewMajorAscendEntry >,                
{
    type Item = ArrayMapping::ViewMinorDescendEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}


//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY)
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain comb) * (the pivot block of the matching array).
#[derive(Copy, Clone)]
pub struct CombCodomainInvTimesMappingMatchedBlock< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where     
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    umatch_ref:     & 'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,  
}

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    
    CombCodomainInvTimesMappingMatchedBlock
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,    
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
        RingOperator:                           Clone,
        OrderComparatorViewMajorAscendEntry:    Clone,
        OrderComparatorViewMinorDescendEntry:   Clone,        
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlock`].
    pub fn new( umatch_ref: &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlock{ umatch_ref }
    }
}


//  INDICES AND COEFFICIENTS
//  -------------------------------------------------------------------------------------------------------------- 


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    IndicesAndCoefficients for

    &'a CombCodomainInvTimesMappingMatchedBlock
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >         

    where
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,            

{   type KeyMin = ArrayMapping::KeyMin; type KeyMaj = ArrayMapping::KeyMaj; type SnzVal = ArrayMapping::SnzVal; } 



//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT OracleMajorAscend FOR < KEYMAJ, KEYMIN >
//  -------------------------------------------------------------------------------------------------------------- 



impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    OracleMajorAscend for

    &'a CombCodomainInvTimesMappingMatchedBlock
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    where   
        ArrayMapping:               OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>> >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<std::slice::Iter<'a, (ArrayMapping::KeyMin, ArrayMapping::SnzVal)>>>,
        // OrderComparatorViewMajorAscendEntry: 'a, RingOperator: 'a, ArrayMapping::SnzVal: 'a, ArrayMapping::KeyMaj: 'a, ArrayMapping::KeyMin: 'a, ArrayMapping::ViewMajorAscend: 'a,            

{
    type ViewMajorAscend            =   LinearCombinationSimplified< 
                                                OnlyIndicesInsideCollection< ArrayMapping::ViewMajorAscendIntoIter, &'a HashMap<ArrayMapping::KeyMin, usize>, ArrayMapping::KeyMin, ArrayMapping::SnzVal>, 
                                                ArrayMapping::KeyMin, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry 
                                            >;
    type ViewMajorAscendEntry       =   ArrayMapping::ViewMajorAscendEntry;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;


    fn view_major_ascend( &self, keymaj: ArrayMapping::KeyMaj ) -> Self::ViewMajorAscend
        {
        
        // define a row vector
        // NB: we have to translate the index `keymaj` which has type `ArrayMapping::KeyMaj` into the index `ordmaj` which has type `usize`, because `self.umatch.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal` is indexed by unsigned integers 
        let ordmaj  =   self.umatch_ref.array_matching.keymaj_to_ordmaj( &keymaj ).unwrap();

        let combining_coefficients 
            = self.umatch_ref.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal().view_major_ascend( ordmaj );     

        // the matched columns of the mapping array
        let matched_cols_of_mapping_array : OnlyKeyMinInsideCollection< OracleRef< 'a, ArrayMapping >, &'a HashMap<ArrayMapping::KeyMin, usize>, >
            =   self.umatch_ref.array_mapping_matched_cols_only();
        // let matched_cols_of_mapping_array // : OnlyKeyMinInsideCollection<'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, HashMap<ArrayMapping::KeyMin, usize>, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal>
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
            =   hit_merge_by_predicate( iter_over_scaled_views, self.umatch_ref.order_comparator_viewmajorascendentry.clone() )
                    .simplify( self.umatch_ref.ring_operator.clone() );                         

        // wrap in a LinearCombinationSimplified struct
        return LinearCombinationSimplified{ linear_combination_simplified: product_vector }

    }
}



//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT OracleMinorDescend FOR < KEYMAJ, KEYMIN >
//  --------------------------------------------------------------------------------------------------------------


/// A wrapper struct for descending minor views of `CombCodomainInvTimesMappingMatchedBlock` 
pub struct  CombCodomainInvTimesMappingMatchedBlockViewMinorDescend
                < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >
    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        ArrayMapping::ViewMajorAscend:          IntoIterator, // this is required by virtually everything that involves a umatch struct
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // this is required by virtually everything that involves a umatch struct        
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal >,
        OrderComparatorViewMinorDescendEntry:                   Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,                
            { 
                view_minor_descend_of_comb_codomain_inv_times_mapping_matched_block:
                    LinearCombinationSimplified<
                            MapByTransform<
                                    Chain<
                                            Once<(usize, ArrayMapping::SnzVal)>, 
                                            VecOfVecSimpleViewMinorDescend<
                                                    'a,
                                                    usize,
                                                    ArrayMapping::SnzVal
                                                >
                                        >, 
                                    ArrayMapping::ViewMinorDescendEntry, 
                                    ReindexEntry<
                                            (usize, ArrayMapping::SnzVal), 
                                            ArrayMapping::ViewMinorDescendEntry, 
                                            usize, 
                                            ArrayMapping::KeyMaj, 
                                            ArrayMapping::SnzVal, 
                                            &'a Vec<ArrayMapping::KeyMaj>
                                        >
                                >, 
                                ArrayMapping::KeyMaj, 
                                ArrayMapping::SnzVal, 
                                RingOperator, 
                                OrderComparatorReverse<OrderComparatorViewMinorDescendEntry>,
                        >,           
            }

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >

    Iterator for 

    CombCodomainInvTimesMappingMatchedBlockViewMinorDescend
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >
    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        ArrayMapping::ViewMajorAscend:          IntoIterator, // this is required by virtually everything that involves a umatch struct
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // this is required by virtually everything that involves a umatch struct        
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal >,
        OrderComparatorViewMinorDescendEntry:                   Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >, 
{
    type Item = ArrayMapping::ViewMinorDescendEntry;

    fn next( &mut self ) -> Option< Self::Item > {
        self.view_minor_descend_of_comb_codomain_inv_times_mapping_matched_block.next()
    }
}


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    OracleMinorDescend for

    &'a CombCodomainInvTimesMappingMatchedBlock
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    where   
        ArrayMapping:                           OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                   Clone,
        ArrayMapping::ViewMinorDescend:         IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:    KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        ArrayMapping::ViewMajorAscend:          IntoIterator, // this is required by virtually everything that involves a umatch struct
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // this is required by virtually everything that involves a umatch struct        
        RingOperator:                           Clone + Semiring< ArrayMapping::SnzVal >,
        OrderComparatorViewMinorDescendEntry:    Clone + StrictlyLess<  ArrayMapping::ViewMinorDescendEntry >,
{

    type ViewMinorDescend           =   CombCodomainInvTimesMappingMatchedBlockViewMinorDescend
                                            < 'a, ArrayMapping, RingOperator, OrderComparatorViewMinorDescendEntry >;
    type ViewMinorDescendEntry      =   ArrayMapping::ViewMinorDescendEntry;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;


    fn   view_minor_descend( &self, keymin: Self::KeyMin ) -> Self::ViewMinorDescend {

        let matched_row_submatrix_of_mapping_array = self.umatch_ref.array_mapping_matched_rows_only();
        let comb_codomain_inv_matched_block = self.umatch_ref.array_comb_codomain_inv_matched_block_indexed_by_keymaj();

        let column      =   matched_row_submatrix_of_mapping_array.view_minor_descend( keymin );
        
        let column_new  = vector_matrix_multiply_minor_descend_simplified(
                    column,
                    comb_codomain_inv_matched_block,
                    self.umatch_ref.ring_operator.clone(),
                    self.umatch_ref.order_comparator_viewminordescendentry.clone(),
                );

        CombCodomainInvTimesMappingMatchedBlockViewMinorDescend{ 
                view_minor_descend_of_comb_codomain_inv_times_mapping_matched_block: column_new 
            }

        // LinearCombinationSimplified<
        //         MapByTransform<
        //                 Chain<
        //                         Once<(usize, ArrayMapping::SnzVal)>, 
        //                         VecOfVecSimpleViewMinorDescend<
        //                                 usize,
        //                                 ArrayMapping::SnzVal
        //                             >
        //                     >, 
        //                 ArrayMapping::ViewMinorDescendEntry, 
        //                 ReindexEntry<
        //                         (usize, SnzVal), 
        //                         ArrayMapping::ViewMinorDescendEntry, 
        //                         usize, 
        //                         ArrayMapping::KeyMaj, 
        //                         ArrayMapping::SnzVal, 
        //                         &Vec<ArrayMapping::KeyMaj>
        //                     >
        //             >, 
        //             ArrayMapping::KeyMaj, 
        //             ArrayMapping::SnzVal, 
        //             RingOperator, 
        //             OrderComparatorReverse<OrderComparatorViewMinorDescendEntry>,
        //     >     
        
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
pub struct CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where     
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:              IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:         KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    umatch_ref:    &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,  
}

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    
    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:              IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:         KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,      
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj`].
    fn new( umatch_ref: &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj{ umatch_ref }
    }
}


//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT ORACLE FOR ROWS INDEXED BY usize, COLUMNS BY KEYMIN
//  --------------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    IndicesAndCoefficients for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >       

    where
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,            

{   type KeyMin = ArrayMapping::KeyMin; type KeyMaj = usize; type SnzVal = ArrayMapping::SnzVal; }  



impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    OracleMajorAscend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    where   
        ArrayMapping:               OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                     Clone,
        ArrayMapping::ViewMajorAscend:            IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        RingOperator:               Clone + Semiring< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:            Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>> >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<std::slice::Iter<'a, (ArrayMapping::KeyMin, ArrayMapping::SnzVal)>>>,
        // OrderComparatorViewMajorAscendEntry: 'a, RingOperator: 'a, ArrayMapping::SnzVal: 'a, ArrayMapping::KeyMaj: 'a, ArrayMapping::KeyMin: 'a, ArrayMapping::ViewMajorAscend: 'a,            

{
    type ViewMajorAscend            =   LinearCombinationSimplified< 
                                                OnlyIndicesInsideCollection< ArrayMapping::ViewMajorAscendIntoIter, &'a HashMap<ArrayMapping::KeyMin, usize>, ArrayMapping::KeyMin, ArrayMapping::SnzVal>, 
                                                ArrayMapping::KeyMin, 
                                                ArrayMapping::SnzVal, 
                                                RingOperator, 
                                                OrderComparatorViewMajorAscendEntry 
                                            >;
    type ViewMajorAscendEntry       =   ArrayMapping::ViewMajorAscendEntry;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, ordmaj: usize ) -> Self::ViewMajorAscend 
        {
        
        // define a row vector
        let combining_coefficients 
            = self.umatch_ref.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal().view_major_ascend( ordmaj );     

        // the matched columns of the mapping array
        let matched_cols_of_mapping_array : OnlyKeyMinInsideCollection< OracleRef< 'a, ArrayMapping >, &'a HashMap<ArrayMapping::KeyMin, usize>, >
            =   self.umatch_ref.array_mapping_matched_cols_only();
        // let matched_cols_of_mapping_array // : OnlyKeyMinInsideCollection<'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, HashMap<ArrayMapping::KeyMin, usize>, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal>
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
            =   hit_merge_by_predicate( iter_over_scaled_views, self.umatch_ref.order_comparator_viewmajorascendentry.clone() )
                    .simplify( self.umatch_ref.ring_operator.clone() );

        // wrap in a LinearCombinationSimplified struct
        return LinearCombinationSimplified{ linear_combination_simplified: product_vector }

    }
}




impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    OracleMinorDescend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    where   
        ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                       Clone,
        ArrayMapping::ViewMajorAscend:              IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:         KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        ArrayMapping::ViewMinorDescend:             IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:        KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal >,
        // OrderComparatorViewMajorAscendEntry:                       Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>> >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<std::slice::Iter<'a, (ArrayMapping::KeyMin, ArrayMapping::SnzVal)>>>,
        // OrderComparatorViewMajorAscendEntry: 'a, RingOperator: 'a, ArrayMapping::SnzVal: 'a, ArrayMapping::KeyMaj: 'a, ArrayMapping::KeyMin: 'a, ArrayMapping::ViewMajorAscend: 'a,            

{
    type ViewMinorDescend           =   LinearCombinationSimplified< 
                                                Chain<
                                                        Once<(usize,ArrayMapping::SnzVal)>,
                                                        VecOfVecSimpleViewMinorDescend<'a, usize, ArrayMapping::SnzVal >,
                                                    >,
                                                usize, 
                                                ArrayMapping::SnzVal, 
                                                RingOperator, 
                                                OrderComparatorReverse< OrderComparatorAutoLtByKey< usize, ArrayMapping::SnzVal, (usize, ArrayMapping::SnzVal)> >, 
                                            >;
    type ViewMinorDescendEntry      =   (usize,ArrayMapping::SnzVal);
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;

    fn view_minor_descend( &self, keymin: ArrayMapping::KeyMin ) -> Self::ViewMinorDescend 
        {

        // the matched rows of the mapping array
        let matched_rows_of_mapping_array 
            =   self.umatch_ref.array_mapping_matched_rows_only();   

        // one of its columns
        let column = matched_rows_of_mapping_array.view_minor_descend( keymin );

        // reindex the column
        let column_reindexed = 
            ChangeIndexSimple::new( 
                    column, 
                    self.umatch_ref.array_matching.bimap_maj_ref().hashmap_val_to_ord() 
                );
        
        // multiply by R^{-1}
        let column_product = 
            vector_matrix_multiply_minor_descend_simplified(
                    column_reindexed,
                    self.umatch_ref.array_comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal(),
                    self.umatch_ref.ring_operator.clone(),
                    OrderComparatorAutoLtByKey::new(),                  
                );

        return column_product
    }
}















//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- ROWS INDEXED BY KEYMIN, COLUMNS BY KEYMIN
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain comb) * (the pivot block of the matching array).
/// 
/// This struct is almost identical to [`CombCodomainInvTimesMappingMatchedBlock`].  The only differences are that 
/// (i) the corresponding matrix oracle has rows indexed by matched minor keys of the  major keys, rather than 
/// by the matched major keys themselves, and (ii) the minor descending views iterate over entries in descending order 
/// of *minor* index, rather than major index.
#[derive(Copy, Clone)]
pub struct CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin< 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where     
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:              IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:         KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
{
    umatch_ref:    &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,  
}

impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    
    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    where   
        ArrayMapping:                               OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:              IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:         KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,      
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin`].
    pub fn new( umatch_ref: &'a UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin{ umatch_ref }
    }
}


//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT ORACLE FOR < KEYMIN, KEYMIN >
//  --------------------------------------------------------------------------------------------------------------


impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >

    IndicesAndCoefficients for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >       

    where
        ArrayMapping:                           OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::ViewMajorAscend:          IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:     KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,            

{   type KeyMin = ArrayMapping::KeyMin; type KeyMaj = ArrayMapping::KeyMin; type SnzVal = ArrayMapping::SnzVal; }  








impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    OracleMajorAscend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    where   
        ArrayMapping:               OracleMajorAscend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                     Clone,
        ArrayMapping::ViewMajorAscend:            IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        RingOperator:               Clone + Semiring< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:            Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
{
    type ViewMajorAscend            =   LinearCombinationSimplified< 
                                                OnlyIndicesInsideCollection< ArrayMapping::ViewMajorAscendIntoIter, &'a HashMap<ArrayMapping::KeyMin, usize>, ArrayMapping::KeyMin, ArrayMapping::SnzVal>, 
                                                ArrayMapping::KeyMin, 
                                                ArrayMapping::SnzVal, 
                                                RingOperator, 
                                                OrderComparatorViewMajorAscendEntry 
                                            >;
    type ViewMajorAscendEntry       =   ArrayMapping::ViewMajorAscendEntry;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymin: ArrayMapping::KeyMin ) -> Self::ViewMajorAscend 
        {
            // This implementation simply invokes a major view of CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
            // This works because the two matrices have the same set of major views -- they simply look up these views via different indices 
            let ordmaj = self.umatch_ref.array_matching.keymin_to_ordmaj( &keymin ).unwrap();
            let proxy 
                = self.umatch_ref.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_ordmaj();
            return proxy.view_major_ascend( ordmaj );
        }
}        









impl < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    OracleMinorDescend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 

    where   
        ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        ArrayMapping::SnzVal:                       Clone,
        ArrayMapping::ViewMajorAscend:              IntoIterator,
        ArrayMapping::ViewMajorAscendEntry:         Clone + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValNew< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; `Clone` is required by 
        ArrayMapping::ViewMinorDescend:             IntoIterator,
        ArrayMapping::ViewMinorDescendEntry:        KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal >,
        OrderComparatorViewMajorAscendEntry:        Clone + StrictlyLess<  ArrayMapping::ViewMajorAscendEntry >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<Iter<'a, (usize, ArrayMapping::SnzVal)>> >,
        // &'a VecOfVecSimple<usize, ArrayMapping::SnzVal>: OracleMajorAscend<usize, Cloned<std::slice::Iter<'a, (ArrayMapping::KeyMin, ArrayMapping::SnzVal)>>>,
        // OrderComparatorViewMajorAscendEntry: 'a, RingOperator: 'a, ArrayMapping::SnzVal: 'a, ArrayMapping::KeyMaj: 'a, ArrayMapping::KeyMin: 'a, ArrayMapping::ViewMajorAscend: 'a,            

{
    type ViewMinorDescend           =   IterWrappedVec< ArrayMapping::ViewMajorAscendEntry >;
    type ViewMinorDescendEntry      =   ArrayMapping::ViewMajorAscendEntry;
    type ViewMinorDescendIntoIter   =   IterWrappedVec< ArrayMapping::ViewMajorAscendEntry >;

    fn view_minor_descend( &self, keymin: ArrayMapping::KeyMin ) -> Self::ViewMinorDescend 
        {

        // a column of the matched block of R^{-1}*D; entries of the column are indexed by integers
        let col_indexed_by_ordmaj   =   self.umatch_ref.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_ordmaj().view_minor_descend( keymin.clone() );

        // reindex the entries with minor keys
        let mut col_reindexed      =   col_indexed_by_ordmaj.map( 
                                                |x| 
                                                ArrayMapping::ViewMajorAscendEntry::new(
                                                        self.umatch_ref.array_matching.ordmaj_to_keymin( x.0 ),
                                                        x.1
                                                    )
                                            )
                                            .collect_vec();
        col_reindexed.shrink_to_fit();

        // now all we have to do is sort the column according to (minor key) index!

        // repackage the order comparator for entries in major ascending views; we do this because the resulting struct has methods to produce std::comp::Ordering, which we will need to sort the column vector
        let order_decider 
                =   InferTotalOrderFromStrictlyLess::new( 
                            OrderComparatorReverse::new( // we have to reverse the order because we want our output to be sorted in *descending* order
                                    self.umatch_ref.order_comparator_viewmajorascendentry.clone() 
                                )
                        );
        // wrap the order_decider in a closure
        let order_decider_function = |x: &ArrayMapping::ViewMajorAscendEntry, y: &ArrayMapping::ViewMajorAscendEntry|  order_decider.decide_cmp(x, y);
        // use the closure to the sort the vector
        col_reindexed.sort_by( order_decider_function );


        return IterWrappedVec::new( col_reindexed )
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
// pub struct CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal< 'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
//     where     
//         ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::ViewMajorAscend:            IntoIterator,
//         ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
// {
//     umatch_ref:     UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > ,  
// }

// impl < 'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
    
//     CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal
//         < 'a, 'b, ArrayMapping, ArrayMapping::ViewMajorAscend, ArrayMapping::KeyMin, ArrayMapping::KeyMaj, ArrayMapping::SnzVal, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry > 
//     where   
//         ArrayMapping::KeyMin:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::KeyMaj:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         ArrayMapping::ViewMajorAscend:            IntoIterator,
//         ArrayMapping::ViewMajorAscendEntry:      KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >,
//         RingOperator:               Clone,
//         OrderComparatorViewMajorAscendEntry:            Clone,
// {
//     // Make a new [`CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal`].
//     fn new( umatch_ref: & UmatchRowMajor< ArrayMapping, RingOperator, OrderComparatorViewMajorAscendEntry, OrderComparatorViewMinorDescendEntry >  ) -> Self {
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

    #[test]
    fn test_umatchrowmajor_comprehensive_small() {

        // import packages
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::matrices::matrix_types::{vec_of_vec::VecOfVecSimple, oracle_ref::OracleRef};
        use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorLtByKey, OrderComparatorAutoLt};
        use crate::matrices::debug::verify_that_product_is_identity;
        use crate::matrices::operations::umatch::row_major::{UmatchRowMajor, new_umatchrowmajor};
        use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;
        use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend};
        use crate::utilities::iterators::is_sorted::IsSortedBy;
        use itertools::Itertools;

        // define the coefficient ring
        let modulus                     =   5;
        let ring_operator                   =   PrimeOrderFieldOperator::new( modulus );        

        // define the matrix we wish to factor
        let num_indices_major           =   2;
        let num_indices_minor           =   3;
        let array_mapping_data              =   VecOfVecSimple::new( 
            vec![   
                                vec![(0,1), (1,2), (2,3)], 
                                vec![              (2,1)]  
                            ] );
        let array_mapping                               =   & array_mapping_data;

        // compute the U-match factorization
        let umatch
            =   new_umatchrowmajor(
                    array_mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator.clone(), 
                    OrderComparatorAutoLt::<usize>::new(), 
                    OrderComparatorAutoLt::<usize>::new(),
                );
        
        // extract R, R^{-1}, C, C^{-1}, and M
        let array_matching = umatch.array_matching_ref();
        let array_comb_codomain = umatch.array_comb_codomain();
        let array_comb_codomain_inv = umatch.array_comb_codomain_inv();        
        let array_comb_domain = umatch.array_comb_domain();        
        let array_comb_domain_inv = umatch.array_comb_domain_inv(); 

        // get references to R, R^{-1}, C, C^{-1}, and M        
        let array_comb_codomain_ref         =   OracleRef::new( & array_comb_codomain );
        let array_comb_codomain_inv_ref         =   OracleRef::new( & array_comb_codomain_inv );
        let array_comb_domain_ref         =   OracleRef::new( & array_comb_domain );
        let array_comb_domain_inv_ref         =   OracleRef::new( & array_comb_domain_inv );   
        
        // compute some products
        let product_domain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_domain_ref, array_comb_domain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );
        let product_codomain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_ref, array_comb_codomain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );        
        let product_codomain_comb_inv_times_mapping = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_inv_ref, array_mapping, ring_operator.clone(), OrderComparatorAutoAnyType );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrixLazyMajorAscendSimplified::new( product_codomain_comb_inv_times_mapping, array_comb_domain_ref, ring_operator.clone(), OrderComparatorAutoAnyType );                


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

    }


}    


//  ---------------------------------------------------------------------
//  Unit tests
//  ---------------------------------------------------------------------

#[cfg(test)]
mod unit_tests {
    use std::{iter::Cloned, array};

    use itertools::{Product, assert_equal, Itertools};

    use crate::{matrices::{display::print_indexed_major_views, operations::{multiply::vector_matrix_multiply_major_ascend_simplified, umatch::row_major::{UmatchRowMajor}}, matrix_types::{oracle_ref::OracleRef, transpose::AntitransposeLazy}, }, utilities::partial_order::is_sorted_strictly};
    use crate::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
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
        use crate::matrices::operations::umatch::row_major::get_codomain_comb_inv_off_diag_pivot_block;
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
                                            // OrderComparatorAutoAnyType,
                                        );       
                                        
        // println!("{:?}", & codomain_comb_inv_pivot_block.0);
        // println!("{:?}", & codomain_comb_inv_pivot_block.1);        
    }


    /// This test targets the initial computation of the pivot block of the inverse of the codomain COMB.    
    #[test]
    fn test_initial_decomposition_another_example() {
        use itertools::Itertools;

        use crate::matrices::operations::umatch::row_major::get_codomain_comb_inv_off_diag_pivot_block;
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

        let array_mapping_transformed =   VecWiseTransformed::new( 
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
            // OrderComparatorAutoAnyType,
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

        use crate::matrices::operations::umatch::row_major::get_codomain_comb_inv_off_diag_pivot_block;
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

        let array_mapping_transformed  =   VecWiseTransformed::new( 
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
            // OrderComparatorAutoAnyType,
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
    //  RECOVERY OF MAJOR VIEWS OF COMBS -- COMPREHENSIVE
    //  ===================================================================================    


    /// Checks that Umatch decomposition is correct (using a small example matrix, D) in the following sense:
    /// R^{-1} * R = I
    /// C^{-1} * C = I
    /// R^{-1} * D * C = M
    /// And the rows of R, R^{-1}, C, and C^{-1} appear in strictly ascending order
    #[test]
    fn test_umatchrowmajor_comprehensive_small() {

        use crate::matrices::random_constructors::random_vec_of_vec_simple;
        use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorLtByKey, OrderComparatorAutoLt};
        use crate::matrices::debug::verify_that_product_is_identity;
        use crate::matrices::operations::umatch::row_major::{new_umatchrowmajor};
        use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;
        use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend};
        use crate::utilities::iterators::is_sorted::IsSortedBy;

        let num_indices_major           =   2;
        let num_indices_minor           =   3;
        // let approximate_density         =   0.3;
        let modulus                     =   3;
        // let allow_nonstructural_zero         =   true;

        let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
        let array_mapping_data       =   VecOfVecSimple::new( 
            vec![   
                                vec![(0,1), (1,2), (2,0)], 
                                vec![              (2,0)]  
                            ] );
        let array_mapping                               =   & array_mapping_data;

        // !!! CONSIDER CHANGING OrderComparatorAutoAnyType TO INCLUDE PHANTOM DATA WHICH FIXES THE TYPE OF THE COMPARED OBJECTS

        
        let umatch: UmatchRowMajor<&VecOfVecSimple<usize, usize>, PrimeOrderFieldOperator, OrderComparatorLtByKey<usize, usize, (usize, usize), OrderComparatorAutoLt<usize>>, OrderComparatorLtByKey<usize, usize, (usize, usize), OrderComparatorAutoLt<usize>>> 
            =   new_umatchrowmajor(
                    array_mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator.clone(), 
                    OrderComparatorAutoLt::<usize>::new(), 
                    OrderComparatorAutoLt::<usize>::new(),
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
        
        // check that the major and minor views of the inverse of the codomain COMB agree
        let comb_codomain_vec_of_vec_simple     
            =   VecOfVecSimple::from_iterable( (0..num_indices_major).map( |k| array_comb_codomain_inv.view_major_ascend(k) ) );
        for row_index in 0..num_indices_major {
            itertools::assert_equal( 
                    array_comb_codomain.view_minor_descend( row_index ),
                    (& comb_codomain_vec_of_vec_simple).view_minor_descend( row_index )
                )
        }

        // check that the major and minor views of `CombCodomainInv` agree
        let comb_codomain_inv_vec_of_vec_simple     
            =   VecOfVecSimple::from_iterable( (0..num_indices_major).map( |k| array_comb_codomain_inv.view_major_ascend(k) ) );
        for row_index in 0..num_indices_major {
            println!("PRINTING HERE: see below");
            println!("ROW INDEX = {:?}", row_index );
            println!("{:?}", array_comb_codomain_inv.view_minor_descend( row_index ).collect_vec());
            println!("{:?}", (& comb_codomain_inv_vec_of_vec_simple).view_minor_descend( row_index ).collect_vec());            
            assert_equal( 
                    array_comb_codomain_inv.view_minor_descend( row_index ).collect_vec(),
                    (& comb_codomain_inv_vec_of_vec_simple).view_minor_descend( row_index ).collect_vec()
                )
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
        use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorAutoLt};
        use crate::matrices::debug::verify_that_product_is_identity;
        use crate::matrices::operations::umatch::row_major::{new_umatchrowmajor};
        use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;
        use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend};
        use crate::utilities::iterators::is_sorted::IsSortedBy;        

        let num_indices_major           =   10;
        let num_indices_minor           =   20;
        let approximate_density           =   0.2;
        let modulus                     =   17;
        let allow_nonstructural_zero     =   true;

        let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
        let array_mapping_data       =   random_vec_of_vec_simple( num_indices_major, num_indices_minor, approximate_density, modulus, allow_nonstructural_zero );
        let array_mapping                               =   & array_mapping_data;

        let umatch 
            =   new_umatchrowmajor( 
                    array_mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator.clone(), 
                    OrderComparatorAutoLt::new(), 
                    OrderComparatorAutoLt::new(),
                );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
        let array_matching = umatch.array_matching_ref();

        let array_comb_codomain = umatch.array_comb_codomain();
        let array_comb_codomain_inv = umatch.array_comb_codomain_inv();        
        let array_comb_domain = umatch.array_comb_domain();        
        let array_comb_domain_inv = umatch.array_comb_domain_inv();  
        let array_comb_codomain_inv_times_mapping_matched_block = umatch.array_comb_codomain_inv_times_mapping_array_matched_block();  
        let array_comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin = umatch.array_comb_codomain_inv_times_mapping_array_matched_block_with_rows_indexed_by_matched_keymin();

        let array_comb_codomain_ref         =   OracleRef::new( & array_comb_codomain );
        let array_comb_codomain_inv_ref         =   OracleRef::new( & array_comb_codomain_inv );
        let array_comb_domain_ref         =   OracleRef::new( & array_comb_domain );
        let array_comb_domain_inv_ref         =   OracleRef::new( & array_comb_domain_inv );            
        let array_comb_codomain_inv_times_mapping_matched_block_ref     =   & array_comb_codomain_inv_times_mapping_matched_block;
        let array_comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref = & array_comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin;
        

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

        // check that the major and minor views of `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` agree (i.e. that, taken all together, they run over the same entries)     
        verify_viewmajorascend_compatible_with_viewminordescend(
                OracleRef::new( array_comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref ),
                umatch.array_matching_ref().bimap_min_ref().vec_ord_to_val().iter().cloned(),
                umatch.array_matching_ref().bimap_min_ref().vec_ord_to_val().iter().cloned(),                
            );
        
        // check that the major and minor views of `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` are sorted
        for keymin in umatch.array_matching_ref().bimap_min_ref().vec_ord_to_val().iter().cloned() {
            assert!(    array_comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_minor_descend( keymin ).is_sorted_by( |x, y| x.0 > y.0 )     );
        }
        for keymin in umatch.array_matching_ref().bimap_min_ref().vec_ord_to_val().iter().cloned() {
            assert!(    array_comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_major_ascend( keymin ).is_sorted_by( |x, y| x.0 < y.0 )     );
        }     



        // ----------------------------------------------------------------------------------------------------------------
        // check that the major and minor views of `CombCodomain` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                OracleRef::new( & array_comb_codomain_ref ),
                0..num_indices_major,
                0..num_indices_major,
            );

        // check that the minor views of `CombCodomain` are sorted
        for keymin in 0..num_indices_major {
            let view_minor_descend = array_comb_codomain_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend.collect_vec(), 
                    & umatch.order_comparator_viewminordescendentry_reverse() 
                ) );
        }

        // check that the major views of `CombCodomain` are sorted
        for keymin in 0..num_indices_major {
            let view_major_ascend = array_comb_codomain_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_comparator_viewmajorascendentry() 
                ) );
        }        


        // ----------------------------------------------------------------------------------------------------------------
        // check that the major and minor views of `CombCodomainInv` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                OracleRef::new( & array_comb_codomain_inv_ref ),
                0..num_indices_major,
                0..num_indices_major,
            );

        // check that the minor views of `CombCodomainInv` are sorted
        for keymin in 0..num_indices_major {
            let view_minor_descend = array_comb_codomain_inv_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend.collect_vec(), 
                    & umatch.order_comparator_viewminordescendentry_reverse() 
                ) );
        }

        // check that the major views of `CombCodomainInv` are sorted
        for keymin in 0..num_indices_major {
            let view_major_ascend = array_comb_codomain_inv_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_comparator_viewmajorascendentry() 
                ) );
        }        

        // ----------------------------------------------------------------------------------------------------------------        
        // check that the major and minor views of `CombDomain` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                OracleRef::new( & array_comb_domain_ref ),
                0..num_indices_minor,
                0..num_indices_minor,
            );

        // check that the minor views of `CombDomain` are sorted
        for keymin in 0..num_indices_minor {
            let view_minor_descend = array_comb_domain_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend.collect_vec(), 
                    & umatch.order_comparator_viewminordescendentry_reverse() 
                ) );
        }  

        // check that the major views of `CombDomain` are sorted
        for keymin in 0..num_indices_minor {
            let view_major_ascend = array_comb_domain_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_comparator_viewmajorascendentry() 
                ) );
        }          
        
        // ----------------------------------------------------------------------------------------------------------------        
        // check that the major and minor views of `CombDomainInv` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                OracleRef::new( & array_comb_domain_inv_ref ),
                0..num_indices_minor,
                0..num_indices_minor,
            );

        // check that the minor views of `CombDomainInv` are sorted
        for keymin in 0..num_indices_minor {
            let view_minor_descend = array_comb_domain_inv_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend, 
                    & umatch.order_comparator_viewminordescendentry_reverse() 
                ) );
        }   
        
        // check that the major views of `CombDomainInv` are sorted
        for keymin in 0..num_indices_minor {
            let view_major_ascend = array_comb_domain_inv_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_comparator_viewmajorascendentry() 
                ) );
        }           






        // check that `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` is upper triangular
        for keymin in umatch.array_matching_ref().bimap_min_ref().vec_ord_to_val().iter().cloned() {
            assert!(    array_comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_minor_descend( keymin.clone() ).next().unwrap().0 == keymin     );
        }        


// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:BEGIN        
        // check that columns are sorted in strictly descending order
        //  NOTE: THIS IS UNNECESSARY FOR THE COMBS, SINCE WE TEST THAT THEIR MINOR VIEWS EQUAL THOSE OF VecOfVecSimple objects, WHOSE MINOR DESCENDING VIEWS ARE *ALWAYS* STRICTLY DECREASING IN INDEX
        for keymaj in 0 .. num_indices_minor { 
            assert!(    array_mapping.view_minor_descend( keymaj.clone()             ).is_sorted_by( |x, y| x.0 > y.0 )     );
            assert!(    array_comb_codomain.view_minor_descend( keymaj.clone()       ).is_sorted_by( |x, y| x.0 > y.0 )     );
            // assert!(    array_comb_codomain_inv.view_minor_descend( keymaj.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            // assert!(    array_comb_domain.view_minor_descend( keymaj.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            // assert!(    array_comb_domain_inv.view_minor_descend( keymaj.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }          
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END  

        
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG: BEGIN   
        // check that the major and minor views of the inverse of the codomain COMB agree
        let comb_codomain_vec_of_vec_simple     
            =   VecOfVecSimple::from_iterable( (0..num_indices_major).map( |k| array_comb_codomain.view_major_ascend(k) ) );
        for keymaj in 0..num_indices_major {
            println!("VIEW MAJOR DESCEND IS STARTING FOR THIS ROUND: keymaj = {:?}", keymaj);
            println!("VIEW MAJOR DESCEND LAZY CONSTRUCTION: {:?}", array_comb_codomain.view_minor_descend( keymaj ).collect_vec());
            println!("VIEW MAJOR DESCEND FROM MAJOR VIEW: {:?}", (& comb_codomain_vec_of_vec_simple).view_minor_descend( keymaj ).collect_vec());            
            println!("VIEW MAJOR DESCEND IS FINISHED FOR THIS ROUND");
            itertools::assert_equal( 
                    array_comb_codomain.view_minor_descend( keymaj ),
                    (& comb_codomain_vec_of_vec_simple).view_minor_descend( keymaj )
                )
        }      
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END      
        
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

        use crate::matrices::operations::umatch::row_major::{new_umatchrowmajor, get_codomain_comb_inv_off_diag_pivot_block, UmatchRowMajor, CombCodomainInvTimesMappingMatchedBlock};
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
        
        let order_comparator_viewmajorascendentry                =   OrderComparatorAutoAnyType;     
        let order_comparator_viewminordescendentry                =   OrderComparatorAutoAnyType;                
        let ring_operator       =   PrimeOrderFieldOperator::new( 13 );

        let array_mapping_ref   =   & array_mapping;
        let umatch  
            =   new_umatchrowmajor( 
                        array_mapping_ref, 
                        (0..3).rev(), 
                        ring_operator.clone(),
                        order_comparator_viewmajorascendentry.clone(),
                        order_comparator_viewminordescendentry.clone(),                        
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

        use crate::matrices::operations::umatch::row_major::{new_umatchrowmajor, get_codomain_comb_inv_off_diag_pivot_block, UmatchRowMajor, CombCodomainInvTimesMappingMatchedBlock};
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
        let order_comparator_viewmajorascendentry = OrderComparatorAutoAnyType;
        let order_comparator_viewminordescendentry = OrderComparatorAutoAnyType;        

        let umatch_root = 
                new_umatchrowmajor( 
                        array_mapping_ref, 
                        (0 .. num_rows).rev(), 
                        ring_operator.clone(),
                        order_comparator_viewmajorascendentry.clone(),
                        order_comparator_viewminordescendentry.clone(),
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

        use crate::matrices::operations::umatch::row_major::{new_umatchrowmajor, get_codomain_comb_inv_off_diag_pivot_block, UmatchRowMajor, CombCodomainInvTimesMappingMatchedBlock};
        use crate::matrices::matrix_types::vec_of_vec::{VecOfVecSimple, VecOfVecSimpleViewMinorDescend};
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::rings::operator_structs::ring_native::RingNative;
        use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorAutoLtByKey};   
        use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;     

        let num_rows = 10; let num_cols = 10; let modulus = 7;
        let array_mapping = random_m_by_n_matrix(num_rows, num_cols, modulus);
        let array_mapping_ref = & array_mapping;
        let ring_operator = PrimeOrderFieldOperator::new( modulus );
        let order_comparator_viewmajorascendentry = OrderComparatorAutoAnyType;
        let order_comparator_viewminordescendentry = OrderComparatorAutoAnyType;        

        let umatch_root = 
                new_umatchrowmajor( 
                        array_mapping_ref, 
                        (0 .. num_rows).rev(), 
                        ring_operator.clone(),
                        order_comparator_viewmajorascendentry.clone(),
                        order_comparator_viewminordescendentry.clone(),                        
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















