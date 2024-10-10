//! Tools used internally to compute a U-match factorization



use super::*;



use crate::algebra::vectors::entries::{KeyValSet, KeyValNew};

use crate::algebra::matrices::types::matching::{GeneralizedMatchingArrayWithMajorOrdinals};


use crate::algebra::matrices::types::vec_of_vec::sorted::{VecOfVec, VecOfVecMatrixColumnReverse};
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients};


use crate::utilities::iterators::general::{SkipUntil, IterWrappedVec};
use crate::utilities::iterators::merge::hit::{HitMerge, hit_bulk_insert, hit_merge_by_predicate};
use crate::algebra::vectors::entries::{KeyValGet};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};

use crate::utilities::order::{JudgePartialOrder, OrderOperatorByKey, OrderOperatorByKeyCutsom, ReverseOrder};


use crate::algebra::vectors::operations::{Scale, VectorOperations, Simplify, OnlyIndicesInsideCollection, ChangeIndexSimple, LinearCombinationSimplified};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::{Peekable, Once, Chain, Rev};




use crate::algebra::matrices::operations::multiply::{vector_matrix_multiply_minor_descend_simplified};
use crate::algebra::matrices::operations::transform_vector_wise::{OnlyKeyMinInsideCollection};




//  HELPER FUNCTION 
//  ----------------------------------------------------------------------------------

/// Helper function for `try_writing_next_viewmaj_of_codomain_comb_to_buffer`.  That function
/// involves multiplying a row of a codomain COMB by a mapping array; sometimes entries
/// at the beginning of the resulting row vector can be deleted without penalty; the 
/// current function performs this deletion on a (scalar multiple) of an individual row of the 
/// mapping array.
fn codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array<
        Mapping,
        RingOperator,
        OrderOperatorRowEntries,
    > 
    (
        codomain_comb_inv_entry:            ( usize, Mapping::Coefficient ),
        scale_factor:                       Mapping::Coefficient,        
        truncation_limit:                   & Mapping::EntryMajor,
        mapping:                            & Mapping,     
        matching:                           & GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,        
        ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
        order_operator_major:               OrderOperatorRowEntries,            
    )
    ->
    Peekable<
            Scale < 
                    Mapping::ViewMajorAscendIntoIter,  // a major view of the mapping array
                    Mapping::ColIndex,
                    RingOperator, 
                    Mapping::Coefficient,
                >, 
        >
    where 
        Mapping:                 ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                 Clone + Hash + PartialEq + std::cmp::Eq,
        Mapping::RowIndex:                 Clone + Hash + PartialEq + std::cmp::Eq,         
        Mapping::Coefficient:                 Clone,               
        RingOperator:           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder <  Mapping::EntryMajor >,  
        Mapping::ViewMajorAscend:        IntoIterator,        
        Mapping::EntryMajor:  KeyValSet < Mapping::ColIndex, Mapping::Coefficient >,  
{

    let ( codomain_comb_inv_entry_index, codomain_comb_inv_entry_coeff ) = codomain_comb_inv_entry;

    // we multiply the entry coefficient `codomain_comb_inv_entry.val()` by the `scale_factor`, 
    let scale_factor_aggregate = ring_operator.multiply( scale_factor, codomain_comb_inv_entry_coeff ); 
    
    // then we multiply `codomain_comb_inv_entry.val() * scale_factor_aggregate` by the view of the **mapping array**
    // (not matching array) that is indexed by `codomain_comb_inv_entry.key()`
    
    mapping // <-- note: mapping, **not matching**
                    .view_major_ascend( 
                            // this picks out the row index corresponding to the entry; this index might *not* be an integer, though `codomain_comb_inv_entry.key()` *is* an integer
                            matching.ord_to_keymaj( codomain_comb_inv_entry_index )  
                        )
                    .into_iter()
                    .scale( scale_factor_aggregate, ring_operator )
                    // making the iterator peekable will help us prune some entries
                    .peekable() 
                    // we discard the entries of this iterator that have index <= `leading_entry_to_eliminate.key()` (because those entries would eventually be cancelled out, anyway -- this saves possibly a lot of heap insert + pop operations)
                    // in doing so, we must take care not to throw out too many entries by accident; that is why peekable is used                    
                    .skip_until( |x| {
                        let this = &order_operator_major; this.judge_partial_cmp(truncation_limit, x) == Some(Ordering::Less) }  ) 
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
/// - the `mapping` field stores matrix `M` in row-major form
/// - the portions of the matching array and inverse-codomain-COMB are stored in the fields
/// `matching` and `array_codomain_comb_inv_off_diag`.
/// - the field `order_operator_major` correctly determines which minor key precedes which 
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
        Mapping,
        RingOperator,
        OrderOperatorRowEntries,
    > 
    (
        codomain_comb_inv_off_diag_view_buffer:      &mut Vec< ( usize, Mapping::Coefficient ) >,          
        entries_to_elim_simplified_heap:    & mut Simplify <
                                                        HitMerge < 
                                                                Peekable<
                                                                        // the HitMerge struct merges a collection of iterators into a single iterator; the type of those iterators is given below (they are essentially scaled vectors)
                                                                        Scale < 
                                                                                Mapping::ViewMajorAscendIntoIter,  // a major view of the mapping array
                                                                                Mapping::ColIndex,
                                                                                RingOperator, // the ring_operator operator
                                                                                Mapping::Coefficient,
                                                                            >, 
                                                                    >,
                                                                // the thing that declares whether one major key comes before of after another    
                                                                OrderOperatorRowEntries                                                                
                                                            >,
                                                        Mapping::ColIndex,
                                                        RingOperator,
                                                        Mapping::Coefficient,
                                                    >,          
        mapping:                      & Mapping,    
        matching:                     & GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
        array_codomain_comb_inv_off_diag:   & Vec< Vec< (usize, Mapping::Coefficient) > >,  // the off-diagonal part of the (block indexed by pivot indices of) inverse of the codomain COMB       
        ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
        order_operator_major:             OrderOperatorRowEntries,            
    ) 
    ->

    Option< Mapping::EntryMajor >

    where 
        Mapping:                 ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                 Clone + Hash + PartialEq + std::cmp::Eq + Debug, // !!! remove debug
        Mapping::RowIndex:                 Clone + Hash + PartialEq + std::cmp::Eq + Debug,      // !!! remove debug   
        Mapping::Coefficient:                 Clone + Debug,                // remove debug
        RingOperator:           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:   Clone + JudgePartialOrder <  Mapping::EntryMajor >,  
        Mapping::ViewMajorAscend:        IntoIterator,        // !!!!! REMOVE THE DEBUG + CLONE
        Mapping::EntryMajor:  KeyValSet < Mapping::ColIndex, Mapping::Coefficient > + Debug,       // !!!!! REMOVE THE DEBUG + CLONE REQUIREMENT WHEN DONE DEBUGGING!
        // HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorRowEntries>: Clone, // !!! remove this later
{
    // println!("initiating construction of next row");  // !!! delete later
    // println!("matching {:?}", matching);  // !!! delete later    

    // we use this while loop because a for-loop results in a borrowing conflict (the conflict arrises because the for loop 
    // iteratos over `entries_to_elim_simplified_heap`, but we modify `entries_to_elim_simplified_heap` within the for-loop)
    while let Some( leading_entry_to_eliminate ) = entries_to_elim_simplified_heap.next() {

        // if entries_to_elim_simplified_heap.unsimplified.len() > 10 { println!("styx[[{:?}]]", entries_to_elim_simplified_heap.unsimplified.len()) }

        // println!("WHILE LOOP CHECKPOINT: leading_entry_to_eliminate: {:?}", &leading_entry_to_eliminate); // !!! REMOVE LATER
        // println!("WHILE LOOP CHECKPOINT: entries_to_elim_simplified_heap: {:?}", entries_to_elim_simplified_heap.clone().collect_vec() ); // !!! REMOVE LATER        
    
        // IF THE MINOR INDEX OF THE ENTRY IS MATCHED, THEN WE CAN ELIMINATE
        if let Some( ordmaj_matched )       =   matching.keymin_to_ord( & leading_entry_to_eliminate.key() ) {

            let scale_factor                =   ring_operator.negate(
                                                                ring_operator.divide(
                                                                        leading_entry_to_eliminate.val(),                                                                            
                                                                        matching.ord_to_snzval( ordmaj_matched ),
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
                                        ( ordmaj_matched, RingOperator::one() ) 
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
                                //                 & leading_entry_to_eliminate, // truncation_limit:                   Mapping::EntryMajor,
                                //                 mapping, //                     & Mapping,     
                                //                 matching, //                    & GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,        
                                //                 ring_operator.clone(), //                    RingOperator, // the operator for the coefficient ring_operator    
                                //                 order_operator_major.clone() //               OrderOperatorRowEntries,                                      
                                //             )
                                //             .collect_vec()
                                //     );

                                codomain_comb_entry_to_scaled_truncated_viewmaj_of_mapping_array(
                                        codomain_comb_inv_entry,
                                        scale_factor.clone(),
                                        & leading_entry_to_eliminate, // truncation_limit:                   Mapping::EntryMajor,
                                        mapping, //                     & Mapping,     
                                        matching, //                    & GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,        
                                        ring_operator.clone(), //                    RingOperator, // the operator for the coefficient ring_operator    
                                        order_operator_major.clone() //               OrderOperatorRowEntries,                                      
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
                                .cloned() // note: we have to do this b/c `codomain_comb_inv_off_diag_view_buffer` is a `VecOfVec`, not a `VecOfVec`
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
    None    
}





//  FUNCTION(S) TO EXTRACT Ripiv (the pivot portion of R^{-1}[pivot_indices, pivot_indices])
//  ------------------------------------------------------------------------------------------

/// Returns the block of the codomain COMB indexed by pivot indices, and a representation of the matching matrix.
/// 
/// For details on this factorization, see [this preprint](https://arxiv.org/pdf/2108.08831.pdf).
/// 
/// This U-match factorization is computed via the standard "cohomology algorithm."
/// 
/// **It is important that `order_operator_major` compares order of two entries based *solely* on
/// the value of the associated indices (not the associated coefficients).  This is easy
/// to get wrong, and it's hard to debug, so we are keeping the function private for now**
pub fn codomain_comb_inv_off_diag_pivot_block_unsafecomparator< Mapping, RingOperator, IterRowIndex, OrderOperatorOfEntries >
    ( 
            mapping:                      &Mapping,  // use of a reference poses no problem for future calls to this function, since output has no references        
            iter_keymaj:                        IterRowIndex,
            ring_operator:                      RingOperator,            
            order_operator_of_entries:        OrderOperatorOfEntries,
            // mut order_operator_of_keymaj:     OrderOperatorColEntries,
    ) 
    -> 
    ( 
        VecOfVec< usize, Mapping::Coefficient >,          
        GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
    )
    where   Mapping:                     ViewRowAscend + IndicesAndCoefficients,
            IterRowIndex:                 Iterator < Item = Mapping::RowIndex >,
            Mapping::ColIndex:                     Clone + Hash + std::cmp::Eq + Debug, 
            Mapping::RowIndex:                     Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            Mapping::Coefficient:                     Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            Mapping::ViewMajorAscend:            IntoIterator + Clone, // !!! remove clone
            Mapping::EntryMajor:      KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderOperatorOfEntries:   Clone + JudgePartialOrder <  Mapping::EntryMajor >, // !!! remove clone
            // OrderOperatorColEntries:       JudgePartialOrder<  Mapping::RowIndex >,
            HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorOfEntries>: Clone // !!!! remove this

{
    
    // Initialize some objects
    let mut entries_to_elim_simplified_heap    =   HitMerge::new( order_operator_of_entries.clone() )
                                                                            .simplify( ring_operator.clone() ); // an empty (simplified) merger of iterators
    let mut matching          =   GeneralizedMatchingArrayWithMajorOrdinals::new(); // an all-zero generalized matching array
    let mut array_codomain_comb_inv_off_diag: Vec< Vec< ( usize, Mapping::Coefficient ) > >  =   Vec::new();    
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
        //     match order_operator_of_keymaj.judge_lt( prior_keyaj, &keymaj ) {
        //         true => {  },
        //         false => { panic!("Umatch factorization cannot be completed because the sequence of major indices is not strictly ascending") }
        //     }               
        // }
        // prior_keymaj_opt.replace( keymaj.clone() );

        // clear the collection of entries to eliminate
        entries_to_elim_simplified_heap.unsimplified.clear();        

        // insert the sequence of entries in row `keymaj`
        entries_to_elim_simplified_heap.unsimplified.insert_one_iter(
                mapping.view_major_ascend( keymaj.clone() )
                    .into_iter()
                    .scale( RingOperator::one(), ring_operator.clone() ) // !!! might be more efficient to use a two-type iter than to perform this multiplication, which only serves to make the iterator compatible with the HitMerge struct
                    .peekable()
        );

        codomain_comb_inv_off_diag_view_buffer.clear();
        let leading_entry_uneliminable_opt =
        try_writing_next_viewmaj_of_codomain_comb_to_buffer(
                & mut codomain_comb_inv_off_diag_view_buffer,
                & mut entries_to_elim_simplified_heap,        
                  mapping,
                & matching,
                & array_codomain_comb_inv_off_diag,
                ring_operator.clone(),
                order_operator_of_entries.clone(),
            );

        match leading_entry_uneliminable_opt {
            // in this case we do not have a pivot
            Some( leading_entry_uneliminable )    => { 
                let pivot_keymin        =   leading_entry_uneliminable.key();
                let pivot_coeff         =   leading_entry_uneliminable.val();
                // println!("leading_entry_uneliminable {:?}", &leading_entry_uneliminable);              
                // update the matching array
                matching.push( keymaj, pivot_keymin, pivot_coeff ); 
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
                
                // update the vector that stores the (off diagonal entries of the) "pivot" part of the inverse of the codomain COMB
                array_codomain_comb_inv_off_diag.push( vec_to_push );
            }, 

            // in this case we don't have a pivot entry, so there is nothing to do
            None   => {}
        }
    }
    
    // if there are no pairs, return empty matrices (otherwise we will need to perform some transformations)
    if array_codomain_comb_inv_off_diag.is_empty() { return ( VecOfVec::new(array_codomain_comb_inv_off_diag), matching ) }

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
    matching.reverse();

    // return the (off-diagonal entries of the pivot block of the) codomain COMB
    ( VecOfVec::new(array_codomain_comb_inv_off_diag), matching )

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
/// **It is important that `order_operator_major` compares order of two entries based *solely* on
/// the value of the associated indices (not the associated coefficients).  This is easy
/// to get wrong, and it's hard to debug, so we are keeping the function private for now**
pub fn codomain_comb_inv_off_diag_pivot_block_unsafecomparator_skipmatched< Mapping, RingOperator, IterRowIndex, KeyBoth, OrderOperatorOfEntries >
    ( 
            mapping:                      &Mapping,  // use of a reference poses no problem for future calls to this function, since output has no references        
            iter_keymaj:                        IterRowIndex,
            ring_operator:                      RingOperator,            
            order_operator_of_entries:        OrderOperatorOfEntries,
            // mut order_operator_of_keymaj:     OrderOperatorColEntries,
    ) 
    -> 
    ( 
        VecOfVec< usize, Mapping::Coefficient >, 
        GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
    )
    where   Mapping:                           ViewRowAscend + IndicesAndCoefficients< RowIndex = KeyBoth, ColIndex = KeyBoth >,
            IterRowIndex:                             Iterator < Item = Mapping::RowIndex >,
            Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq + Debug, 
            Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            Mapping::Coefficient:                   Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            Mapping::ViewMajorAscend:          IntoIterator + ParetoShortCircuit<Mapping::EntryMajor>, // !!! remove clone
            Mapping::EntryMajor:     KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + Debug,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderOperatorOfEntries:               Clone + JudgePartialOrder <  Mapping::EntryMajor >, // !!! remove clone
            // OrderOperatorColEntries:       JudgePartialOrder<  Mapping::RowIndex >,
            // HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorOfEntries>: Clone // !!!! remove this

{
    
    // Initialize some objects
    let mut entries_to_elim_simplified_heap    =   HitMerge::new( order_operator_of_entries.clone() )
                                                                            .simplify( ring_operator.clone() ); // an empty (simplified) merger of iterators
    let mut matching          =   GeneralizedMatchingArrayWithMajorOrdinals::new(); // an all-zero generalized matching array
    let mut array_codomain_comb_inv_off_diag: Vec< Vec< ( usize, Mapping::Coefficient ) > >  =   Vec::new();    
    let mut codomain_comb_inv_off_diag_view_buffer   =   Vec::new();
    let mut codomain_comb_inv_off_diag_view_buffer_simplified   =   Vec::new();    

    // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
    // set up a variable to track the last major index produced by `iter_keymaj`; we will use this to ensure that the major keys returned by `iter_keymaj` appear in strictly ascending order
    // let mut prior_keymaj_opt = None;

    // let mut counter = 0;
    // let pb = indicatif::ProgressBar::new(100000);
    let mut sc_counter = 0;

    // build the (pivot block of the) codomain COMB row by row
    for keymaj in iter_keymaj {


        // print!("row={}",counter);
        // pb.println(format!("row = {}", counter));
        // counter +=1;        
        // print!("out");        
        // print!("simplex={:?}",keymaj.clone());
        // use std::{thread, time};
        // thread::sleep(time::Duration::from_millis(10));
    
        if matching.contains_keymin( & keymaj ) { 
            continue }      
          

        // UNCOMMENT THIS TO RECOVER A PORTION OF THE CODE THAT CHECKS THAT THE keymaj APPEAR IN STRICTLY ASCENDING ORDER
        // check that this major key is strictly greater than the last
        // if let Some( ref prior_keyaj ) = prior_keymaj_opt {
        //     match order_operator_of_keymaj.judge_lt( prior_keyaj, &keymaj ) {
        //         true => {  },
        //         false => { panic!("Umatch factorization cannot be completed because the sequence of major indices is not strictly ascending") }
        //     }               
        // }
        // prior_keymaj_opt.replace( keymaj.clone() );

        // clear the collection of entries to eliminate
        entries_to_elim_simplified_heap.unsimplified.clear();        

        // if possible, short circuit
        let next_iter = mapping.view_major_ascend( keymaj.clone() );
        if let Some( pivot_entry ) = next_iter.pareto_short_circuit() {
            let pivot_keymin = pivot_entry.key();   let pivot_coeff = pivot_entry.val();
            matching.push( keymaj, pivot_keymin, pivot_coeff );
            array_codomain_comb_inv_off_diag.push( vec![] );
            // print!("sc: {}", sc_counter );
            sc_counter += 1;
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
                    mapping,
                    & matching,
                    & array_codomain_comb_inv_off_diag,
                    ring_operator.clone(),
                    order_operator_of_entries.clone(),
                );


        match leading_entry_uneliminable_opt {
            // in this case we do not have a pivot
            Some( leading_entry_uneliminable )    => { 
                let pivot_keymin        =   leading_entry_uneliminable.key();
                let pivot_coeff         =   leading_entry_uneliminable.val();
                // println!("leading_entry_uneliminable {:?}", &leading_entry_uneliminable);              
                // update the matching array
                matching.push( keymaj, pivot_keymin, pivot_coeff ); 
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
                
                // update the vector that stores the (off diagonal entries of the) "pivot" part of the inverse of the codomain COMB
                array_codomain_comb_inv_off_diag.push( vec_to_push );
            }, 

            // in this case we don't have a pivot entry, so there is nothing to do
            None   => {}
        }       
    }
    
    // if there are no pairs, return empty matrices (otherwise we will need to perform some transformations)
    if array_codomain_comb_inv_off_diag.is_empty() { return ( VecOfVec::new(array_codomain_comb_inv_off_diag), matching ) }

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
    matching.reverse();

    // println!("pareto pairs: {:?}", sc_counter );
    // println!("number of matched pairs: {:?}", matching.num_pairs() );    

    // return the (off-diagonal entries of the pivot block of the) codomain COMB
    ( VecOfVec::new(array_codomain_comb_inv_off_diag), matching )

}











/// Calculates two quantities: (1) the square submatrix of the inverse of the codomain COMB 
/// indexed by row pivot indices, and (2) the matching array.
/// 
/// # Design notes
/// 
/// We use `VecOfVec< usize, Mapping::Coefficient >` to store the pivot block of the inverse of 
/// the codomain COMB rather than a `VecOfVec< ... >` for the following reasons: (1) we typically
/// need 
/// 
pub fn get_codomain_comb_inv_off_diag_pivot_block< Mapping, RingOperator, IterRowIndex, OrderOperatorRowEntries, > // OrderOperatorColEntries >
    ( 
            mapping:                      &Mapping,             
            iter_keymaj:                        IterRowIndex,
            ring_operator:                      RingOperator,            
            order_operator_for_keymin:        OrderOperatorRowEntries,
            // order_operator_for_keymaj:        OrderOperatorColEntries,            
    ) 
    -> 
    ( 
        VecOfVec< usize, Mapping::Coefficient >,
        GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
    )
    where   Mapping:                 ViewRowAscend + IndicesAndCoefficients,
            IterRowIndex:             Iterator < Item = Mapping::RowIndex >,
            Mapping::ColIndex:                 Clone + Hash + std::cmp::Eq + Debug, 
            Mapping::RowIndex:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            Mapping::Coefficient:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            Mapping::ViewMajorAscend:        IntoIterator + Clone, // !!! remove clone
            Mapping::EntryMajor:  KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderOperatorRowEntries:   Clone + JudgePartialOrder <  Mapping::ColIndex >, // !!! remove clone
            // OrderOperatorColEntries:   JudgePartialOrder <  Mapping::RowIndex >, // !!! remove clone            
            HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorRowEntries>>: Clone // !!!! remove this

{
    // let order_operator_for_entries_with_minor_keys = 
    //         |x: &Mapping::EntryMajor, y: &Mapping::EntryMajor | 
    //             order_operator_for_keymin.judge_lt( &x.key(), &y.key() );

    let order_operator_for_entries_with_minor_keys = OrderOperatorByKeyCutsom::new( order_operator_for_keymin );
    // let order_operator_for_entries_with_major_keys = OrderOperatorByKeyCutsom::new( order_operator_for_keymaj );
    
    let ( array_codomain_comb_inv_off_diag_pivot_block, matching ) = 
    codomain_comb_inv_off_diag_pivot_block_unsafecomparator(
            mapping,
            iter_keymaj,
            ring_operator,
            order_operator_for_entries_with_minor_keys,
            // order_operator_for_keymaj,            
    );

    ( array_codomain_comb_inv_off_diag_pivot_block, matching )
}



/// Sames as [get_codomain_comb_inv_off_diag_pivot_block], but applies the clearing optimization.
pub fn get_codomain_comb_inv_off_diag_pivot_block_with_clearing< Mapping, RingOperator, IterRowIndex, KeyBoth, OrderOperatorRowEntries, > // OrderOperatorColEntries >
    ( 
            mapping:                      &Mapping,             
            iter_keymaj:                        IterRowIndex,
            ring_operator:                      RingOperator,            
            order_operator_for_keymin:        OrderOperatorRowEntries,
            // order_operator_for_keymaj:        OrderOperatorColEntries,            
    ) 
    -> 
    ( 
        VecOfVec< usize, Mapping::Coefficient >,
        GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
    )
    where   Mapping:                           ViewRowAscend + IndicesAndCoefficients<ColIndex = KeyBoth, RowIndex=KeyBoth >,
            IterRowIndex:                             Iterator < Item = Mapping::RowIndex >,
            Mapping::ColIndex:                 Clone + Hash + std::cmp::Eq + Debug, 
            Mapping::RowIndex:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            Mapping::Coefficient:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            Mapping::ViewMajorAscend:        IntoIterator + ParetoShortCircuit<Mapping::EntryMajor>, // !!! remove clone
            Mapping::EntryMajor:  KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!
            OrderOperatorRowEntries:   Clone + JudgePartialOrder <  Mapping::ColIndex >, // !!! remove clone
            // OrderOperatorColEntries:   JudgePartialOrder <  Mapping::RowIndex >, // !!! remove clone            
            // HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorRowEntries>>: Clone // !!!! remove this

{
    // let order_operator_for_entries_with_minor_keys = 
    //         |x: &Mapping::EntryMajor, y: &Mapping::EntryMajor | 
    //             order_operator_for_keymin.judge_lt( &x.key(), &y.key() );

    let order_operator_for_entries_with_minor_keys = OrderOperatorByKeyCutsom::new( order_operator_for_keymin );
    // let order_operator_for_entries_with_major_keys = OrderOperatorByKeyCutsom::new( order_operator_for_keymaj );
    
    let ( array_codomain_comb_inv_off_diag_pivot_block, matching ) = 
    codomain_comb_inv_off_diag_pivot_block_unsafecomparator_skipmatched(
            mapping,
            iter_keymaj,
            ring_operator,
            order_operator_for_entries_with_minor_keys,
            // order_operator_for_keymaj,            
    );

    ( array_codomain_comb_inv_off_diag_pivot_block, matching )
}







//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- ROWS INDEXED BY KEYMIN, COLUMNS BY KEYMIN
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain COMB) * (the pivot block of the matching array).
/// 
/// This struct is almost identical to [`CombCodomainInvTimesMappingMatchedBlock`].  The only differences are that 
/// (i) the corresponding matrix oracle has rows indexed by matched minor keys of the  major keys, rather than 
/// by the matched major keys themselves, and (ii) the minor descending views iterate over entries in descending order 
/// of *minor* index, rather than major index.
#[derive(Copy, Clone)]
pub struct CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where     
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:              IntoIterator,
        Mapping::EntryMajor:         KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    pub umatch_ref:    &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,  
}

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    
    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:              IntoIterator,
        Mapping::EntryMajor:         KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,      
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin`].
    pub fn new( umatch_ref: &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin{ umatch_ref }
    }
}


//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT ORACLE FOR < KEYMIN, KEYMIN >
//  --------------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    IndicesAndCoefficients for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >       

    where
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,            

{   
    type EntryMajor = Mapping::EntryMajor;
    type EntryMinor = Mapping::EntryMajor;    
    type ColIndex = Mapping::ColIndex; 
    type RowIndex = Mapping::ColIndex; 
    type Coefficient = Mapping::Coefficient; 
}  








impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    ViewRowAscend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    where   
        Mapping:               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                     Clone,
        Mapping::ViewMajorAscend:            IntoIterator,
        Mapping::EntryMajor:      KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        RingOperator:               Clone + Semiring< Mapping::Coefficient >,
        OrderOperatorRowEntries:            Clone + JudgePartialOrder<  Mapping::EntryMajor >,
{
    type ViewMajorAscend            =   LinearCombinationSimplified< 
                                                OnlyIndicesInsideCollection< Mapping::ViewMajorAscendIntoIter, &'a HashMap<Mapping::ColIndex, usize>, Mapping::ColIndex, Mapping::Coefficient>, 
                                                Mapping::ColIndex, 
                                                Mapping::Coefficient, 
                                                RingOperator, 
                                                OrderOperatorRowEntries 
                                            >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymin: Mapping::ColIndex ) -> Self::ViewMajorAscend 
        {
            // This implementation simply invokes a major view of CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
            // This works because the two matrices have the same set of major views -- they simply look up these views via different indices 
            let ordmaj = self.umatch_ref.matching.keymin_to_ord( &keymin ).unwrap();
            let proxy 
                = self.umatch_ref.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_ordmaj();
            proxy.view_major_ascend( ordmaj )
        }
}        









impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    ViewColDescend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    where   
        Mapping:                               ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                       Clone,
        Mapping::ViewMajorAscend:              IntoIterator,
        Mapping::EntryMajor:         Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; `Clone` is required by 
        Mapping::ViewMinorDescend:             IntoIterator,
        Mapping::EntryMinor:        KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        RingOperator:                               Clone + Semiring< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        // &'a VecOfVec<usize, Mapping::Coefficient>: ViewRowAscend<usize, Cloned<Iter<'a, (usize, Mapping::Coefficient)>> >,
        // &'a VecOfVec<usize, Mapping::Coefficient>: ViewRowAscend<usize, Cloned<std::slice::Iter<'a, (Mapping::ColIndex, Mapping::Coefficient)>>>,
        // OrderOperatorRowEntries: 'a, RingOperator: 'a, Mapping::Coefficient: 'a, Mapping::RowIndex: 'a, Mapping::ColIndex: 'a, Mapping::ViewMajorAscend: 'a,            

{
    type ViewMinorDescend           =   IterWrappedVec< Mapping::EntryMajor >;
    type ViewMinorDescendIntoIter   =   IterWrappedVec< Mapping::EntryMajor >;

    fn view_minor_descend( &self, keymin: Mapping::ColIndex ) -> Self::ViewMinorDescend 
        {

        // a column of the matched block of R^{-1}*D; entries of the column are indexed by integers
        let col_indexed_by_ordmaj   =   self.umatch_ref.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_ordmaj().view_minor_descend( keymin );

        // reindex the entries with minor keys
        let mut col_reindexed      =   col_indexed_by_ordmaj.map( 
                                                |x| 
                                                Mapping::EntryMajor::new(
                                                        self.umatch_ref.matching.ord_to_keymin( x.0 ),
                                                        x.1
                                                    )
                                            )
                                            .collect_vec();
        col_reindexed.shrink_to_fit();

        // now all we have to do is sort the column according to (minor key) index!

        // repackage the order comparator for entries in major ascending views; we do this because the resulting struct has methods to produce std::comp::Ordering, which we will need to sort the column vector
        let order_decider =
                //     InferTotalOrderFromJudgePartialOrder::new( 
                            ReverseOrder::new( // we have to reverse the order because we want our output to be sorted in *descending* order
                                    self.umatch_ref.order_operator_major.clone() 
                                );
                        // );
        // wrap the order_decider in a closure
        let order_decider_function = |x: &Mapping::EntryMajor, y: &Mapping::EntryMajor|  order_decider.judge_partial_cmp(x, y).unwrap();
        // use the closure to the sort the vector
        col_reindexed.sort_by( order_decider_function );


        IterWrappedVec::new( col_reindexed )
    }

}









//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- INDEXED BY MINOR KEY ORDINAL INTEGERS
//  =========================================================================================================


//  !!! UNDER CONSTRUCTION
//  !!! AN IMPORTANT CHALLENGE IS THAT MATCHING MATRICES ARE NOT CURRENTLY DESIGNED TO STORE MINOR KEY ORDINALS



// /// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
// /// Concretely, this matrix equals (the pivot block of the inverse of the codomain COMB) * (the pivot block of the matching array).
// /// 
// /// This marix is indexed by **integers** (its major and minor keys are `usize`).  Specifically, row `i`
// /// corresponds to `\rho_i` and column `j` corresponds to `\kappa_j`, in the notation of "U-match factorization, ...".
// pub struct CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal< 'a, 'b, Mapping, Mapping::ViewMajorAscend, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
//     where     
//         Mapping::ColIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         Mapping::RowIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         Mapping::ViewMajorAscend:            IntoIterator,
//         Mapping::EntryMajor:      KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
// {
//     umatch_ref:     Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,  
// }

// impl < 'a, 'b, Mapping, Mapping::ViewMajorAscend, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    
//     CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal
//         < 'a, 'b, Mapping, Mapping::ViewMajorAscend, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
//     where   
//         Mapping::ColIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         Mapping::RowIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         Mapping::ViewMajorAscend:            IntoIterator,
//         Mapping::EntryMajor:      KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
//         RingOperator:               Clone,
//         OrderOperatorRowEntries:            Clone,
// {
//     // Make a new [`CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal`].
//     fn new( umatch_ref: & Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  ) -> Self {
//         CombCodomainInvTimesMappingMatchedBlockIndexedByColumnOrdinal{ umatch_ref: (*umatch_ref).clone() }
//     }
// }

    


//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- INDEXED BY MAJOR KEY ORDINAL INTEGERS
//  =========================================================================================================


//  !!! UNDER CONSTRUCTION
//
//  !!! ISSUE: THE OBJECT mapping DOES NOT RETURN ENTRIES IN ASCENDING ORDER OF MAJOR KEY ORDINAL; WOULD HAVE TO GENERATE ALL ENTRIES AND RETURN AS A VECTOR




//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- ROWS INDEXED BY ORDMAJ
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain COMB) * (the pivot block of the matching array).
/// 
/// This struct is almost identical to [`CombCodomainInvTimesMappingMatchedBlock`].  The only difference is 
/// that the corresponding matrix oracle has rows indexed by the integer ordinals of the matched major keys, rather than 
/// by the matched major keys themselves.
#[derive(Copy)]
pub struct CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where     
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:              IntoIterator,
        Mapping::EntryMajor:         KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    pub umatch_ref:    &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,  
}

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    
    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:              IntoIterator,
        Mapping::EntryMajor:         KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,      
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj`].
    pub fn new( umatch_ref: &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj{ umatch_ref }
    }
}

//  IMPLEMENT CLONE
//  --------------------------------------------------------------------------------------------------------------

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    Clone for 

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where     
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:              IntoIterator,
        Mapping::EntryMajor:         KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    fn clone(&self) -> Self { Self{ umatch_ref: self.umatch_ref } } 
}


//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT ORACLE FOR ROWS INDEXED BY usize, COLUMNS BY KEYMIN
//  --------------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    IndicesAndCoefficients for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >       

    where
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,            

{   
    type EntryMajor = Mapping::EntryMajor;
    type EntryMinor = (usize,Mapping::Coefficient);
    type ColIndex = Mapping::ColIndex; 
    type RowIndex = usize; 
    type Coefficient = Mapping::Coefficient; 
}  



impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    ViewRowAscend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    where   
        Mapping:               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                     Clone,
        Mapping::ViewMajorAscend:            IntoIterator,
        Mapping::EntryMajor:      KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        RingOperator:               Clone + Semiring< Mapping::Coefficient >,
        OrderOperatorRowEntries:            Clone + JudgePartialOrder<  Mapping::EntryMajor >,

{
    type ViewMajorAscend            =   LinearCombinationSimplified< 
                                                OnlyIndicesInsideCollection< Mapping::ViewMajorAscendIntoIter, &'a HashMap<Mapping::ColIndex, usize>, Mapping::ColIndex, Mapping::Coefficient>, 
                                                Mapping::ColIndex, 
                                                Mapping::Coefficient, 
                                                RingOperator, 
                                                OrderOperatorRowEntries 
                                            >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, ordmaj: usize ) -> Self::ViewMajorAscend 
        {
        
        // define a row vector
        let combining_coefficients 
            = self.umatch_ref.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal().view_major_ascend( ordmaj );     

        // the matched columns of the mapping array
        let matched_cols_of_mapping_array : OnlyKeyMinInsideCollection< &'a Mapping, &'a HashMap<Mapping::ColIndex, usize>, >
            =   self.umatch_ref.mapping_matched_cols_only();
        // let matched_cols_of_mapping_array // : OnlyKeyMinInsideCollection<'a, 'b, Mapping, Mapping::ViewMajorAscend, HashMap<Mapping::ColIndex, usize>, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient>
        //     =   self.umatch.mapping_ref();            

        // the terms whose sum equals the product of the vector with the matrix
        let iter_over_scaled_views = 
            combining_coefficients
                    // .iter()
                    .map(   
                            |( ordmaj, snzval)|  
                            matched_cols_of_mapping_array.view_major_ascend( 
                                    self.umatch_ref.matching.ord_to_keymaj( ordmaj ) 
                                )
                                .scale( snzval, self.umatch_ref.ring_operator.clone() )                                
                        );

        // sum the terms
        hit_merge_by_predicate( iter_over_scaled_views, self.umatch_ref.order_operator_major.clone() )
                    .simplify( self.umatch_ref.ring_operator.clone() )
    }
}




impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    ViewColDescend for

    CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    where   
        Mapping:                               ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                       Clone,
        Mapping::ViewMajorAscend:              IntoIterator,
        Mapping::EntryMajor:         KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        Mapping::ViewMinorDescend:             IntoIterator,
        Mapping::EntryMinor:        KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        RingOperator:                               Clone + Semiring< Mapping::Coefficient >,

{
    type ViewMinorDescend           =   LinearCombinationSimplified< 
                                                Chain<
                                                        Once<(usize,Mapping::Coefficient)>,
                                                        //VecOfVecMatrixColumnReverse<'a, usize, Mapping::Coefficient >,
                                                        Cloned<Rev<std::slice::Iter<'a, (usize, <Mapping as IndicesAndCoefficients>::Coefficient)>>>,
                                                    >,
                                                usize, 
                                                Mapping::Coefficient, 
                                                RingOperator, 
                                                ReverseOrder< OrderOperatorByKey< usize, Mapping::Coefficient, (usize, Mapping::Coefficient)> >, 
                                            >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;

    fn view_minor_descend( &self, keymin: Mapping::ColIndex ) -> Self::ViewMinorDescend 
        {

        // the matched rows of the mapping array
        let matched_rows_of_mapping_array 
            =   self.umatch_ref.mapping_matched_rows_only();   

        // one of its columns
        let column = matched_rows_of_mapping_array.view_minor_descend( keymin );

        // reindex the column
        let column_reindexed = 
            ChangeIndexSimple::new( 
                    column, 
                    self.umatch_ref.matching.bimap_maj_ref().val_to_ord_hashmap() 
                );
        
        // multiply by R^{-1}
        let column_product = 
            vector_matrix_multiply_minor_descend_simplified(
                    column_reindexed,
                    self.umatch_ref.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal(),
                    self.umatch_ref.ring_operator.clone(),
                    OrderOperatorByKey::new(),                  
                );

        column_product
    }
}

