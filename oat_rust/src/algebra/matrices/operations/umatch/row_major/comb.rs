//! Columar ordered matching bases
//! 
//! The notion of a columnar ordered matching basis (COMB) was introduced in [U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co)homology](https://arxiv.org/abs/2108.08831), by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.
//! 
//! This module defines structs for the domain and codomain COMBs, as well as their inverses.
//! 
//! # Design notes
//! 
//! We also define wrapper structs for views of the COMBs.  We use structs and not type
//! aliases because structs allow one to take advantage of type bounds on type parameters
//! (specifically, bounds that ensure existence of associated types), whereas parameters on type
//! aliases can't currently be constrained with `where` clauses.


// use debugit::DebugIt as D;
// use debugit::debugit;

//  =========================================================================================================
//  COMB'S (COLUMNAR ORDERED MATCHING BASES -- UNCOMPRESSED)
//  =========================================================================================================

use super::Umatch;
use super::construction::*;

use crate::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMinorKeys};
use crate::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified, vector_matrix_multiply_minor_descend_simplified};
use crate::algebra::matrices::operations::transform_vector_wise::{OnlyKeyMinInsideCollection};
use crate::algebra::matrices::operations::solve::triangle::{TriangularSolverMinorDescend};




use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVecMatrixColumnReverse;
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients};
use crate::algebra::vectors::entries::{ReindexEntry};
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet, KeyValNew};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};


use crate::algebra::vectors::operations::RemapEntryTuple;
use crate::utilities::functions::evaluate::{ IdentityFunction, EvaluateFunctionFnMutWrapper };
use crate::utilities::iterators::general::{IterTwoType, IterWrappedVec, OncePeekable, MapByTransform};
use crate::utilities::iterators::merge::hit::{HitMerge, hit_merge_by_predicate};
use crate::utilities::iterators::merge::two_type::MergeTwoItersByOrderOperator;
use crate::utilities::order::{JudgePartialOrder, OrderOperatorByKey, ReverseOrder};


use crate::algebra::vectors::operations::{Scale, VectorOperations, Simplify, OnlyIndicesInsideCollection, OnlyIndicesOutsideCollection, ChangeIndexSimple, ChangeEntryType, LinearCombinationSimplified};

use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::Rev;
use std::iter::{Cloned, Peekable, Once, Chain};
use std::marker::PhantomData;
use std::slice::{Iter};

use derive_getters::Dissolve;
use derive_new::new;
use itertools::Itertools;





// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the codomain COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct CombCodomain< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:            IntoIterator,
        Mapping::EntryMajor:      KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    pub umatch:     &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the inverse of the codomain COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct CombCodomainInv< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,    
        Mapping::ColIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:            IntoIterator,
        Mapping::EntryMajor:      KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    pub umatch:     &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the domain COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct CombDomain< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,    
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    pub umatch:     &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,
}

// <IGNORED> /// Uses the "inner identities" defined in [`Umatch factorization`](https://arxiv.org/pdf/2108.08831.pdf) 
// <IGNORED> /// to construct rows of the inverse of the domain COMB in a lazy fashion.
#[derive(Copy, Clone, new, Dissolve)]
pub struct CombDomainInv< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                           ViewRowAscend + IndicesAndCoefficients, 
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    pub umatch:     &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,
}




//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: DOMAIN
//  ---------------------------------------------------------------------------------------------------------



//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    IndicesAndCoefficients for 

    CombDomain
            < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    where   
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,    
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,            
            
{   
    type EntryMajor = Mapping::EntryMajor;
    type EntryMinor = Mapping::EntryMajor;    
    type RowIndex = Mapping::ColIndex; 
    type ColIndex = Mapping::ColIndex; 
    type Coefficient = Mapping::Coefficient;  
}            



//  ORACLE MAJOR ASCEND
//  ---------------------------------------------------------------------------------------------------------


//  FOR UNIT TESTS, RUN A TEXT SEARCH FOR <#CombDomainImplementViewRowAscend>

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    ViewRowAscend for 

    CombDomain< 
            'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries,
        >

    where 
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::Coefficient:                       Clone,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,      
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValNew is used (i) to create the diagonal elements of the array, and (ii) to create the 'problem vector' b in a linear problem xA = b whose solution is used in our constructor
        Mapping:                               ViewRowAscend + IndicesAndCoefficients
{
    type ViewMajorAscend            =   CombDomainViewMajorAscend< 
                                                'a, Mapping, RingOperator, OrderOperatorRowEntries, 
                                            >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymin: Self::RowIndex ) -> Self::ViewMajorAscend       

    {
        match self.umatch.matching.contains_keymin( &keymin ) { 
            false => { 

                CombDomainViewMajorAscend{
                                iter_unwrapped:         IterTwoType::Iter1(
                                                                std::iter::once(  Mapping::EntryMajor::new( keymin, RingOperator::one() )  ) 
                                                            ),
                                phantom_arraymapping:   PhantomData,
                            }
            }
            true => { 

                // The matrix A from Hang et al., "Umatch factorization ...", with rows indexed by major key ordinals
                // Note: this struct only contains a reference
                let seed = self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_ordmaj();


                // Struct that encodes the matching from minor keys to major key ordinals
                let matching_ref = self.umatch.matching_ref();
                let keymin_to_ordmaj_wrapped = matching_ref.bimap_min_ref().val_to_ord_hashmap();

                // obtain a row of A^{-1}
                let seed_inv_row_vec = 
                    EchelonSolverMajorAscendWithMinorKeys::solve(
                            std::iter::once( Mapping::EntryMajor::new( keymin.clone(), RingOperator::one() ) ), // the standard unit vector supported on `keymin`
                            seed.clone(), // matrix A
                            keymin_to_ordmaj_wrapped.clone(),
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_operator_major(),
                        )
                        .solution().unwrap();

                // obtain a copy of R_{\rho \rho}
                let comb_codomain_inv_matched_block = self.umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal();

                // obtain a copy of the column submatrix of the mapping array indexed by unmatched column indices
                let mapping_npcols = self.umatch.mapping_matchless_cols_only();

                // y = - seed_inv_row_vec * R_{\rho \rho}^{-1} * D_{ \rho \bar \kappa}
                let x 
                    =   vector_matrix_multiply_major_ascend_simplified(
                                seed_inv_row_vec
                                    .iter()
                                    .map(   |x| 
                                            (   
                                                matching_ref.keymin_to_ord(& x.key() ).unwrap(),  // replace integer indices (major ordinals) with major key indices
                                                self.umatch.ring_operator.negate( x.val() ) // multiply by -1 (recally that we are computing *MINUS* seed_inv_row_vec * R_{\rho \rho}^{-1} * D_{ \rho \bar \kappa})
                                            ) 
                                        ),
                                comb_codomain_inv_matched_block,
                                self.umatch.ring_operator.clone(),
                                OrderOperatorByKey::new(),
                            );

                let y
                    =   vector_matrix_multiply_major_ascend_simplified(
                                x.map( 
                                        | ( ordmaj, snzval ) |
                                        ( 
                                                matching_ref.ord_to_keymaj( ordmaj ),
                                                snzval
                                            )
                                    ),
                                mapping_npcols,
                                self.umatch.ring_operator.clone(),
                                self.umatch.order_operator_major.clone(),                                
                            );

                // rescale entries of seed_inv_row_vec, transforming it into row `keymin` of `A^{-1} M_{\rho \kappa}`
                // NB: this has to occur AFTER we've used `seed_inv_row_vec` to construct `y`
                let mut seed_inv_row_vec_times_matching = seed_inv_row_vec;
                for entry in seed_inv_row_vec_times_matching.iter_mut() {
                    entry.set_val( 
                            self.umatch.ring_operator.multiply( 
                                    entry.val(),  
                                    self.umatch.matching_ref().keymin_to_snzval( &entry.key() )
                                )                         
                        )
                }

                // merge seed_inv_row_vec
                let z = MergeTwoItersByOrderOperator::new( 
                                                            IterWrappedVec::new( seed_inv_row_vec_times_matching ).peekable(), 
                                                            y.peekable(), 
                                                            self.umatch.order_operator_major.clone()
                                                        );


                // println!("MADE IT THROUGH END OF CombDomain.view_major_ascend (MATCHED INDEX)");                
                CombDomainViewMajorAscend{
                                iter_unwrapped:         IterTwoType::Iter2( z ),
                                phantom_arraymapping:   PhantomData,
                            }

            }
        }
    }
}


pub struct CombDomainViewMajorAscend< 
                    'a, Mapping, RingOperator, OrderOperatorRowEntries, 
                >

    where 
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::Coefficient:                       Clone,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValNew is used (i) to create the diagonal elements of the array, and (ii) to create the 'problem vector' b in a linear problem xA = b whose solution is used in our constructor
        Mapping:                               ViewRowAscend + IndicesAndCoefficients                

{
    iter_unwrapped:     IterTwoType< 
                                Once< Mapping::EntryMajor >,
                                MergeTwoItersByOrderOperator<
                                        Peekable< IterWrappedVec< Mapping::EntryMajor > >, 
                                        Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection< Mapping::ViewMajorAscendIntoIter, &'a HashMap<Mapping::ColIndex, usize>, Mapping::ColIndex, Mapping::Coefficient>, 
                                                        Mapping::ColIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries
                                                    >, 
                                            >,
                                        OrderOperatorRowEntries
                                    >,
                            >,
    phantom_arraymapping:   PhantomData< Mapping >
}  


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, >

    Iterator for

    CombDomainViewMajorAscend< 
            'a, Mapping, RingOperator, OrderOperatorRowEntries, 
        >   

    where 
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::Coefficient:                       Clone,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValNew is used (i) to create the diagonal elements of the array, and (ii) to create the 'problem vector' b in a linear problem xA = b whose solution is used in our constructor
        Mapping:                               ViewRowAscend + IndicesAndCoefficients                        

{
    type Item = Mapping::EntryMajor;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}         



//  ORACLE MINOR DESCEND
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    ViewColDescend for 

    CombDomain< 
            'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries,
        >
    where
        Mapping:                               ViewRowAscend + ViewColDescend + IndicesAndCoefficients,        
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::Coefficient:                       Clone,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor  >,
        OrderOperatorColEntries:       Clone + JudgePartialOrder<  Mapping::EntryMinor >,                
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
{
    type ViewMinorDescend           =   CombDomainViewMinorDescend< 
                                                'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries,
                                            >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;
    
    fn  view_minor_descend( &self, index: Mapping::ColIndex ) -> Self::ViewMinorDescend {

        match self.umatch.matching.bimap_min_ref().ord( & index ) {

            Some( ordmaj ) => {                     
                let matched_snzval  =   self.umatch.matching.ord_to_snzval(ordmaj);
                // let matched_keymaj = self.umatch.matching.ordmaj_to_keymaj(ordmaj);
                let unit_vec = OncePeekable::new( Mapping::EntryMajor::new( index, matched_snzval ) );
                let col = 
                    TriangularSolverMinorDescend::solve(
                        unit_vec,
                        self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin(), // matrix A
                        self.umatch.ring_operator.clone(),
                        self.umatch.order_operator_major.clone(), // the constructor function `EchelonSolverMajorAscendWithMinorKeys::solve(` takes care of reversing the order comparator, so we don't have to                        
                    );
                    // EchelonSolverMinorDescendWithMinorKeys::solve(
                    //         unit_vec,
                    //         self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin(), // matrix A
                    //         IdentityFunction{},
                    //         self.umatch.ring_operator.clone(),
                    //         self.umatch.order_operator_major.clone(), // the constructor function `EchelonSolverMajorAscendWithMinorKeys::solve(` takes care of reversing the order comparator, so we don't have to
                    //     )
                    //     .quotient();   
                // let col = ChangeEntryType::new( col ); // change entries from tuples to ArrayMappingEntryMajor
                CombDomainViewMinorDescend{ 
                                iter_unwrapped:     IterTwoType::Iter1( col ) 
                            }
            }
            None => {
                // a column c of D[ matched_rows, : ] 
                let col_matched_rows_only 
                    =   self.umatch.mapping_matched_rows_only().view_minor_descend( index.clone() )
                            .negate( self.umatch.ring_operator() ); // this multiplies the column by -1
                            // .map( | entry | (self.umatch.matching.bimap_maj.ord( entry.key() ), entry.val()) ); // reindex so that entries are indexed by ordinals of matched major keys
                // get an object that represents the REVERSED total order on minor keys
                let order_operator_major_reverse_total = self.umatch.order_operator_major_reverse();
                // a linear combination v of columns of R_{\rho \rho}^{-1}; we reindex by minor key, then sort by minor key
                let mut lin_comb_r 
                    =   vector_matrix_multiply_minor_descend_simplified(
                                col_matched_rows_only, 
                                self.umatch.comb_codomain_inv_matched_block_indexed_by_keymaj(), 
                                self.umatch.ring_operator(), 
                                self.umatch.order_operator_minor()
                            )
                            .map(   |x| //|(key,val)| 
                                    Mapping::EntryMajor::new(
                                            self.umatch.matching.keymaj_to_keymin(&x.key() ).unwrap(),
                                            x.val(),
                                        )
                                )
                            .collect_vec(); // collect into a vector
                lin_comb_r.sort_by( |a,b| order_operator_major_reverse_total.judge_partial_cmp( a, b ).unwrap() );  // sort the vector according to minor key
                // let lin_comb_r = lin_comb_r.into_iter();
                // A^{-1} v  <-- this is equal to the matched part of the column we want to construct
                let matched_part = 
                    TriangularSolverMinorDescend::solve(
                        lin_comb_r, // the vector b in "solve Ax = b"
                        self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin(), // matrix A
                        self.umatch.ring_operator(),
                        self.umatch.order_operator_major(),
                    );
                    // EchelonSolverMinorDescendWithMinorKeys::solve(
                    //         IterWrappedVec::new( lin_comb_r ),
                    //         self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin(), // matrix A
                    //         IdentityFunction{},
                    //         self.umatch.ring_operator.clone(),
                    //         self.umatch.order_operator_major.clone(), // the constructor function `EchelonSolverMajorAscendWithMinorKeys::solve(` takes care of reversing the order comparator, so we don't have to
                    //     )
                    //     .quotient();
                // change the entry type of this iterator from (RowIndex, Coefficient) to Mapping::EntryMinor
                let matched_part =   ChangeEntryType::new( matched_part, ); // rust infers the new entry type that we want, automatically
                let col =   MergeTwoItersByOrderOperator::new(
                                    matched_part.peekable(),
                                    OncePeekable::new( Mapping::EntryMajor::new( index, RingOperator::one() ) ),
                                    self.umatch.order_operator_major_reverse(), // recall that we want entries returned in *descending* order
                                );                             
                CombDomainViewMinorDescend{ 
                                    iter_unwrapped:     IterTwoType::Iter2( col )
                                }                                
            }
        }
    }
}



pub struct CombDomainViewMinorDescend< 
                    'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries,
                >
    where
        Mapping:                               ViewRowAscend + ViewColDescend + IndicesAndCoefficients,        
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::Coefficient:                       Clone,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor  >,
        OrderOperatorColEntries:       Clone + JudgePartialOrder<  Mapping::EntryMinor >,                
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,                        
{
    iter_unwrapped:    IterTwoType<     
                                TriangularSolverMinorDescend<
                                        OncePeekable< Mapping::EntryMajor >, 
                                        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
                                            <'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries>,
                                        Mapping::ColIndex,
                                        RingOperator,
                                        OrderOperatorRowEntries,                                                                                                                                
                                    >,                                 
                                // ChangeEntryType<
                                //         // EchelonSolverMinorDescendWithMinorKeys<
                                //         //         OncePeekable< Mapping::EntryMajor >, 
                                //         //         IdentityFunction, 
                                //         //         CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
                                //         //             <'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries>, 
                                //         //         RingOperator, 
                                //         //         OrderOperatorRowEntries
                                //         //     >,                                           
                                //         Mapping::EntryMajor,
                                //         Mapping::ColIndex,
                                //         Mapping::Coefficient,                              
                                //     >,                    
                                MergeTwoItersByOrderOperator<
                                        Peekable< 
                                                ChangeEntryType<
                                                        TriangularSolverMinorDescend<
                                                                Vec< Mapping::EntryMajor >,
                                                                CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
                                                                    <'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries>,
                                                                    Mapping::ColIndex,
                                                                    RingOperator,
                                                                    OrderOperatorRowEntries,                                                                                                                                
                                                            >,
                                                        // EchelonSolverMinorDescendWithMinorKeys<
                                                        //         IterWrappedVec< Mapping::EntryMajor >, 
                                                        //         IdentityFunction, 
                                                        //         CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin
                                                        //             <'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries>, 
                                                        //         RingOperator, 
                                                        //         OrderOperatorRowEntries
                                                        //     >, 
                                                        Mapping::EntryMajor,
                                                        Mapping::ColIndex,
                                                        Mapping::Coefficient,                                    
                                                    >,
                                            >,
                                        OncePeekable<Mapping::EntryMajor>, 
                                        ReverseOrder< OrderOperatorRowEntries >
                                    >,        
                            >,
}

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    Iterator for 

    CombDomainViewMinorDescend< 
            'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries,
        >

    where
        Mapping:                               ViewRowAscend + ViewColDescend + IndicesAndCoefficients,        
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array
        Mapping::Coefficient:                       Clone,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor  >,
        OrderOperatorColEntries:       Clone + JudgePartialOrder<  Mapping::EntryMinor >,                
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,                        
{
    type Item = Mapping::EntryMajor;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}            










//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: DOMAIN INVERSE
//  ---------------------------------------------------------------------------------------------------------

//  FOR UNIT TESTS, RUN A TEXT SEARCH FOR <#CombDomainInvImplementViewRowAscend>

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    IndicesAndCoefficients for 

    CombDomainInv
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >   
        
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


//  ORACLE MAJOR ASCEND
//  ---------------------------------------------------------------------------------------------------------

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    ViewRowAscend for 

    CombDomainInv
            < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    where 
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        Mapping::Coefficient:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,
{
    type ViewMajorAscend            =   CombDomainInvViewMajorAscend
                                            < Mapping, RingOperator, OrderOperatorRowEntries, >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymin: Mapping::ColIndex )  // this function signature looks strange b/c the major keys of the domain COMB have type Mapping::ColIndex!  
        -> Self::ViewMajorAscend      

    {
        match self.umatch.matching.keymin_to_ord( &keymin ) { // this function call looks strange b/c the major keys of the domain COMB have type Mapping::ColIndex!  
            None => { 
                // this is a non-pivot row, so return a unit vector
                CombDomainInvViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter1(
                                                            std::iter::once(  Mapping::EntryMajor::new( keymin, RingOperator::one() )  ) 
                                                        )
                    }
            }
            Some( ordmaj ) => {  
                let scalar  =   self.umatch.matching.ord_to_snzval( ordmaj );
                let scalar_inv      =   self.umatch.ring_operator.invert( scalar );
                // let key_to_codomain_comb    =   self.umatch.matching.ordmaj_to_keymaj(ordmaj);

                // this equals row `keymaj` of M_{rk}^{-1} * R^{-1}_{\rho \rho}
                let lefthand_factor_vec     =   
                        self.umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal()                                                    
                            .view_major_ascend( ordmaj )
                            .scale( scalar_inv, self.umatch.ring_operator.clone() )
                            .map(   |(x,y)|   // re-index the entries, so that indices align with the row indices of the mapping array
                                    (self.umatch.matching.ord_to_keymaj(x), y) 
                                );

                let iter = vector_matrix_multiply_major_ascend_simplified ( 
                            lefthand_factor_vec,
                            & self.umatch.mapping,
                            self.umatch.ring_operator.clone(),
                            self.umatch.order_operator_major.clone(),
                        );

                CombDomainInvViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter2(  iter  )
                            }                        
            }
        }
    }
}


pub struct  CombDomainInvViewMajorAscend
                < Mapping, RingOperator, OrderOperatorRowEntries, >
    where 
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        Mapping::Coefficient:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,                
{
    iter_unwrapped:     IterTwoType< 
                                Once
                                    < Mapping::EntryMajor >,
                                LinearCombinationSimplified
                                    < Mapping::ViewMajorAscendIntoIter, Mapping::ColIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries >,
                            >,
}

impl < Mapping, RingOperator, OrderOperatorRowEntries, >

    Iterator for

    CombDomainInvViewMajorAscend
        < Mapping, RingOperator, OrderOperatorRowEntries, >
        
    where 
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        Mapping::Coefficient:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping:                               ViewRowAscend + IndicesAndCoefficients,                
        
{
    type Item = Mapping::EntryMajor;
    
    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}



//  ORACLE MINOR DESCEND
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    ViewColDescend for 

    CombDomainInv
            < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    where 
        Mapping:                               ViewRowAscend + ViewColDescend + IndicesAndCoefficients,    
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        Mapping::Coefficient:                       Clone + Debug,  // !!! DELETE THE DEBUG REQUIREMENT LATER; IT'S JUST USED FOR DEBUGGING
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,        
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor >,        
        OrderOperatorColEntries:       Clone + JudgePartialOrder<  Mapping::EntryMinor >,                
{
    type ViewMinorDescend           =   Vec< Mapping::EntryMajor >;
    type ViewMinorDescendIntoIter   =   < Self::ViewMinorDescend as IntoIterator >::IntoIter;    

    fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend {
        let col     =   self.umatch.mapping_matched_rows_only().view_minor_descend( index.clone() );
        let mut solution    =   
            vector_matrix_multiply_minor_descend_simplified(
                    col, 
                    self.umatch.comb_codomain_inv_matched_block_indexed_by_keymaj(), 
                       self.umatch.ring_operator(), 
                    self.umatch.order_operator_minor(),
                )
                .map(
                    |x|
                    {
                        let ordmaj = self.umatch.matching.keymaj_to_ord( &x.key() ).unwrap();
                        let matched_keymin = self.umatch.matching.ord_to_keymin( ordmaj );
                        let matched_snzval = self.umatch.matching.ord_to_snzval( ordmaj );
                        let coefficient = self.umatch.ring_operator.divide( x.val(), matched_snzval );
                        Mapping::EntryMajor::new( matched_keymin, coefficient )
                    }
                )
                .collect_vec();
        // if `index` is not a matched minor key, then append an entry of form (index, 1)
        if ! self.umatch.matching.contains_keymin( & index ) {
            solution.push( Mapping::EntryMajor::new( index, RingOperator::one() ) );
        }

        // sort the entries in the solution in descending order
        let order_operator_major_reverse_extended = self.umatch.order_operator_major_reverse(); //_extended_functionality();
        solution.sort_by( |a,b| order_operator_major_reverse_extended.judge_partial_cmp(a, b).unwrap() ); // sort the entries of the vector
        
        solution

    }
}


//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: CODOMAIN
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    IndicesAndCoefficients for

    CombCodomain
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    where
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,            

{   
    type EntryMajor = Mapping::EntryMinor;
    type EntryMinor = Mapping::EntryMinor;    
    type ColIndex = Mapping::RowIndex; 
    type RowIndex = Mapping::RowIndex; 
    type Coefficient = Mapping::Coefficient; 
}        



//  IMPLEMENT ORACLE MAJOR ASCEND FOR COMB: CODOMAIN
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    ViewRowAscend for 

    CombCodomain
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    where 
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::Coefficient:                   Clone,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:    Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        OrderOperatorColEntries:   Clone + JudgePartialOrder<  Mapping::EntryMinor >,
        Mapping::ViewMajorAscend:          IntoIterator,        
        Mapping::EntryMajor:            Clone + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping::EntryMinor:            Clone + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,        
        Mapping:                        ViewRowAscend + IndicesAndCoefficients
{
    type ViewMajorAscend            =   CombCodomainViewMajorAscend< Mapping::EntryMinor >;
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
            keymaj: Mapping::RowIndex 
        ) 
        -> Self::ViewMajorAscend
    {

        // define the matrix A that will fit into an equation xA = b
        let seed_with_integer_indexed_rows = self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_ordmaj();

        // define the problem vector b that will fit into an equation xA = b
        let mapping_matched_cols_only = self.umatch.mapping_matched_cols_only();
        let problem_vector = mapping_matched_cols_only.view_major_ascend( keymaj.clone() );

        // Struct that encodes the matching from minor keys of A to major key ordinals of A
        let matching_ref = self.umatch.matching_ref();
        let keymin_to_ordmaj = | keymin: Mapping::ColIndex | -> Option< usize > { matching_ref.keymin_to_ord( &keymin ) };
        let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

        // Solve xA = b for x
        let mut solution_vec = 
        EchelonSolverMajorAscendWithMinorKeys::solve( // solution to xA = b
                problem_vector, // mabrix b
                seed_with_integer_indexed_rows, // matrix A
                keymin_to_ordmaj_wrapped,
                self.umatch.ring_operator.clone(),
                self.umatch.order_operator_major.clone(),
            )
            .solution().unwrap();

        // Sort solution vector x
        // println!("REDO THIS SECTION OF CODE AFTER REFACTORING <CompareOrder>");
        let order_operator_major_clone = self.umatch.order_operator_major.clone(); // we have to clone the order comparator in order to compare order, since the `lt` method takes `& mut self` as an argument -- but to get a mutable reference in the present context we would have to borrow the umatch struct as mutable
        let compare_order_for_vector_sort 
            = | x: &Mapping::EntryMajor, y: &Mapping::EntryMajor |-> Ordering { // we have "repackage" order_operator_minor so that it returns an std::cmp::Ordering rather than a bool
            match order_operator_major_clone.judge_lt( x, y ) {
                true => { Ordering::Less },
                false => { 
                    match order_operator_major_clone.judge_lt( y, x ) {
                        true => { Ordering::Greater },
                        false => { Ordering:: Equal }
                    }
                }
            }
        };
        solution_vec.sort_by( compare_order_for_vector_sort );

        // Reindex solution vector x
        let mut solution_vec_reindexed = 
            solution_vec.into_iter()
                .map(   | item |
                        Mapping::EntryMinor::new( 
                            self.umatch.matching.keymin_to_keymaj( &item.key() ).unwrap(),
                            item.val(),
                        )
                    )
                .collect_vec();

        // a function that sends (a,b) to `Less` if a < b and to `Greater` otherwise
        let order_operator_minor = self.umatch.order_operator_minor.clone(); // we have to make a clone because, as currently written, `self` is behind a mutable reference and `order_operator_minor` may mutate itself when it compares two objects
        solution_vec_reindexed.sort_by( 
                |a,b| {  
                    if order_operator_minor.judge_lt( a, b ) { Ordering::Less }
                    else { Ordering::Greater }
                    }
            ); 
        let solution_vec_reindexed_and_sorted = IterWrappedVec::new( solution_vec_reindexed );                            

        // Possibly append an entry with coefficient 1
        match self.umatch.matching_ref().contains_keymaj( & keymaj ){
            false =>    { // if keymaj is UNMATCHED then merge in a copy of ( keymaj, 1 )
                CombCodomainViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter1(
                                                            std::iter::once( Mapping::EntryMinor::new(keymaj, RingOperator::one() ) )
                                                                .chain( solution_vec_reindexed_and_sorted ) 
                                                        )
                            }

            }
            true =>     { // otherwise just return the reindexed, sorted solution vector
                CombCodomainViewMajorAscend{
                                iter_unwrapped:     IterTwoType::Iter2(
                                                            solution_vec_reindexed_and_sorted
                                                        )
                            }
            }
        }
    }   
}


pub struct CombCodomainViewMajorAscend< T: Clone >
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

    CombCodomainViewMajorAscend< T >
{
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}






//  IMPLEMENT ORACLE MINOR DESCEND FOR COMB: CODOMAIN
//  ---------------------------------------------------------------------------------------------------------

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    ViewColDescend for

    CombCodomain
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector      
        Mapping::EntryMinor:    Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector     
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorColEntries:   Clone + JudgePartialOrder<  Mapping::EntryMinor >,
        OrderOperatorRowEntries:    Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        Mapping::ColIndex:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        Mapping::Coefficient:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING        
{
    type ViewMinorDescend           =   CombCodomainViewMinorDescend
                                            < Mapping, RingOperator, OrderOperatorColEntries >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;


    fn   view_minor_descend( &self, keymaj: Mapping::RowIndex ) -> Self::ViewMinorDescend {


        match self.umatch.matching_ref().keymaj_to_keymin( & keymaj ) {

            // follow this branch if keymaj is matched to a minor key, keymin
            Some( keymin ) => {

                // 
                let unit_vector_keymin = std::iter::once( Mapping::EntryMajor::new( keymin, RingOperator::one() ) );   // the standard unit vector supported on `keymaj`                                

                // the matched block of R^{-1} D, indexed by minor keys along rows and columns
                let seed = self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin();

                // DEPRECATED (OK TO DELETE)
                // the matching relation from rows of A to columns of A
                // let keymaj_to_keymin = self.umatch.matching_ref().bijection_keymaj_to_keymin();
                
                // a column c of A^{-1}
                let column_of_a_inv = TriangularSolverMinorDescend::solve(
                        unit_vector_keymin,
                        seed, // matrix A
                        self.umatch.ring_operator(),
                        self.umatch.order_operator_major(), // the constructor function `EchelonSolverMajorAscendWithMinorKeys::solve(` takes care of reversing the order comparator, so we don't have to
                    );
                
// ------------------ diagnostic: begin

                // // MAJOR KEYS

                // let unit_vector_keymin_copy = std::iter::once( Mapping::EntryMajor::new( keymin.clone(), RingOperator::one() ) );   // the standard unit vector supported on `keymaj`                                                
                // let seed_copy = self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin();
                // let mut column_of_A_inv_copy = EchelonSolverMinorDescendWithMinorKeys::solve(
                //     unit_vector_keymin_copy,
                //     seed_copy, // matrix A
                //     EvaluateFunctionFnMutWrapper::new( |x| x ),
                //     self.umatch.ring_operator.clone(),
                //     self.umatch.order_operator_major.clone(), // the constructor function `EchelonSolverMajorAscendWithMinorKeys::solve(` takes care of reversing the order comparator, so we don't have to
                // );

                // for p in 0 .. 3 { println!("column_of_A_inv_copy_majorkeys.next() = {:?}", column_of_A_inv_copy.next()) }            
                // println!("LENGTH( column_of_A_inv_copy_majorkeys ) = (next line)");
                // println!("{:?}", column_of_A_inv_copy.count() );


                // // MINOR KEYS
                // // let keymaj_to_keymin_copy = self.umatch.matching_ref().bijection_keymaj_to_keymin();
                // // let unit_vector_copy = std::iter::once( Mapping::EntryMinor::new( keymaj.clone(), RingOperator::one() ) );   // the standard unit vector supported on `keymaj`
                // // let mut column_of_A_inv_copy = EchelonSolverMinorDescendWithMinorKeys::solve(
                // //         unit_vector_copy,
                // //         & seed, // matrix A
                // //         keymaj_to_keymin_copy,
                // //         self.umatch.ring_operator.clone(),
                // //         self.umatch.order_operator_minor.clone(), // the constructor function `EchelonSolverMajorAscendWithMinorKeys::solve(` takes care of reversing the order comparator, so we don't have to
                // //     );    
                // // for p in 0 .. 3 { println!("column_of_A_inv_copy.next() = {:?}", column_of_A_inv_copy.next()) }            
                // // println!("LENGTH( column_of_a_inv ) = (next line)");
                // // println!("{:?}", column_of_A_inv_copy.count() );

// ------------------ diagnostic: end


                // the product vector D_{m k} * c  (refer to the matching identities in the Umatch facotrization paper for explanation)
                let column_of_comb  =   vector_matrix_multiply_minor_descend_simplified(
                                    column_of_a_inv,
                                                self.umatch.mapping_ref(),
                                                self.umatch.ring_operator.clone(),
                                                self.umatch.order_operator_minor.clone(), // the constructor handles reversing the order for us, so we don't have to
                                            );                                          
                
                // wrap the product vector in an enum
                CombCodomainViewMinorDescend{
                                iter_unwrapped:     IterTwoType::Iter1( column_of_comb )
                            }

            }

            // follow this branch if keymaj is not matched to a minor key
            None => {

                let unit_vector_keymaj = std::iter::once( Mapping::EntryMinor::new( keymaj.clone(), RingOperator::one() ) );   // the standard unit vector supported on `keymaj`                
                CombCodomainViewMinorDescend{
                                iter_unwrapped:     IterTwoType::Iter2( unit_vector_keymaj )
                            }                                
            }
        }


        
        
    }
        
}




pub struct CombCodomainViewMinorDescend
                < Mapping, RingOperator, OrderOperatorColEntries >

    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector      
        Mapping::EntryMinor:    Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector     
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorColEntries:   Clone + JudgePartialOrder<  Mapping::EntryMinor >,
        Mapping::ColIndex:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        Mapping::Coefficient:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING                        
{
    iter_unwrapped:     IterTwoType<
                                LinearCombinationSimplified
                                    < Mapping::ViewMinorDescendIntoIter, Mapping::RowIndex, Mapping::Coefficient, RingOperator, ReverseOrder< OrderOperatorColEntries > >,
                                std::iter::Once< Mapping::EntryMinor >,
                            >
}             

impl    < Mapping, RingOperator, OrderOperatorColEntries >

        Iterator for 

        CombCodomainViewMinorDescend
                < Mapping, RingOperator, OrderOperatorColEntries >

    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector      
        Mapping::EntryMinor:    Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct; KeyValNew is used to contruct a unit vector     
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorColEntries:   Clone + JudgePartialOrder<  Mapping::EntryMinor >,
        Mapping::ColIndex:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING
        Mapping::Coefficient:                   Debug, // !!!! DELETE THIS, ONLY USED FOR DEBUGGING                        
{
    type Item = Mapping::EntryMinor;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}





//  ---------------------------------------------------------------------------------------------------------
//  IMPLEMENT ORACLE FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    IndicesAndCoefficients for

    CombCodomainInv
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >    

    where
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,            

{   
    type EntryMajor = Mapping::EntryMinor;
    type EntryMinor = Mapping::EntryMinor;    
    type ColIndex = Mapping::RowIndex; 
    type RowIndex = Mapping::RowIndex; 
    type Coefficient = Mapping::Coefficient; 
}  



//  IMPLEMENT ORACLE MAJOR ASCEND FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------



impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    ViewRowAscend for 

    CombCodomainInv 
            < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >

    where 
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::Coefficient:                   Clone,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:    Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        OrderOperatorColEntries:   Clone + JudgePartialOrder<  Mapping::EntryMinor >,
        Mapping::ViewMajorAscend:          IntoIterator,        
        Mapping::EntryMajor:        Clone + 
                                    KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping::EntryMinor:        Clone + 
                                    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
{

    type ViewMajorAscend            =   // Vec< Mapping::EntryMinor >;
                                        CombCodomainInvViewMajorAscend
                                           < 'a, Mapping, RingOperator, >;
                                        // ChangeEntryType<   
                                        //         IterTwoType< 
                                        //                 MergeTwoItersByOrderOperator<
                                        //                         Peekable<
                                        //                                 ChangeIndexSimple<
                                        //                                         LinearCombinationSimplified
                                        //                                             <Chain<Once<(usize, Mapping::Coefficient)>, Cloned<Iter<'a, (usize, Mapping::Coefficient)>>>, usize, Mapping::Coefficient, RingOperator, OrderOperatorByKey<usize, Mapping::Coefficient, (usize, Mapping::Coefficient)>>, 
                                        //                                         &'a Vec<Mapping::RowIndex>, 
                                        //                                         usize, 
                                        //                                         Mapping::RowIndex, 
                                        //                                         Mapping::Coefficient
                                        //                                     >
                                        //                             >, 
                                        //                         OncePeekable<(Mapping::RowIndex, Mapping::Coefficient)>, 
                                        //                         OrderOperatorColEntries
                                        //                     >,                        
                                        //                 ChangeIndexSimple<
                                        //                         Chain<
                                        //                                 Once<(usize, Mapping::Coefficient)>, 
                                        //                                 Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
                                        //                             >, 
                                        //                         &'a Vec<Mapping::RowIndex>, usize, Mapping::RowIndex, Mapping::Coefficient
                                        //                     >
                                        //             >, 
                                        //         Mapping::EntryMinor,
                                        //         Mapping::RowIndex, 
                                        //         Mapping::Coefficient,    
                                        //     >;
                                //         IterTwoType<
                                //         std::iter::Chain<
                                //                 OncePeekable<<Mapping as IndicesAndCoefficients>::EntryMinor>, 
                                //                 MapByTransform<
                                //                         Simplify<
                                //                                 HitMerge<
                                //                                         Scale<
                                //                                                 std::iter::Chain<
                                //                                                         std::iter::Once<(usize, <Mapping as IndicesAndCoefficients>::Coefficient)>, 
                                //                                                         Cloned<std::slice::Iter<'a, (usize, <Mapping as IndicesAndCoefficients>::Coefficient)>>
                                //                                                     >, 
                                //                                                 usize, 
                                //                                                 RingOperator, 
                                //                                                 <Mapping as IndicesAndCoefficients>::Coefficient
                                //                                             >, 
                                //                                             OrderOperatorByKey<
                                //                                                     usize, 
                                //                                                     <Mapping as IndicesAndCoefficients>::Coefficient, 
                                //                                                     (usize, <Mapping as IndicesAndCoefficients>::Coefficient)
                                //                                                 >
                                //                                     >, 
                                //                                 usize, 
                                //                                 RingOperator, 
                                //                                 <Mapping as IndicesAndCoefficients>::Coefficient
                                //                             >, 
                                //                         <Mapping as IndicesAndCoefficients>::EntryMinor, 
                                //                         RemapEntryTuple<
                                //                                 usize, 
                                //                                 <Mapping as IndicesAndCoefficients>::RowIndex, 
                                //                                 &'a Vec<<Mapping as IndicesAndCoefficients>::RowIndex>, 
                                //                                 <Mapping as IndicesAndCoefficients>::Coefficient, 
                                //                                 <Mapping as IndicesAndCoefficients>::Coefficient, 
                                //                                 IdentityFunction, 
                                //                                 <Mapping as IndicesAndCoefficients>::EntryMinor
                                //                             >
                                //                     >
                                //             >, 
                                //         MapByTransform<
                                //                 Chain<
                                //                         Once<(usize, Mapping::Coefficient)>, 
                                //                         Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
                                //                     >,
                                //                 Mapping::EntryMinor,
                                //                 RemapEntryTuple<
                                //                         usize, 
                                //                         <Mapping as IndicesAndCoefficients>::RowIndex, 
                                //                         &'a Vec<<Mapping as IndicesAndCoefficients>::RowIndex>, 
                                //                         <Mapping as IndicesAndCoefficients>::Coefficient, 
                                //                         <Mapping as IndicesAndCoefficients>::Coefficient, 
                                //                         IdentityFunction, 
                                //                         <Mapping as IndicesAndCoefficients>::EntryMinor
                                //                     >                      
                                //             >          
                                // >;                                        
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;//< Vec< Mapping::EntryMinor > as IntoIterator>::IntoIter;


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
            keymaj: Mapping::RowIndex 
        ) 
        -> 
        Self::ViewMajorAscend
    {
        // a struct that realizes the function sending (ordmaj: usize, coefficient) --> t: Mapping::EntryMinor,
        // where t has 
        //   - key equal to the ordmaj'th major key in the matching matrix
        //   - val equal to coefficient
        let entry_changer = RemapEntryTuple::new(
                                    self.umatch.matching_ref().bimap_maj_ref().ord_to_val_vec(),
                                    IdentityFunction{},
                                );

        match self.umatch.matching_ref().keymaj_to_ord( & keymaj ) {
            // unmatched row index
            None => {
                // define the matrix A that will fit into an equation xA = b
                let seed_with_integer_indexed_rows = self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_ordmaj();
                // let seed_with_integer_indexed_rows_ref = & seed_with_integer_indexed_rows;

                // define the problem vector b that will fit into an equation xA = b
                let mapping_matched_cols_only = self.umatch.mapping_matched_cols_only();
                let problem_vector 
                    = mapping_matched_cols_only
                        .view_major_ascend( keymaj.clone() )
                        .negate( self.umatch.ring_operator() );  // recall that there is a MINUS SIGN in the formula in Theorem 6 of Hang et al., "Umatch factorization: ..."

                // Struct that encodes the matching from minor keys of A to major key ordinals of A
                let matching_ref = self.umatch.matching_ref();
                let keymin_to_ordmaj = | keymin: Mapping::ColIndex | -> Option< usize > { matching_ref.keymin_to_ord( &keymin ) };
                let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

                // Solve xA = b for x
                let solution_vec = 
                EchelonSolverMajorAscendWithMinorKeys::solve( // solution to xA = b
                        problem_vector, // mabrix b
                        seed_with_integer_indexed_rows, // matrix A
                        keymin_to_ordmaj_wrapped,
                        self.umatch.ring_operator.clone(),
                        self.umatch.order_operator_major.clone(),
                    )
                    .quotient();
                let solution_vec_integer_indexed 
                    = ChangeIndexSimple::new( solution_vec, self.umatch.matching_ref().bimap_min_ref().val_to_ord_hashmap() );
                    
                // Multiply the solution x by matrix R_{\rho \rho}^{-1}
                let comb_codomain_inv_matched_block = self.umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal();

                // Compute the portion of desired matrix view indexed by matched major keys -- however, this vector will actually have usize indices which we must change next
                let vec_with_matched_indices 
                    = vector_matrix_multiply_major_ascend_simplified( 
                            solution_vec_integer_indexed, 
                            comb_codomain_inv_matched_block, 
                            self.umatch.ring_operator(),
                            OrderOperatorByKey::new(),
                        );

                // reindex the preceding vector so that it is indeed indexed by major keys
                let vec_with_matched_indices_reindexed 
                    = MapByTransform::new( 
                            vec_with_matched_indices,
                            entry_changer.clone(),             
                        );
                
                // Compute the portion of desired matrix view indexed by unmatched major keys -- this portion of the vector has exactly one nonzero entry, and it is the first entry in the vector, because it's the diagonal element of a row of a triangular matrix
                let vec_with_unmatched_index = OncePeekable::new(  Mapping::EntryMinor::new(keymaj, RingOperator::one() )  ); 

                // Merge the two portions of the vector together, to form a whole
                let merged  
                    =   vec_with_unmatched_index.chain( vec_with_matched_indices_reindexed );                  
                
                return  CombCodomainInvViewMajorAscend::new( IterTwoType::Iter1( merged ) ) //CombCodomainInvViewMajorAscend{ iter_unwrapped: ChangeEntryType::new( IterTwoType::Iter1( merged ) ) }
            }
            // matched row index
            Some( ordmaj ) => {
                let reindexed_iter = 
                    MapByTransform::new( 
                            self.umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal().view_major_ascend( ordmaj ),
                            entry_changer,
                        );
                return  CombCodomainInvViewMajorAscend::new( IterTwoType::Iter2( reindexed_iter ) ) // CombCodomainInvViewMajorAscend{ iter_unwrapped: ChangeEntryType::new( IterTwoType::Iter2( reindexed_iter ) ) }
            }
        }
        




        // ON 2022/04/29 THIS CODE SEEMS TO BE UNCESSSARY/DEPRECATED; CONSIDER DELETING
        // // define the matrix A that will fit into an equation xA = b
        // let seed_with_integer_indexed_rows = self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_ordmaj();

        // // define the problem vector b that will fit into an equation xA = b
        // let mapping_matched_cols_only = self.umatch.mapping_matched_cols_only();
        // let problem_vector = mapping_matched_cols_only.view_major_ascend( keymaj.clone() );

        // // Struct that encodes the matching from minor keys of A to major key ordinals of A
        // let matching_ref = self.umatch.matching_ref();
        // let keymin_to_ordmaj = | keymin: Mapping::ColIndex | -> usize { matching_ref.keymin_to_ordmaj( &keymin ).unwrap() };
        // let keymin_to_ordmaj_wrapped = EvaluateFunctionFnMutWrapper::new( keymin_to_ordmaj );        

        // // Solve xA = b for x
        // let mut solution_vec = 
        // EchelonSolverMajorAscendWithMinorKeys::solve( // solution to xA = b
        //         problem_vector, // mabrix b
        //         && seed_with_integer_indexed_rows, // matrix A
        //         keymin_to_ordmaj_wrapped,
        //         self.umatch.ring_operator.clone(),
        //         self.umatch.order_operator_major.clone(),
        //     )
        //     .collect_vec();

        // // Sort solution vector x
        // let mut order_operator_major_clone = self.umatch.order_operator_major.clone(); // we have to clone the order comparator in order to compare order, since the `lt` method takes `& mut self` as an argument -- but to get a mutable reference in the present context we would have to borrow the umatch struct as mutable
        // let compare_order_for_vector_sort = | x: &Mapping::EntryMajor, y: &Mapping::EntryMajor |-> Ordering {
        //     match order_operator_major_clone.judge_lt( x, y ) {
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
        //                     self.umatch.matching.keymin_to_keymaj( &item.key() ).unwrap(),
        //                     item.val(),
        //                 )
        //             )
        //         .collect_vec();
        // let solution_vec_reindexed_and_sorted = IterWrappedVec::new( solution_vec_reindexed );                            

        // // Possibly append an entry with coefficient 1
        // match self.umatch.matching_ref().contains_keymaj( & keymaj ){
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




#[derive(Clone, Dissolve, new)]
pub struct CombCodomainInvViewMajorAscend
            < 'a, Mapping, RingOperator, >

    where 
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::Coefficient:                   Clone,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        // OrderOperatorRowEntries:    Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        // OrderOperatorColEntries:   Clone + JudgePartialOrder< Mapping::EntryMinor >,
        Mapping::ViewMajorAscend:          IntoIterator,        
        Mapping::EntryMajor:     KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,

{
    iter_unwrapped:
    IterTwoType<
            std::iter::Chain<
                    OncePeekable<<Mapping as IndicesAndCoefficients>::EntryMinor>, 
                    MapByTransform<
                            Simplify<
                                    HitMerge<
                                            Scale<
                                                    std::iter::Chain<
                                                            std::iter::Once<(usize, <Mapping as IndicesAndCoefficients>::Coefficient)>, 
                                                            Cloned<std::slice::Iter<'a, (usize, <Mapping as IndicesAndCoefficients>::Coefficient)>>
                                                        >, 
                                                    usize, 
                                                    RingOperator, 
                                                    <Mapping as IndicesAndCoefficients>::Coefficient
                                                >, 
                                                OrderOperatorByKey<
                                                        usize, 
                                                        <Mapping as IndicesAndCoefficients>::Coefficient, 
                                                        (usize, <Mapping as IndicesAndCoefficients>::Coefficient)
                                                    >
                                        >, 
                                    usize, 
                                    RingOperator, 
                                    <Mapping as IndicesAndCoefficients>::Coefficient
                                >, 
                            <Mapping as IndicesAndCoefficients>::EntryMinor, 
                            RemapEntryTuple<
                                    usize, 
                                    <Mapping as IndicesAndCoefficients>::RowIndex, 
                                    &'a Vec<<Mapping as IndicesAndCoefficients>::RowIndex>, 
                                    <Mapping as IndicesAndCoefficients>::Coefficient, 
                                    <Mapping as IndicesAndCoefficients>::Coefficient, 
                                    IdentityFunction, 
                                    <Mapping as IndicesAndCoefficients>::EntryMinor
                                >
                        >
                >, 
            MapByTransform<
                    Chain<
                            Once<(usize, Mapping::Coefficient)>, 
                            Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
                        >,
                    Mapping::EntryMinor,
                    RemapEntryTuple<
                            usize, 
                            <Mapping as IndicesAndCoefficients>::RowIndex, 
                            &'a Vec<<Mapping as IndicesAndCoefficients>::RowIndex>, 
                            <Mapping as IndicesAndCoefficients>::Coefficient, 
                            <Mapping as IndicesAndCoefficients>::Coefficient, 
                            IdentityFunction, 
                            <Mapping as IndicesAndCoefficients>::EntryMinor
                        >                      
                >          
    >
}    

//     iter_unwrapped:
//     IterTwoType<
//         std::iter::Chain<
//             OncePeekable<
//                 <Mapping as IndicesAndCoefficients>::EntryMinor
//             >,
//             MapByTransform<
//                 Simplify<
//                     HitMerge<
//                         Scale<
//                             std::iter::Chain<
//                                 std::iter::Once<
//                                     (usize, <Mapping as IndicesAndCoefficients>::Coefficient)
//                                 >,
//                                 Cloned<
//                                     std::slice::Iter<'_, (usize, <Mapping as IndicesAndCoefficients>::Coefficient)>
//                                 >
//                             >,
//                             usize,
//                             RingOperator,
//                             <Mapping as IndicesAndCoefficients>::Coefficient
//                         >,
//                         OrderOperatorByKey<
//                             usize,
//                             <Mapping as IndicesAndCoefficients>::Coefficient,
//                             (usize, <Mapping as IndicesAndCoefficients>::Coefficient)
//                         >
//                     >,
//                     usize,
//                     RingOperator,
//                     <Mapping as IndicesAndCoefficients>::Coefficient
//                 >
//             >,
//             <Mapping as IndicesAndCoefficients>::EntryMinor,
//             RemapEntryTuple<
//                 usize,
//                 <Mapping as IndicesAndCoefficients>::RowIndex,
//                 &Vec<<Mapping as IndicesAndCoefficients>::RowIndex>,
//                 <Mapping as IndicesAndCoefficients>::Coefficient,
//                 <Mapping as IndicesAndCoefficients>::Coefficient,
//                 IdentityFunction,
//                 <Mapping as IndicesAndCoefficients>::EntryMinor
//             >
//         >,
//         _
//     >    
        // IterTwoType<
        // Chain<
        //         OncePeekable<<Mapping as IndicesAndCoefficients>::EntryMinor>, 
        //         MapByTransform<
        //                 Simplify<
        //                         HitMerge<
        //                                 Scale<
        //                                         std::iter::Chain<
        //                                                 std::iter::Once<(usize, <Mapping as IndicesAndCoefficients>::Coefficient)>, 
        //                                                 Cloned<
        //                                                         std::slice::Iter<'_, (usize, <Mapping as IndicesAndCoefficients>::Coefficient)>
        //                                                     >
        //                                             >, 
        //                                         usize, 
        //                                         RingOperator, 
        //                                         <Mapping as IndicesAndCoefficients>::Coefficient
        //                                     >, 
        //                                 OrderOperatorByKey<
        //                                         usize, 
        //                                         <Mapping as IndicesAndCoefficients>::Coefficient, 
        //                                         (usize, <Mapping as IndicesAndCoefficients>::Coefficient)
        //                                     >
        //                             >, 
        //                             usize, 
        //                             RingOperator, 
        //                             <Mapping as IndicesAndCoefficients>::Coefficient
        //                         >, 
        //                     <Mapping as IndicesAndCoefficients>::EntryMinor, 
        //                     RemapEntryTuple<
        //                             usize, 
        //                             <Mapping as IndicesAndCoefficients>::RowIndex, 
        //                             &Vec<<Mapping as IndicesAndCoefficients>::RowIndex
        //                         >, <Mapping as IndicesAndCoefficients>::Coefficient, <Mapping as IndicesAndCoefficients>::Coefficient, IdentityFunction, <Mapping as IndicesAndCoefficients>::EntryMinor>>>, _>
// }
//     iter_unwrapped:     ChangeEntryType<   
//                                 IterTwoType< 
//                                         MergeTwoItersByOrderOperator<
//                                                 Peekable<
//                                                         ChangeIndexSimple<
//                                                                 LinearCombinationSimplified
//                                                                     <Chain<Once<(usize, Mapping::Coefficient)>, Cloned<Iter<'a, (usize, Mapping::Coefficient)>>>, usize, Mapping::Coefficient, RingOperator, OrderOperatorByKey<usize, Mapping::Coefficient, (usize, Mapping::Coefficient)>>, 
//                                                                 &'a Vec<Mapping::RowIndex>, 
//                                                                 usize, 
//                                                                 Mapping::RowIndex, 
//                                                                 Mapping::Coefficient
//                                                             >
//                                                     >, 
//                                                 OncePeekable<(Mapping::RowIndex, Mapping::Coefficient)>, 
//                                                 OrderOperatorColEntries
//                                             >,                        
//                                         ChangeIndexSimple<
//                                                 Chain<
//                                                         Once<(usize, Mapping::Coefficient)>, 
//                                                         Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
//                                                     >, 
//                                                 &'a Vec<Mapping::RowIndex>, usize, Mapping::RowIndex, Mapping::Coefficient
//                                             >
//                                     >, 
//                                 Mapping::EntryMinor,
//                                 Mapping::RowIndex, 
//                                 Mapping::Coefficient,    
//                             >,
// }

// impl < 'a, Mapping, RingOperator, OrderOperatorColEntries >

//     CombCodomainInvViewMajorAscend
//             < 'a, Mapping, RingOperator, OrderOperatorColEntries >

//     where 
//         Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq,
//         Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq,
//         Mapping::Coefficient:                   Clone,
//         RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
//         // OrderOperatorRowEntries:    Clone + JudgePartialOrder<  Mapping::EntryMajor >,
//         OrderOperatorColEntries:   Clone + JudgePartialOrder<  ( Mapping::RowIndex, Mapping::Coefficient ) >,
//         Mapping::ViewMajorAscend:          IntoIterator,        
//         Mapping::EntryMajor:     KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
//         Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,
//         Mapping:                           ViewRowAscend + IndicesAndCoefficients,
// {
//     /// Returns the wrapped iterator, and consumes the wrapper.
//     pub fn unwrap_iter( self ) ->     ChangeEntryType<   
//                                     IterTwoType< 
//                                             MergeTwoItersByOrderOperator<
//                                                     Peekable<
//                                                             ChangeIndexSimple<
//                                                                     LinearCombinationSimplified<   
//                                                                             Chain<
//                                                                                     Once<(usize, Mapping::Coefficient)>, 
//                                                                                     Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
//                                                                                 >, 
//                                                                             usize, 
//                                                                             Mapping::Coefficient, 
//                                                                             RingOperator, 
//                                                                             OrderOperatorByKey<
//                                                                                     usize, 
//                                                                                     Mapping::Coefficient, 
//                                                                                     (usize, Mapping::Coefficient)
//                                                                                 >
//                                                                         >, 
//                                                                     &'a Vec<Mapping::RowIndex>, 
//                                                                     usize, 
//                                                                     Mapping::RowIndex, 
//                                                                     Mapping::Coefficient
//                                                                 >
//                                                         >, 
//                                                     OncePeekable<(Mapping::RowIndex, Mapping::Coefficient)>, 
//                                                     OrderOperatorColEntries
//                                                 >,                        
//                                             ChangeIndexSimple<
//                                                     Chain<
//                                                             Once<(usize, Mapping::Coefficient)>, 
//                                                             Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
//                                                         >, 
//                                                     &'a Vec<Mapping::RowIndex>, usize, Mapping::RowIndex, Mapping::Coefficient
//                                                 >
//                                         >, 
//                                     Mapping::EntryMinor,
//                                     Mapping::RowIndex, 
//                                     Mapping::Coefficient,    
//                                 > 
//     {
//         self.iter_unwrapped
//     }
// }


impl < 'a, Mapping, RingOperator, >

    Iterator for

    CombCodomainInvViewMajorAscend
            < 'a, Mapping, RingOperator, >

    where 
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq,
        Mapping::Coefficient:                   Clone,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        // OrderOperatorRowEntries:    Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        // OrderOperatorColEntries:   Clone + JudgePartialOrder<  ( Mapping::RowIndex, Mapping::Coefficient ) >,
        Mapping::ViewMajorAscend:          IntoIterator,        
        Mapping::EntryMajor:     KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,            
{
    type Item = Mapping::EntryMinor;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}       


// type CombCodomainInvViewMajorAscend
//         < 'a, Mapping, RingOperator, OrderOperatorColEntries >
//     =
//         IterTwoType< 
//                 Chain<
//                         Once< Mapping::EntryMinor >, 
//                         RemapEntryTuple<
//                                 LinearCombinationSimplified<   
//                                         Chain<
//                                                 Once< Mapping::EntryMinor >, 
//                                                 MapByTransform<
//                                                         Cloned<Iter<'a, (usize, Mapping::Coefficient)>>,
//                                                         Mapping::EntryMinor,
//                                                         RemapEntryTuple<
//                                                                 usize,
//                                                                 Mapping::RowIndex,
//                                                                 &'a Vec< Matrix::KeyMaj >,                                                                        
//                                                                 Mapping::Coefficient,
//                                                                 Mapping::Coefficient,
//                                                                 IdentityFunction,
//                                                                 Mapping::EntryMinor,
//                                                             >
//                                                     >
                                            
//                                             >, 
//                                             Mapping::IndexMajor, 
//                                         Mapping::Coefficient, 
//                                         RingOperator, 
//                                         OrderOperatorByKey<
//                                                 usize, 
//                                                 Mapping::Coefficient, 
//                                                 (usize, Mapping::Coefficient)
//                                             >
//                                     >, 
//                                 &'a Vec<Mapping::RowIndex>, 
//                                 usize, 
//                                 Mapping::RowIndex, 
//                                 Mapping::Coefficient
//                             >
//                         OncePeekable<(Mapping::RowIndex, Mapping::Coefficient)>, 
//                         OrderOperatorColEntries
//                     >,                        
//                 ChangeIndexSimple<
//                         Chain<
//                                 Once<(usize, Mapping::Coefficient)>, 
//                                 Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
//                             >, 
//                         &'a Vec<Mapping::RowIndex>, usize, Mapping::RowIndex, Mapping::Coefficient
//                     >
//             >;

            // ChangeEntryType<   
            //                                     IterTwoType< 
            //                                             MergeTwoItersByOrderOperator<
            //                                                     Peekable<
            //                                                             ChangeIndexSimple<
            //                                                                     LinearCombinationSimplified<   
            //                                                                             Chain<
            //                                                                                     Once<(usize, Mapping::Coefficient)>, 
            //                                                                                     Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
            //                                                                                 >, 
            //                                                                             usize, 
            //                                                                             Mapping::Coefficient, 
            //                                                                             RingOperator, 
            //                                                                             OrderOperatorByKey<
            //                                                                                     usize, 
            //                                                                                     Mapping::Coefficient, 
            //                                                                                     (usize, Mapping::Coefficient)
            //                                                                                 >
            //                                                                         >, 
            //                                                                     &'a Vec<Mapping::RowIndex>, 
            //                                                                     usize, 
            //                                                                     Mapping::RowIndex, 
            //                                                                     Mapping::Coefficient
            //                                                                 >
            //                                                         >, 
            //                                                     OncePeekable<(Mapping::RowIndex, Mapping::Coefficient)>, 
            //                                                     OrderOperatorColEntries
            //                                                 >,                        
            //                                             ChangeIndexSimple<
            //                                                     Chain<
            //                                                             Once<(usize, Mapping::Coefficient)>, 
            //                                                             Cloned<Iter<'a, (usize, Mapping::Coefficient)>>
            //                                                         >, 
            //                                                     &'a Vec<Mapping::RowIndex>, usize, Mapping::RowIndex, Mapping::Coefficient
            //                                                 >
            //                                         >, 
            //                                     Mapping::EntryMinor,
            //                                     Mapping::RowIndex, 
            //                                     Mapping::Coefficient,    
            //                                 >             



//  IMPLEMENT ORACLE MINOR DESCEND FOR COMB: CODOMAIN INV
//  ---------------------------------------------------------------------------------------------------------



impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    ViewColDescend for

    CombCodomainInv
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    
    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,    
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::ViewMinorDescendIntoIter: Iterator,        
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient> + KeyValSet< Mapping::RowIndex, Mapping::Coefficient> + KeyValNew< Mapping::RowIndex, Mapping::Coefficient>,
        OrderOperatorColEntries:   Clone + JudgePartialOrder< Mapping::EntryMinor >,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries:    Clone + JudgePartialOrder< Mapping::EntryMajor >,
        
{
    type ViewMinorDescend                   =   CombCodomainInvViewMinorDescend
                                                    < 'a, Mapping, RingOperator, OrderOperatorColEntries, >;
    type ViewMinorDescendIntoIter           =   Self::ViewMinorDescend;

    fn view_minor_descend
    ( 
        &self, 
        keymaj: Mapping::RowIndex
    ) 
    -> 
    Self::ViewMinorDescend
    {
        match self.umatch.matching.contains_keymaj( &keymaj ) {
            true => {
                // get the ordinal of the matched major key
                let ordmaj = self.umatch.matching.keymaj_to_ord( &keymaj ).unwrap();
                // a column, c, of (R_{\rho \rho}^{-1})  [THIS NOTATION COMES FROM INNER IDENTITIES, C.F. UMATCH FACTORIZATION]
                // we have to reindex this column with minor keys, then sort it according to the order on minor keys, because our next step will be to run a triangular solve, and the triangular matrix is only triangular with respect to the order on minor keys                                
                let mut comb_codomain_inv_matched_block_col     = 
                    self.umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal()
                        .view_minor_descend( ordmaj )
                        .map(   |(key,val)|   // here we convert indices to minor keys
                                Mapping::EntryMajor::new(
                                        self.umatch.matching.ord_to_keymin(key), 
                                        val
                                    )  
                            )
                        .collect_vec(); // collect entries into a vector to simplify sorting
                // sort the entries of the column
                // let cmp_style_comparator = InferTotalOrderFromJudgePartialOrder::new( // we have to construct an order comparator that spits out values of type Ord rather than type bool
                //         self.umatch.order_operator_major() 
                //     );    
                comb_codomain_inv_matched_block_col.sort_by( |x,y| self.umatch.order_operator_major.judge_partial_cmp( x, y ).unwrap() );
                let A = self.umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin();
                // A^{-1} * c
                let echelon_solution    =   TriangularSolverMinorDescend::solve(
                                                    comb_codomain_inv_matched_block_col.into_iter(),
                                                    A,
                                                    self.umatch.ring_operator(),
                                                    self.umatch.order_operator_major(),
                                                );
                let echelon_solution_minus = echelon_solution.negate( self.umatch.ring_operator() );
                let unmatched_part_of_solution = vector_matrix_multiply_minor_descend_simplified(
                        echelon_solution_minus,
                        self.umatch.mapping_matchless_rows_only(),
                        self.umatch.ring_operator(),
                        self.umatch.order_operator_minor(),
                    );

                // this equals the variable comb_codomain_inv_matched_block_col defined above; we need a second copy, and calling the constructor a second time avoids the need to impose "Clone" on the struct
                let matched_part_of_solution     = 
                    self.umatch.comb_codomain_inv_matched_block_indexed_by_keymaj()
                        .view_minor_descend( keymaj );
                      
                // let   matched_part_of_solution = IterTwoType::Iter1( matched_part_of_solution );
                // let unmatched_part_of_solution = IterTwoType::Iter2( unmatched_part_of_solution );

                let solution = MergeTwoItersByOrderOperator::new(
                            matched_part_of_solution.peekable(),
                            unmatched_part_of_solution.peekable(),
                            self.umatch.order_operator_minor_reverse(),
                );
                let solution = ChangeEntryType::new( IterTwoType::Iter1( solution ) );

                // wrap the output in a struct that has a simpler type signature
                CombCodomainInvViewMinorDescend{ iter_unwrapped: solution }
                
            }
            false => {
                let solution = IterTwoType::Iter2( 
                        std::iter::once( 
                                Self::EntryMinor::new( keymaj, RingOperator::one()  ) 
                            ),
                    );
                let solution = ChangeEntryType::new( solution );
                // wrap the output in a struct that has a simpler type signature                    
                CombCodomainInvViewMinorDescend{ iter_unwrapped: solution }
            }
        }
    }
}


pub struct CombCodomainInvViewMinorDescend
                < 'a, Mapping, RingOperator, OrderOperatorColEntries, >
    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,    
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::ViewMinorDescendIntoIter: Iterator,        
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient> + KeyValSet< Mapping::RowIndex, Mapping::Coefficient> + KeyValNew< Mapping::RowIndex, Mapping::Coefficient>,
        OrderOperatorColEntries:   Clone + JudgePartialOrder< Mapping::EntryMinor >,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        // OrderOperatorRowEntries:    Clone + JudgePartialOrder< Mapping::EntryMajor >,                
{
    iter_unwrapped: 
                    // <   EntryIter, EntryNew, Index, RingElement, >
                    ChangeEntryType<
                            IterTwoType<
                                    MergeTwoItersByOrderOperator<
                                            Peekable<
                                                    MapByTransform<
                                                            Chain<
                                                                    Once<(usize, Mapping::Coefficient)>, 
                                                                    // VecOfVecMatrixColumnReverse<'a, usize, Mapping::Coefficient>
                                                                    Cloned<Rev<std::slice::Iter<'a, (usize, <Mapping as IndicesAndCoefficients>::Coefficient)>>>,
                                                                >, 
                                                            Mapping::EntryMinor, 
                                                            ReindexEntry<(usize, Mapping::Coefficient), Mapping::EntryMinor, usize, Mapping::RowIndex, Mapping::Coefficient, &'a Vec<Mapping::RowIndex>>
                                                        >,                                                                
                                                >,
                                            Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection<
                                                                Mapping::ViewMinorDescendIntoIter, 
                                                                &'a HashMap< Mapping::RowIndex, usize>, 
                                                                Mapping::RowIndex, 
                                                                Mapping::Coefficient, 
                                                            >,
                                                            Mapping::RowIndex, 
                                                            Mapping::Coefficient, 
                                                            RingOperator, 
                                                            ReverseOrder<OrderOperatorColEntries>
                                                    >,                                                                                                                             
                                                >,                                                        
                                            ReverseOrder< 
                                                    OrderOperatorColEntries
                                                >,                                                                
                                        >,
                                    std::iter::Once< Mapping::EntryMinor >,
                                >,
                            Mapping::EntryMinor,
                            Mapping::RowIndex,
                            Mapping::Coefficient,
                        >,
}

impl < 'a, Mapping, RingOperator, OrderOperatorColEntries, >

    CombCodomainInvViewMinorDescend
                < 'a, Mapping, RingOperator, OrderOperatorColEntries, >
    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,    
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::ViewMinorDescendIntoIter: Iterator,        
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient> + KeyValSet< Mapping::RowIndex, Mapping::Coefficient> + KeyValNew< Mapping::RowIndex, Mapping::Coefficient>,
        OrderOperatorColEntries:   Clone + JudgePartialOrder< Mapping::EntryMinor >,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        // OrderOperatorRowEntries:    Clone + JudgePartialOrder< Mapping::EntryMajor >,                
{
    pub fn unwrap_iter( self ) -> 
                    ChangeEntryType<
                            IterTwoType<
                                    MergeTwoItersByOrderOperator<
                                            Peekable<
                                                    MapByTransform<
                                                            Chain<
                                                                    Once<(usize, Mapping::Coefficient)>, 
                                                                    // VecOfVecMatrixColumnReverse<'a, usize, Mapping::Coefficient>
                                                                    Cloned<Rev<std::slice::Iter<'a, (usize, <Mapping as IndicesAndCoefficients>::Coefficient)>>>,
                                                                >, 
                                                            Mapping::EntryMinor, 
                                                            ReindexEntry<(usize, Mapping::Coefficient), Mapping::EntryMinor, usize, Mapping::RowIndex, Mapping::Coefficient, &'a Vec<Mapping::RowIndex>>
                                                        >,                                                                
                                                >,
                                            Peekable< 
                                                LinearCombinationSimplified<
                                                        OnlyIndicesOutsideCollection<
                                                                Mapping::ViewMinorDescendIntoIter, 
                                                                &'a HashMap< Mapping::RowIndex, usize>, 
                                                                Mapping::RowIndex, 
                                                                Mapping::Coefficient, 
                                                            >,
                                                            Mapping::RowIndex, 
                                                            Mapping::Coefficient, 
                                                            RingOperator, 
                                                            ReverseOrder<OrderOperatorColEntries>
                                                    >,                                                                                                                             
                                                >,                                                        
                                            ReverseOrder< 
                                                    OrderOperatorColEntries
                                                >,                                                                
                                        >,
                                    std::iter::Once< Mapping::EntryMinor >,
                                >,
                            Mapping::EntryMinor,
                            Mapping::RowIndex,
                            Mapping::Coefficient,
                        >
    {
        self.iter_unwrapped
    }    
}



impl < 'a, Mapping, RingOperator, OrderOperatorColEntries, >

    Iterator for

    CombCodomainInvViewMinorDescend
                < 'a, Mapping, RingOperator, OrderOperatorColEntries, >
    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,    
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     Clone + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + KeyValNew< Mapping::ColIndex, Mapping::Coefficient >, // KeyValNew is required because at one point we have to solver a triangualr system of equations where indices are minor keys
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::ViewMinorDescendIntoIter: Iterator,        
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient> + KeyValSet< Mapping::RowIndex, Mapping::Coefficient> + KeyValNew< Mapping::RowIndex, Mapping::Coefficient>,
        OrderOperatorColEntries:   Clone + JudgePartialOrder< Mapping::EntryMinor >,
        RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        // OrderOperatorRowEntries:    Clone + JudgePartialOrder< Mapping::EntryMajor >,                
{
    type Item = Mapping::EntryMinor;

    fn next( &mut self ) -> Option< Self::Item > {
        self.iter_unwrapped.next()
    }
}


//  =========================================================================================================
//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY)
//  =========================================================================================================





/// Represents the matrix A defined on p22 of "U-match factorization, ..." by Hang, Ziegelmeier, Giusti, and Henselman-Petrusek.
/// Concretely, this matrix equals (the pivot block of the inverse of the codomain COMB) * (the pivot block of the matching array).
#[derive(Copy, Clone, Dissolve)]
pub struct CombCodomainInvTimesMappingMatchedBlock< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where     
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    pub umatch_ref:     & 'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > ,  
}

impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    
    CombCodomainInvTimesMappingMatchedBlock
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,    
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
        RingOperator:                           Clone,
        OrderOperatorRowEntries:    Clone,
        OrderOperatorColEntries:   Clone,        
{
    // Make a new [`CombCodomainInvTimesMappingMatchedBlock`].
    pub fn new( umatch_ref: &'a Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  ) -> Self {
        CombCodomainInvTimesMappingMatchedBlock{ umatch_ref }
    }
}


//  INDICES AND COEFFICIENTS
//  -------------------------------------------------------------------------------------------------------------- 


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >

    IndicesAndCoefficients for

    &'a CombCodomainInvTimesMappingMatchedBlock
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >         

    where
        Mapping:                           ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,            

{   
    type EntryMajor = Mapping::EntryMajor;
    type EntryMinor = Mapping::EntryMinor;
    type ColIndex = Mapping::ColIndex; 
    type RowIndex = Mapping::RowIndex; 
    type Coefficient = Mapping::Coefficient; 
} 



//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT ViewRowAscend FOR < KEYMAJ, KEYMIN >
//  -------------------------------------------------------------------------------------------------------------- 



impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    ViewRowAscend for

    &'a CombCodomainInvTimesMappingMatchedBlock
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    where   
        Mapping:               ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        RingOperator:                           Clone + Semiring< Mapping::Coefficient >,
        OrderOperatorRowEntries:    Clone + JudgePartialOrder<  Mapping::EntryMajor >,
        // &'a VecOfVec<usize, Mapping::Coefficient>: ViewRowAscend<usize, Cloned<Iter<'a, (usize, Mapping::Coefficient)>> >,
        // &'a VecOfVec<usize, Mapping::Coefficient>: ViewRowAscend<usize, Cloned<std::slice::Iter<'a, (Mapping::ColIndex, Mapping::Coefficient)>>>,
        // OrderOperatorRowEntries: 'a, RingOperator: 'a, Mapping::Coefficient: 'a, Mapping::RowIndex: 'a, Mapping::ColIndex: 'a, Mapping::ViewMajorAscend: 'a,            

{
    type ViewMajorAscend            =   LinearCombinationSimplified< 
                                                OnlyIndicesInsideCollection< Mapping::ViewMajorAscendIntoIter, &'a HashMap<Mapping::ColIndex, usize>, Mapping::ColIndex, Mapping::Coefficient>, 
                                                Mapping::ColIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries 
                                            >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;


    fn view_major_ascend( &self, keymaj: Mapping::RowIndex ) -> Self::ViewMajorAscend
        {
        
        // define a row vector
        // NB: we have to translate the index `keymaj` which has type `Mapping::RowIndex` into the index `ordmaj` which has type `usize`, because `self.umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal` is indexed by unsigned integers 
        let ordmaj  =   self.umatch_ref.matching.keymaj_to_ord( &keymaj ).unwrap();

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



//  MATCHED BLOCK OF (CODOMAIN COMB INV * MAPPING ARRAY) -- IMPLEMENT ViewColDescend FOR < KEYMAJ, KEYMIN >
//  --------------------------------------------------------------------------------------------------------------


/// A wrapper struct for descending minor views of `CombCodomainInvTimesMappingMatchedBlock` 
pub struct  CombCodomainInvTimesMappingMatchedBlockViewMinorDescend
                < 'a, Mapping, RingOperator, OrderOperatorColEntries >
    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        Mapping::ViewMajorAscend:          IntoIterator, // this is required by virtually everything that involves a umatch struct
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >, // this is required by virtually everything that involves a umatch struct        
        RingOperator:                           Clone + Semiring< Mapping::Coefficient >,
        OrderOperatorColEntries:                   Clone + JudgePartialOrder<  Mapping::EntryMinor >,                
            { 
                view_minor_descend_of_comb_codomain_inv_times_mapping_matched_block:
                    LinearCombinationSimplified<
                            MapByTransform<
                                    Chain<
                                            Once<(usize, Mapping::Coefficient)>, 
                                            // VecOfVecMatrixColumnReverse<
                                            //         'a,
                                            //         usize,
                                            //         Mapping::Coefficient
                                            //     >
                                            Cloned<
                                                    Rev<
                                                            std::slice::Iter<
                                                                    'a, 
                                                                    (usize, <Mapping as IndicesAndCoefficients>::Coefficient)
                                                                >
                                                        >
                                                >,                                            
                                        >, 
                                    Mapping::EntryMinor, 
                                    ReindexEntry<
                                            (usize, Mapping::Coefficient), 
                                            Mapping::EntryMinor, 
                                            usize, 
                                            Mapping::RowIndex, 
                                            Mapping::Coefficient, 
                                            &'a Vec<Mapping::RowIndex>
                                        >
                                >, 
                                Mapping::RowIndex, 
                                Mapping::Coefficient, 
                                RingOperator, 
                                ReverseOrder<OrderOperatorColEntries>,
                        >,           
            }

impl < 'a, Mapping, RingOperator, OrderOperatorColEntries >

    Iterator for 

    CombCodomainInvTimesMappingMatchedBlockViewMinorDescend
        < 'a, Mapping, RingOperator, OrderOperatorColEntries >
    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        Mapping::ViewMajorAscend:          IntoIterator, // this is required by virtually everything that involves a umatch struct
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >, // this is required by virtually everything that involves a umatch struct        
        RingOperator:                           Clone + Semiring< Mapping::Coefficient >,
        OrderOperatorColEntries:                   Clone + JudgePartialOrder<  Mapping::EntryMinor >, 
{
    type Item = Mapping::EntryMinor;

    fn next( &mut self ) -> Option< Self::Item > {
        self.view_minor_descend_of_comb_codomain_inv_times_mapping_matched_block.next()
    }
}


impl < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    ViewColDescend for

    &'a CombCodomainInvTimesMappingMatchedBlock
        < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 

    where   
        Mapping:                           ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::Coefficient:                   Clone,
        Mapping::ViewMinorDescend:         IntoIterator,
        Mapping::EntryMinor:    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct        
        Mapping::ViewMajorAscend:          IntoIterator, // this is required by virtually everything that involves a umatch struct
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >, // this is required by virtually everything that involves a umatch struct        
        RingOperator:                           Clone + Semiring< Mapping::Coefficient >,
        OrderOperatorColEntries:    Clone + JudgePartialOrder<  Mapping::EntryMinor >,
{

    type ViewMinorDescend           =   CombCodomainInvTimesMappingMatchedBlockViewMinorDescend
                                            < 'a, Mapping, RingOperator, OrderOperatorColEntries >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;


    fn   view_minor_descend( &self, keymin: Self::ColIndex ) -> Self::ViewMinorDescend {

        let matched_row_submatrix_of_mapping_array = self.umatch_ref.mapping_matched_rows_only();
        let comb_codomain_inv_matched_block = self.umatch_ref.comb_codomain_inv_matched_block_indexed_by_keymaj();

        let column      =   matched_row_submatrix_of_mapping_array.view_minor_descend( keymin );
        
        let column_new  = vector_matrix_multiply_minor_descend_simplified(
                    column,
                    comb_codomain_inv_matched_block,
                    self.umatch_ref.ring_operator.clone(),
                    self.umatch_ref.order_operator_minor.clone(),
                );

        CombCodomainInvTimesMappingMatchedBlockViewMinorDescend{ 
                view_minor_descend_of_comb_codomain_inv_times_mapping_matched_block: column_new 
            }

        // LinearCombinationSimplified<
        //         MapByTransform<
        //                 Chain<
        //                         Once<(usize, Mapping::Coefficient)>, 
        //                         VecOfVecMatrixColumnReverse<
        //                                 usize,
        //                                 Mapping::Coefficient
        //                             >
        //                     >, 
        //                 Mapping::EntryMinor, 
        //                 ReindexEntry<
        //                         (usize, Coefficient), 
        //                         Mapping::EntryMinor, 
        //                         usize, 
        //                         Mapping::RowIndex, 
        //                         Mapping::Coefficient, 
        //                         &Vec<Mapping::RowIndex>
        //                     >
        //             >, 
        //             Mapping::RowIndex, 
        //             Mapping::Coefficient, 
        //             RingOperator, 
        //             ReverseOrder<OrderOperatorColEntries>,
        //     >     
        
    }
        
}

