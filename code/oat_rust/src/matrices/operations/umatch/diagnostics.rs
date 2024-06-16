use itertools::Itertools;

use crate::{matrices::{matrix_types::{oracle_ref::OracleRef, vec_of_vec::VecOfVecSimple}, display::print_indexed_major_views, debug::verify_viewmajorascend_compatible_with_viewminordescend, matrix_oracle_traits::IndicesAndCoefficients}, utilities::{partial_order::{is_sorted_strictly, StrictlyLess, OrderComparatorAutoLtByKey}, iterators::merge::heap_of_iterators::HitMerge}, };
use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
use crate::rings::operator_traits::{Semiring,Ring,DivisionRing};
use crate::matrices::random_constructors::random_vec_of_vec_simple;
use crate::utilities::partial_order::{OrderComparatorAutoAnyType, OrderComparatorAutoLt};
use crate::matrices::debug::verify_that_product_is_identity;
use crate::matrices::operations::umatch::row_major::{new_umatchrowmajor};
use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;
use crate::matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend};
use crate::utilities::iterators::is_sorted::IsSortedBy;  
use std::{hash::Hash, iter::Peekable};
use std::fmt::Debug;
use crate::vectors::operations::Scale;
use crate::entries::{KeyValGet, KeyValSet, KeyValNew};
use crate::utilities::partial_order::OrderComparatorLtByKey;


fn test_umatchrowmajor_comprehensive< ArrayMapping, IterKeyMin, IterKeyMaj, RingOperator, OrderComparatorKeyMin, OrderComparatorKeyMaj >
        ( 
            array_mapping:              ArrayMapping,
            iter_keymin:                IterKeyMin,
            iter_keymaj:                IterKeyMaj,
            ring_operator:              RingOperator,
            order_comparator_keymin:    OrderComparatorKeyMin,
            order_comparator_keymaj:    OrderComparatorKeyMaj,
        ) 

where   ArrayMapping:                               OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients, // we require OracleMinorDescend for the function `new(..)`, but not for the `UmatchRowMajor` struct itself.
        IterKeyMin:                                 Clone + Iterator < Item = ArrayMapping::KeyMin >,        
        IterKeyMaj:                                 Clone + Iterator < Item = ArrayMapping::KeyMaj >,
        ArrayMapping::KeyMin:                       Clone + Hash + std::cmp::Eq + Debug, 
        ArrayMapping::KeyMaj:                       Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        ArrayMapping::SnzVal:                       Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:                               Clone + Semiring< ArrayMapping::SnzVal > + Ring< ArrayMapping::SnzVal > + DivisionRing< ArrayMapping::SnzVal >,
        ArrayMapping::ViewMajorAscend:              Clone + IntoIterator, // !!! remove clone        
        ArrayMapping::ViewMajorAscendEntry:         Clone + Debug + KeyValGet< ArrayMapping::KeyMin, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMin, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        ArrayMapping::ViewMinorDescend:             IntoIterator,        
        ArrayMapping::ViewMinorDescendEntry:        Clone + KeyValGet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal > + KeyValSet< ArrayMapping::KeyMaj, ArrayMapping::SnzVal >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct                
        OrderComparatorKeyMin:                      Clone + StrictlyLess <  ArrayMapping::KeyMin >, // !!! remove clone
        OrderComparatorKeyMaj:                      Clone + StrictlyLess <  ArrayMapping::KeyMaj >,
        HitMerge<Peekable<Scale<<ArrayMapping::ViewMajorAscend as IntoIterator>::IntoIter, ArrayMapping::KeyMin, RingOperator, ArrayMapping::SnzVal>>, OrderComparatorLtByKey< ArrayMapping::KeyMin, ArrayMapping::SnzVal, ArrayMapping::ViewMajorAscendEntry, OrderComparatorKeyMin>>: Clone, // !!!! remove this                 
{

      

    // let num_indices_major           =   10;
    // let num_indices_minor           =   20;
    // let approximate_density           =   0.2;
    // let modulus                     =   17;
    // let allow_nonstructural_zero     =   true;

    // let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
    // let array_mapping_data       =   random_vec_of_vec_simple( num_indices_major, num_indices_minor, approximate_density, modulus, allow_nonstructural_zero );
    // let array_mapping                               =   & array_mapping_data;

    let umatch 
        =   new_umatchrowmajor( 
                    array_mapping, //:              ArrayMapping, 
                    iter_keymaj,//:                IterKeyMaj,
                    ring_operator,//:              RingOperator,
                    order_comparator_keymin,//:    OrderComparatorKeyMin,
                    order_comparator_keymaj,//:    OrderComparatorKeyMaj,                
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
    

    let product_domain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_domain_ref, array_comb_domain_inv_ref, ring_operator.clone(), OrderComparatorAutoLtByKey::new() );
    let product_codomain = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_ref, array_comb_codomain_inv_ref, ring_operator.clone(), OrderComparatorAutoAnyType );        
    let product_codomain_comb_inv_times_mapping = ProductMatrixLazyMajorAscendSimplified::new( array_comb_codomain_inv_ref, array_mapping, ring_operator.clone(), OrderComparatorAutoAnyType );      
    let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrixLazyMajorAscendSimplified::new( product_codomain_comb_inv_times_mapping, array_comb_domain_ref, ring_operator.clone(), OrderComparatorAutoAnyType );                        


    println!("array_mapping:");
    print_indexed_major_views( & array_mapping, iter_keymaj.clone() );
    println!("array_matching:");
    print_indexed_major_views( & array_matching, iter_keymaj.clone() );        
    println!("comb_domain:");
    print_indexed_major_views( & array_comb_domain, iter_keymin.clone() );        
    println!("comb_domain_inv:");
    print_indexed_major_views( & array_comb_domain_inv, iter_keymin.clone() );     
    println!("comb_codomain:");
    print_indexed_major_views( & array_comb_codomain, iter_keymaj.clone() );        
    println!("comb_codomain_inv:");
    print_indexed_major_views( & array_comb_codomain_inv, iter_keymaj.clone() );                        
    println!("comb_codomain_inv * mapping * comb_domain:");
    print_indexed_major_views( & product_codomain_comb_inv_times_mapping_times_domain_comb, iter_keymaj.clone() );                                        
    for column_index in iter_keymin.clone() {
        println!("{:?}", product_domain.view_major_ascend( column_index ).collect_vec() );
        itertools::assert_equal( product_domain.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
    }



    // check that the product of the domain comb with its inverse is identity: C * C^{-1} = I
    for keymin in iter_keymin.clone() { 
        assert_eq!(
            product_domain.view_major_ascend( keymin ).collect_vec(),
            vec![ (keymin, 1) ]
        ) 
    }

    // check that the product of the codomain comb with its inverse is identity R * R^{-1} = I
    for keymaj in iter_keymaj.clone() { 
        assert_eq!(
            product_codomain.view_major_ascend( keymaj ).collect_vec(),
            vec![ (keymaj, 1) ]
        ) 
    }    
    
    // check the factorization R^{-1} * D * C = M
    for keymaj in iter_keymaj.clone() { 
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
            iter_keymaj.clone(),
            iter_keymaj.clone(),
        );

    // check that the minor views of `CombCodomain` are sorted
    for keymin in iter_keymaj.clone() {
        let view_minor_descend = array_comb_codomain_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend.collect_vec(), 
                & umatch.order_comparator_viewminordescendentry_reverse() 
            ) );
    }

    // check that the major views of `CombCodomain` are sorted
    for keymin in iter_keymaj.clone() {
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
            iter_keymaj.clone(),
            iter_keymaj.clone(),
        );

    // check that the minor views of `CombCodomainInv` are sorted
    for keymin in iter_keymaj.clone() {
        let view_minor_descend = array_comb_codomain_inv_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend.collect_vec(), 
                & umatch.order_comparator_viewminordescendentry_reverse() 
            ) );
    }

    // check that the major views of `CombCodomainInv` are sorted
    for keymin in iter_keymaj.clone() {
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
            iter_keymin.clone(),
            iter_keymin.clone(),
        );

    // check that the minor views of `CombDomain` are sorted
    for keymin in iter_keymin.clone() {
        let view_minor_descend = array_comb_domain_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend.collect_vec(), 
                & umatch.order_comparator_viewminordescendentry_reverse() 
            ) );
    }  

    // check that the major views of `CombDomain` are sorted
    for keymin in iter_keymin.clone() {
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
            iter_keymin.clone(),
            iter_keymin.clone(),
        );

    // check that the minor views of `CombDomainInv` are sorted
    for keymin in iter_keymin.clone() {
        let view_minor_descend = array_comb_domain_inv_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend, 
                & umatch.order_comparator_viewminordescendentry_reverse() 
            ) );
    }   
    
    // check that the major views of `CombDomainInv` are sorted
    for keymin in iter_keymin.clone() {
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
    for keymaj in iter_keymin.clone() { 
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
        =   VecOfVecSimple::from_iterable( (iter_keymaj.clone()).map( |k| array_comb_codomain.view_major_ascend(k) ) );
    for keymaj in iter_keymaj.clone() {
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
    for keymaj in iter_keymaj.clone() { 
        assert!(    array_mapping.view_major_ascend( keymaj.clone()             ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    array_comb_codomain.view_major_ascend( keymaj.clone()       ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    array_comb_codomain_inv.view_major_ascend( keymaj.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    array_comb_domain.view_major_ascend( keymaj.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    array_comb_domain_inv.view_major_ascend( keymaj.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
    }  
            

}

