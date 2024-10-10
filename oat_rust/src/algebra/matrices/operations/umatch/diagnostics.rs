use itertools::Itertools;

use crate::{algebra::matrices::{types::{oracle_ref::OracleRef, vec_of_vec::sorted::VecOfVec}, display::print_indexed_major_views, debug::verify_viewmajorascend_compatible_with_viewminordescend, query::IndicesAndCoefficients}, utilities::{order::{is_sorted_strictly, JudgePartialOrder, OrderOperatorByKey}, iterators::merge::hit::HitMerge}, };
use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
use crate::algebra::rings::operator_traits::{Semiring,Ring,DivisionRing};
use crate::algebra::matrices::random_constructors::random_mod_p_with_density;
use crate::utilities::order::{OrderOperatorAuto, };
use crate::algebra::matrices::debug::verify_that_product_is_identity;
use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
use crate::algebra::matrices::operations::multiply::ProductMatrix;
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend};
use crate::utilities::iterators::is_sorted::IsSortedBy;  
use std::{hash::Hash, iter::Peekable};
use std::fmt::Debug;
use crate::algebra::vectors::operations::Scale;
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet, KeyValNew};
use crate::utilities::order::OrderOperatorByKeyCutsom;


fn test_umatchrowmajor_comprehensive< Mapping, IterKeyMin, IterKeyMaj, RingOperator, OrderOperatorKeyMin, OrderOperatorKeyMaj >
        ( 
            mapping:              Mapping,
            iter_keymin:                IterKeyMin,
            iter_keymaj:                IterKeyMaj,
            ring_operator:              RingOperator,
            order_operator_keymin:    OrderOperatorKeyMin,
            order_operator_keymaj:    OrderOperatorKeyMaj,
        ) 

where   Mapping:                               ViewRowAscend + ViewColDescend + IndicesAndCoefficients, // we require ViewColDescend for the function `new(..)`, but not for the `Umatch` struct itself.
        IterKeyMin:                                 Clone + Iterator < Item = Mapping::ColIndex >,        
        IterKeyMaj:                                 Clone + Iterator < Item = Mapping::RowIndex >,
        Mapping::ColIndex:                       Clone + Hash + std::cmp::Eq + Debug, 
        Mapping::RowIndex:                       Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        Mapping::Coefficient:                       Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        Mapping::ViewMajorAscend:              Clone + IntoIterator, // !!! remove clone        
        Mapping::EntryMajor:         Clone + Debug + KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + KeyValSet< Mapping::ColIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + KeyValSet< Mapping::RowIndex, Mapping::Coefficient >, // KeyValSet is required in order to construct the `Simplify` struct which is wrapped in the LinearCombinationSimplified struct                
        OrderOperatorKeyMin:                      Clone + JudgePartialOrder <  Mapping::ColIndex >, // !!! remove clone
        OrderOperatorKeyMaj:                      Clone + JudgePartialOrder <  Mapping::RowIndex >,
        HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorKeyMin>>: Clone, // !!!! remove this                 
{

      

    // let num_indices_major           =   10;
    // let num_indices_minor           =   20;
    // let approximate_density           =   0.2;
    // let modulus                     =   17;
    // let allow_nonstructural_zero     =   true;

    // let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
    // let mapping_data       =   random_mod_p_with_density( num_indices_major, num_indices_minor, approximate_density, modulus, allow_nonstructural_zero );
    // let mapping                               =   & mapping_data;

    let umatch 
        =   Umatch( 
                    mapping, //:              Mapping, 
                    iter_keymaj,//:                IterKeyMaj,
                    ring_operator,//:              RingOperator,
                    order_operator_keymin,//:    OrderOperatorKeyMin,
                    order_operator_keymaj,//:    OrderOperatorKeyMaj,                
                );
    // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
    let matching = umatch.matching_ref();

    let comb_codomain = umatch.comb_codomain();
    let comb_codomain_inv = umatch.comb_codomain_inv();        
    let comb_domain = umatch.comb_domain();        
    let comb_domain_inv = umatch.comb_domain_inv();  
    let comb_codomain_inv_times_mapping_matched_block = umatch.comb_codomain_inv_times_mapping_matched_block();  
    let comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin = umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin();

    let comb_codomain_ref         =   & comb_codomain;
    let comb_codomain_inv_ref         =   OracleRef::new( & comb_codomain_inv );
    let comb_domain_ref         =   OracleRef::new( & comb_domain );
    let comb_domain_inv_ref         =   OracleRef::new( & comb_domain_inv );            
    let comb_codomain_inv_times_mapping_matched_block_ref     =   & comb_codomain_inv_times_mapping_matched_block;
    let comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref = & comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin;
    

    let product_domain = ProductMatrix::new( comb_domain_ref, comb_domain_inv_ref, ring_operator.clone(), OrderOperatorByKey::new() );
    let product_codomain = ProductMatrix::new( comb_codomain_ref, comb_codomain_inv_ref, ring_operator.clone(), OrderOperatorAuto );        
    let product_codomain_comb_inv_times_mapping = ProductMatrix::new( comb_codomain_inv_ref, mapping, ring_operator.clone(), OrderOperatorAuto );      
    let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrix::new( product_codomain_comb_inv_times_mapping, comb_domain_ref, ring_operator.clone(), OrderOperatorAuto );                        


    println!("mapping:");
    print_indexed_major_views( & mapping, iter_keymaj.clone() );
    println!("matching:");
    print_indexed_major_views( & matching, iter_keymaj.clone() );        
    println!("comb_domain:");
    print_indexed_major_views( & comb_domain, iter_keymin.clone() );        
    println!("comb_domain_inv:");
    print_indexed_major_views( & comb_domain_inv, iter_keymin.clone() );     
    println!("comb_codomain:");
    print_indexed_major_views( & comb_codomain, iter_keymaj.clone() );        
    println!("comb_codomain_inv:");
    print_indexed_major_views( & comb_codomain_inv, iter_keymaj.clone() );                        
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
            matching.view_major_ascend( keymaj ).collect_vec()
        ) 
    }    

    // check that the major and minor views of `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` agree (i.e. that, taken all together, they run over the same entries)     
    verify_viewmajorascend_compatible_with_viewminordescend(
            OracleRef::new( comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref ),
            umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned(),
            umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned(),                
        );
    
    // check that the major and minor views of `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` are sorted
    for keymin in umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned() {
        assert!(    comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_minor_descend( keymin ).is_sorted_by( |x, y| x.0 > y.0 )     );
    }
    for keymin in umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned() {
        assert!(    comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_major_ascend( keymin ).is_sorted_by( |x, y| x.0 < y.0 )     );
    }     



    // ----------------------------------------------------------------------------------------------------------------
    // check that the major and minor views of `CombCodomain` agree
    verify_viewmajorascend_compatible_with_viewminordescend(
            OracleRef::new( & comb_codomain_ref ),
            iter_keymaj.clone(),
            iter_keymaj.clone(),
        );

    // check that the minor views of `CombCodomain` are sorted
    for keymin in iter_keymaj.clone() {
        let view_minor_descend = comb_codomain_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend.collect_vec(), 
                & umatch.order_operator_minor_reverse() 
            ) );
    }

    // check that the major views of `CombCodomain` are sorted
    for keymin in iter_keymaj.clone() {
        let view_major_ascend = comb_codomain_ref.view_major_ascend( keymin );
        assert!( is_sorted_strictly(  
                & view_major_ascend.collect_vec(), 
                & umatch.order_operator_major() 
            ) );
    }        


    // ----------------------------------------------------------------------------------------------------------------
    // check that the major and minor views of `CombCodomainInv` agree
    verify_viewmajorascend_compatible_with_viewminordescend(
            OracleRef::new( & comb_codomain_inv_ref ),
            iter_keymaj.clone(),
            iter_keymaj.clone(),
        );

    // check that the minor views of `CombCodomainInv` are sorted
    for keymin in iter_keymaj.clone() {
        let view_minor_descend = comb_codomain_inv_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend.collect_vec(), 
                & umatch.order_operator_minor_reverse() 
            ) );
    }

    // check that the major views of `CombCodomainInv` are sorted
    for keymin in iter_keymaj.clone() {
        let view_major_ascend = comb_codomain_inv_ref.view_major_ascend( keymin );
        assert!( is_sorted_strictly(  
                & view_major_ascend.collect_vec(), 
                & umatch.order_operator_major() 
            ) );
    }        

    // ----------------------------------------------------------------------------------------------------------------        
    // check that the major and minor views of `CombDomain` agree
    verify_viewmajorascend_compatible_with_viewminordescend(
            OracleRef::new( & comb_domain_ref ),
            iter_keymin.clone(),
            iter_keymin.clone(),
        );

    // check that the minor views of `CombDomain` are sorted
    for keymin in iter_keymin.clone() {
        let view_minor_descend = comb_domain_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend.collect_vec(), 
                & umatch.order_operator_minor_reverse() 
            ) );
    }  

    // check that the major views of `CombDomain` are sorted
    for keymin in iter_keymin.clone() {
        let view_major_ascend = comb_domain_ref.view_major_ascend( keymin );
        assert!( is_sorted_strictly(  
                & view_major_ascend.collect_vec(), 
                & umatch.order_operator_major() 
            ) );
    }          
    
    // ----------------------------------------------------------------------------------------------------------------        
    // check that the major and minor views of `CombDomainInv` agree
    verify_viewmajorascend_compatible_with_viewminordescend(
            OracleRef::new( & comb_domain_inv_ref ),
            iter_keymin.clone(),
            iter_keymin.clone(),
        );

    // check that the minor views of `CombDomainInv` are sorted
    for keymin in iter_keymin.clone() {
        let view_minor_descend = comb_domain_inv_ref.view_minor_descend( keymin );
        assert!( is_sorted_strictly(  
                & view_minor_descend, 
                & umatch.order_operator_minor_reverse() 
            ) );
    }   
    
    // check that the major views of `CombDomainInv` are sorted
    for keymin in iter_keymin.clone() {
        let view_major_ascend = comb_domain_inv_ref.view_major_ascend( keymin );
        assert!( is_sorted_strictly(  
                & view_major_ascend.collect_vec(), 
                & umatch.order_operator_major() 
            ) );
    }           






    // check that `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` is upper triangular
    for keymin in umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned() {
        assert!(    comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_minor_descend( keymin.clone() ).next().unwrap().0 == keymin     );
    }        


// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:BEGIN        
    // check that columns are sorted in strictly descending order
    //  NOTE: THIS IS UNNECESSARY FOR THE COMBS, SINCE WE TEST THAT THEIR MINOR VIEWS EQUAL THOSE OF VecOfVec objects, WHOSE MINOR DESCENDING VIEWS ARE *ALWAYS* STRICTLY DECREASING IN INDEX
    for keymaj in iter_keymin.clone() { 
        assert!(    mapping.view_minor_descend( keymaj.clone()             ).is_sorted_by( |x, y| x.0 > y.0 )     );
        assert!(    comb_codomain.view_minor_descend( keymaj.clone()       ).is_sorted_by( |x, y| x.0 > y.0 )     );
        // assert!(    comb_codomain_inv.view_minor_descend( keymaj.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
        // assert!(    comb_domain.view_minor_descend( keymaj.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
        // assert!(    comb_domain_inv.view_minor_descend( keymaj.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
    }          
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END  

    
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG: BEGIN   
    // check that the major and minor views of the inverse of the codomain COMB agree
    let comb_codomain_vec_of_vec_simple     
        =   VecOfVec::from_iterable( (iter_keymaj.clone()).map( |k| comb_codomain.view_major_ascend(k) ) );
    for keymaj in iter_keymaj.clone() {
        println!("VIEW MAJOR DESCEND IS STARTING FOR THIS ROUND: keymaj = {:?}", keymaj);
        println!("VIEW MAJOR DESCEND LAZY CONSTRUCTION: {:?}", comb_codomain.view_minor_descend( keymaj ).collect_vec());
        println!("VIEW MAJOR DESCEND FROM MAJOR VIEW: {:?}", (& comb_codomain_vec_of_vec_simple).view_minor_descend( keymaj ).collect_vec());            
        println!("VIEW MAJOR DESCEND IS FINISHED FOR THIS ROUND");
        itertools::assert_equal( 
                comb_codomain.view_minor_descend( keymaj ),
                (& comb_codomain_vec_of_vec_simple).view_minor_descend( keymaj )
            )
    }      
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END      
    
    // check that rows are sorted in strictly ascending order
    for keymaj in iter_keymaj.clone() { 
        assert!(    mapping.view_major_ascend( keymaj.clone()             ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    comb_codomain.view_major_ascend( keymaj.clone()       ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    comb_codomain_inv.view_major_ascend( keymaj.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    comb_domain.view_major_ascend( keymaj.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
        assert!(    comb_domain_inv.view_major_ascend( keymaj.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
    }  
            

}

