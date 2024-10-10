use oat_rust::{topology::simplicial::{simplices::filtered::SimplexFiltered, from::graph_weighted::{ChainComplexVrFiltered}}, utilities::order::OrderOperatorAuto};        
use oat_rust::topology::point_cloud::unit_circle;    

use oat_rust::algebra::vectors::entries::KeyValGet;
use oat_rust::algebra::rings::operator_structs::ring_native::{FieldRational64};    
use oat_rust::algebra::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
use oat_rust::algebra::matrices::query::{ViewColDescend, ViewRowAscend};
use oat_rust::algebra::matrices::operations::umatch::row_major::{Umatch};    
use oat_rust::algebra::matrices::types::third_party::IntoCSR;

use oat_rust::utilities::order::{ is_sorted_strictly, OrderOperatorByKeyCutsom, };
use oat_rust::utilities::distances::{rowwise_distances};
use oat_rust::utilities::iterators::general::minmax;

use ordered_float::OrderedFloat; 
use itertools::Itertools;  
use std::sync::Arc;


fn main() {

      
    let npoints = 200;
    let dimension_max = 2;

    let pcloud = unit_circle( npoints, Some(-1.0 .. 1.0));

    let dissimilarity_matrix_data
        = rowwise_distances(pcloud)
            .into_iter()
            .map(|x| x.into_iter().enumerate().collect_vec() )
            .collect_vec()
            .into_csr( npoints, npoints );
    let dissimilarity_matrix = & dissimilarity_matrix_data;

    
    let dissimilarity_value_min = OrderedFloat(0.0);
    let dissimilarity_value_max = 
        minmax( 
                (0..npoints).map(
                        |x| 
                        dissimilarity_matrix.view_major_ascend(x).into_iter().map(
                                |x| 
                                x.val()
                            ) 
                    ) 
            ).unwrap_or( dissimilarity_value_min.clone() ); 

    let ring_operator = FieldRational64::new();
    let chain_complex_data = ChainComplexVrFiltered::new( & dissimilarity_matrix, npoints, dissimilarity_value_max, dissimilarity_value_min, ring_operator );
    // let chain_complex_ref = & chain_complex; 
    // let chain_complex = ChainComplexVrFilteredArc{ arc: Arc::new(chain_complex_data) };   
    let chain_complex = Arc::new(chain_complex_data);           
    let keymaj_vec = chain_complex.cliques_in_order(dimension_max);
    let keymin_vec = chain_complex.cliques_in_order(dimension_max+1);


    // verify_viewmajorascend_compatible_with_viewminordescend(
    //         chain_complex.clone(),
    //         keymin_vec.iter().cloned(),
    //         keymaj_vec.iter().cloned(),
    //     );        

    let iter_keymaj = keymaj_vec.iter().cloned();    

    // println!("check that oracle has strictly sorted rows");
    // // print_indexed_major_views( & chain_complex_ref, iter_keymaj.clone() );  // print the major views       
    // for keymaj in iter_keymaj.clone() {        
    //     assert!( is_sorted_strictly( 
    //                                     & chain_complex.view_major_ascend(keymaj.clone()).collect_vec() , 
    //                                     & OrderOperatorByKeyCutsom::new( OrderOperatorAuto ) 
    //                                 ) );
    // }

    // println!("press enter to continue");
    // let mut guess = String::new();
    // io::stdin()
    //     .read_line(&mut guess)
    //     .expect("Failed to read line");        

    // println!("check that oracle has strictly sorted columns");
    // // print_indexed_minor_views( & chain_complex_ref, iter_keymaj.clone() );  // print the major views        
    // for keymaj in iter_keymaj.clone() {
    //     assert!( is_sorted_strictly(    & chain_complex.view_minor_descend(keymaj).iter().cloned().collect_vec() , 
    //                                     & OrderOperatorByKeyCutsom::new( OrderOperatorAutoReverse::new() )  // NOTE THAT HERE WE USE GT 
    //                                 ) );
    // }    

    // println!("press enter to continue");
    // let mut guess = String::new();
    // io::stdin()
    //     .read_line(&mut guess)
    //     .expect("Failed to read line");       

    println!("starting umatch");
    let umatch = Umatch::factor_with_clearing(
            chain_complex.clone(), 
            iter_keymaj.clone(), 
            ring_operator, 
            OrderOperatorAuto, 
            OrderOperatorAuto, 
        );      
          
}
