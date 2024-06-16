use itertools::Itertools;
use sohar_rust::simplicial::from::graph_weighted::data_structures::CliqueBoundaryMatrix;
use sohar_rust::simplicial::from::graph_weighted::user_interface::{CliqueParams, cliques_in_order};
use oat_rust::matrices::matrix_oracle_traits::{OracleMinorDescend, OracleMajorAscend};
use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor_with_clearing, ParetoShortCircuit};
use oat_rust::rings::operator_structs::ring_native::{DivisionRingNative, field_rational_i64};
use sohar_rust::point_cloud::unit_circle;
use oat_rust::utilities::distances::{rowwise_distances, minmax};
use oat_rust::utilities::partial_order::OrderComparatorAutoLt;
use std::f64::INFINITY;
use std::time::{Instant};
use ordered_float::OrderedFloat;


fn main() {

    // MODIFY REDUCTION TO ACCEPT A CLOSURE THAT WILL PRINT ERROR INFO

    // PARAMETERS

    let maxdim = 1;
    let m = 2000; // numpoints
    let n = 2; // dimension
    let ring_operator = field_rational_i64();  
    // let ring_operator = BooleanFieldOperator::new();        

    // POINT CLOUD

    // let pcloud = random_matrix(10, 50);
    let pcloud = unit_circle( m, Some(-1.0 .. 1.0));

    let dismat = rowwise_distances(pcloud);

    let params = CliqueParams::new( &dismat, maxdim, None, ring_operator.clone() );

    let boundary_matrix = CliqueBoundaryMatrix::new( & params );
    let boundary_matrix_ref = & boundary_matrix;

    let now = Instant::now();    
    let keymaj_vec = cliques_in_order( & params );
    println!("time to collect row indices into vec:              {}", now.elapsed().as_millis());
    let iter_keymaj = keymaj_vec.iter().cloned();

    let now = Instant::now();    
    for keymaj in iter_keymaj.clone() { let _ = boundary_matrix_ref.view_major_ascend(keymaj); }
    println!("time to iterate rows:                              {}", now.elapsed().as_millis());    

    // let now = Instant::now();    
    // for keymaj in iter_keymaj.clone() { let _ = boundary_matrix_ref.view_major_ascend(keymaj).collect_vec(); }
    // println!("time to expand all rows:                           {}", now.elapsed().as_millis());  
    
    let now = Instant::now();   
    let mut num_short_circuits = 0; 
    for keymaj in iter_keymaj.clone() { 
        let iter = boundary_matrix_ref.view_major_ascend(keymaj);
        if let Some(t) = iter.pareto_short_circuit() { num_short_circuits += 1; }
        else { let _ = iter.collect_vec(); } 
    }
    println!("time to expand all rows except short circuited:    {}", now.elapsed().as_millis());  

    let now = Instant::now();        
    let umatch = new_umatchrowmajor_with_clearing(
        boundary_matrix_ref, 
        iter_keymaj.clone(), 
        ring_operator.clone(), 
        OrderComparatorAutoLt::new(), 
        OrderComparatorAutoLt::new(), 
    );  
    println!("time to umatch:                                    {}", now.elapsed().as_millis());  
    println!("num rows:           {}", keymaj_vec.len() );  
    println!("num short circuits: {}", num_short_circuits);


    
    // Get the matching array (cf the  paper on umatch factorization)
    let array_matching = umatch.array_matching_ref();
    
    // The set {columns of the domain comb that are not matched upward or downward} 
    // forms a basis for homology
    let mut bettis = vec![0; maxdim+1];

    let mut dim;
    for simplex in iter_keymaj.clone() {        
        if array_matching.contains_keymaj( &simplex ) || array_matching.contains_keymin( &simplex ) 
            { continue }
        dim = simplex.dim();
        bettis[dim] += 1;
    }

    let mut barcode = Vec::new();
    for keymaj in iter_keymaj.clone() {        
        if let Some(keymin) = array_matching.keymaj_to_keymin( &keymaj )  {
            if keymaj.filvalue != keymin.filvalue {
                barcode.push( ( keymaj.dim(), f64::from(keymaj.filvalue), f64::from(keymin.filvalue) )  )
            }
        } else if ! array_matching.contains_keymin( &keymaj ) {
            barcode.push( ( keymaj.dim(), f64::from(keymaj.filvalue), INFINITY )  )
        }
    }    

    let mut bar_lengths = vec![ Vec::new(); maxdim+1];
    for bar in barcode.iter() { bar_lengths[bar.0].push( OrderedFloat(bar.2-bar.1) ) }
    for lengths in bar_lengths.iter_mut() { lengths.sort() }

    println!("bettis: {:?}", &bettis);
    println!("bar lengths dim 1: {:?}",  bar_lengths[1])
    // for dim in 0 .. maxdim+1 {
    //     println!("bar lengths dim {:?}: {:?}", dim, &bar_lengths[dim]);
    // }
    
// ;

    




//     let maxdim = 2;
//     let dowker_simplices: Vec<Vec<usize>> = vec![ (0..100).collect(), vec![0,151], vec![1,151] ];

//     println!("with clearning:");
//     let now = Instant::now();
//     let basis = homology_basis_from_dowker( &dowker_simplices, maxdim, UseClearing::Yes );
//     println!("{}", now.elapsed().as_millis());
//     let bettis: Vec<usize> = (0usize..maxdim+1).map(|x| basis[x].len() ).collect();    
//     println!("{:?}",bettis );

//     // println!("without clearning:");   
//     // let now = Instant::now();    
//     // let basis = homology_basis_from_dowker( &dowker_simplices, maxdim, UseClearing::No );     
//     // println!("{}", now.elapsed().as_millis());    
//     // // println!("{:?}", basis);    
//     // let bettis: Vec<usize> = (0usize..maxdim+1).map(|x| basis[x].len() ).collect();    
//     // println!("{:?}",bettis );    

//     println!("{:?}", random_matrix(2,2));
}    
