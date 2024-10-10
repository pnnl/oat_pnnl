// use oat_rust::utilities::homology::dowker::{homology_basis_from_dowker, UseClearing};
// use oat_rust::topology::simplicial::from::relation::{homology_basis_from_dowker, UseClearing};
// use oat_rust::algebra::rings::operator_structs::ring_native::FieldRational64;
// use std::time::{Instant};


fn main() {

    // let maxdim = 2;
    // let dowker_simplices: Vec<Vec<usize>> = vec![ (0..200).collect() ]; //, vec![0,151], vec![1,151] ];
    // let ring_operator = FieldRational64::new();

    // println!("with clearning:");
    // let now = Instant::now();
    // let basis = homology_basis_from_dowker( &dowker_simplices, maxdim, ring_operator, UseClearing::Yes );
    // println!("{}", now.elapsed().as_millis());
    // let bettis: Vec<usize> = ( 0usize .. maxdim as usize + 1 ).map(|x| basis[x].len() ).collect();    
    // println!("{:?}",bettis );

    // println!("without clearning:");   
    // let now = Instant::now();    
    // let basis = homology_basis_from_dowker( &dowker_simplices, maxdim, UseClearing::No );     
    // println!("{}", now.elapsed().as_millis());    
    // // println!("{:?}", basis);    
    // let bettis: Vec<usize> = (0usize..maxdim+1).map(|x| basis[x].len() ).collect();    
    // println!("{:?}",bettis );    

}    


// to profile: 


// use std::thread::sleep;

// fn main() {
//    let now = Instant::now();

//    // we sleep for 2 seconds
//    sleep(Duration::new(2, 0));
//    // it prints '2'
//    println!("{}", now.elapsed().as_secs());
// }