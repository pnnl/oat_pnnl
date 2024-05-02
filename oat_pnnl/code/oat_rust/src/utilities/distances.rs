use crate::rings::operator_traits::{Semiring, Ring};

use super::{partial_order::StrictlyLess, iterators::general::PeekUnqualified};




/// Return the Euclidean pairwise distances between rows, in
/// a square vec-of-vec matrix.
pub fn rowwise_distances( mat: Vec<Vec<f64>> ) -> Vec< Vec< f64 > > {
    
    let nrows = mat.len();
    let mut dist = vec![ vec![0f64;nrows]; nrows ];
    let l: f64;
    
    for rowind in 0 .. nrows {

        let rowvec = &mat[rowind];
        for colind in rowind .. nrows {
            let colvec  =   &mat[colind];
            let l   =   rowvec
                            .iter().zip( colvec.iter() )
                            .map(|(x,y)| (x-y)*(x-y) )
                            .sum::<f64>()
                            .sqrt();
            dist[ rowind ][ colind ] = l;
            dist[ colind ][ rowind ] = l;            
        }        
    }
    return dist
}


/// Return the minimum of the maximum of the row vectors
pub fn minmax<T: Ord>( mat: & Vec< Vec< T > > ) -> &T {
    return mat.iter().map(|x| x.iter().max().unwrap() ).min().unwrap()
}





