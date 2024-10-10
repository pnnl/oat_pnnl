//! Calculating Euclidean distances

use ordered_float::OrderedFloat;








/// Return the Euclidean pairwise distances between rows, in
/// a square vec-of-vec matrix.
pub fn rowwise_distances( mat: Vec<Vec<f64>> ) -> Vec< Vec< OrderedFloat<f64> > > {
    
    let nrows = mat.len();
    let mut dist = vec![ vec![ OrderedFloat(0f64);nrows]; nrows ];
    let _l: f64;
    
    for rowind in 0 .. nrows {

        let rowvec = &mat[rowind];
        for colind in rowind .. nrows {
            let colvec  =   &mat[colind];
            let l   =   rowvec
                            .iter().zip( colvec.iter() )
                            .map(|(x,y)| (x-y)*(x-y) )
                            .sum::<f64>()
                            .sqrt();
            let l       =   OrderedFloat(l);
            dist[ rowind ][ colind ] = l;
            dist[ colind ][ rowind ] = l;            
        }        
    }
    dist
}







