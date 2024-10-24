//! Random object generators


use itertools::Itertools;
use ordered_float::OrderedFloat;
use rand::Rng;



/// A sequence of combinations sampled with given probability.
/// 
/// For `k` in `cardinalities` we generate the sequence of all cardinality-`k` subsets of  `{0.. size_ambient_set}`, in lexicographic order.
/// Each element of the outer sequence is pushed to a storage vector with probability `probability`.  The storage vector is returned, when done.
pub fn rand_sequences< T: Iterator<Item=usize> >( size_ambient_set: usize, cardinalities: T, probability: f64 ) -> Vec< Vec< usize > > 
{
    use rand::distributions::{Bernoulli, Distribution};

    let d = Bernoulli::new(  probability ).unwrap();
    let mut rng = rand::thread_rng();    
    
    let mut combinations    =   vec![];
    for cardinality in cardinalities {
        for combo in (0 .. size_ambient_set ).combinations( cardinality ) {
            if d.sample(&mut rng) {
                combinations.push( combo )
            }
        }
    }  

    combinations

}



/// Generates an m x n matrix with entries sampled from the uniform 
/// distribution on [0,1]
pub fn random_matrix( m: usize, n:usize ) -> Vec< Vec< OrderedFloat<f64> > > {
    let mut rng = rand::thread_rng();

    let mut matrix = vec![ Vec::with_capacity(n); m ];
    for rowind in 0..m {
        for _ in 0..n {
            matrix[rowind].push( OrderedFloat(rng.gen_range(0.0 .. 1.0)) )
        }
    }
    matrix
}

/// Generates an m x n symmetric matrix of form (D+D')/2, where D is
/// a random matrix with entries sampled from the uniform  distribution on [0,1]
pub fn random_symmetric_matrix( m: usize ) -> Vec< Vec< OrderedFloat<f64> > > {

    let mut matrix = random_matrix(m,m);
    for p in 0..m {
        for q in p..m {
            let t = (matrix[p][q] + matrix[q][p]) / 2.0 ;
            matrix[p][q] = t;
            matrix[q][p] = t;
        }        
    }
    matrix
}

/// Generates an m x n symmetric matrix of form (D+D')/2, where D is
/// a random matrix with entries sampled from the uniform  distribution on [0,1]
pub fn random_symmetric_matrix_zero_diag( m: usize ) -> Vec< Vec< OrderedFloat<f64> > > {

    let mut matrix = random_matrix(m,m);
    for p in 0..m {
        for q in p..m {
            let t = (matrix[p][q] + matrix[q][p]) / 2.0 ;
            matrix[p][q] = t;
            matrix[q][p] = t;
        }        
    }
    for p in 0 .. m { matrix[p][p] = OrderedFloat(0.0) }
    matrix
}