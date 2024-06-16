


use itertools::Itertools;
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

    return combinations

}



/// Generates an m x n matrix with entries sampled from the uniform 
/// distribution on [0,1]
pub fn random_matrix( m: usize, n:usize ) -> Vec< Vec< f64 > > {
    let mut rng = rand::thread_rng();

    let mut matrix = vec![ Vec::with_capacity(n); m ];
    for rowind in 0..m {
        matrix[rowind].push( rng.gen_range(0.0 .. 1.0) )
    }
    return matrix
}