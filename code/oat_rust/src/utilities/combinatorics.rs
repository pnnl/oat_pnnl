//! Mathematical combinatorics.


// use crate::utilities::index::EndIndex;
use std::iter::FromIterator;



/// SEE BELOW FOR A TEST OF THIS FUNCTION
/// Returns a vector that runs over all sequences with a given sum that respect the
/// given capacity vector.
pub fn  fixed_sum_sequences(
        caps:           & Vec< usize >,
        target_sum:     usize
        )
        ->
        Vec< Vec< usize> >
{
    let cap_aggregate       =   caps.iter().sum();

    // case 0: problem is insoluble
    if target_sum > cap_aggregate { return Vec::with_capacity(0) }

    // case 1: prolbem is trivial because our sequence must have length 0
    else if caps.is_empty() { return vec![ vec![] ] }

    // all remaining cases:
    // we will recursively solve the problem for sequences of 1-shorter length; to do so we
    // make a truncated sequence of caps and several alternative target sums for that 
    // truncated sequence.
    let trunc_caps              =   Vec::from_iter( caps.iter().cloned().take( caps.len() -1 ) );
    let trunc_cap_agg: usize    =   trunc_caps.iter().sum();
    
    // calculate the max and min possible values in the deleted end slot
    let last_min        =   
        match trunc_cap_agg < target_sum {
            true    =>  target_sum - trunc_cap_agg,
            false   =>  0
        };
   
    let mut last_cap    =   caps.last().unwrap().clone();
    if target_sum < last_cap { last_cap = target_sum.clone() }
    
    // make a container to store results
    let mut sequences   =   Vec::new();

    // we must recursively solve the problem for each possible value we will place in the deleted end slot
    for last_val    in   last_min .. last_cap + 1 {
        
        let trunc_target_sum        =   target_sum - last_val;
        let mut trunc_sequences     =   fixed_sum_sequences(
                                            & trunc_caps,
                                            trunc_target_sum,
                                        );
        for trunc_seq in trunc_sequences.iter_mut() { trunc_seq.push( last_val.clone() ) }  // complete each truncated sequence to a full-length sequence
        sequences.append( &mut trunc_sequences ) // collect the results into our growing pool
    }

    sequences

}

/// STARTED BUT NOT FINISHED
/// Counts the number of minimal elements of the iterator.
// pub fn  count_minimal_elements< I, T >( 
//             iter: I 
//         ) 
//         -> 
//         Option< (T, usize) >
//     where   I:  IntoIterator< Item=T >,
//             T:  Ord + Eq + Clone,
// {
//     let mut min_val_opt             =   None;
//     let mut count                   =   0;
//     for item in iter {
//         if min_val_opt.is_none() { min_val_opt = Some( item.clone() ); count = 1 }
//         else {
//             let min_val             =   min_val_opt.unwrap();
//             match item.cmp( & min_val ) {
//                 Less    =>  { min_val_opt = Some( item.clone() ); count = 1 },
//                 Equal   =>  { count += 1 },
//                 Greater =>  {}
//             }               
//         }
//     }

//     if min_val_opt.is_none() { return None }
//     else { return Some( (min_val_opt.unwrap(), count ) ) }
// }            



#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use rand::Rng;
    use itertools::Itertools;


    #[test]
    fn test_fixed_sum_sequences() {

        let max_seq_length              =   4;
        let max_capacity_param_ceiling  =   4;
    
        let mut rng = rand::thread_rng();
        for seq_len in 0 .. max_seq_length {

            for max_capacity in 0 .. max_capacity_param_ceiling{
                let caps        =   Vec::from_iter( (0..seq_len).map(|_x| rng.gen_range( 0 .. max_capacity + 1 )) ); // random sequence of capacitites
                for target_sum in 0 .. caps.iter().sum() {
                    let mut a   =   Vec::from_iter(  
                                        caps.iter()
                                            .map(|x| 0 .. x+1)
                                            .multi_cartesian_product()
                                            .filter(|x| x.iter().sum::<usize>() == target_sum )
                                    );
                    let mut b   =   fixed_sum_sequences( & caps, target_sum );
    
                    a.sort();
                    b.sort();
    
                    assert_eq!( a, b );
                }
            }
        }
    } 


}    