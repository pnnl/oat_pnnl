//! Binary search on sorted lists
//! 
//! Most operations in this module are unsafe, because they assume (without checking) that inputs are sorted.

use std::cmp::Ordering;




/// Find a value `p` in `{ min, .., max }` that meets a given criterion, or return `None`.
/// 
/// The returned value is **not** gauranteed to be the first one that satisfies the criterion.
/// 
/// We implicitly assume that there exists an ambient sequence of form
///  
/// `v = [ false, false, .., true, true, ..., false, false, ..]`
/// 
/// where at least one  `true` value exists, and all `true` values are contiguous
/// (however the `true` values may lie outside the interval `{ min, .., max }`).
/// We assume that `search_direction(p)` is equal to `Ordering::Equal` if `v[p]==true`,
/// is equal to `Ordering::Less` if the first true value occurs below `p`, and 
/// is equal to `Ordering::Greater` if the first true value occurs above `p`
/// 
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::binary_search::find_sorted_binary_oracle;
/// 
/// let v = vec![ (3,1.0), (5,1.0) ];
/// 
/// assert_eq!( find_sorted_binary_oracle(0,1, |p| 1.cmp( & v[p as usize].0 ) ), None    ); // search for an entry with index 1
/// assert_eq!( find_sorted_binary_oracle(0,1, |p| 3.cmp( & v[p as usize].0 ) ), Some(0) ); // search for an entry with index 3
/// assert_eq!( find_sorted_binary_oracle(0,1, |p| 5.cmp( & v[p as usize].0 ) ), Some(1) ); // search for an entry with index 5
/// assert_eq!( find_sorted_binary_oracle(0,1, |p| 7.cmp( & v[p as usize].0 ) ), None    ); // search for an entry with index 7
/// ```
/// 
/// # Notes
/// This code is unit-tested on all 0-1 sparse vectors of length < 8; see source code for details.
pub fn find_sorted_binary_oracle< F >( 
        mut min:            isize, 
        mut max:            isize, 
        search_direction:   F 
    ) -> 
    Option< isize > 
    where
        F: Fn(isize) -> Ordering

{  
let mut mid;
if max < min { return None }    
while min <= max {
    mid = (min + max)/2;
    match search_direction(mid) {
        Ordering::Equal     =>  { return Some(mid) },
        Ordering::Less      =>  { max = mid - 1; },
        Ordering::Greater   =>  { min = mid + 1; },
    }
}
None
}


/// Find an entry with index `n` in a *strictly sorted* sequence of index-value pairs.
/// 
/// The input, `sparsevec`, is a Rust vector of form `[ (i0,v0), (i1,v1), .. ]` where 
/// `i0 ≤ i1 ≤ ..`.
/// 
/// The output is either
/// - an index `a` such that `sparsevec[a][0] = `n`, or
/// - `None`, if no such index exists
/// 
/// The value `a` is **not** guaranteed to be the first index where `n` appears.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::binary_search::find_sorted_binary_tuple;
/// 
/// let sparsevec = vec![ (3,1.0), (5,1.0) ];
/// 
/// assert_eq!( find_sorted_binary_tuple(&sparsevec, 1), None      );
/// assert_eq!( find_sorted_binary_tuple(&sparsevec, 3), Some( 0 ) );
/// assert_eq!( find_sorted_binary_tuple(&sparsevec, 5), Some( 1 ) );
/// assert_eq!( find_sorted_binary_tuple(&sparsevec, 7), None      );
/// ```
/// 
/// # Notes
/// This code is unit-tested on all 0-1 sparse vectors of length < 8; see source code for details.
pub fn find_sorted_binary_tuple<T>( 
    sparsevec:      & Vec< (usize,T) >, 
    n:              usize 
) -> Option< usize > 
{  

    let search_direction = |p: isize| n.cmp( & sparsevec[p as usize].0 );
    let a = 0;
    let b = (sparsevec.len() as isize) -1 ;
    find_sorted_binary_oracle(a, b, search_direction).map(|x| x as usize)

    // let mut low = 0;
    // let mut hig = sparsevec.len();
    // let mut mid = hig / 2;
    // let mut hit = sparsevec[mid].0;
    // loop {
    //     match n.cmp( &hit ) {
    //         Ordering::Equal     => {  },
    //         Ordering::Less      => {
    //             hig = (low + hig) / 2;
    //             mid = (low + mid) / 2;
    //         },
    //         Ordering::Greater  => {
    //             low = (low + hig) / 2;
    //             mid = (mid + hig) / 2;
    //         },            
    //     }
    //     hit = sparsevec[mid].0;
    //     if hit == n { return Some( &sparsevec[mid] ) }        
    //     if low == mid { return None }
    // }
}



/// Find an entry with value `n` in a *(not necessarily strictly) sorted* `sequence`.
/// 
/// Look-up is performed by a binary search.
/// 
/// The output is either
/// - an index `a` such that `sparsevec[a] = `n`, or
/// - `None`, if no such index exists
/// 
/// The index `a` is **not** guaranteed to be the first index where `n` appears.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::binary_search::find_in_sorted_sequence;
/// 
/// let sequence = vec![ 3, 5 ];
/// 
/// assert_eq!( find_in_sorted_sequence(&sequence, & 1), None      );
/// assert_eq!( find_in_sorted_sequence(&sequence, & 3), Some( 0 ) );
/// assert_eq!( find_in_sorted_sequence(&sequence, & 5), Some( 1 ) );
/// assert_eq!( find_in_sorted_sequence(&sequence, & 7), None      );
/// ```
/// 
/// # Notes
/// This code is unit-tested on all 0-1 sparse vectors of length < 8; see source code for details.
pub fn find_in_sorted_sequence< T: Ord >( 
            sequence:       & Vec< T >, 
            n:              & T
        ) -> 
        Option< usize > 
{  
    let search_direction = |p: isize| n.cmp( & sequence[p as usize] );
    let a = 0;
    let b = (sequence.len() as isize) -1 ;
    find_sorted_binary_oracle(a, b, search_direction).map(|x| x as usize)
}


/// Find an entry with value `n` in a user-specified subsequence of a *(not necessarily strictly) sorted* `sequence`.
/// 
/// Look-up is performed by a binary search.
/// 
/// The output is either
/// - an index `a` such that `lower ≤ a < upper` and `sparsevec[a] = `n`, or
/// - `None`, if no such index exists
/// 
/// The `lower` bound defaults to 0 and the `upper` bound defaults to `sequence.len()`.
/// 
/// The index `a` is **not** guaranteed to be the first index where `n` appears.
/// 
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::binary_search::find_in_sorted_sequence_within_bounds;
/// 
/// let sequence    =   vec![ 3, 5 ];
/// let lower       =   Some(1);
/// let upper       =   Some(1);
/// 
/// assert_eq!( find_in_sorted_sequence_within_bounds(&sequence, & 3, lower, None  ), None      );
/// assert_eq!( find_in_sorted_sequence_within_bounds(&sequence, & 5, lower, None  ), Some( 1 ) );
/// assert_eq!( find_in_sorted_sequence_within_bounds(&sequence, & 3, None,  upper ), Some( 0 ) );
/// assert_eq!( find_in_sorted_sequence_within_bounds(&sequence, & 5, None,  upper ), None      );
/// ```
/// 
/// # Notes
/// This code is unit-tested on all 0-1 sparse vectors of length < 8; see source code for details.
pub fn find_in_sorted_sequence_within_bounds< T: Ord >( 
            sequence:       & Vec< T >, 
            n:              & T,
            lower:          Option< usize >,     
            upper:          Option< usize >,
        ) -> 
        Option< usize > 
{      
    let search_direction = |p: isize| n.cmp( & sequence[p as usize] );
    let a = lower.unwrap_or(0) as isize;
    let b = upper.unwrap_or(sequence.len()) as isize - 1;
    find_sorted_binary_oracle(a, b, search_direction).map(|x| x as usize)
}



/// Determine set containment.
/// 
/// *Both input vectors must be sored in (not necessarily strictly) ascending oder*.
/// 
/// ```
/// use oat_rust::utilities::binary_search::contains_subset;
/// 
/// let a: Vec<usize>   =   vec![];
/// let b: Vec<usize>   =   vec![0,1];
/// let c: Vec<usize>   =   vec![1,2];
/// let d: Vec<usize>   =   vec![0,1,2];
/// 
/// assert!(   contains_subset( &b, &a )  );
/// assert!(   contains_subset( &d, &b )  );
/// assert!(   contains_subset( &d, &d )  );
/// assert!( ! contains_subset( &a, &d )  );
/// assert!( ! contains_subset( &b, &c )  );
/// ```
pub fn contains_subset< T: Ord >( 
                superset:   & Vec< T >, 
                subset:     & Vec< T > 
        ) -> bool {
    let mut lower = Some( 0 ); // lower bound for a search
    let upper = Some( superset.len() );
    for element in subset {
        match find_in_sorted_sequence_within_bounds( superset, element, lower, upper) {
            None                    =>  { return false }
            Some( index )       =>  { lower.replace(index); }
        }
    }
    true
}




/// Given an increasing sequence `left_limits` and an integer `pigeon`, find
/// an index `p` such that `left_limits[p] <= pigoen < left_limits[p+1]`.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::binary_search::find_window;
/// 
/// let left_limits = vec![0,2,2,3];
/// 
/// let u: Vec< Option<usize> > = (0..5).map( |x| find_window( &left_limits, x ) ).collect();
/// let v = vec![Some(0),Some(0),Some(2),None,None];
/// 
/// assert_eq!( u, v );
/// ```
pub fn find_window( left_limits: &Vec<usize>, pigeon: usize ) -> Option< usize > {
    let search_direction = |p: isize| -> Ordering {
            let min = left_limits[p as usize];
            let max = left_limits[p as usize +1 ];
            if max <= pigeon { return Ordering::Greater }
            if min >  pigeon { return Ordering::Less }
            Ordering::Equal
        };
    let a = 0;
    let b = ( left_limits.len() as isize ) - 2; // the value couldn't be greater than 2, as we need to index "up" by one to get the upper limit
    find_sorted_binary_oracle(a, b, search_direction).map(|x| x as usize)
}