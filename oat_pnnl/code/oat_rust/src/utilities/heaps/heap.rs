//! Min-heaps with user-defined comparison functions.
//!
//! In this context a **heap** means a binary heap realized as a 0-indexed vector.
//!
//! There are already nice rust packages for working with heaps.  This package
//! is different because
//!
//! * it lets the user define their own comparison function `less_than(x, y)`
//! * the heaps are 0-indexed
//!     * this is helpful when your heap contains complex objects, so 
//!     filling in the 0-position with a "filler" object is hard.
//!
//! 




use num_traits::identities::Zero;
use std::iter::IntoIterator;




//  -----------------------------------------------------------------------------
//  UTILITY FUNCTIONS 
//  -----------------------------------------------------------------------------


/// The index of the first node in a binary tree that shares the same row as
/// the input node.
///
/// # Examples
///
/// ```
/// // The binary tree: 
/// //                     0
/// //             1               2
/// //         3       4       5       6
/// //       7  8    9  10   11 12   13 14  
/// use oat_rust::utilities::heaps::heap::{min_node_same_row};
/// assert_eq!( 1, min_node_same_row( &2  ) );
/// assert_eq!( 3, min_node_same_row( &4  ) );
/// assert_eq!( 7, min_node_same_row( &11 ) );
/// ```
pub fn min_node_same_row( m : &usize ) -> usize {
    let n = m.clone();
    let mut k = 0;
    let mut exp = 1;

    while n >= k + exp {
        k += exp;
        exp *= 2;
    }
    k
}

/// The parent node in a 0-indexed heap; returns 0 if the input is 0.
///
/// This is used for bulk updates.
///
/// # Examples
///
/// ```
/// // The binary tree: 
/// //                     0
/// //             1               2
/// //         3       4       5       6
/// //       7  8    9  10   11 12   13 14  
///
/// use oat_rust::utilities::heaps::heap::parent_or_0;
///
/// assert_eq!( 1, parent_or_0( &4 ) );
/// assert_eq!( 0, parent_or_0( &0 ) );
/// ```
pub fn parent_or_0( m : &usize ) -> usize {
    let n = m.clone();
    match n.is_zero() {
        true  => 0,
        false => (n-1)/2 
    }
}

/// The parent node in a 0-indexed heap.  Returns `None` if input is 0.
///
/// # Examples
///
/// ```rust
/// // The binary tree: 
/// //                     0
/// //             1               2
/// //         3       4       5       6
/// //       7  8    9  10   11 12   13 14  
///
/// use oat_rust::utilities::heaps::heap::parent;
///
/// assert_eq!( Some(1), parent( &4 ) );
/// assert_eq!( None,    parent( &0 ) );
/// ```
pub fn parent( m : &usize ) -> Option<usize> {
    let n = m.clone();
    match n.is_zero() {
        true  => None,
        false => Some( (n-1)/2 )
    }
}

/// Index of the first child node in a 0-indexed binary heap.
///
/// # Examples
///
/// ```rust
/// // The binary tree: 
/// //                     0
/// //             1               2
/// //         3       4       5       6
/// //       7  8    9  10   11 12   13 14  
///
/// use oat_rust::utilities::heaps::heap::child_a;
///
/// assert_eq!( 11, child_a( &5 ) );
/// assert_eq!( 13, child_a( &6 ) );
/// ```
pub fn child_a( m : &usize ) -> usize { 2 * m + 1 }

/// Index of the second child node in a 0-indexed binary heap.
///
/// # Examples
///
/// ```rust
/// // The binary tree: 
/// //                     0
/// //             1               2
/// //         3       4       5       6
/// //       7  8    9  10   11 12   13 14  
///
/// use oat_rust::utilities::heaps::heap::child_b;
///
/// assert_eq!( 4, child_b( &1 ) );
/// assert_eq!( 8, child_b( &3 ) );
/// ```
pub fn child_b( m : &usize ) -> usize { 2 * m + 2 }


//  -----------------------------------------------------------------------------
//  SIFT
//  -----------------------------------------------------------------------------


//  SIFT DOWN
//
/// Sift down element at `index` (`heap` is a min-heap wrt the ordering)
pub fn sift_down<T, S>(heap: &mut [T], index: usize, mut less_than: S)
    where S: FnMut(&T, &T) -> bool
{
    debug_assert!(index <= heap.len());
    let mut pos = index;
    let mut child = child_a( &pos );
    // the `pos` conditional is to avoid a bounds check
    while pos < heap.len() && child < heap.len() {
        let right = child + 1;

        // pick the smaller of the two children
        if right < heap.len() && less_than(&heap[right], &heap[child]) {
            child = right;
        }

        // sift down is done if we are already in order
        if !less_than(&heap[child], &heap[pos]) {
            return;
        }
        heap.swap(pos, child);
        pos = child;
        child = child_a( &pos );
    }
}


//  -----------------------------------------------------------------------------
//  HEAPIFY (ALL / TAIL) 
//  -----------------------------------------------------------------------------


//  HEAPIFY WHOLE VECTOR 
//
/// Heapify (min-heap) with respect to `less_than`.
pub fn heapify<T, S>(data: &mut [T], mut less_than: S)
    where S: FnMut(&T, &T) -> bool
{
    if data.is_empty() {
        return; // there is nothing to do
    }
    for i in (0..data.len() / 2).rev() {
        sift_down(data, i, &mut less_than);
    }
}

//  HEAPIFY TAIL 
//
/// Heapify (min-heap) with respect to `less_than`, assuming that 
/// `data[0..tail_base]` (excluding `tail_base`) is alraedy a heap.
/// 
/// This method is based off Figure 3 of Elmasry + Katajainen: Towards ultimate
/// binary heaps.    
pub fn heapify_tail<T, S>(data: &mut [T], mut less_than: S, tail_base: &usize)
    where S: FnMut(&T, &T) -> bool
{
    if tail_base >= &data.len() {
        // in this case the vector is already heapified; nothing to do
        return;
    }
    if data.len() < 2 {
        // in this case the vector is already heapified; nothing to do
        return;
    }

    let mut right = data.len()-1;
    let mut left = std::cmp::max( tail_base.clone(), parent( &right ).unwrap() );
    
    while right != 0 {
        left  = parent_or_0( &left  ); // leftmost parent
        right = parent_or_0( &right ); // rightmost parent
        for i in ( left..( right + 1 ) ).rev() {
            sift_down( data, i, &mut less_than );
        }
    }
}

//  ---------------------------------------------------------------------------
//  INSERT
//  ---------------------------------------------------------------------------

/// Add new elements to a heapified vector, and heapify.
///
/// Heapification is performed not with [`heapify`] but with [`heapify_tail`],
/// which is based on off Figure 3 of Elmasry + Katajainen: Towards ultimate 
/// binary heaps.
///
pub fn bulk_insert< I, F >  (   heap: &mut Vec< <I as IntoIterator>::Item >, 
                                less_than: F, 
                                iter: I   
                            ) 
    where   I: IntoIterator,
            F: FnMut(&<I as IntoIterator>::Item,
                     &<I as IntoIterator>::Item) -> bool
{
    let base = heap.len();
    heap.extend( iter.into_iter() );
    heapify_tail( heap, less_than, &base);
}


//  ---------------------------------------------------------------------------
//  POP
//  ---------------------------------------------------------------------------


/// Pop the top (meaning smallest) element from the heap; replace this with
/// the last element and sift down, producing a new heap with 1 less element.
///
/// # Examples
///
/// ```
/// use oat_rust::utilities::heaps::heap::pop;
///
/// let mut vec = vec![ 0, 3, 5, 4 ];
/// let top = pop( &mut vec, |p, q| &p < &q ).unwrap(); // recall that <unwrap> is a function to extract x from Some(x)
/// assert_eq!( 0,   top );
/// assert_eq!( vec, vec![ 3, 4, 5] );
///
/// ```
pub fn pop< T, F> ( heap: &mut Vec <T>, less_than: F ) -> Option<T>
    where F: FnMut( &T, &T) -> bool
{
    if heap.len().is_zero() { return None }

    let val = heap.swap_remove( 0 );

    if heap.len() > 0 {
        sift_down( heap, 0, less_than );
    }

    return Some( val )
}


//  ---------------------------------------------------------------------------
//  HEAP ITERATOR 
//  ---------------------------------------------------------------------------


/// A struct to iteratively pop elements off a heap.
///
/// Assuming `heap` is a min-heap with respect to `less_than`, this iterator
/// will return elements in (ascending, wrt `less_than`) sorted order.
struct HeapIterator< T, F >
    where F: FnMut( &T, &T) -> bool
{
    heap: Vec< T >,
    less_than: F
}

impl< T, F > Iterator for HeapIterator< T, F >
    where F: FnMut( &T, &T) -> bool
{
   type Item = T;

   fn next( &mut self) -> Option< T > {
        let less_than = &mut self.less_than;
        pop( &mut self.heap, |p, q| less_than( &p, &q) )
   }
}


//  -----------------------------------------------------------------------------
//  TESTS 
//  -----------------------------------------------------------------------------

use rand::{Rng};

/// Generate a vector of length `n` entries drawn from the uniform distribution
/// on 0, ..., k-1.
pub fn randgen_n_of_k( n: usize, k: usize) -> Vec<usize> {
    let mut rng = rand::thread_rng();
    let v : Vec<usize> = (0..n).map(|_| rng.gen_range(0..k)).collect();
    return v
}

/// True iff the input vector is already heapified with respect to `less_than`.
pub fn is_heapified< T, F >( vec: Vec<T>, mut less_than: F ) -> bool 
    where   T: PartialOrd + Clone + std::fmt::Debug, 
            F: FnMut( &T, &T) -> bool
{
    let mut heapified = true;
    for pos in (1..vec.len()).rev() {
        let par = parent( & pos ).unwrap();
        if less_than( & vec[pos], & vec[par] ) {
            heapified = false;
            println!("{:?}", vec.clone());
            println!("{:?}", pos.clone());
            println!("{:?}", par.clone());
            break;
        }
    }
    heapified
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heap_functions() {

        let n = 10;
        for _ in 0..n {
                    
            let vec = randgen_n_of_k( n, n/2);

            let precedes = |p: &usize, q: &usize| p < q;

            for a in 0..n {
                for b in a..n {
                    // excise two parts of a vector
                    let veca : Vec<usize> = ( 0..a ).map(|x| vec[x].clone()).collect();
                    let vecb : Vec<usize> = ( a..b ).map(|x| vec[x].clone()).collect();

                    // heapify the fist part
                    // ----------------------
                    
                    let mut heap_a = veca.clone();
                    heapify( &mut heap_a, precedes );
                    assert!( is_heapified( heap_a.clone(), precedes ) );

                    // append the second, and heapify
                    // -------------------------------
                    
                    let mut heap_b = heap_a.clone();
                    heap_b.extend( vecb.clone() );
                    heapify_tail( &mut heap_b, precedes, &a );
                    assert!( is_heapified( heap_b.clone(), precedes ) );
                    
                    // alternatively, bulk-insert the tail
                    // -----------------------------------
                    
                    let mut heap_c = heap_a.clone();
                    bulk_insert( &mut heap_c, precedes, vecb.clone() );
                    assert!( is_heapified( heap_c.clone(), precedes ) );


                    // check that popping out elements returns a sorted vector
                    // -------------------------------------------------------

                    let sorted_a : Vec<usize> = HeapIterator{ heap: heap_a.clone(), less_than: |p, q| &p < &q}.collect();
                    let sorted_b : Vec<usize> = HeapIterator{ heap: heap_b.clone(), less_than: |p, q| &p < &q}.collect();
                    let sorted_c : Vec<usize> = HeapIterator{ heap: heap_c.clone(), less_than: |p, q| &p < &q}.collect();
                    
                    let mut sorted_true_a : Vec<usize> = (0..a).map(|x| vec[x].clone()).collect();
                    let mut sorted_true_b : Vec<usize> = (0..b).map(|x| vec[x].clone()).collect();
                    sorted_true_a.sort();
                    sorted_true_b.sort();

                    assert_eq!( sorted_a, sorted_true_a.clone() );
                    assert_eq!( sorted_b, sorted_true_b.clone() );
                    assert_eq!( sorted_c, sorted_true_b.clone() );
                }
            }
        }
    }

}
