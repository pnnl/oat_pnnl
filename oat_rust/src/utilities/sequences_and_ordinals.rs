//! Sequences and ordinals.
//! 
//! A sequence of elements in a set `S` is a function `f{0, .., N} -> S`.  A bijective
//! sequence is called an "enumeration".  The "ordinal" of an element `s` in `S` is the
//! integer `n` such that `f(n) = s`.  This module offers tools to work with
//! sequences, enumerations, and ordinals.

use std::{collections::HashMap, marker::PhantomData};
use std::hash::Hash;
use std::iter::FromIterator;

use serde::{Deserialize, Serialize};

use super::binary_search::{find_in_sorted_sequence, find_in_sorted_sequence_within_bounds, contains_subset};
use super::functions::evaluate::{IdentityFunction, EvaluateFunction, ReferencedEvaluator};
use super::iterators::general::MapByTransformTrait;
use super::order::{is_sorted_strictly, OrderOperatorAuto};

use std::ops::Index;

//  ---------------------------------------------------------------------------
//  STRICTLY SORTED SEQUENCES
//  ---------------------------------------------------------------------------


/// Represents a strictly-increasing sequence of elements.
/// 
/// This object is just a wrapper for a vector; it's primary purpose
/// is to provide a "certificate" that this vector has entries 
/// sorted in strictly asceding order.
#[derive(Clone,Debug,PartialEq,PartialOrd,Eq,Hash,Ord)]
pub struct SortedVec< T: Ord > {
    sequence: Vec< T >
}
 impl < T: Ord > 
    
    SortedVec
        < T > 
{

    /// Construct a new `SortedVec` by wrapping a vector.
    /// 
    /// Panics if the input sequence is not sorted in strictly ascending order.
    pub fn new( sequence: Vec< T > ) -> Result< Self, Vec< T > >  {
        if ! is_sorted_strictly( &sequence, &OrderOperatorAuto ){

            return Err::< Self, Vec< T > >( sequence );
        }
        Ok( SortedVec { sequence } )
    }

    /// Construct a new `SortedVec` from an iterator.
    /// 
    /// Panics if the input sequence is not sorted in strictly ascending order.
    pub fn from_iter< I: IntoIterator< Item=T > >( sequence: I ) -> Self {
        SortedVec { sequence: sequence.into_iter().collect() }
    }

    /// Return the inner sequence, consuming the wrapper.
    pub fn into_vec( self ) -> Vec< T > { self.sequence }

    /// Immutable reference to the internally stored sorted vector
    pub fn vec( &self ) -> & Vec< T > { & self.sequence }

    /// Return the index where `element` appears in the sequence (or `None`, if it does not appear)
    pub fn find( &self, element: &T ) -> Option< usize > {
        find_in_sorted_sequence( self.vec(), element )
    }

    /// Return the index where `element` appears in the sequence (or `None`,
    /// if it does not appear in the index range `[lower, upper)` )
    pub fn find_in_range( &self, element: &T, lower: Option<usize>, upper: Option<usize> ) -> Option< usize > {
        find_in_sorted_sequence_within_bounds( self.vec(), element, lower, upper )
    }    

    /// Length of the sequence.
    pub fn len( &self ) -> usize { self.sequence.len() }

    /// Return the `true` iff the sequence contains `element`.
    pub fn contains( &self, element: &T ) -> bool {
        find_in_sorted_sequence( self.vec(), element ).is_some()
    }    

    /// Return the `true` iff `subset` is a subset of `self`.
    pub fn contains_subset( &self, subset: & SortedVec<T> ) -> bool {
        contains_subset( self.vec(), subset.vec() )
    }  

    /// Return the `true` iff `superset` is a superset of `self`.
    pub fn includes_into( &self, superset: & SortedVec<T> ) -> bool {
        contains_subset( superset.vec(), self.vec() )
    }            

}


impl < T >

    Index< usize > for 
    
    SortedVec
        < T >

    where
        T:  Ord
{
    type Output = T;
    fn index( & self, index: usize ) -> & Self::Output { &self.sequence[ index ] }
}        


//  ---------------------------------------------------------------------------
//  A COUNTERPART TO ITERTOOLS "COMBINATIONS", BUT REVERSE LEXICOGRAPHIC ORDER
//  ---------------------------------------------------------------------------



/// Behaves like `Itertools.Combinations`, but returns combinations in *reverse order*.
/// 
/// See [`CombinationsReverse::new`], [`CombinationsReverse::remap_elements`], [`CombinationsReverse::from_iterable`] for different behaviors.
#[derive(Clone,Debug)]
pub struct CombinationsReverse< T, F > {
    combination:        Vec< usize >,
    num_elements:       usize,
    is_empty:           bool,
    remap:              F,
    phantom_element:    PhantomData<T>,
}


impl CombinationsReverse
        < usize, IdentityFunction > {

    /// Iterates over all sequences of length `sequence_length` sorted in strictly ascending order,
    /// whose elements lie in { 0 .. num_elements } in **reverse lexicographic order**.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::sequences_and_ordinals::CombinationsReverse;
    /// use oat_rust::utilities::functions::evaluate::IdentityFunction;
    /// 
    /// let num_elements        =   4;
    /// let sequence_length     =   2;
    /// let combinations        =   CombinationsReverse::new( num_elements, sequence_length );
    /// 
    /// itertools::assert_equal( combinations, vec![ vec![2,3], vec![1,3], vec![1,2], vec![0,3], vec![0,2], vec![0,1] ] );
    /// ```
    pub fn new( num_elements: usize, sequence_length: usize, ) -> Self {
        let remap   =   IdentityFunction::new();
        if sequence_length <= num_elements {
            let combination: Vec<usize> = ( num_elements - sequence_length .. num_elements ).collect();
            let is_empty = false;
            CombinationsReverse { combination, num_elements, is_empty, remap, phantom_element:PhantomData }
        } else {
            let combination: Vec<usize> = vec![];
            let is_empty = true;
            CombinationsReverse { combination, num_elements, is_empty, remap, phantom_element:PhantomData }            
        }
    }
}


impl < 'a, T >

    CombinationsReverse
        < T, &'a Vec<T> > 

    where
        T:  Clone        
{

    /// Iterate over combinations drawn from a vector, in the reverse of the order returned by `Itertools.Combinations`.
    /// 
    /// Syntax sugar for `remap_elements( num_elements, sequence_length, remap )`, where
    /// 
    /// - `sequence_length` is the length of `v`
    /// - `remap` is the function sending `i` to `v[i]`
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::sequences_and_ordinals::CombinationsReverse;
    /// use oat_rust::utilities::functions::evaluate::IdentityFunction;
    /// 
    /// let num_elements        =   4;
    /// let sequence_length     =   2;
    /// let elements: Vec<i32>  =   vec![0,-1,-2,-3];
    /// let combinations: CombinationsReverse< i32, &Vec<i32>>         
    ///     =   CombinationsReverse::from_vec( sequence_length, &elements );
    /// let ground_truth: Vec< Vec< i32 > >     =   
    ///         vec![ vec![-2i32,-3], vec![-1,-3], vec![-1,-2], vec![-0,-3], vec![-0,-2], vec![-0,-1] ];
    /// itertools::assert_equal( combinations, ground_truth );
    /// 
    /// let elements: Vec<i32>  =   vec![0,0,0,0];
    /// let combinations: CombinationsReverse< i32, &Vec<i32>>         
    ///     =   CombinationsReverse::from_vec( sequence_length, &elements );
    /// let ground_truth: Vec< Vec< i32 > >     =   
    ///         vec![ vec![0,0], vec![0,0], vec![0,0], vec![0,0], vec![0,0], vec![0,0] ];
    /// itertools::assert_equal( combinations, ground_truth );
    /// ```
    pub fn from_vec( sequence_length: usize, elements: &'a Vec< T > ) -> Self 
    {
        let num_elements = elements.len();
        CombinationsReverse::remap_elements( num_elements, sequence_length, elements )
    }  
}



impl < T, F >

    CombinationsReverse
        < T, F >
{
    /// Remaps elements with `remap`.
    /// 
    /// Returns the same sequence of elements as `CombinationsReverse( num_elements, sequence_length )`,
    /// except that elements are remapped using the function `remap`.  Note: this means that 
    /// individual sequences might have repeated elements.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::sequences_and_ordinals::CombinationsReverse;
    /// use oat_rust::utilities::functions::evaluate::IdentityFunction;
    /// 
    /// let num_elements        =   4;
    /// let sequence_length     =   2;
    /// let remap               =   vec![0,-1,-2,-3]; // vectors implement EvaluateFunction automatically
    /// let combinations: CombinationsReverse< i32, &Vec<i32>>       =   
    ///     CombinationsReverse::remap_elements( num_elements, sequence_length, &remap );
    /// let ground_truth: Vec< Vec< i32 > >     =   
    ///         vec![ vec![-2i32,-3], vec![-1,-3], vec![-1,-2], vec![-0,-3], vec![-0,-2], vec![-0,-1] ];
    /// itertools::assert_equal( combinations, ground_truth );
    /// 
    /// let remap: Vec< i32 >   =   vec![0,0,0,0]; // vectors implement EvaluateFunction automatically
    /// let combinations: CombinationsReverse< i32, &Vec<i32>>       =   
    ///     CombinationsReverse::remap_elements( num_elements, sequence_length, &remap );
    /// let ground_truth: Vec< Vec< i32 > >     =   
    ///         vec![ vec![0,0], vec![0,0], vec![0,0], vec![0,0], vec![0,0], vec![0,0] ];
    /// itertools::assert_equal( combinations, ground_truth );
    /// ```
    pub fn remap_elements( num_elements: usize, sequence_length: usize, remap: F ) -> Self {
        if sequence_length <= num_elements {
            let combination: Vec<usize> = ( num_elements - sequence_length .. num_elements ).collect();
            let is_empty = false;
            CombinationsReverse { combination, num_elements, is_empty, remap, phantom_element: PhantomData }
        } else {
            let combination: Vec<usize> = vec![];
            let is_empty = true;
            CombinationsReverse { combination, num_elements, is_empty, remap, phantom_element: PhantomData }            
        }
    }    


      

    /// Return `true` iff we can decrement `self.combination[index]` and still have a strictly ascending sequence.
    fn can_decrement( &self, index: usize ) -> bool {
        if index == 0 { self.combination[0] > 0 }
        else { self.combination[index-1] + 1 < self.combination[index] }       
    }
}




impl < T, F >

    Iterator for 
    
    CombinationsReverse
        < T, F > 

    where
        F:      Clone + EvaluateFunction< usize, T >
{
    type Item = Vec< T >;

    fn next( &mut self ) -> Option< Self::Item > {  
        if self.is_empty { return None }
        let return_val = self.combination.iter().cloned()
            .map_by_transform( ReferencedEvaluator::new( &self.remap ) )
            .collect();
        let m = self.combination.len();

        for k in ( 0 .. m ).rev() {
            // working from the end of the vector to the front, for each element

            if self.can_decrement(k) {
                // if we can decrement, then decrement
                self.combination[k] -= 1;
                // and replace all following elements with their maximum possible values
                let offset = self.num_elements - m;
                for l in k+1 .. m { self.combination[l] = offset + l }
                // then return 
                return Some( return_val )                  
            }
        }

        // if we make it here, then we cannot generate any new sequences
        self.is_empty = true;
        Some( return_val )        
    }
}








//  ---------------------------------------------------------------------------
//  BIMAPS WITH {0, .., N}
//  ---------------------------------------------------------------------------

/// Represents a surjective map {0,..,N} -> S and another map S -> {0, .., N}; if
/// one of these maps is a bijection, then the other should be its inverse.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BiMapSequential< T >
    where T : Hash + Eq
{ 
    ord_to_val: Vec< T >, 
    val_to_ord: HashMap< T, usize > 
}

impl < T > BiMapSequential < T > 
    where T: Clone + Hash + Eq
{
    // /// Evaluate the function {0, ..., N} -> S
    // pub fn elt( &self, ord: usize ) -> Option< T >   { if ord < self.ord_to_val.len()self.ord_to_val[ ord ].clone() }

    // /// Evaluate the function S -> {0, ..., N}
    // pub fn ord( &self, elt_ref: &T     ) -> usize { self.val_to_ord.get( elt_ref ).unwrap().clone() }   

    /// Returns an immutable reference to the private field `self.val_to_ord`
    pub fn val_to_ord_hashmap( &self ) -> & HashMap< T, usize > { & self.val_to_ord }

    /// Returns an immutable reference to the private field `self.ord_to_val`.
    pub fn ord_to_val_vec( &self ) -> & Vec < T > { & self.ord_to_val }

    /// Returns the number of elements in the map.
    pub fn len( &self ) -> usize { self.ord_to_val.len() }

    /// Returns `true` if the sequence has length 0.
    pub fn is_empty( &self ) -> bool { self.len() == 0 }    

    /// Returns `true` if the given key is present in the set of matched keys.
    pub fn contains_key( &self, key_ref: &T ) -> bool { self.val_to_ord.contains_key( key_ref ) }

    /// Evaluate the function {0, ..., N} -> S
    pub fn ord( &self, a: &T ) -> Option< usize > { 
        self.val_to_ord.get( a ).copied() 
    }

    /// Returns the `a`th element of the sequence.
    pub fn val( &self, a: usize ) -> T { 
        self.ord_to_val[ a ].clone()
        // if a < self.ord_to_val.len() { Some( self.ord_to_val[ a ].clone() ) } else { None }
    }      

    /// Reverses the order of elements; if the `BiMapSequential` stores a bijeciton `phi: X -> {0.. N}` then
    /// the resulting bijection `psi` sends `x` to `N - 1 - phi(x)`
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::sequences_and_ordinals::BiMapSequential;
    /// use itertools::Itertools;
    /// 
    /// let mut a = BiMapSequential::from_vec( vec![ 1, 2, 3 ] );
    /// 
    /// let internal_ord_to_val = a.ord_to_val_vec().iter().cloned();
    /// itertools::assert_equal( internal_ord_to_val , vec![ 1, 2, 3] );
    /// 
    /// let mut internal_val_ord_pairs = a.val_to_ord_hashmap().iter().map(|(&x,&y)| (x.clone(), y.clone())).collect_vec();
    /// internal_val_ord_pairs.sort();
    /// itertools::assert_equal( internal_val_ord_pairs, vec![ (1,0), (2,1), (3,2) ] );  
    /// 
    /// a.reverse();
    /// 
    /// let internal_ord_to_val = a.ord_to_val_vec().iter().cloned();
    /// itertools::assert_equal( internal_ord_to_val , vec![ 3, 2, 1] );
    /// 
    /// let mut internal_val_ord_pairs = a.val_to_ord_hashmap().iter().map(|(&x,&y)| (x.clone(), y.clone())).collect_vec();
    /// internal_val_ord_pairs.sort();
    /// itertools::assert_equal( internal_val_ord_pairs, vec![ (1,2), (2,1), (3,0) ] ); 
    /// ```
    pub fn reverse( &mut self ) {
        // reverse the order of entries in the vector `ord_to_val`
        self.ord_to_val.reverse();
        // reverse the ordinals mapped to by the hashmap `val_to_ord`
        let num_pairs_minus_one = self.ord_to_val.len() - 1;
        for val in self.val_to_ord.values_mut() {
            *val = num_pairs_minus_one - *val;
        }
    }

    /// Create a new, empty sequential bimap.
    pub fn new() -> BiMapSequential< T > {
        BiMapSequential{ ord_to_val: Vec::new(), val_to_ord: HashMap::new() }
    }
    
    /// Create sequential bimap
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::sequences_and_ordinals::BiMapSequential;        
    /// use assert_panic::assert_panic;
    /// 
    /// // Construct a BiMapSequential
    /// let bimap   =   BiMapSequential::from_vec( vec![ 'a', 'b' ] );
    /// assert_eq!( 0, bimap.ord( &'a' ).unwrap() );
    /// assert_eq!( 'b', bimap.val( 1 ) );        
    /// 
    /// // Affirm that the compiler panics if we try to construct a BiMapSequential from a vector with duplicate elements
    /// assert_panic!( 
    ///         let _bimap = BiMapSequential::from_vec( vec![ 4, 3, 3, ] )  
    ///     );
    /// ```
    pub fn from_vec( vec: Vec< T > ) -> BiMapSequential< T >
    {
        let hash    =   HashMap::from_iter(
                            vec.iter().cloned().enumerate().map(|x| (x.1, x.0) )
                        );
        if hash.len() < vec.len() { panic!("Attempt to create `BiMapSequential` from a vector with two or more equal entries.") }
        BiMapSequential{ ord_to_val: vec, val_to_ord: hash}
    }

    /// Push a new element to the end of the sequence.  Panics if the element is already present in the sequence.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::sequences_and_ordinals::BiMapSequential;        
    /// 
    /// // Construct a BiMapSequential
    /// let mut bimap   =   BiMapSequential::from_vec( vec![ 'a', 'b' ] );
    /// bimap.push( 'c' );
    /// assert_eq!( 2, bimap.ord( &'c' ).unwrap() );
    /// assert_eq!( 'c', bimap.val( 2 ) );       
    /// 
    /// // The following induces a panic, because pushing an element which already exists in the sequence results in a panic.
    /// // let mut bimap   =   BiMapSequential::from_vec( vec![ 1, 2 ] );
    /// // let _ = bimap.push( 2 );
    /// ```
    pub fn push( &mut self, new_element: T ) {
        // Push the new key-value pair to the internal hashmap
        let prior_value = self.val_to_ord.insert( 
                    new_element.clone(), 
                    self.len() 
                );
        // Panic if the element already exists in the sequence
        if prior_value.is_some() { panic!( "Attempted to push a value to a BiMapSequential, but the value is already present in the sequence." ) }
        // Append the element to the end of the sequence vector
        self.ord_to_val.push( new_element );
    }
}

impl    < T > 
        FromIterator< T > 
        for 
        BiMapSequential < T > 
    where   T:  Clone + Hash + std::cmp::Eq
{
    fn from_iter< I: IntoIterator<Item=T>>(iter: I) -> Self {

        let vec     =   Vec::from_iter( iter );

        BiMapSequential::from_vec( vec )
    }
}


//  ---------------------------------------------------------------------------
//  PRIMITIVE ORDINALS
//  ---------------------------------------------------------------------------


// #[derive(Clone, Debug, PartialEq)]
// pub struct OrdinalData < T : Ord + Eq + PartialOrd + PartialEq + Hash > {
//     pub ord_to_val:  Vec< T >,
//     pub val_to_ord:  HashMap< T, usize >
// }

// impl    < T >
//         OrdinalData
//         < T >
//         where T : Ord + Eq + PartialOrd + PartialEq + Hash + Clone
// {
//     /// The ordinal of the raw filtration value
//     pub fn ord( &self, a: &T ) -> Option< usize > { 
//         self.val_to_ord.get( a ).map(|x| x.clone()) 
//     }
//     /// The raw filtration value of the ordinal
//     pub fn val( &self, a: usize ) -> Option< T > { 
//         if a < self.ord_to_val.len() { Some( self.ord_to_val[ a ].clone() ) } else { None }
//     }    
// }


/// Given a vector of elements of a poset, first sort the vector and delete 
/// duplicate entries; the resulting vector represents a bijection from 
/// {0, .., n} to the set of unique values in the vector.  Store this new vector
/// in an OrdinalData struct together with a hashmap representing the inverse bijection
pub fn ordinate_unique_vals < FilRaw > ( v: & Vec< FilRaw > ) -> BiMapSequential< FilRaw > 
    where FilRaw: Ord + Hash + Clone
{
    let mut a       =   v.clone();
    let mut b       =   HashMap::new();
    a.sort();       // sort entries
    a.dedup();      // remove duplicates

    for (i, t) in a.iter().enumerate() {
        b.insert( t.clone(), i );
    }

    BiMapSequential { ord_to_val: a, val_to_ord: b }
}

/// If we view a vector as a surjective function `f: { 0 .. n } ->> X`,
/// then `reverse_hash_sequential` returns a function `X -> { 0 .. n }`
/// which is inverse to `f` if and only if `f` is bijective.
pub fn  reverse_hash_sequential< T: Hash + std::cmp::Eq + Clone >( 
            vec: & Vec< T >
        ) 
        -> 
        HashMap< T, usize >
{
    let mut rev_hash    =   HashMap::new();

    for (i, t) in vec.iter().enumerate() {
        rev_hash.insert( t.clone(), i );
    }

    rev_hash
}







//  TESTS
//  =========================================================================================================


//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    

    #[test] 
    fn test_bimap_from_vec() {
        use crate::utilities::sequences_and_ordinals::BiMapSequential;        
        use assert_panic::assert_panic;

        // Construct a BiMapSequential
        let bimap   =   BiMapSequential::from_vec( vec![ 'a', 'b' ] );
        assert_eq!( 0, bimap.ord( &'a' ).unwrap() );
        assert_eq!( 'b', bimap.val( 1 ) );       

        // Affirm that the compiler panics if we try to construct a BiMapSequential from a vector with duplicate elements
        assert_panic!( 
                let _bimap = BiMapSequential::from_vec( vec![ 4, 3, 3, ] )  
            );
    }  
    

    #[test] 
    fn test_bimap_push() {
        use crate::utilities::sequences_and_ordinals::BiMapSequential;        

        // Construct a BiMapSequential
        let mut bimap   =   BiMapSequential::from_vec( vec![ 'a', 'b' ] );
        bimap.push( 'c' );
        assert_eq!( 2, bimap.ord( &'c' ).unwrap() );
        assert_eq!( 'c', bimap.val( 2 ) );       

        // The following induces a panic, because pushing an element which already exists in the sequence results in a panic.
        // let mut bimap   =   BiMapSequential::from_vec( vec![ 1, 2 ] );
        // let _ = bimap.push( 2 );
    }      


    #[test]
    fn test_bimap_reverse() {
        use crate::utilities::sequences_and_ordinals::BiMapSequential;
        use itertools::Itertools;

        let mut a = BiMapSequential::from_vec( vec![ 1, 2, 3 ] );

        let internal_ord_to_val = a.ord_to_val_vec().iter().cloned();
        itertools::assert_equal( internal_ord_to_val , vec![ 1, 2, 3] );

        let mut internal_val_ord_pairs = a.val_to_ord_hashmap().iter().map(|(&x,&y)| (x, y)).collect_vec();
        internal_val_ord_pairs.sort();
        itertools::assert_equal( internal_val_ord_pairs, vec![ (1,0), (2,1), (3,2) ] );  
        
        a.reverse();

        let internal_ord_to_val = a.ord_to_val_vec().iter().cloned();
        itertools::assert_equal( internal_ord_to_val , vec![ 3, 2, 1] );

        let mut internal_val_ord_pairs = a.val_to_ord_hashmap().iter().map(|(&x,&y)| (x, y)).collect_vec();
        internal_val_ord_pairs.sort();
        itertools::assert_equal( internal_val_ord_pairs, vec![ (1,2), (2,1), (3,0) ] );                 
    }


    // fn test_iter() {
    //     let v = Arc::new( vec![ vec![0,1], vec![0,1] ] );
    //     let u = v.iter();
    
    //     let w = vec![ vec![0,1], vec![0,1] ];
    //     let z = w.iter();        
    // }
}