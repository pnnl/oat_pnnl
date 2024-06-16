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

use crate::entries::KeyValGet;
use crate::utilities::partial_order::StrictlyLess;








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
    pub fn hashmap_val_to_ord( &self ) -> & HashMap< T, usize > { & self.val_to_ord }

    /// Returns an immutable reference to the private field `self.ord_to_val`.
    pub fn vec_ord_to_val( &self ) -> & Vec < T > { & self.ord_to_val }

    /// Returns the number of elements in the map.
    pub fn len( &self ) -> usize { self.ord_to_val.len() }

    /// Returns `true` if the sequence has length 0.
    pub fn is_empty( &self ) -> bool { self.len() == 0 }    

    /// Returns `true` if the given key is present in the set of matched keys.
    pub fn contains_key( &self, key_ref: &T ) -> bool { self.val_to_ord.contains_key( key_ref ) }

    /// Evaluate the function {0, ..., N} -> S
    pub fn ord( &self, a: &T ) -> Option< usize > { 
        self.val_to_ord.get( a ).map(|x| x.clone()) 
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
    /// let internal_ord_to_val = a.vec_ord_to_val().iter().cloned();
    /// itertools::assert_equal( internal_ord_to_val , vec![ 1, 2, 3] );
    /// 
    /// let mut internal_val_ord_pairs = a.hashmap_val_to_ord().iter().map(|(&x,&y)| (x.clone(), y.clone())).collect_vec();
    /// internal_val_ord_pairs.sort();
    /// itertools::assert_equal( internal_val_ord_pairs, vec![ (1,0), (2,1), (3,2) ] );  
    /// 
    /// a.reverse();
    /// 
    /// let internal_ord_to_val = a.vec_ord_to_val().iter().cloned();
    /// itertools::assert_equal( internal_ord_to_val , vec![ 3, 2, 1] );
    /// 
    /// let mut internal_val_ord_pairs = a.hashmap_val_to_ord().iter().map(|(&x,&y)| (x.clone(), y.clone())).collect_vec();
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
        if ! prior_value.is_none() { panic!( "Attempted to push a value to a BiMapSequential, but the value is already present in the sequence." ) }
        // Append the element to the end of the sequence vector
        self.ord_to_val.push( new_element.clone() );
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
        b.insert( t.clone(), i.clone() );
    }

    BiMapSequential { ord_to_val: a, val_to_ord: b }
}


pub fn  reverse_hash_sequential< T: Hash + std::cmp::Eq + Clone >( 
            vec: & Vec< T >
        ) 
        -> 
        HashMap< T, usize >
{
    let mut rev_hash    =   HashMap::new();

    for (i, t) in vec.iter().enumerate() {
        rev_hash.insert( t.clone(), i.clone() );
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

        let internal_ord_to_val = a.vec_ord_to_val().iter().cloned();
        itertools::assert_equal( internal_ord_to_val , vec![ 1, 2, 3] );

        let mut internal_val_ord_pairs = a.hashmap_val_to_ord().iter().map(|(&x,&y)| (x.clone(), y.clone())).collect_vec();
        internal_val_ord_pairs.sort();
        itertools::assert_equal( internal_val_ord_pairs, vec![ (1,0), (2,1), (3,2) ] );  
        
        a.reverse();

        let internal_ord_to_val = a.vec_ord_to_val().iter().cloned();
        itertools::assert_equal( internal_ord_to_val , vec![ 3, 2, 1] );

        let mut internal_val_ord_pairs = a.hashmap_val_to_ord().iter().map(|(&x,&y)| (x.clone(), y.clone())).collect_vec();
        internal_val_ord_pairs.sort();
        itertools::assert_equal( internal_val_ord_pairs, vec![ (1,2), (2,1), (3,0) ] );                 
    }
}