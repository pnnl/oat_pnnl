//! Provides functionality for sets, e.g. trait for determining whether a collection (e.g. a hashmap or a sequence) contains a key or an element.
//! 
//! See also the [oat_rust module for general iterators](crate::utilities::iterators::general), for set operations on ordered iterators.

use std::collections::{HashMap, BTreeMap, HashSet, BTreeSet};
use std::hash::Hash;

use itertools::{Merge, Itertools};

use super::iterators::general::OnlyDuplicates;


/// Contains one method, [`map_has_key_or_sequence_has_element`], which determines whether a collection contains an element.
pub trait MapHasKeyOrSequenceHasElement< Key > {
    /// Returns `true` if the hash object contains the given key.
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool;
}

impl < Key: Hash + std::cmp::Eq, Val> MapHasKeyOrSequenceHasElement< Key > for HashMap< Key, Val > {
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool { self.contains_key( key ) }    
}

impl < Key: Ord, Val> MapHasKeyOrSequenceHasElement< Key > for BTreeMap< Key, Val > {
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool { self.contains_key( key ) }    
}

impl < Key: Hash + std::cmp::Eq > MapHasKeyOrSequenceHasElement< Key > for HashSet< Key > {
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool { self.contains( key ) }    
}

impl < Key: Ord > MapHasKeyOrSequenceHasElement< Key > for BTreeSet< Key > {
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool { self.contains( key ) }    
}

impl < Key: std::cmp::PartialEq > MapHasKeyOrSequenceHasElement< Key > for Vec< Key > {
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool { self.contains( key ) }    
}

// Whenever `T` implements the trait, so does `&'a T`
impl < 'a, Key, T: MapHasKeyOrSequenceHasElement< Key > > MapHasKeyOrSequenceHasElement< Key > for &'a T {
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool { (*self).map_has_key_or_sequence_has_element( key ) }    
}

/// Implements the trait `MapHasKeyOrSequenceHasElement` on a reference.
/// 
/// This object solves the following problem: suppose you have a type `T` that implements the trait `MapHasKeyOrSequenceHasElement`,
/// but you need to pass this object to a function that consumes its arguments.  Rather than losing `T`, you can wrap a reference
/// to `T` in [`MapHasKeyOrSequenceHasElementRefWrapper`].  The wrapper will also implement the look-up method [`map_has_key_or_sequence_has_element`].
#[derive(Copy, Clone)]
pub struct MapHasKeyOrSequenceHasElementRefWrapper< 'a, T > { 
    pub map_has_key_or_sequence_has_element_ref: &'a T 
}

impl < 'a, T > MapHasKeyOrSequenceHasElementRefWrapper< 'a, T > {
    pub fn new( map_has_key_or_sequence_has_element_ref: &'a T ) -> MapHasKeyOrSequenceHasElementRefWrapper< 'a, T > { MapHasKeyOrSequenceHasElementRefWrapper{ map_has_key_or_sequence_has_element_ref } }
}

impl < 'a, Key, T: MapHasKeyOrSequenceHasElement< Key > > MapHasKeyOrSequenceHasElement< Key > for MapHasKeyOrSequenceHasElementRefWrapper< 'a, T > {
    fn map_has_key_or_sequence_has_element( &self, key: &Key ) -> bool { self.map_has_key_or_sequence_has_element( key ) }    
}