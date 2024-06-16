//! A hashmap of vectors.  The keys of the hashmap represent major keys of the matrix.

use std::collections::HashMap;
use std::hash::Hash;
use std::marker::PhantomData;

use crate::matrices::matrix_oracle_traits::{OracleMajor, IndicesAndCoefficients};
use crate::utilities::partial_order::StrictlyLess;




//  & HashMap< Vec < Entry > >
//  --------------------------------------------------------------------------------


//  THERE IS NOTHING WRONG WITH THIS IMPLEMENTATION PER SE, BUT IF WE WANT IT TO IMPLEMENT THE ORACLE TRAIT WE WILL HAVE TO ADD SOME PHANTOM DATA TO HashMap OR FIND SOME OTHER WAY TO ASSOCIATE A KeyMin and SnzVal TYPES

//
// impl < 'a, IndexCoeffPair, KeyMaj >
    
//     OracleMajor for 

//     &'a HashMap
//         < KeyMaj, Vec< IndexCoeffPair > >

// where   //IndexCoeffPair:       KeyValGet < KeyMin, SnzVal >,
//         // OrderComparator:     StrictlyLess<  IndexCoeffPair >,   
//         KeyMaj:                 std::cmp::Eq + Hash,
// {
//     type ViewMajor          =   Vec< IndexCoeffPair >;
//     type ViewMajorEntry     =   < Self::ViewMajor as IntoIterator >::Item;
//     type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;
    
//     fn view_major( & self, index: Self::KeyMaj ) -> Self::ViewMajor {
//         return self.get( &index ).unwrap().as_slice()
//     } 
// }



//  HashOfVec
//  --------------------------------------------------------------------------------
/// A hashmap of vectors.  The keys of the hashmap represent major keys.
pub struct HashOfVec< KeyMin, KeyMaj, SnzVal, IndexCoeffPair, OrderComparator > { 
    hash_of_vec:        HashMap< KeyMaj, Vec< IndexCoeffPair > >,
    order_comparator:   OrderComparator,
    phantom_keymin:     PhantomData< KeyMin >,
    phantom_snzval:     PhantomData< SnzVal >,
}


impl < 'a, KeyMin, KeyMaj, SnzVal, IndexCoeffPair, OrderComparator >

    IndicesAndCoefficients for

    &'a HashOfVec
        < KeyMin, KeyMaj, SnzVal, IndexCoeffPair, OrderComparator >

{ type KeyMin = KeyMin; type KeyMaj = KeyMaj; type SnzVal = SnzVal; }            


impl < 'a, KeyMin, KeyMaj, SnzVal, IndexCoeffPair, OrderComparator >
    
    OracleMajor  for 

    &'a HashOfVec
        < KeyMin, KeyMaj, SnzVal, IndexCoeffPair, OrderComparator >

where   //IndexCoeffPair:       KeyValGet < KeyMin, SnzVal >,
        // OrderComparator:     StrictlyLess<  IndexCoeffPair >,   
        KeyMaj:                 std::cmp::Eq + Hash,
{
    type ViewMajor          =   &'a[IndexCoeffPair];
    type ViewMajorEntry     =   < Self::ViewMajor as IntoIterator >::Item;
    type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;

    fn view_major( & self, index: Self::KeyMaj ) -> Self::ViewMajor {
        return self.hash_of_vec.get( &index ).unwrap().as_slice()
    } 
}