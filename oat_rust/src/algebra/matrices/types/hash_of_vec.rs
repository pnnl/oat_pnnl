//! A hashmap of vectors.  The keys of the hashmap represent major keys of the matrix.

use std::collections::HashMap;
use std::hash::Hash;
use std::marker::PhantomData;

use crate::algebra::matrices::query::{ViewRow, IndicesAndCoefficients};





//  & HashMap< Vec < Entry > >
//  --------------------------------------------------------------------------------


//  THERE IS NOTHING WRONG WITH THIS IMPLEMENTATION PER SE, BUT IF WE WANT IT TO IMPLEMENT THE ORACLE TRAIT WE WILL HAVE TO ADD SOME PHANTOM DATA TO HashMap OR FIND SOME OTHER WAY TO ASSOCIATE A ColIndex and Coefficient TYPES

//
// impl < 'a, IndexCoeffPair, RowIndex >
    
//     ViewRow for 

//     &'a HashMap
//         < RowIndex, Vec< IndexCoeffPair > >

// where   //IndexCoeffPair:       KeyValGet < ColIndex, Coefficient >,
//         // OrderOperator:     JudgePartialOrder<  IndexCoeffPair >,   
//         RowIndex:                 std::cmp::Eq + Hash,
// {
//     type ViewMajor          =   Vec< IndexCoeffPair >;
//     type EntryMajor = < Self::ViewMajor as IntoIterator >::Item;
//     type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;
    
//     fn view_major( & self, index: Self::RowIndex ) -> Self::ViewMajor {
//         return self.get( &index ).unwrap().as_slice()
//     } 
// }



//  HashOfVec
//  --------------------------------------------------------------------------------
/// A hashmap of vectors.  The keys of the hashmap represent major keys.
pub struct HashOfVec< ColIndex, RowIndex, Coefficient, IndexCoeffPair, > { 
    hash_of_vec:        HashMap< RowIndex, Vec< IndexCoeffPair > >,
    phantom_keymin:     PhantomData< ColIndex >,
    phantom_snzval:     PhantomData< Coefficient >,
}


impl < 'a, ColIndex, RowIndex, Coefficient, IndexCoeffPair, >

    IndicesAndCoefficients for

    &'a HashOfVec
        < ColIndex, RowIndex, Coefficient, IndexCoeffPair, >

{ 
    type ColIndex = ColIndex; 
    type RowIndex = RowIndex; 
    type Coefficient = Coefficient; 
    type EntryMajor = &'a IndexCoeffPair; 
    type EntryMinor = (RowIndex,Coefficient);
}            


impl < 'a, ColIndex, RowIndex, Coefficient, IndexCoeffPair, >
    
    ViewRow  for 

    &'a HashOfVec
        < ColIndex, RowIndex, Coefficient, IndexCoeffPair, >

where   //IndexCoeffPair:       KeyValGet < ColIndex, Coefficient >,
        // OrderOperator:     JudgePartialOrder<  IndexCoeffPair >,   
        RowIndex:                 std::cmp::Eq + Hash,
{
    type ViewMajor          =   &'a[IndexCoeffPair];
    type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;

    fn view_major( & self, index: Self::RowIndex ) -> Self::ViewMajor {
        return self.hash_of_vec.get( &index ).unwrap().as_slice()
    } 
}