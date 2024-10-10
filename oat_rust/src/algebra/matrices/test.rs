use std::marker::PhantomData;

use std::slice::Iter;
use std::iter::Cloned;

use super::query::{IndicesAndCoefficients, ViewRow};






pub struct Vv< 'a >{
    vv:                 Vec< Vec< (usize, usize) > >,
    phantom_lifetime:   PhantomData< &'a usize >
}

impl < 'a > Vv< 'a > {
    pub fn new( vv: Vec< Vec< (usize, usize) > > ) -> Self { Vv{ vv, phantom_lifetime: PhantomData } }
}


//  IndicesAndCoefficients
impl < 'a > 

    IndicesAndCoefficients for 
    Vv< 'a >

{ type ColIndex = usize; type RowIndex = usize; type Coefficient = usize; }   



impl < 'a, > 

    ViewRow for 
    
    Vv< 'a >

    where

    {
        type ViewMajor            =   Cloned<Iter<'a, (usize, usize) >>;
        type ViewMajorIntoIter    =   Cloned<Iter<'a, (usize, usize) >>;
        type EntryMajor =   (usize, usize);                

        fn   view_major( &self, index: usize ) -> Self::ViewMajor { self.vv[index].iter().cloned() }
    }