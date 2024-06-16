
//!  General utilities for work with iterators.

use std::{iter::{Iterator, Peekable, Once}, marker::PhantomData, collections::HashMap};
use std::hash::Hash;

use itertools::{Merge, Itertools};

use crate::utilities::functions::evaluate::EvaluateFunction;


//  ---------------------------------------------------------------------------
//  FIND MINIMUM INDEX

/// Returns the index of the minimum value of the iterator.
/// 
/// If several elements are equally minimum, the index of the first element is returned. If the iterator is empty, None is returned.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::find_min;
/// 
/// let a               =   vec![ 1, 3, 0, 2, 0 ];
/// let b: Vec<usize>   =   vec![ ];
/// 
/// assert_eq!( find_min(a), Some(2) );
/// assert_eq!( find_min(b), None    );
/// ```
pub fn find_min< I >( iter: I ) -> Option< usize > 
    where 
        I:          IntoIterator,
        I::Item:    Ord,
{
    iter.into_iter()
        .enumerate()
        .min_by(|x,y| x.1.cmp(&y.1) )
        .map(|x| x.0 )
}


//  ---------------------------------------------------------------------------
//  PEEKING


/// Similar to itertools::PeakingNext, but without the `F` parameter
/// (we hypothesize that omitting this closure will help with type
/// inference)
pub trait PeekUnqualified : Iterator {
    fn peek_unqualified( &mut self ) -> Option < & Self::Item >;
}

impl < I : Iterator > PeekUnqualified for Peekable< I >
{
    fn peek_unqualified( &mut self ) -> Option < &<Self as Iterator>::Item > { self.peek() }
}


//  ---------------------------------------------------------------------------
//  SKIP UNTIL


/// Skips each item returned by the iterator, until some item `i` satisfies 
/// `predicate(&i)=true`.  The resulting iterator returns `i` and every item 
/// thereafter.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::SkipUntil;
/// 
/// let vec = vec![ 0, 1, 2, 1, 0];
/// let iter = vec.iter().cloned().peekable().skip_until( |x| x > &1);
/// itertools::assert_equal( iter, vec![ 2, 1, 0 ])  
/// ```
pub trait SkipUntil : Iterator {
    fn skip_until< P >( self, predicate: P ) -> Self
    where
        P: FnMut(&Self::Item) -> bool;
}

impl < I : Iterator + PeekUnqualified > SkipUntil for I 
{
    fn skip_until< P >( mut self, mut predicate: P ) -> Self 
    where
        P: FnMut(&Self::Item) -> bool    
    
    {  
        while let Some( peek ) = self.peek_unqualified() {   if predicate( peek ) { break } else { let _ = self.next(); } }
        return self
    }
}






//  ---------------------------------------------------------------------------
//  TRANSFORM ENTRY-WISE


/// Similar to `std::iter::Map` but instead of taking a transformation that implements `FnMut`, 
/// this object takes a transformation that implements `EvaluateFunction`.
/// 
/// See also [`FilterMapByTransform`].
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::MapByTransform;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
/// 
/// let entrywise_transform_unwrapped = |x : usize| -> i32 { x as i32 };
/// let entrywise_transform_wrapped = EvaluateFunctionFnMutWrapper::new(entrywise_transform_unwrapped);
/// let vec = vec![ 1, 2, 3];
/// let iter = vec.iter().cloned();
/// 
/// let transformed_entrywise = MapByTransform::new( iter, entrywise_transform_wrapped );
/// itertools::assert_equal( transformed_entrywise, vec![ 1i32, 2, 3] );  
/// ```
pub struct MapByTransform
                < Iter, ItemNew, ItemwiseTransform >
    where
        Iter:               Iterator,
        ItemwiseTransform:  EvaluateFunction< Iter::Item, ItemNew >
{
    iter:                   Iter,
    itemwise_transform:     ItemwiseTransform,
    phantom_itemnew:        PhantomData< ItemNew >,
}                

impl < Iter, ItemNew, ItemwiseTransform >

    MapByTransform
        < Iter, ItemNew, ItemwiseTransform >

    where
        Iter:               Iterator,
        ItemwiseTransform:  EvaluateFunction< Iter::Item, ItemNew > 
{
    pub fn new( iter: Iter, itemwise_transform: ItemwiseTransform ) -> MapByTransform < Iter, ItemNew, ItemwiseTransform > {  
        MapByTransform{ iter, itemwise_transform, phantom_itemnew: PhantomData }
    }
}          


impl < Iter, ItemNew, ItemwiseTransform >

    Iterator for

    MapByTransform
        < Iter, ItemNew, ItemwiseTransform >

    where
        Iter:               Iterator,
        ItemwiseTransform:  EvaluateFunction< Iter::Item, ItemNew >        
{
    type Item = ItemNew;
    fn next( &mut self ) -> Option< Self::Item > {
        match self.iter.next() {
            Some( item_old ) => { return Some( self.itemwise_transform.evaluate_function( item_old ) ) },
            None => { return None }
        }
    }
}


impl < Iter, ItemNew, ItemwiseTransform >

    Clone for

    MapByTransform
        < Iter, ItemNew, ItemwiseTransform >

    where
        Iter:               Clone + Iterator,
        ItemwiseTransform:  Clone + EvaluateFunction< Iter::Item, ItemNew >        
{
    fn clone( & self ) -> Self { MapByTransform::new( self.iter.clone(), self.itemwise_transform.clone() )}
}


pub trait MapByTransformTrait< ItemwiseTransform > {
    fn map_by_transform< ItemNew >( self, itemwise_transform: ItemwiseTransform ) 
        -> 
        MapByTransform
            < Self, ItemNew, ItemwiseTransform >
    where
        Self:               Iterator + Sized,
        ItemwiseTransform:  EvaluateFunction< Self::Item, ItemNew >
    {
        MapByTransform { iter: self, itemwise_transform, phantom_itemnew: PhantomData }
    }
}

impl < ItemwiseTransform, Iter >

    MapByTransformTrait
        < ItemwiseTransform > for
    
    Iter

    where 
        Iter:   Iterator
{}        


//  ---------------------------------------------------------------------------
//  FILTER MAP -- FILTERTING WITH OBJECTS THAT IMPLMENT EvaluateFunction INSTEAD OF FnMut




/// Similar to `std::iter::FilterMap`, but takes an object that satisfies `FilterMapTransform< Input, Output >` rather than an object
/// that implements `FnMut( Input ) -> Output`.
/// 
/// The `std::iter::filter_map` method for iterators is highly versatile, however it requires an
/// argument that implements `FnMut`.  It is often hard to provide such an argument without
/// using closure operators which have `opaque type`.  As an alternative, the oat_rust library
/// offers this struct, which takes any object which imlements `FilterMapTransform`.
/// 
/// See also [`MapByTransform`].
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::{FilterMapByTransform};
/// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
/// 
/// // define a struct to execute the filter-map function
/// struct ShiftDown;
/// 
/// // implement the `EvaluateFunction` trait
/// impl EvaluateFunction< usize, Option<i32> > for ShiftDown{
///     fn evaluate_function( &self, input: usize ) -> Option< i32 > {
///         match input > 3 {
///             true    => { return None },
///             false   => { return Some( input as i32 - 3 ) }
///         }
///     }
/// }
/// 
/// // define an iterator
/// let iter = 0..9;
/// 
/// // create a `filter-mapped` object
/// let iter_transformed = FilterMapByTransform::new( iter, ShiftDown );
/// 
/// // verify this object is correct        
/// itertools::assert_equal( iter_transformed, vec![ -3, -2, -1, 0 ]) 
/// ```
pub struct FilterMapByTransform< I, FilteringObject, ItemNew > 
    where
        I:                  Iterator,
        FilteringObject:    EvaluateFunction< I::Item, Option< ItemNew > >,
{
    iter:               I,
    filtering_object:   FilteringObject,
    phantom_itemnew:    PhantomData< ItemNew >
}

impl < I, FilteringObject, ItemNew >  

    FilterMapByTransform
        < I, FilteringObject, ItemNew > 

    where
        I:                  Iterator,
        FilteringObject:    EvaluateFunction< I::Item, Option< ItemNew > >,
{
    pub fn new( iter: I, filtering_object: FilteringObject ) -> FilterMapByTransform< I, FilteringObject, ItemNew > 
        { FilterMapByTransform{ iter, filtering_object, phantom_itemnew: PhantomData } }
}        


impl < I, FilteringObject, ItemNew > 

    Iterator for 
    
    FilterMapByTransform< I, FilteringObject, ItemNew > 

    where
        I:                  Iterator,
        FilteringObject:    EvaluateFunction< I::Item, Option< ItemNew > >

{
    type Item = ItemNew;

    fn next( &mut self ) -> Option< Self::Item >  {
        while let Some( next_item ) = self.iter.next() {
            match  self.filtering_object.evaluate_function( next_item ) {
                None    => continue,
                Some( new_item ) => { return Some( new_item ) }
            }
        }
        return None
    }
}

impl < I, FilteringObject, ItemNew > 

    Clone for 
    
    FilterMapByTransform< I, FilteringObject, ItemNew > 

    where
        I:                  Clone + Iterator,
        FilteringObject:    Clone + EvaluateFunction< I::Item, Option< ItemNew > >    
{
    fn clone( & self ) -> Self { FilterMapByTransform::new( self.iter.clone(), self.filtering_object.clone() )}
}


// Converts 
pub trait FilterMapByTransformTrait< ItemNew >
    where
        Self:               Iterator + Sized,
        // FilteringObject:    FilterMapTransform< Self::Item, ItemNew >
{
    /// Wraps an iterator in a [`FilterMapByTransform`] struct.  
    /// 
    /// **Note** this method is convenient for chaining transformations of an iterator, however it requires the
    /// iterator to implement `Sized`.  To create a `FilterMapTransform` when the iterator does not implement
    /// `Sized`, one can use [`FilterMapByTransform::new`]
    fn filter_map_by_transform< 
                FilteringObject: EvaluateFunction< Self::Item, Option< ItemNew > > 
            >
        ( self, filtering_object: FilteringObject ) 
        -> 
        FilterMapByTransform< Self, FilteringObject, ItemNew > ;
}

impl < Iter, ItemNew > 

    FilterMapByTransformTrait
        < ItemNew > for
    
    Iter

    where
        Iter:       Iterator
{
    fn filter_map_by_transform< FilteringObject: EvaluateFunction< Self::Item, Option< ItemNew > > >
        ( self, filtering_object: FilteringObject ) -> FilterMapByTransform< Self, FilteringObject, ItemNew > {
        FilterMapByTransform::new(self, filtering_object)
    }
}

// 2-TYPE ITERATOR
// ---------------------------------------------------------------------------

/// An enum that combines two different iterator types into a single iterator type.
///
/// There are many operations available in Rust for combining iterators (chaining, merging, etc.).
/// Quite often, the function that performs the combination requires every iterator to have the same
/// type.  However, sometimes we need to combine iterators with different types.
/// This enum allows us to do so.  
/// Given an iterator `iter1` of type `I1` and an iterator `iter2` of type `I2`, we can create two
/// new iterators of the same type: `IterTwoType::Iter1( iter 1 )` and `IterTwoType::Iter2( iter 2 )`.
/// The new wrappers will iterator over exactly the same items as the iterators they contain.
pub enum IterTwoType < I1, I2 > 
where   I1:  Iterator,
        I2:  Iterator < Item = I1::Item >,
{
    Iter1( I1 ),
    Iter2( I2 ),
}

impl < I1, I2 > 
        
    Iterator       
    for IterTwoType < I1, I2 > 
    where   I1:  Iterator,
            I2:  Iterator < Item = I1::Item >,    
    
{
    type Item = I1::Item;
    
    fn next( & mut self ) -> Option< Self::Item > {
        use IterTwoType::*; // not sure why we have to include it but this seems to be necessary, c.f. https://stackoverflow.com/questions/33925232/how-to-match-over-self-in-an-enum

        match self {
            Iter1( iter1 ) =>   {
                iter1.next()
            },
            Iter2( iter2 ) =>   {
                iter2.next()
            }
        }
    }
}

impl < I1, I2 > 
        
    Clone       
    for IterTwoType < I1, I2 > 
    where   I1:  Clone + Iterator,
            I2:  Clone + Iterator < Item = I1::Item >,    
    
{
    fn clone(&self) -> Self {
        match self {
            Self::Iter1( iter1 ) => { Self::Iter1( iter1.clone() ) },
            Self::Iter2( iter2 ) => { Self::Iter2( iter2.clone() ) },            
        }
    }
}


// OnlyDuplicates
// ---------------------------------------------------------------------------


/// Returns one copy of each item that repeats one or more times, and skips non-repeated items.
/// 
/// This struct wraps around another iterator.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::OnlyDuplicates;
/// 
/// let a   =   vec![ 1, 2, 2, 3, 3, 3, 4, 5, 5 ];
/// let b   =   vec![ 2, 3, 5 ];
/// 
/// let only_dups       =   OnlyDuplicates::new( a.iter() );
/// 
/// // Check that `only_dups` returns the same sequence of items as `b.iter()`
/// assert!(    only_dups.eq( b.iter() )    );
/// 
/// ```
pub struct OnlyDuplicates< I: Iterator >{
    iter_unwrapped: I,
    next_candidate: Option< I::Item >,
}

impl < I > 

    OnlyDuplicates< I > 

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,
{
    pub fn new( mut iter_unwrapped: I ) -> OnlyDuplicates< I > { 
        let next_candidate  =   iter_unwrapped.next();
        OnlyDuplicates { iter_unwrapped, next_candidate: next_candidate } 
    }
}

impl    < I > 
        
    Iterator for

    OnlyDuplicates< I >

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,        
{
    type Item = I::Item;

    fn next( &mut self ) -> Option< Self::Item > {

        // initialize two unasigned variables
        let mut is_duplicate; 
        let mut following_candidate;

        loop {
            if self.next_candidate.is_none() { return None }

            following_candidate     =   self.iter_unwrapped.next();
            is_duplicate            =   self.next_candidate == following_candidate;

            //  iterate over all duplicates
            while following_candidate == self.next_candidate {
                following_candidate = self.iter_unwrapped.next();
            }

            if is_duplicate{
                // swap bindings of our two variables
                std::mem::swap( &mut self.next_candidate, &mut following_candidate );
                // return the old candidate (whose value is now bound to `following_candidate`)
                return following_candidate
            } else {
                self.next_candidate     =   following_candidate;
            }
        }
    }
}    


// SkipDuplicates
// ---------------------------------------------------------------------------


/// Skips items that repeat any number of times.
/// 
/// This struct wraps around another iterator.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::SkipDuplicates;
/// use itertools;
/// 
/// // EXAMPLE 1
/// let a   =   vec![ 1, 2, 2, 3, 3, 3, 4, 5 ];
/// let b   =   vec![ 1, 4, 5 ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// // EXAMPLE 2
/// let a   =   vec![ 1, 2, 2, 3, 3, 3, 4, 5, 5 ];
/// let b   =   vec![ 1, 4 ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// // EXAMPLE 3
/// let a: Vec<usize>   =   vec![  ];
/// let b: Vec<usize>   =   vec![  ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// // EXAMPLE 4
/// let a: Vec<usize>   =   vec![ 1, 1, 2, 2 ];
/// let b: Vec<usize>   =   vec![  ];
/// 
/// let skip_dups       =   SkipDuplicates::new( a.iter() );
/// 
/// // Check that `skip_dups` returns the same sequence of items as `b.iter()`
/// itertools::assert_equal(    skip_dups,  b.iter()     );
/// 
/// ```
pub struct SkipDuplicates< I: Iterator >{
    iter_unwrapped: I,
    next_candidate: Option< I::Item >,
}

impl < I > 

    SkipDuplicates< I > 

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,
{
    pub fn new( mut iter_unwrapped: I ) -> SkipDuplicates< I > { 
        let next_candidate  =   iter_unwrapped.next();
        SkipDuplicates { iter_unwrapped, next_candidate: next_candidate } 
    }
}

impl    < I > 
        
    Iterator for

    SkipDuplicates< I >

    where 
        I:          Iterator,
        I::Item:    std::cmp::PartialEq,        
{
    type Item = I::Item;

    fn next( &mut self ) -> Option< Self::Item > {

        // initialize two unasigned variables
        let mut is_duplicate; 
        let mut following_candidate;

        loop {
            if self.next_candidate.is_none() { return None }

            following_candidate     =   self.iter_unwrapped.next();
            is_duplicate            =   self.next_candidate == following_candidate;

            //  iterate over all duplicates
            while following_candidate == self.next_candidate {
                following_candidate = self.iter_unwrapped.next();
            }

            if is_duplicate{
                self.next_candidate     =   following_candidate;
            } else {
                // swap bindings of our two variables
                std::mem::swap( &mut self.next_candidate, &mut following_candidate );
                // return the old candidate (whose value is now bound to `following_candidate`)
                return following_candidate
            }
        }
    }
}    


//  SET OPERTIONS ON ORDERED ITERATORS
//  ==================================================================

/// Determine if one ordered iterator contains another.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::ordered_iter_contains;
/// 
/// let a   =   Vec::<usize>::new();
/// let b   =   vec![ 0, 1,     ];
/// let c   =   vec![    1, 2,  ];
/// let d   =   vec![ 0, 1, 2   ];
/// 
/// assert!(   ordered_iter_contains( d.iter(), a.iter() ) );
/// assert!(   ordered_iter_contains( d.iter(), b.iter() ) );
/// assert!(   ordered_iter_contains( d.iter(), c.iter() ) );
/// assert!( ! ordered_iter_contains( c.iter(), b.iter() ) );
/// assert!( ! ordered_iter_contains( b.iter(), c.iter() ) );
/// assert!( ! ordered_iter_contains( a.iter(), d.iter() ) );
/// ```
pub fn ordered_iter_contains<I>( mut iter_large: I, iter_small: I ) -> bool 
    where 
        I:          Iterator,
        I::Item:    Ord + PartialEq + std::fmt::Debug,
{

    // use size_hints;
    for item in iter_small{
        loop{
            match iter_large.next(){
                None => { return false }
                Some( x ) =>{
                    println!("{:?}", (&x, &item));
                    match x.cmp( &item ) {
                        std::cmp::Ordering::Less    =>  { continue }
                        std::cmp::Ordering::Equal   =>  { break },
                        std::cmp::Ordering::Greater =>  { return false }
                    }
                }
            }   
        }
    }
    return true
}

/// Iterates over the intersection of `n` ordered iterators
/// 
/// Each iterator must return elements in strictly ascending order.
/// 
/// # Examples
/// ```
/// use oat_rust::utilities::iterators::general::IntersectOrderedIteratorsUnsafe;
/// 
/// //  TEST 1: intersection of 3 iterators
/// 
/// let ordered_sets    =   vec![
///                             vec![ 0, 1, 2, 3, 4,    ],
///                             vec![       2, 3, 4,    ],
///                             vec![ 0,    2,    4, 5, ],
///                         ];
/// let intersection    =   IntersectOrderedIteratorsUnsafe::new(
///                             ordered_sets.iter().map(|x| x.iter() )
///                         );
/// let v               =   vec![ 2, 4, ];
/// 
/// itertools::assert_equal( v.iter(), intersection );
/// 
/// //  TEST 2: intersection of 0 iterators
/// 
/// let ordered_sets: Vec<Vec<usize>>    =   vec![];
/// let intersection    =   IntersectOrderedIteratorsUnsafe::new(
///                             ordered_sets.iter().map(|x| x.iter() )
///                         );
/// let v: Vec<usize>   =   vec![];
/// 
/// itertools::assert_equal( v.iter(), intersection );
/// ```
pub struct IntersectOrderedIteratorsUnsafe< I > 
    where
        I:          Iterator,
        I::Item:    PartialEq + Ord,
{
    iterators:  Vec< I >,
    next_elts:  Vec< Option< I::Item > >
}


impl < I > 

    IntersectOrderedIteratorsUnsafe< I >

    where
        I:          Iterator,
        I::Item:    Clone + PartialEq + Ord,
{
    pub fn new< T: IntoIterator<Item=I> >( iterators: T ) -> IntersectOrderedIteratorsUnsafe<I> {
        let mut iterators: Vec<_>   =   iterators.into_iter().collect();
        iterators.shrink_to_fit();
        let mut next_elts = Vec::with_capacity( iterators.len() );
        for elt in iterators.iter_mut().map(|x| x.next() ) {
            next_elts.push( elt );
        }
        IntersectOrderedIteratorsUnsafe{ iterators, next_elts }
    }
}        



impl < I > 

    Iterator for 

    IntersectOrderedIteratorsUnsafe< I >

    where
        I:          Iterator,
        I::Item:    PartialEq + Ord,
{
    type Item       =   I::Item;

    fn next( &mut self ) -> Option< Self::Item > {

        // The intersection of 0 iterators is emtpy
        if self.iterators.is_empty() { return None }

        let mut min_index;
        let mut min_val;        
        let mut all_min;
        let mut indices_to_check;

        loop {
            // find the index of the minimum value
            min_index   =   find_min( self.next_elts.iter() ).unwrap();

            // pull out the next element from the iterator at index min_index, and swap its bound value with self.next_elts[min_index]
            min_val     =   self.iterators[ min_index ].next(); // this will become the true min val on the next line
            std::mem::swap(&mut min_val, &mut self.next_elts[min_index]); // now min_val equals the minimum value of next_elts, before we modified next_elts

            // If one of the iterators is empty then (i) min_val = None, and (ii) there are no more elements in the intersection
            if min_val.is_none() { return None }
    
            // Check whether every element equals min, and replace the min vals
            all_min     =   true;
            indices_to_check = (0..min_index).chain( min_index+1 .. self.next_elts.len()); // iterate over every index except min_index
            for p in indices_to_check {
                match self.next_elts[p] == min_val {
                    true    =>  { self.next_elts[p] = self.iterators[p].next()  },
                    false   =>  { all_min = false }
                }
            }

            if all_min { return min_val }
        }
    }

}


/// Returns an iterator that runs over elments in the inersection of `iter1` and `iter2`.
/// 
/// Each iterator must be sorted in **strictly** ascending order.
/// 
/// # Examples
/// 
/// # Opportunities for improvement
/// 
/// terminate search after one iterator terminates
fn intersect_ordered_iterators< I >( iter1: I, iter2: I )
    ->
    OnlyDuplicates< Merge< I, I > >

    where
        I:          Iterator,
        I::Item:    PartialOrd,
{
    OnlyDuplicates::new( iter1.merge( iter2 ))
}

/// Returns an iterator that runs over elments in the symetric difference of `iter1` and `iter2`
/// 
/// Each iterator must be sorted in **strictly** ascending order.
fn symmetric_difference_of_ordered_iterators< I >( iter1: I, iter2: I )
    ->
    OnlyDuplicates< Merge< I, I > >

    where
        I:          Iterator,
        I::Item:    PartialOrd,
{
    OnlyDuplicates::new( iter1.merge( iter2 ))
}



// fn union_of_ordered_iterators()

fn set_difference_of_ordered_iterators<I>( iter1: I, iter2: I ){
    // let iter1 = iter1.map(|x| (x,1));
    // let iter2 = iter2.map(|x| (x,2));    
    // let symmetric_difference = symmetric_difference_of_ordered_iterators_by_cmp(iter1, iter2, |x,y| x.0 == y.0);
    // return symmetric_difference.filter_map(|x| x.1 == 1, |x| x.0)
}





// OncePeekable
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct OncePeekable< T > { item_opt: Option< T > }

impl < T > OncePeekable< T > {
    pub fn new( item: T ) -> OncePeekable< T > { OncePeekable{ item_opt: Some(item) } }
}

impl < T > Iterator for OncePeekable< T > {
    type Item = T;
    fn next( &mut self ) -> Option< Self::Item > {
        std::mem::replace( &mut self.item_opt, None ) // replace the internally stored item with None, and return whatever was stored there previously
    }
}

impl < T > 

    PeekUnqualified for 

    OncePeekable< T >
{
    fn peek_unqualified( &mut self ) -> Option < & Self::Item > { self.item_opt.as_ref() }
}


// VECTOR WRAPPED IN AN ITERATOR
// ---------------------------------------------------------------------------

/// Wrapper containing a `Vec` and a pointer `usize`; iterates over the elements of
/// the vector, incrementing the pointer as it goes.
/// 
/// This struct is essentially a vector that implements `Iterator`.  It differs from 
/// `Vec<T>` because `Vec<T>` does not implement `Iterator`.  Moreover, it differs 
/// from the iterator we usually associate with a vector, `Itera< '_, T >` in that 
/// `IterWrappedVec< T >` owns its own data.  Thus, for example, we can write
/// a function which takes no arguments, but which returns a nonempty 
/// `IterWrappedVec< T >`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::iterators::general::IterWrappedVec;
/// 
/// // error[E0106]: missing lifetime specifier
/// // fn iter_from_nothing() -> std::slice::Iter<'_, usize> { ( vec![ 0, 0, 0 ] ).iter() } 
/// 
/// // no error
/// fn vec_wrapped_in_iter_from_nothing() -> IterWrappedVec< usize > { IterWrappedVec::new( vec![ 0, 0, 0 ] ) } 
/// 
/// let iter = vec_wrapped_in_iter_from_nothing();
/// itertools::assert_equal(iter, vec![0,0,0])
/// ```
pub struct IterWrappedVec< T > { wrapped_vec: Vec< T >, pointer_to_next_val: usize }

impl < T > IterWrappedVec < T > {
    /// Create a new `IterWrappedVec< T >` from a `Vec< T >`.  
    pub fn new( wrapped_vec: Vec< T > ) -> IterWrappedVec < T > { IterWrappedVec{ wrapped_vec, pointer_to_next_val: 0  } }
}

impl < T: Clone > Iterator for IterWrappedVec< T > { 
    type Item = T;

    fn next( &mut self ) -> Option< Self::Item > {
        match self.pointer_to_next_val < self.wrapped_vec.len() {
            true    =>  { 
                let a = self.wrapped_vec[ self.pointer_to_next_val ].clone(); // retreive the next item
                self.pointer_to_next_val += 1; // increment the pointer
                return Some( a ) // return the item
            }
            false   =>  { return None }
        }
    }
 }






//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    use crate::utilities::iterators::general::FilterMapByTransform;
    use crate::utilities::functions::evaluate::EvaluateFunction;


    #[test]
    fn test_skip_until() {
        use crate::utilities::iterators::general::SkipUntil;

        let vec = vec![ 0, 1, 2, 1, 0];
        let iter = vec.iter().cloned().peekable().skip_until( |x| x > &1);
        itertools::assert_equal( iter, vec![ 2, 1, 0 ])        
    }


    #[test]
    fn test_filter_map_by_transform() {
        use crate::utilities::iterators::general::FilterMapByTransform;

        // define a struct to execute the filter-map function
        struct ShiftDown;

        // implement the `FilterMapTransform` trait
        impl EvaluateFunction< usize, Option< i32 > > for ShiftDown{
            fn evaluate_function( &self, input: usize ) -> Option< i32 > {
                match input > 3 {
                    true    => { return None },
                    false   => { return Some( input as i32 - 3 ) }
                }
            }
        }

        // define an iterator
        let iter = 0..9;

        // create a `filter-mapped` object
        let iter_transformed = FilterMapByTransform::new( iter, ShiftDown );

        // verify this object is correct        
        itertools::assert_equal( iter_transformed, vec![ -3, -2, -1, 0 ])        
    }    



    #[test]
    fn test_filter_map_by_transform_2() {
        use crate::utilities::iterators::general::{FilterMapByTransformTrait};
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        hash.insert( 2usize, -2);
        
        let iter_data = vec![0, 1, 2, 3];
        let iter = iter_data.iter().cloned();
        
        let iter_filtermapped = iter.filter_map_by_transform( & hash );
        itertools::assert_equal( iter_filtermapped.cloned(), vec![-1, -2] );  
    }   

    #[test]
    fn test_implementation_of_filtermapobjectmethod_on_hashmap() {
        use crate::utilities::iterators::general::{FilterMapByTransformTrait};
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        hash.insert( 2usize, -2);
        
        let iter_data = vec![0, 1, 2, 3];
        let iter = iter_data.iter().cloned();
        
        let iter_filtermapped = iter.filter_map_by_transform( & hash );
        itertools::assert_equal( iter_filtermapped.cloned(), vec![-1, -2] );       
    }

    #[test]
    fn test_map_by_transform() {
        use crate::utilities::iterators::general::MapByTransform;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;

        let entrywise_transform_unwrapped = |x : usize| -> i32 { x as i32 };
        let entrywise_transform_wrapped = EvaluateFunctionFnMutWrapper::new(entrywise_transform_unwrapped);
        let vec = vec![ 1, 2, 3];
        let iter = vec.iter().cloned();

        let transformed_entrywise = MapByTransform::new( iter, entrywise_transform_wrapped );
        itertools::assert_equal( transformed_entrywise, vec![ 1i32, 2, 3] );        
    }
}    


//  ---------------------------------------------------------------------
//  Unit tests
//  ---------------------------------------------------------------------

#[cfg(test)]
mod unit_tests {
}