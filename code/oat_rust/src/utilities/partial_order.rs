//! Order plays an important role in the oat_rust for several reasons.  Some examples include
//!
//!   - it is easier to to write a lazy iterator for the sum of two sparse vectors (represented as lists of nonzero entires) if the lists are already sorted.
//!   - min-heaps rely on order comparisons to function
//!
//! Often we find that a single application requires us to order a set of elements in one way at one time, and another way at another time.
//! To achieve this without changing the type of elements being ordered, we use a "helper" object that implements traits for
//! comparing two elements -- traits such as [`StrictlyLess`], [`ComparePartialOrder`], and [`CompareOrder`].


use std::cmp::Ordering;


//  ---------------------------------------------------------------------------
//  SORTING
//  ---------------------------------------------------------------------------


/// Returns true iff `order_comparator.strictly_less( a, b ) == true` for every consequtive pair of
/// elements `a`, `b`, in the slice.  Returns `true` by default for slices of length `<= 1`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::partial_order::{ is_sorted_strictly, OrderComparatorAutoAnyType };
/// 
/// assert!( is_sorted_strictly(   &vec![ 0, 1, 2], & OrderComparatorAutoAnyType ) );
/// assert!( ! is_sorted_strictly( &vec![ 0, 1, 1], & OrderComparatorAutoAnyType ) );
/// ```
pub fn is_sorted_strictly < T, OrderComparator, > ( data: &[T], order_comparator: & OrderComparator ) -> bool
    where OrderComparator: StrictlyLess<  T >,
{
    data.windows(2).all(|w| order_comparator.strictly_less( &w[0], &w[1]) )
}



//  ---------------------------------------------------------------------------
//  ORDER COMPARATOR TRAITS
//  ---------------------------------------------------------------------------

//  StrictlyLess
//  ---------------------------------------------------------------------------

//  REMARK: This was copied from the itertools library.  The exhact developers
//  suspect that the trait was introduced in order to supply type information 
//  that ordinary closures might miss.  At least, we tried removing the trait
//  from the HitMerge code we adapted from Itertools, and ran into some type 
//  erros at compilation.

/// A trait with two methods, `strictly_less` and `strictly_greater`, which compare pairs of
/// elements.
pub trait StrictlyLess< T > {
    fn strictly_less(& self, a: &T, b: &T) -> bool;
    fn strictly_greater(& self, a: &T, b: &T) -> bool { self.strictly_less(a, b) }
}


//  DecidePartialOrder + DecideOrder
//  ---------------------------------------------------------------------------


/// Modeled off the Rust trait `std::cmp::PartialOrd`.
pub trait DecidePartialOrder< T > {
    fn decide_partial_cmp(&self,lhs: &T, rhs: &T) -> Option<Ordering> {
        let is_le = self.decide_le(lhs, rhs);
        let is_ge = self.decide_ge(lhs, rhs);
        if is_le {
            if is_ge { Some( Ordering::Equal ) }
            else { Some( Ordering::Less ) }
        }
        else {
            if is_ge { Some( Ordering::Greater ) }
            else { None }
        }
    }

    fn decide_lt(&self, lhs: &T, rhs: &T ) -> bool { self.decide_le( lhs, rhs ) && ! self.decide_ge( lhs, rhs ) }
    fn decide_le(&self,lhs: &T, rhs: &T) -> bool;
    fn decide_gt(&self,lhs: &T, rhs: &T) -> bool { self.decide_lt( rhs, lhs ) }
    fn decide_ge(&self,lhs: &T, rhs: &T) -> bool { self.decide_le( rhs, lhs ) }
}

/// Modeled off the Rust trait `std::cmp::Ord`.
pub trait DecideOrder< T > : DecidePartialOrder< T > 
{
    fn decide_cmp(&self, lhs: &T, rhs: &T ) -> Ordering {
        // match self.decide_lt(&lhs, &rhs) {
        //     true => { Ordering::Less },
        //     false => {  
        //         match self.decide_gt(lhs, rhs) {
        //             true => { Ordering::Greater },
        //             false => { Ordering::Equal }
        //         }
        //     }
        // }
        self.decide_partial_cmp(lhs, rhs).unwrap()
    }

    fn decide_max(&self, lhs: T, rhs: T) -> T { match self.decide_lt(&lhs, &rhs) { true => { rhs }, false => { lhs } } }
    fn decide_min(&self, lhs: T, rhs: T) -> T { match self.decide_lt(&lhs, &rhs) { true => { lhs }, false => { rhs } } }
    fn decide_clamp(&self, clampee: T, min: T, max: T) -> T { 
        if self.decide_lt( &max, &min ) { panic!("Cannot call `self.decide_clamp( clampee, min, max )` when max < min") };

        if self.decide_lt( &clampee, &min ) { return min }
        else if self.decide_lt( &max, &clampee ) { return max }
        else { return clampee }
    }
}


//  ---------------------------------------------------------------------------
//  IMPLEMENTORS OF `DecideOrder`
//  ---------------------------------------------------------------------------

/// Wrapper around any struct `ComparatorStrictlyLess` that implements the trait `StrictlyLess`; the wrapper then implements `DecidePartialOrder` and `DecideOrder`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::partial_order::{InferTotalOrderFromStrictlyLess, OrderComparatorAutoAnyType, DecidePartialOrder, DecideOrder};
/// 
/// let order_comparator      =   InferTotalOrderFromStrictlyLess::new( OrderComparatorAutoAnyType );
/// assert_eq!( order_comparator.decide_le( &2, &2), true );
/// assert_eq!( order_comparator.decide_lt( &2, &2), false);
/// assert_eq!( order_comparator.decide_le( &1, &2), true );
/// assert_eq!( order_comparator.decide_lt( &1, &2), true );
/// assert_eq!( order_comparator.decide_le( &2, &1), false);
/// assert_eq!( order_comparator.decide_lt( &2, &1), false);
/// ```
pub struct InferTotalOrderFromStrictlyLess< ComparatorStrictlyLess >{ comparator_strictly_less: ComparatorStrictlyLess }

impl < ComparatorStrictlyLess > InferTotalOrderFromStrictlyLess< ComparatorStrictlyLess > {
    pub fn new( comparator_strictly_less: ComparatorStrictlyLess ) -> InferTotalOrderFromStrictlyLess< ComparatorStrictlyLess > { InferTotalOrderFromStrictlyLess{ comparator_strictly_less  } }
}

impl < ComparatorStrictlyLess, T > 

    DecidePartialOrder< T > for
    InferTotalOrderFromStrictlyLess< ComparatorStrictlyLess >

    where 
        ComparatorStrictlyLess:     StrictlyLess< T >,
{
    fn decide_le(&self,lhs: &T, rhs: &T) -> bool { ! self.comparator_strictly_less.strictly_less(rhs, lhs) }
    fn decide_lt(&self,lhs: &T, rhs: &T) -> bool { self.comparator_strictly_less.strictly_less(lhs, rhs) }  
}

impl < ComparatorStrictlyLess, T > 

    DecideOrder< T > for
    InferTotalOrderFromStrictlyLess< ComparatorStrictlyLess >

    where 
        ComparatorStrictlyLess:     StrictlyLess< T >,
{}


//  ---------------------------------------------------------------------------
//  IMPLEMENTORS OF `StrictlyLess`
//  ---------------------------------------------------------------------------


//  Reverse order
//  -------------

/// Represents the reverse of a total order
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::partial_order::{ OrderComparatorReverse, OrderComparatorAutoAnyType, StrictlyLess };
/// 
/// let mut order_comaprator_reverse = OrderComparatorReverse::new( OrderComparatorAutoAnyType );
/// assert!( order_comaprator_reverse.strictly_less( &2, &1) )
/// ```
#[derive(Clone, Debug)]
pub struct OrderComparatorReverse< UnreversedOrderComparator > { unreversed_order_comparator: UnreversedOrderComparator }

impl < UnreversedOrderComparator > OrderComparatorReverse< UnreversedOrderComparator > { 
    pub fn new( unreversed_order_comparator: UnreversedOrderComparator ) ->  OrderComparatorReverse< UnreversedOrderComparator > { OrderComparatorReverse{ unreversed_order_comparator } }
}

impl < T, UnreversedOrderComparator: StrictlyLess< T > > 

    StrictlyLess
        < T > for 

    OrderComparatorReverse
        < UnreversedOrderComparator >

{
    fn strictly_less( & self, a: &T, b: &T) -> bool { self.unreversed_order_comparator.strictly_less(b, a) }
}



//  Less-than (auto implementor)
//  ----------------------------

/// Zero-memory struct that arbitrates the (partial) order on structs that implement [PartialOrd](` std::cmp::PartialOrd`)
/// or [Ord](` std::cmp::Ord`); 
#[derive(Clone, Debug)]
pub struct OrderComparatorAutoLt< T: PartialOrd >{ phantom_compared_type: PhantomData< T > }

impl < T: PartialOrd > OrderComparatorAutoLt< T > {
    pub fn new() -> OrderComparatorAutoLt< T > { OrderComparatorAutoLt{ phantom_compared_type: PhantomData } }
}

impl<T: PartialOrd> StrictlyLess< T > for OrderComparatorAutoLt< T > {
    fn strictly_less( & self, a: &T, b: &T) -> bool { a < b }
}

impl<T: PartialOrd> DecidePartialOrder< T > for OrderComparatorAutoLt< T > {
    fn decide_le(& self, lhs: &T, rhs: &T) -> bool { lhs <= rhs }
    fn decide_lt(& self, lhs: &T, rhs: &T) -> bool { lhs <  rhs }    
}

impl<T: PartialOrd + Ord> DecideOrder< T > for OrderComparatorAutoLt< T > {} // all methods in this trait are auto-implemented


//  Less-than (auto implementor, any type)
//  ----------------------------

/// Zero-memory struct that arbitrates the (partial) order on structs that implement [PartialOrd](` std::cmp::PartialOrd`)
/// or [Ord](` std::cmp::Ord`); 
#[derive(Clone, Debug)]
pub struct OrderComparatorAutoAnyType;

impl OrderComparatorAutoAnyType {
    pub fn new() -> OrderComparatorAutoAnyType { OrderComparatorAutoAnyType }
}

impl<T: PartialOrd> StrictlyLess< T > for OrderComparatorAutoAnyType {
    fn strictly_less( & self, a: &T, b: &T) -> bool { a < b }
}

impl<T: PartialOrd> DecidePartialOrder< T > for OrderComparatorAutoAnyType {
    fn decide_le(& self, lhs: &T, rhs: &T) -> bool { lhs <= rhs }
    fn decide_lt(& self, lhs: &T, rhs: &T) -> bool { lhs <  rhs }    
}

impl<T: PartialOrd + Ord> DecideOrder< T > for OrderComparatorAutoAnyType {} // all methods in this trait are auto-implemented



//  Greater-than (auto implementor)
//  -------------------------------

/// Zero-memory struct that arbitrates the OPPOSITE (partial) order on structs that implement [PartialOrd](` std::cmp::PartialOrd`)
/// or [Ord](` std::cmp::Ord`); 
#[derive(Clone, Debug)]
pub struct OrderComparatorAutoGt;

impl OrderComparatorAutoGt {
    pub fn new() -> OrderComparatorAutoGt { OrderComparatorAutoGt }
}

impl<T: PartialOrd> StrictlyLess< T> for OrderComparatorAutoGt {
    fn strictly_less( & self, a: &T, b: &T) -> bool {
        a > b
    }
}

impl<T: PartialOrd> DecidePartialOrder< T > for OrderComparatorAutoGt {
    fn decide_le(& self, lhs: &T, rhs: &T) -> bool { lhs >= rhs }
    fn decide_lt(& self, lhs: &T, rhs: &T) -> bool { lhs >  rhs }    
}

impl<T: PartialOrd + Ord> DecideOrder< T > for OrderComparatorAutoGt {} // all methods in this trait are auto-implemented



//  Less-than-by-key (auto implementor)
//  -----------------------------------

/// Zero-memory struct representing the "strictly less than" relation for objects that implement [`KeyValGet<Key, Val>`](`KeyValGet`).
/// Calling `strictly_less( &x, &y, )` returns true iff `x.key() < y.key()`.  The `x.key()` struct must implement
/// [PartialOrd](`PartialOrd`).  
/// 

#[derive(Clone, Debug)]
pub struct OrderComparatorAutoLtByKey< Key, Val, KeyValPair >{
    phantom_key:    PhantomData < Key >,
    phantom_val:    PhantomData < Val >,
    phantom_pair:   PhantomData < KeyValPair >,        
}

// Constructor
impl< Key, Val, KeyValPair > OrderComparatorAutoLtByKey< Key, Val, KeyValPair >{
    pub fn new() -> OrderComparatorAutoLtByKey< Key, Val, KeyValPair > { 
        OrderComparatorAutoLtByKey{
                phantom_key:    PhantomData,
                phantom_val:    PhantomData,
                phantom_pair:   PhantomData, 
            }         
    }
}

// StrictlyLess trait implementation
impl< Key, Val, KeyValPair >

    StrictlyLess 
        < KeyValPair > for 
    
    OrderComparatorAutoLtByKey 
        < Key, Val, KeyValPair >

    where   Key:            PartialOrd,
            KeyValPair:     KeyValGet < Key, Val, >
    
    {

    fn strictly_less( & self, a: & KeyValPair, b: & KeyValPair) -> bool {
        a.key() < b.key()
    }
}



//  Less-than-by-key (auto implementor for 2-tuples)
//  -----------------------------------

/// Zero-memory struct representing the "strictly less than" relation for tuples of length 2.
/// Calling `strictly_less( &x, &y, )` returns true iff `x.0 < y.0`.  The `x.0` struct must implement
/// [PartialOrd](`PartialOrd`).  
/// 
/// This struct is strictly less general than `OrderComparatorAutoLtByKey`, but it seems to
/// give the compiler fewer worries over trait implementation requirements.
#[derive(Clone, Debug)]
pub struct OrderComparatorAutoLtByFirstTupleEntry{}

// Constructor
impl OrderComparatorAutoLtByFirstTupleEntry{
    pub fn new() -> OrderComparatorAutoLtByFirstTupleEntry{    OrderComparatorAutoLtByFirstTupleEntry{}    }
}

// StrictlyLess trait implementation
impl< Key, Val >

    StrictlyLess 
        < (Key, Val) > for 
    
    OrderComparatorAutoLtByFirstTupleEntry //< Key, Val >

    where   
        Key:            PartialOrd,
    {
    fn strictly_less( & self, a: & (Key, Val), b: & (Key, Val) ) -> bool {    a.0 < b.0    }
}


//  Greater-than-by-key (auto implementor)
//  --------------------------------------

/// Zero-memory struct representing the "greater than" relation for objects that implement [`KeyValGet<Key, Val>`](`KeyValGet`).
/// Calling `strictly_less( &x, &y, )` returns true `x.key() > y.key()`.  The `x.key()` struct must implement
/// [PartialOrd](`PartialOrd`).  
#[derive(Clone, Debug)]
pub struct OrderComparatorAutoGtByKey < Key, Val, KeyValPair > {
    phantom_key:    PhantomData < Key >,
    phantom_val:    PhantomData < Val >,
    phantom_pair:   PhantomData < KeyValPair >,        
}

// Constructor
impl< Key, Val, KeyValPair > OrderComparatorAutoGtByKey< Key, Val, KeyValPair >{
    pub fn new() -> OrderComparatorAutoGtByKey< Key, Val, KeyValPair > { 
        OrderComparatorAutoGtByKey{
                phantom_key:    PhantomData,
                phantom_val:    PhantomData,
                phantom_pair:   PhantomData, 
            }         
    }
}

// StrictlyLess trait implementation
impl< Key, Val, KeyValPair >

    StrictlyLess 
        < KeyValPair >
    
    for OrderComparatorAutoGtByKey 
        < Key, Val, KeyValPair >

    where   Key:            PartialOrd,
            KeyValPair:     KeyValGet < Key, Val, >
    
    {

    fn strictly_less( & self, a: & KeyValPair, b: & KeyValPair) -> bool {
        a.key() > b.key()
    }
}



//  Less-than-by-key-WITH-custom-comparator (auto implementor)
//  --------------------------------------

/// Zero-memory struct representing the "strictly less than" relation for objects that implement [`KeyValGet<Key, Val>`](`KeyValGet`).
/// Calling `strictly_less( &x, &y, )` returns true if `x.key() < y.key()`.
#[derive(Clone, Debug)]
pub struct OrderComparatorLtByKey < Key, Val, KeyValPair, KeyComparator > 
{
    key_comparator: KeyComparator,
    phantom_key:    PhantomData < Key >,
    phantom_val:    PhantomData < Val >,
    phantom_pair:   PhantomData < KeyValPair >,        
}

// Constructor
impl< Key, Val, KeyValPair, KeyComparator > OrderComparatorLtByKey < Key, Val, KeyValPair, KeyComparator >  {
    pub fn new( key_comparator: KeyComparator ) -> OrderComparatorLtByKey < Key, Val, KeyValPair, KeyComparator >  {
        OrderComparatorLtByKey{ key_comparator: key_comparator, phantom_key: PhantomData, phantom_pair: PhantomData, phantom_val: PhantomData }
    }
}

// StrictlyLess trait implementation
impl< Key, Val, KeyValPair, KeyComparator >

    StrictlyLess 
        < KeyValPair > for 
    
    OrderComparatorLtByKey 
        < Key, Val, KeyValPair, KeyComparator > where

    KeyValPair:     KeyValGet < Key, Val, >,
    KeyComparator:  StrictlyLess <  Key >,
    
    {
    fn strictly_less( & self, a: & KeyValPair, b: & KeyValPair) -> bool {
        self.key_comparator.strictly_less(&a.key(),  &b.key() )
    }
}



//  Mutable closure 
//  ---------------------------------------------

use std::marker::PhantomData;

use crate::entries::KeyValGet;


/// A wrapper for objects that implement `FnMut(&T, &T)->bool`, which implements `StrictlyLess<  T >`.
#[derive(Clone, Debug)]
pub struct OrderComparatorFnWrapper
                < F: Fn(&T, &T)->bool, T > 
    { order_comparator_unwrapped: F, phantom_comparatee: PhantomData< T > }

impl < F: Fn(&T, &T)->bool, T > OrderComparatorFnWrapper< F, T > { 
    pub fn new( order_comparator_unwrapped: F ) -> OrderComparatorFnWrapper< F, T > 
    { OrderComparatorFnWrapper{ order_comparator_unwrapped, phantom_comparatee: PhantomData } } 
}

impl<T, F: Fn(&T, &T)->bool> 

    StrictlyLess
        <T> for 
    OrderComparatorFnWrapper
        < F, T > 
{
    fn strictly_less( & self, a: &T, b: &T) -> bool {
        (self.order_comparator_unwrapped)(a, b)
    }
}






