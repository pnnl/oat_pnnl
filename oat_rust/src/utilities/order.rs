//! Partial and total orders
//! 
//! See also [`indexing_and_bijection`](crate::utilities::indexing_and_bijection).
//! 
//! Rust offers two traits to compare the order of elements: [PartialOrd] and [Ord]. For example, if you want to know if
//! `a` is less than `b`, you can run
//! 
//! ```
//! let a = 0; 
//! let b = 1;
//! 
//! // This will assign a value of `true` to `is_less`.
//! let is_less =   a.lt(&b);
//! ```
//! 
//! However, [PartialOrd] and [Ord] offer limited flexibility,
//! especially in situations where multiple different orders might apply to the same type of object (for example, we might want to
//! sort simplices according to dimension, diameter, or lexicographic order, or some combination of the three).
//! 
//! To provide greater flexibility, OAT provides two similar but different traits, [`JudgePartialOrder`] and [`JudgeOrder`].
//! The functions defined by these traits are similar to [PartialOrd] and [Ord], with one important difference: rather than
//! compare to objects directly, you define a third object which performs the comparison for you.
//! 
//! ```
//! use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorAutoReverse};
//! use crate::oat_rust::utilities::order::JudgePartialOrder;
//! 
//! // Define a and b
//! let a = 0;
//! let b = 1;
//! 
//! // Define an order operator:
//! // this one impelements the default order on objects
//! let order_operator = OrderOperatorAuto;
//! 
//! // Use the order operator to compare a and b
//! assert!( order_operator.judge_lt( &a, &b ) );
//! 
//! // Define a DIFFERENT order operator
//! let order_operator_reverse = OrderOperatorAutoReverse;
//! 
//! // Use the new operator to compare a and b with a DIFFERENT order:
//! // this particular operator happens to reverse the default order on objects,
//! // so it will assert that b is less than a
//! assert!( order_operator_reverse.judge_lt( &b, &a ) );
//! ```
//! 
//! This gives greater
//! flexibility, since you can define many different order operators to compare elements of the same type.
//! 
//! # Terminology
//! 
//! We call an object that implements [JudgePartialOrder] and/or [JudgeOrder] an **order operator**.
//! It is analogous to the idea of a [RingOperator](crate::algebra::rings), in the sense that you use
//! order operators to get order information about other objects, just like you use ring operators to get
//! algebraic information about other structures.
//! 
//! # Design notes
//! 
//! In many situations [OrderOperatorAuto] is this simplest to use.  [OrderOperatorByKey] is also common for comparing vector entries; it has a smaller type signature than [OrderOperatorByKeyCutsom], but for many operations, e.g. [U-match factorization](crate::algebra::matrices::operations::umatch), the larger signature allows users to pass customized order operators.
// //! object that "externalizes" that implements traits for
// //! comparing two elements -- traits such as [`JudgePartialOrder`], [`ComparePartialOrder`], and [`CompareOrder`].
// //! 
// //! Order plays an important role in the OAT for several reasons.  Some examples include
// //!
// //!   - it is easier to to write a lazy iterator for the sum of two sparse vectors (represented as lists of nonzero entires) if the lists are already sorted.
// //!   - min-heaps rely on order comparisons to function
// //!
// //! Often we find that a single application requires us to order a set of elements in one way at one time, and another way at another time.
// //! 
// //! To achieve this without changing the type of elements being ordered, we use a "helper" object that implements traits for
// //! comparing two elements -- traits such as [`JudgePartialOrder`], [`ComparePartialOrder`], and [`CompareOrder`].


use std::cmp::Ordering;


//  ---------------------------------------------------------------------------
//  SORTING
//  ---------------------------------------------------------------------------

//  SEE ALSO: the sequences_and_ordinals module
//  -------------------------------------------


/// Returns true iff `order_operator.judge_lt( a, b ) == true` for every consequtive pair of
/// elements `a`, `b`, in the slice.  Returns `true` by default for slices of length `<= 1`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::order::{ is_sorted_strictly, OrderOperatorAuto };
/// 
/// assert!( is_sorted_strictly(   &vec![ 0, 1, 2], & OrderOperatorAuto ) );
/// assert!( ! is_sorted_strictly( &vec![ 0, 1, 1], & OrderOperatorAuto ) );
/// ```
pub fn is_sorted_strictly < T, OrderOperator, > ( data: &[T], order_operator: & OrderOperator ) -> bool
    where OrderOperator: JudgePartialOrder<  T >,
{
    data.windows(2).all(|w| order_operator.judge_lt( &w[0], &w[1]) )
}



//  ---------------------------------------------------------------------------
//  ORDER OPERATOR TRAITS
//  ---------------------------------------------------------------------------


//  JudgePartialOrder + JudgeOrder
//  ---------------------------------------------------------------------------


/// Modeled off the Rust trait `std::cmp::PartialOrd`.
/// 
/// It is commonly the case that multiple *different* partial orders are
/// relevant to struct.  For example, one might want to order simplices in
/// ascending/descending lexicographic order, or first by dimension and second
/// by lexicographic order, etc.  The Rust traits for `Ord` and `PartialOrd`
/// only accomodate a single order per struct.  Therefore we create third-party
/// objects or "operators" to operationalize different orders.
/// 
///  REMARK: This trait is similar to (and extends) one copied from the itertools 
///  library for `kmerge`.  The OAT developers
///  suspect that the trait was introduced in order to supply type information 
///  that ordinary closures might miss.  At least, we tried removing the trait
///  from the HitMerge code we adapted from Itertools, and ran into some type 
///  erros at compilation.
pub trait JudgePartialOrder< T > {
    /// Similar to the `.cmp` method in `std::cmp::PartialOrd`
    fn judge_partial_cmp
        ( & self, lhs: & T, rhs: & T ) 
        -> 
        Option<Ordering>;    
    // {
    //     let is_le = self.le(lhs, rhs);
    //     let is_ge = self.ge(lhs, rhs);
    //     if is_le {
    //         if is_ge { Some( Ordering::Equal ) }
    //         else { Some( Ordering::Less ) }
    //     }
    //     else {
    //         if is_ge { Some( Ordering::Greater ) }
    //         else { None }
    //     }
    // }
    /// Returns true iff lhs < rhs
    fn judge_lt(&self, lhs: &T, rhs: &T ) -> bool { self.judge_partial_cmp(lhs, rhs) == Some(Ordering::Less) }// { self.le( lhs, rhs ) && ! self.ge( lhs, rhs ) }
    /// Returns true iff lhs ≤ rhs    
    fn judge_le(&self,lhs: &T, rhs: &T) -> bool { match self.judge_partial_cmp(lhs,rhs){ Some(Ordering::Less) => true, Some(Ordering::Equal) => true, _ => false } }
    /// Returns true iff lhs > rhs    
    fn judge_gt(&self,lhs: &T, rhs: &T) -> bool { self.judge_partial_cmp(lhs, rhs) == Some(Ordering::Greater) } // { self.judge_lt( rhs, lhs ) }
    /// Returns true iff lhs ≥ rhs    
    fn judge_ge(&self,lhs: &T, rhs: &T) -> bool { match self.judge_partial_cmp(lhs,rhs){ Some(Ordering::Greater) => true, Some(Ordering::Equal) => true, _ => false } }
}

/// Modeled off the Rust trait `std::cmp::Ord`.
/// 
/// It is commonly the case that multiple *different* partial orders are
/// relevant to struct.  For example, one might want to order simplices in
/// ascending/descending lexicographic order, or first by dimension and second
/// by lexicographic order, etc.  The Rust traits for `Ord` and `PartialOrd`
/// only accomodate a single order per struct.  Therefore we create third-party
/// objects or "operators" to operationalize different orders.
pub trait JudgeOrder< T > : JudgePartialOrder< T > 
{
    /// Similar to the `.cmp` method in `std::cmp::Ord`
    fn judge_cmp(&self, lhs: &T, rhs: &T ) -> Ordering {
        self.judge_partial_cmp(lhs, rhs).unwrap()
    }
    /// Return the greater of lhs, rhs; if they are equal, reutrh lhs
    fn judge_max(&self, lhs: T, rhs: T) -> T { match self.judge_lt(&lhs, &rhs) { true => { rhs }, false => { lhs } } }
    /// Return the lesser of lhs, rhx; if they are equal, return lhs
    fn judge_min(&self, lhs: T, rhs: T) -> T { match self.judge_le(&lhs, &rhs) { true => { lhs }, false => { rhs } } }
    /// Similar to Rust's `num::clamp`    
    fn judge_clamp(&self, clampee: T, min: T, max: T) -> T { 
        if self.judge_lt( &max, &min ) { panic!("Cannot call `self.clamp( clampee, min, max )` when max < min") };

        if self.judge_lt( &clampee, &min ) { min }
        else if self.judge_lt( &max, &clampee ) { return max }
        else { return clampee }
    }
}


//  ---------------------------------------------------------------------------
//  NONE GREATER
//  ---------------------------------------------------------------------------


/// Wrapper for `Option<T>`; implements `Ord` and `PartialOrd`, where `None` is the greatest (as opposed to the least) element
/// 
/// Compare this with `Option<T>`, where `None` is the least element
/// 
/// # Example
/// 
/// ```
/// use oat_rust::utilities::order::NoneGreater;
/// use std::cmp::Ordering;
/// 
/// let a   =   NoneGreater::from_val(0usize);
/// let b   =   NoneGreater::from_val(1usize);
/// let c   =   NoneGreater::from_opt(None);
/// 
/// assert!( a <  b );
/// assert!( a == a );
/// assert!( a <  c );
/// assert!( c == c );
/// 
/// 
/// assert_eq!( a.partial_cmp(&a), Some(Ordering::Equal) );
/// assert_eq!( a.partial_cmp(&b), Some(Ordering::Less ) );
/// assert_eq!( a.partial_cmp(&c), Some(Ordering::Less ) );
/// assert_eq!( c.partial_cmp(&c), Some(Ordering::Equal) );
/// 
/// ```
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct NoneGreater< T > { pub opt: Option< T > }

impl < T > NoneGreater< T > {
    pub fn new( opt: Option<T> ) -> Self { NoneGreater{opt} }    
    pub fn from_opt( opt: Option<T> ) -> Self { NoneGreater{opt} }
    pub fn from_val( val: T ) -> Self{ NoneGreater{ opt: Some(val) } }
    pub fn val( &self ) -> T where T: Clone { self.opt.clone().unwrap() }
}

impl < T: Ord + PartialOrd + Eq > Ord for NoneGreater< T > {
    

    fn cmp(&self, other: &Self) -> Ordering {
        match (&self.opt, &other.opt) {
            (Some(a),Some(b)) => { a.cmp(b) }            
            (None,   None   ) => { Ordering::Equal   }
            (None,   Some(_)) => { Ordering::Greater }
            (Some(_),None   ) => { Ordering::Less    }             
        }
        // match self.opt {
        //     Some( self_val ) => {
        //         match other {
        //             Some( other_val ) => {

        //             }
        //         }
        //     }
        // }
    }
}

impl < T: Ord + PartialOrd > PartialOrd for NoneGreater< T > {
    
    /// ```
    /// use oat_rust::utilities::order::NoneGreater;
    /// 
    /// let a   =   NoneGreater::from_val(0usize);
    /// let b   =   NoneGreater::from_val(1usize);
    /// let c   =   NoneGreater::from_opt(None);
    /// 
    /// assert!( a <  b );
    /// assert!( a == a );
    /// assert!( a <  c );
    /// assert!( c == c );
    /// 
    /// ```
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (&self.opt, &other.opt) {
            (Some(a),Some(b)) => { a.partial_cmp(b) }            
            (None,   None   ) => { Some(Ordering::Equal  ) }
            (None,   Some(_)) => { Some(Ordering::Greater) }
            (Some(_),None   ) => { Some(Ordering::Less   ) }            
        }        
    }
}


//  ---------------------------------------------------------------------------
//  IMPLEMENTORS OF `JudgeOrder`
//  ---------------------------------------------------------------------------

//  WORKS FINE, BUT WE HAVE RESTURCTURED THE ORDER OPERATOR TRAITS, AND THIS STRUCT NO LONGER PROVIDES MUCH FUNCTIONALITY; OK TO DELETE

// /// Wrapper around any struct `ComparatorJudgePartialOrder` that implements the trait `JudgePartialOrder`; the wrapper then implements `JudgePartialOrder` and `JudgeOrder`.
// /// 
// /// # Examples
// /// 
// /// ```
// /// use oat_rust::utilities::order::{InferTotalOrderFromJudgePartialOrder, OrderOperatorAuto, JudgePartialOrder, JudgeOrder};
// /// 
// /// let order_operator      =   InferTotalOrderFromJudgePartialOrder::new( OrderOperatorAuto );
// /// assert_eq!( order_operator.le( &2, &2), true );
// /// assert_eq!( order_operator.judge_lt( &2, &2), false);
// /// assert_eq!( order_operator.le( &1, &2), true );
// /// assert_eq!( order_operator.judge_lt( &1, &2), true );
// /// assert_eq!( order_operator.le( &2, &1), false);
// /// assert_eq!( order_operator.judge_lt( &2, &1), false);
// /// ```
// #[derive(Clone, Copy, Debug, Eq, PartialEq)]
// pub struct InferTotalOrderFromJudgePartialOrder
//             < ComparatorJudgePartialOrder >
//             { comparator_lt: ComparatorJudgePartialOrder }

// impl < ComparatorJudgePartialOrder > 

//     InferTotalOrderFromJudgePartialOrder
//         < ComparatorJudgePartialOrder > 
// {
//     pub fn new( comparator_lt: ComparatorJudgePartialOrder ) -> InferTotalOrderFromJudgePartialOrder< ComparatorJudgePartialOrder > 
//         { InferTotalOrderFromJudgePartialOrder{ comparator_lt  } }
// }

// impl < ComparatorJudgePartialOrder, T > 

//     JudgePartialOrder
//         < T > for
//     InferTotalOrderFromJudgePartialOrder
//         < ComparatorJudgePartialOrder >

//     where 
//         ComparatorJudgePartialOrder:     JudgePartialOrder< T >,
// {
//     fn le(&self,lhs: &T, rhs: &T) -> bool { ! self.comparator_lt.judge_lt(rhs, lhs) }
//     fn lt(&self,lhs: &T, rhs: &T) -> bool { self.comparator_lt.judge_lt(lhs, rhs) }  
// }

// impl < ComparatorJudgePartialOrder, T > 

//     JudgeOrder
//         < T > for
//     InferTotalOrderFromJudgePartialOrder
//         < ComparatorJudgePartialOrder >

//     where 
//         ComparatorJudgePartialOrder:     JudgePartialOrder< T >,
// {}


//  ---------------------------------------------------------------------------
//  IMPLEMENTORS OF `JudgePartialOrder`
//  ---------------------------------------------------------------------------


//  Reverse order
//  -------------

/// Represents the reverse of a total order
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::order::{ ReverseOrder, OrderOperatorAuto, JudgePartialOrder };
/// 
/// let mut order_comaprator_reverse = ReverseOrder::new( OrderOperatorAuto );
/// assert!( order_comaprator_reverse.judge_lt( &2, &1) )
/// ```
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReverseOrder
            < UnreversedOrderOperator > 
{ unreversed_order_operator: UnreversedOrderOperator }

impl < UnreversedOrderOperator > 

    ReverseOrder
        < UnreversedOrderOperator > 
{ 
    pub fn new( unreversed_order_operator: UnreversedOrderOperator ) ->  ReverseOrder< UnreversedOrderOperator > { ReverseOrder{ unreversed_order_operator } }
}

impl < T, UnreversedOrderOperator: JudgePartialOrder< T > > 

    JudgePartialOrder
        < T > for 

    ReverseOrder
        < UnreversedOrderOperator >

{
    fn judge_partial_cmp( & self, lhs: & T, rhs: & T ) -> Option<Ordering> { self.unreversed_order_operator.judge_partial_cmp( rhs, lhs ) }
}

impl < T, UnreversedOrderOperator: JudgePartialOrder< T > + JudgeOrder< T > > 

    JudgeOrder
        < T > for 

    ReverseOrder
        < UnreversedOrderOperator >

{
    fn judge_cmp( & self, lhs: & T, rhs: & T ) -> Ordering { self.unreversed_order_operator.judge_cmp( rhs, lhs ) }
}

/// A convenience trait which allows the user to reverse an order operator via `operator.reverse_order()`
pub trait IntoReverseOrder 
    where Self: Sized
{
    /// Returns an inverted order operator, consuming `self`
    fn into_reverse_order( self ) -> ReverseOrder< Self > { ReverseOrder::new(self) } // for reasons I don't understand, this auto-implementation doesn't seem to work, but the implementation directly below seems to do the trick
}

impl < T: Sized> IntoReverseOrder for T 
{
    /// Returns an inverted order operator, consuming `self`
    fn into_reverse_order( self ) -> ReverseOrder< Self > { ReverseOrder::new(self) }  
}



// //  Less-than (auto implementor)
// //  ----------------------------

// /// Zero-memory struct that arbitrates the (partial) order on structs that implement [PartialOrd]
// /// or [Ord]; 
// #[derive(Debug, Eq, PartialEq)]
// pub struct OrderOperatorAutoLt< T: PartialOrd >{ 
//     phantom_compared_type: PhantomData< T > 
// }

// impl < T: PartialOrd > OrderOperatorAutoLt< T > {
//     pub fn new() -> OrderOperatorAutoLt< T > { OrderOperatorAutoLt{ phantom_compared_type: PhantomData } }
// }

// // Clone
// impl <T> 

//     Clone for 

//     OrderOperatorAutoLt
//         < T > where 

//     T:  PartialOrd
    
// {
//     fn clone(&self) -> Self { Self::new() }
// } 

// // Copy
// impl <T> 

//     Copy for 

//     OrderOperatorAutoLt
//         < T > where 
    
//     T:  PartialOrd        
// {} 

// impl<T: PartialOrd> JudgePartialOrder< T > for OrderOperatorAutoLt< T > {
//     fn judge_partial_cmp( & self, a: &T, b: &T) -> Option<Ordering> { a.partial_cmp(b) }
// }

// // impl<T: PartialOrd> JudgePartialOrder< T > for OrderOperatorAutoLt< T > {
// //     fn le(& self, lhs: &T, rhs: &T) -> bool { lhs <= rhs }
// //     fn lt(& self, lhs: &T, rhs: &T) -> bool { lhs <  rhs }    
// // }

// impl<T: PartialOrd + Ord> JudgeOrder< T > for OrderOperatorAutoLt< T > {} // all methods in this trait are auto-implemented


//  Less-than (auto implementor, any type)
//  ----------------------------

/// Applies the default order to any object that has a default order
/// 
/// That is, you can use this struct to compare the order of `a` and `b`, whenever
/// `a` and `b` implement [PartialOrd]
/// or [Ord]; 
/// 
/// This struct uses zero memory.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OrderOperatorAuto;

impl OrderOperatorAuto {
    pub fn new() -> OrderOperatorAuto { OrderOperatorAuto }
}

impl<T: PartialOrd> JudgePartialOrder< T > for OrderOperatorAuto {
    
    fn judge_partial_cmp
            ( & self, lhs: & T, rhs: & T ) 
            -> 
            Option<Ordering> {
        lhs.partial_cmp( rhs )
    }
}

impl<T: PartialOrd + Ord> 

    JudgeOrder
        < T > for 
        
    OrderOperatorAuto {} // all methods in this trait are auto-implemented



//  Greater-than (auto implementor)
//  -------------------------------

/// Applies the REVERSE of the default order to any object that has a default order
/// 
/// If `a < b`, then this struct will say that `b < a`.
/// 
/// You can use this struct to compare the order of `a` and `b`, whenever
/// `a` and `b` implement [PartialOrd]
/// or [Ord].
/// 
/// This struct uses zero memory.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OrderOperatorAutoReverse;

impl OrderOperatorAutoReverse {
    pub fn new() -> OrderOperatorAutoReverse { OrderOperatorAutoReverse }
}

impl<T: PartialOrd> JudgePartialOrder< T > for OrderOperatorAutoReverse {
    fn judge_partial_cmp
            ( & self, lhs: & T, rhs: & T ) 
            -> 
            Option<Ordering> {
        rhs.partial_cmp(lhs)
    }
}

impl<T: PartialOrd + Ord> JudgeOrder< T > for OrderOperatorAutoReverse {} // all methods in this trait are auto-implemented



//  Less-than-by-key (auto implementor)
//  -----------------------------------

/// Determines the order of key-value pairs according to their keys
/// 
/// You can use this struct to compare elements of type `T`, whenever 
/// - `T` implements [`KeyValGet<Key, Val>`](`KeyValGet`), and
/// - `Key` implements [PartialOrd] or [Ord]
/// 
/// This struct will indicate that `x < y` whenever `x.key() < y.key()`
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq)]
pub struct OrderOperatorByKey< Key, Val, KeyValPair >{
    phantom_key:    PhantomData < Key >,
    phantom_val:    PhantomData < Val >,
    phantom_pair:   PhantomData < KeyValPair >,        
}

// Constructor
impl< Key, Val, KeyValPair > OrderOperatorByKey< Key, Val, KeyValPair >{
    pub fn new() -> OrderOperatorByKey< Key, Val, KeyValPair > { 
        OrderOperatorByKey{
                phantom_key:    PhantomData,
                phantom_val:    PhantomData,
                phantom_pair:   PhantomData, 
            }         
    }
}

// Copy
impl< Key, Val, KeyValPair >

    Clone for 
    
    OrderOperatorByKey 
        < Key, Val, KeyValPair >
{
    fn clone(&self) -> Self { Self::new() }
} 

// Copy
impl< Key, Val, KeyValPair >

    Copy for 
    
    OrderOperatorByKey 
        < Key, Val, KeyValPair >
{}        

// JudgePartialOrder trait implementation
impl< Key, Val, KeyValPair >

    JudgePartialOrder 
        < KeyValPair > for 
    
    OrderOperatorByKey 
        < Key, Val, KeyValPair >

    where   Key:            PartialOrd,
            KeyValPair:     KeyValGet < Key, Val, >
    
    {

    fn judge_partial_cmp
            ( & self, lhs: & KeyValPair, rhs: & KeyValPair ) 
            -> 
            Option<Ordering> {
        ( lhs.key() ).partial_cmp( &rhs.key() )
    }
}

// JudgePartialOrder trait implementation
impl< Key, Val, KeyValPair >

    JudgeOrder 
        < KeyValPair > for 
    
    OrderOperatorByKey 
        < Key, Val, KeyValPair >

    where   Key:            Ord + PartialOrd,
            KeyValPair:     KeyValGet < Key, Val, >
    
    {

    fn judge_cmp
            ( & self, lhs: & KeyValPair, rhs: & KeyValPair ) 
            -> 
            Ordering {
        ( lhs.key() ).cmp( &rhs.key() )
    }
}


//  Less-than-by-key (auto implementor for 2-tuples)
//  -----------------------------------


/// Determines the order of key-value pairs of type `(Key,Val)` according to `Key`
/// 
/// You can use this struct to compare elements of type `(Key,Val)`, whenever `Key` implements
/// [PartialOrd] or [Ord].
/// 
/// This struct will indicate that `(a,x) < (b,y)` whenever `a < b`. It is a little less general
/// than `OrderOperatorByKey`, but it may be more efficient. This is because `OrderOperatorByKey`
/// will clone the first entry of a tuple each time it makes a comparison, whereas `OrderOperatorByKeyTuple`
/// will only make a reference.
/// 
/// This struct uses zero memory.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OrderOperatorByKeyTuple{}

// Constructor
impl OrderOperatorByKeyTuple{
    pub fn new() -> OrderOperatorByKeyTuple{    OrderOperatorByKeyTuple{}    }
}

// JudgePartialOrder trait implementation
impl< Key, Val >

    JudgePartialOrder 
        < (Key, Val) > for 
    
    OrderOperatorByKeyTuple //< Key, Val >

    where   
        Key:            PartialOrd,
    {
    fn judge_partial_cmp
            ( & self, lhs: & (Key, Val), rhs: & (Key, Val) ) 
            -> 
            Option<Ordering> {
        lhs.0.partial_cmp( &rhs.0)
    }
}


//  Greater-than-by-key (auto implementor)
//  --------------------------------------


/// Determines the order of key-value pairs according to their keys (but reversed!)
/// 
/// You can use this struct to compare elements of type `T`, whenever 
/// - `T` implements [`KeyValGet<Key, Val>`](`KeyValGet`), and
/// - `Key` implements [PartialOrd] or [Ord]
/// 
/// This struct will indicate that `x < y` whenever `x.key() > y.key()`
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq)]
pub struct OrderOperatorByKeyReverse < Key, Val, KeyValPair > {
    phantom_key:    PhantomData < Key >,
    phantom_val:    PhantomData < Val >,
    phantom_pair:   PhantomData < KeyValPair >,        
}

// Copy
impl< Key, Val, KeyValPair >

    Clone for 
    
    OrderOperatorByKeyReverse 
        < Key, Val, KeyValPair >
{
    fn clone(&self) -> Self { Self::new() }
} 

// Copy
impl< Key, Val, KeyValPair >

    Copy for 
    
    OrderOperatorByKeyReverse 
        < Key, Val, KeyValPair >
{}   

// Constructor
impl< Key, Val, KeyValPair > OrderOperatorByKeyReverse< Key, Val, KeyValPair >{
    pub fn new() -> OrderOperatorByKeyReverse< Key, Val, KeyValPair > { 
        OrderOperatorByKeyReverse{
                phantom_key:    PhantomData,
                phantom_val:    PhantomData,
                phantom_pair:   PhantomData, 
            }         
    }
}

// JudgePartialOrder trait implementation
impl< Key, Val, KeyValPair >

    JudgePartialOrder 
        < KeyValPair >
    
    for OrderOperatorByKeyReverse 
        < Key, Val, KeyValPair >

    where   Key:            PartialOrd,
            KeyValPair:     KeyValGet < Key, Val, >
    
    {

    fn judge_partial_cmp
            ( & self, lhs: & KeyValPair, rhs: & KeyValPair ) 
            -> 
            Option<Ordering> {
        rhs.key().partial_cmp( &lhs.key() )
    }
}



//  Less-than-by-key-WITH-custom-comparator (auto implementor)
//  --------------------------------------


/// Determines the order of key-value pairs by applying a custom order operator to keys.
/// 
/// This struct is a wrapper around a customized order operator, which the user can provide.
/// The struct will indicate that `(a,x) < (b,y)` whenever the order operator says that `a < b`.
/// 
/// Any object can be treated as a key-value pair as long as it implements [`KeyValGet<Key, Val>`](`KeyValGet`).
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq)]
pub struct OrderOperatorByKeyCutsom < Key, Val, KeyValPair, KeyComparator > 
{
    key_comparator: KeyComparator,
    phantom_key:    PhantomData < Key >,
    phantom_val:    PhantomData < Val >,
    phantom_pair:   PhantomData < KeyValPair >,        
}

// Constructor
impl< Key, Val, KeyValPair, KeyComparator > OrderOperatorByKeyCutsom < Key, Val, KeyValPair, KeyComparator >  {
    pub fn new( key_comparator: KeyComparator ) -> OrderOperatorByKeyCutsom < Key, Val, KeyValPair, KeyComparator >  {
        OrderOperatorByKeyCutsom{ key_comparator, phantom_key: PhantomData, phantom_pair: PhantomData, phantom_val: PhantomData }
    }
}

// Copy
impl< Key, Val, KeyValPair, KeyComparator >

    Clone for 
    
    OrderOperatorByKeyCutsom 
        < Key, Val, KeyValPair, KeyComparator > 

    where 
        KeyComparator:  Clone,
{
    fn clone(&self) -> Self { OrderOperatorByKeyCutsom::new( self.key_comparator.clone() ) }
} 

// Copy
impl< Key, Val, KeyValPair, KeyComparator > 

    Copy for 
    
    OrderOperatorByKeyCutsom 
        < Key, Val, KeyValPair, KeyComparator > 

    where 
        KeyComparator:  Clone + Copy,        
{ }  

// JudgePartialOrder trait implementation
impl< Key, Val, KeyValPair, KeyComparator >

    JudgePartialOrder 
        < KeyValPair > for 
    
    OrderOperatorByKeyCutsom 
        < Key, Val, KeyValPair, KeyComparator > where

    KeyValPair:     KeyValGet < Key, Val, >,
    KeyComparator:  JudgePartialOrder <  Key >,
    
    {
    fn judge_partial_cmp
            ( & self, lhs: & KeyValPair, rhs: & KeyValPair ) 
            -> 
            Option<Ordering> {
        self.key_comparator.judge_partial_cmp( & lhs.key(), & rhs.key())
    }
}

// JudgePartialOrder trait implementation
impl< Key, Val, KeyValPair, KeyComparator >

    JudgeOrder 
        < KeyValPair > for 
    
    OrderOperatorByKeyCutsom 
        < Key, Val, KeyValPair, KeyComparator > where

    KeyValPair:     KeyValGet < Key, Val, >,
    KeyComparator:  JudgePartialOrder <  Key >,
    
    {}


//  Mutable closure 
//  ---------------------------------------------

use std::marker::PhantomData;

use crate::algebra::vectors::entries::KeyValGet;


/// Determines order based on a user-supplied funcion that returns `Option<Ordering>`
/// 
/// Concretely, this is a wrapper for objects that implement `FnMut(&T, &T)->Option<Ordering>`.
/// If `x` is such an object, then we say `a < b` if `x( &a, &b ) = Some(Ordering::Less)`.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq)]
pub struct OrderOperatorByFn
                < F: Fn(&T, &T)->Option<Ordering>, T > 
{ 
    order_operator_unwrapped:         F, 
    phantom_comparatee:               PhantomData< T > 
}


// Constructor
impl < F: Fn(&T, &T)->Option<Ordering>, T > OrderOperatorByFn< F, T > { 
    pub fn new( order_operator_unwrapped: F ) 
            -> 
            OrderOperatorByFn< F, T > 
    { OrderOperatorByFn{ order_operator_unwrapped, phantom_comparatee: PhantomData } } 
}

// Copy
impl <T, F: Fn(&T, &T)->Option<Ordering> > 

    Clone for 
    
    OrderOperatorByFn
        < F, T > 

    where 
        F:  Clone,
{
    fn clone(&self) -> Self { OrderOperatorByFn::new( self.order_operator_unwrapped.clone() ) }
} 

// Copy
impl <T, F: Fn(&T, &T)->Option<Ordering> > 

    Copy for 
    
    OrderOperatorByFn
        < F, T > 

    where 
        F:  Clone + Copy,      
{ }  


// partial order
impl <T, F: Fn(&T, &T)->Option<Ordering> > 

    JudgePartialOrder
        <T> for 
    OrderOperatorByFn
        < F, T > 
{
    // fn lt( & self, a: &T, b: &T) -> bool {
    //     (self.order_operator_unwrapped)(a, b)
    // }
    fn judge_partial_cmp
            ( & self, lhs: & T, rhs: & T ) 
            -> 
            Option<Ordering> {
        (self.order_operator_unwrapped)( lhs, rhs )
    }
}



//  Mutable closure -- generatoed from a "less than" comparison
//  ---------------------------------------------

/// Determines order based on a user-supplied funcion that returns `bool`
/// 
/// Concretely, this is a wrapper for objects that implement `FnMut(&T, &T)->bool`.
/// If `x` is such an object, then we say `a < b` if `x( &a, &b ) = true`.
/// 
/// # See also
/// 
/// See See the documentation for the [order](crate::utilities::order) module for an overview on order operators.
#[derive(Debug, Eq, PartialEq)]
pub struct OrderOperatorByLessThan
                < F: Fn(&T, &T)->bool, T > 
{ 
    order_operator_unwrapped:         F, 
    phantom_comparatee:                 PhantomData< T > 
}

impl < F: Fn(&T, &T)->bool, T > OrderOperatorByLessThan< F, T > { 
    pub fn new( order_operator_unwrapped: F ) 
            -> 
            OrderOperatorByLessThan< F, T > 
    { OrderOperatorByLessThan{ order_operator_unwrapped, phantom_comparatee: PhantomData } } 
}


// Clone
impl <T, F: Fn(&T, &T)->bool > 

    Clone for 
    
    OrderOperatorByLessThan
        < F, T > 

    where 
        F:  Clone,
{
    fn clone(&self) -> Self { OrderOperatorByLessThan::new( self.order_operator_unwrapped.clone() ) }
} 

// Copy
impl <T, F: Fn(&T, &T)->bool > 

    Copy for 
    
    OrderOperatorByLessThan
        < F, T > 

    where 
        F:  Clone + Copy,      
{ }  


// Judge partial order
impl <T, F: Fn(&T, &T)->bool > 

    JudgePartialOrder
        <T> for 

    OrderOperatorByLessThan
        < F, T > 
{
    fn judge_partial_cmp
            ( & self, lhs: & T, rhs: & T ) 
            -> 
            Option<Ordering> {
        if (self.order_operator_unwrapped)( lhs, rhs ) { return Some(Ordering::Less) }
        if (self.order_operator_unwrapped)( rhs, lhs ) { return Some(Ordering::Greater) }
        Some(Ordering::Equal)
    }

    fn judge_lt( &self, lhs: & T, rhs: & T ) -> bool { (self.order_operator_unwrapped)( lhs, rhs ) }
}

impl <T, F: Fn(&T, &T)->bool > 

    JudgeOrder
        <T> for 

    OrderOperatorByLessThan
        < F, T > 
{
    fn judge_cmp
            ( & self, lhs: & T, rhs: & T ) 
            -> 
            Ordering {
        if (self.order_operator_unwrapped)( lhs, rhs ) { return Ordering::Less }
        if (self.order_operator_unwrapped)( rhs, lhs ) { return Ordering::Greater }
        Ordering::Equal
    }
}





