//! Ring operators for rings that already exist in Rust; packaged in convenient, zero-memory wrappers.
//! 
//! To define a semiring (or ring, or division ring) in oat_rust, you define an 
//! object `ring_operator` that implements the semiring (and possibly ring and division ring) trait.  You can then use `ring_operator` to 
//! perform basic operations on the elements of the ring (addition, multiplication, etc.)
//! 
//! Rust already has a number of rings "built in."  The current module provides a 
//! convenient way to generate a ring operator object for one of these built-in rings, using ['SemiringNative`], ['RingNative`], or ['DivisionRingNative`].
//! The these structs use zero memory!
//! `


use num::rational::Ratio;

use crate::rings::operator_traits::{Semiring, Ring, DivisionRing};
use std::marker::PhantomData;

//----------------------------------------------------------
//  SEMIRINGS NATIVE TO RUST
//----------------------------------------------------------


/// Zero-memory struct encoding structure of native Rust semirings.
///
/// # Examples
///
/// ```
/// use oat_rust::rings::operator_structs::ring_native::SemiringNative;
/// use oat_rust::rings::operator_traits::{Semiring};
///
/// let ring  =  < SemiringNative::<usize> >::new();
///
/// assert_eq!( 3, ring.add( 1, 2 ) ); 
/// assert_eq!( 2, ring.multiply( 1, 2 ) );
/// assert_eq!( 0, SemiringNative::<usize>::zero() );
/// assert_eq!( 1, SemiringNative::<usize>::one()  );
/// assert!( ! ring.is_0( 1 ) );
/// assert!(   ring.is_0( 0 ) );
/// assert!(   ring.is_1( 1 ) );
/// assert!( ! ring.is_1( 0 ) );
/// ```
#[derive(Debug, Clone, Copy)]
pub struct SemiringNative< Element >
    where 
        Element:    //num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    //core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    //core::ops::Div < Output = Element > +
                    //std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{ 
    // This phantom field uses zero memory; it is here only 
    // because rust otherwise complains that `Element` is
    // unused.  See the documentation on `PhantomData` for
    // more details.  **Note** that `*const` appears because
    // there is no relevant lifetime parameter for the 
    // struct.  Again, see the docs for `PhantomData`.
    phantom: PhantomData<*const Element> 
}

impl    < Element >
        SemiringNative 
        < Element > 
    where 
        Element:    //num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    //core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    //core::ops::Div < Output = Element > +
                    //std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    // Generate a `SemiringNative`.
    pub fn new( ) -> Self  
    {
        SemiringNative { phantom: PhantomData }
    }
}


impl    < Element > 
        Semiring < Element > for SemiringNative 
        < Element >  
    where 
        Element:    //num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    //core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    //core::ops::Div < Output = Element > +
                    //std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    /// Identity elements
    fn is_0( &self, x: Element ) -> bool { x.is_zero() }
    fn is_1( &self, x: Element ) -> bool { x.is_one() }
    fn zero() -> Element { Element::zero() }
    fn one()  -> Element { Element::one() }

    /// Add
    fn add( &self, x: Element, y: Element ) -> Element { x + y }

    /// Multiply
    fn multiply( &self, x: Element, y: Element ) -> Element { x * y }
}



//----------------------------------------------------------
//  RINGS NATIVE TO RUST
//----------------------------------------------------------

/// Zero-memory struct encoding structure of native Rust rings.
///
/// # Examples
///
/// ```
/// use oat_rust::rings::operator_structs::ring_native::RingNative;
/// use oat_rust::rings::operator_traits::{Semiring, Ring};
///
/// let ring = RingNative::<i64>::new();
/// let a : i64 = 1;
/// let b : i64 = 2;
///
/// assert_eq!( -1, ring.subtract( a, b ) );
/// assert_eq!( -1, ring.negate( a ) );
/// ```
#[derive(Debug, Clone, Copy)]
pub struct RingNative< Element >
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{ 
    // This phantom field uses zero memory; it is here only 
    // because rust otherwise complains that `Element` is
    // unused.  See the documentation on `PhantomData` for
    // more details.  **Note** that `*const` appears because
    // there is no relevant lifetime parameter for the 
    // struct.  Again, see the docs for `PhantomData`.
    phantom: PhantomData<*const Element> 
}

impl    < Element >
        RingNative 
        < Element > 
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    // Generate a `RingNative`.
    pub fn new( ) -> Self  
    {
        RingNative { phantom: PhantomData }
    }
}


impl    < Element > 
        Semiring < Element > for RingNative 
        < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    /// Identity elements
    fn is_0( &self, x: Element ) -> bool { x.is_zero() }
    fn is_1( &self, x: Element ) -> bool { x.is_one() }
    fn zero() -> Element { Element::zero() }
    fn one()  -> Element { Element::one() }

    /// Add
    fn add( &self, x: Element, y: Element ) -> Element { x + y }

    /// Multiply
    fn multiply( &self, x: Element, y: Element ) -> Element { x * y }
}

impl    < Element > 
        Ring < Element > for RingNative 
        < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    /// Subtract `x-y`.
    fn subtract( &self, x: Element, y: Element ) -> Element { x - y }

    /// Additive inverse `-x`. 
    fn negate( &self, x: Element ) -> Element { - x }
}


//----------------------------------------------------------
//  DIVISION RINGS NATIVE TO RUST
//----------------------------------------------------------

/// Zero-memory struct encoding structure of native Rust division rings.
///
/// # Examples
///
/// ```
/// use oat_rust::rings::operator_structs::ring_native::DivisionRingNative;
/// use oat_rust::rings::operator_traits::{Semiring, Ring, DivisionRing};
/// use num::rational::Ratio;
///
///
/// // The `< .. >` brackets around `DivisionRingNative::<Ratio<i64>>`
/// // are used to disambiguate which `new` function should be used
/// // (Rust throws an error if these aren't used)
/// let ring  =     < 
///                     DivisionRingNative::<Ratio<i64>> 
///                 >
///                 ::new();
/// let a = Ratio::new( 2, 3 );
/// let b = Ratio::new( 3, 1 );
/// let c = Ratio::new( 2, 9 );
/// let d = Ratio::new( 3, 2 );
///
/// assert_eq!( c, ring.divide( a, b ) );
/// assert_eq!( d, ring.invert( a ) );
/// ```
///
/// ```
/// use oat_rust::rings::operator_structs::ring_native::DivisionRingNative;
/// use oat_rust::rings::operator_traits::{Semiring, Ring, DivisionRing};
/// use num::rational::Ratio;
/// 
/// 
/// let ring : DivisionRingNative< f64 > =  DivisionRingNative::<f64>::new();
/// let a = 2.0 ;
/// let b = 4.0 ;
///
/// assert_eq!( 0.5, ring.divide( a, b ) );
/// assert_eq!( 0.5, ring.invert( a ) );
///
/// ```
#[derive(Debug, Clone, Copy)]
pub struct DivisionRingNative< Element >
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{ 
    // This phantom field uses zero memory; it is here only 
    // because rust otherwise complains that `Element` is
    // unused.  See the documentation on `PhantomData` for
    // more details.  **Note** that `*const` appears because
    // there is no relevant lifetime parameter for the 
    // struct.  Again, see the docs for `PhantomData`.
    phantom: PhantomData<*const Element> 
}
//{
//    zero: Element, // keep this on hand so it never has to be (de)allocated
//    one: Element, // keep this on hand so it never has to be (de)allocated
//}

impl    < Element >
        DivisionRingNative 
        < Element > 
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone

{
    // Generate a `DivisionRingNative`.
    pub fn new( ) -> Self  
    {
        DivisionRingNative { phantom: PhantomData }
    }
}


impl    < Element > 
        Semiring < Element > for DivisionRingNative 
        < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    /// Identity elements
    fn is_0( &self, x: Element ) -> bool { x.is_zero() }
    fn is_1( &self, x: Element ) -> bool { x.is_one() }
    fn zero() -> Element { Element::zero() }
    fn one()  -> Element { Element::one() }

    /// Add
    fn add( &self, x: Element, y: Element ) -> Element { x + y }

    /// Multiply
    fn multiply( &self, x: Element, y: Element ) -> Element { x * y }
}

impl    < Element > 
        Ring < Element > for DivisionRingNative 
        < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    /// Subtract y from x.
    fn subtract( &self, x: Element, y: Element ) -> Element { x - y }

    /// Additive inverse `-x`. 
    fn negate( &self, x: Element ) -> Element { - x }
}

impl    < Element > 
        DivisionRing < Element > for DivisionRingNative 
        < Element >  
    where 
        Element:    num::traits::Num + 
                    num::traits::Zero +
                    num::traits::One +
                    core::ops::Add < Output = Element >  +
                    core::ops::Sub < Output = Element > +
                    core::ops::Mul < Output = Element >  +
                    core::ops::Div < Output = Element > +
                    std::ops::Neg  < Output = Element > +
                    std::cmp::PartialEq +
                    std::clone::Clone
{
    /// `x/y` if `y` is nonzero.  
    fn divide( &self, x: Element, y: Element ) -> Element { x / y }

    /// `1/x` if `x` is nonzero.  
    fn invert( &self, x: Element ) -> Element { Element::one() / x }
}


//----------------------------------------------------------
//  CREATORS
//----------------------------------------------------------

pub fn field_f64() -> DivisionRingNative < f64 > 
    { DivisionRingNative::new() }

pub fn field_rational_i64() ->  DivisionRingNative < Ratio< i64 > > 
{ DivisionRingNative::new() }
