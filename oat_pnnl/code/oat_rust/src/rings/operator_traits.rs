//! Traits that define what an operator on a semiring, ring, or division ring can do.
//!
//!
//! TO-DO LIST FOR DEVELOPERS:
//! * MAKE ALL FUNCTIONS TAKE INPUTS BY VAL, NOT BY REFERENCE
//! (POINTERS CAN TAKE MORE MEMORY THAN INPUTS)
//! * MAKE AUTO-TRAIT IMPLEMENTATIONS TO HANDLE
//! INPUTS BY REFERENCE


//  ---------------------------------------------------------------------------
//  DESIGN NOTES
//  ---------------------------------------------------------------------------

//  * Advantage of this nested structure: makes it straightforward to define matrix multipication
//  over semirings.
//
//  * Reason for deprecating the function "field name" that tells you the underlying mathematical
//  field: 
//  in general, you always know what struct you're working with; so it suffices to describe the
//  mathematical object underlying the struct in the struct's documentation

use auto_impl::auto_impl;





//  ---------------------------------------------------------------------------
//  THE SEMIRING TRAIT
//  ---------------------------------------------------------------------------

/// Basic operations for semirings.
#[auto_impl(&)] // auto-implement this trait on references to objects that implement the trait
pub trait Semiring < Element > {


    // IDENTITY ELEMENTS

    /// Return the additive identity.
    fn is_0( &self, x : Element ) -> bool;

    /// Return the multiplicative identity.
    fn is_1( &self, x : Element ) -> bool;

    /// Return the additive identity.
    fn zero() -> Element;

    /// Return the multiplicative identity.
    fn one() -> Element;


    // OPERATIONS 
    
    // DESIGN NOTE: if we changed these functions to take
    // non-references as input, we would have to think about
    // wether to require ring elements to implement the
    // copy trait.

    /// Add
    fn add( &self, x : Element, y : Element ) -> Element;

    /// Multiply
    fn multiply( &self, x : Element, y: Element ) -> Element;

}


//  ---------------------------------------------------------------------------
//  THE RING TRAIT
//  ---------------------------------------------------------------------------


/// Basic operations for **unital** rings.
pub trait Ring <Element> : Semiring < Element > {

    /// Subtract y from x.
    fn subtract( &self, x : Element, y: Element ) -> Element;

    /// Reverse the sign of x.
    fn negate( &self, x : Element ) -> Element;

    /// Return `-1`
    fn minus_one( &self ) -> Element {
        self.negate( Self::one() )
    }    

    /// Return `(-1)^k`
    fn minus_one_to_power( &self, k: usize ) -> Element {
        if k.is_even() { Self::one() }
        else { self.minus_one() }
    }    

}


//----------------------------------------------------------
//  THE DIVISION RING TRAIT 
//----------------------------------------------------------

/// Basic operations for division rings.
pub trait DivisionRing <Element> : Ring < Element > {
    
    /// Divide 
    fn divide( &self, x : Element, y: Element ) -> Element;

    /// Invert 
    fn invert( &self, x : Element ) -> Element;

}







//----------------------------------------------------------
//  TRAITS FO CONVENIENCE
//----------------------------------------------------------

use num::integer::Integer;

pub trait MinusOne< RingElement > {
    fn minus_one( &self ) -> RingElement;
}

impl    < RingOperator, RingElement > 
        MinusOne
        < RingElement > 
        for 
        RingOperator
        where   RingOperator: Semiring< RingElement > + Ring< RingElement >
{
    fn minus_one( &self ) -> RingElement {
        self.negate( RingOperator::one() )
    }
}  

pub trait MinusOneToPower< RingElement > {
    fn minus_one_to_power( &self, k: usize ) -> RingElement;
}

impl    < RingOperator, RingElement > 
        MinusOneToPower
        < RingElement > 
        for 
        RingOperator
        where   RingOperator: Semiring< RingElement > + Ring< RingElement >
{
    fn minus_one_to_power( &self, k: usize ) -> RingElement {
        if k.is_even() { RingOperator::one() }
        else { self.minus_one() }
    }
}    