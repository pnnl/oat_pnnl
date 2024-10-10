//! Rings, semirings, and division rings, including finite fields, rational numbers, etc.
//! 
//! 
//! 
//! Many programming lanuages have special object types used to work with rings.  For example, there 
//! might be a `float` type, an `integer` type, and a `bool` type.  
//! 
//! OAT uses two different types
//! to define each ring: one for the elements of the ring, 
//! and one for an object that performs the ring operations (addition, multiplicaton, etc.).
//! In order to be used as a ring operator, an object has to implement one or more of the ring operator traits:
//! [rings](`crate::algebra::rings::operator_traits::Ring`), [semirings](`crate::algebra::rings::operator_traits::Semiring`), and [division rings](``crate::algebra::rings::operator_traits::DivisionRing``).
//! 
//! # Terminology
//! 
//! We call [Ring](`crate::algebra::rings::operator_traits::Ring`), [Semiring](`crate::algebra::rings::operator_traits::Semiring`),
//! and [DivisionRing](``crate::algebra::rings::operator_traits::DivisionRing``) the **ring operator traits**.
//! We call any object that implements one of these traits a **ring operator**. 
//! 
//! # Example
//! 
//! Let's try an example with the two-element (Boolean) field.  
//! 
//! ```
//! // Import the packages that define the ring traits
//! use oat_rust::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing}; 
//! 
//! // Import a package that defines a ring operator for boolean field
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator; 
//! 
//! // Define some elements of the ring.
//! let x = true; // the type of elements is `bool`
//! let y = false;
//! 
//! // Define the ring operator.
//! let ring_operator = BooleanFieldOperator::new(); // the type of the ring operator is `BooleanFieldOperator`
//! 
//! // Use the ring operator to perform arithmetic operations:
//! assert_eq!( true,   ring_operator.add( x, y ) );
//! assert_eq!( true,   ring_operator.subtract( x, y ) );
//! assert_eq!( false,  ring_operator.multiply( x, y ) );
//! assert_eq!( false,  ring_operator.divide( y, x ) );
//! ```
//!    
//! # Predefined rings
//! 
//! OAT offers several pre-defined ring operators.  Here are their constructors
//! 
//! 
//! ```
//! // The ring of integers
//! let ring_operator   =   oat_rust::algebra::rings::operator_structs::ring_native::RingNative::< i64 >::new();
//! 
//! // The two element field (elements are true/false)
//! let ring_operator   =   oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator::new();
//! 
//! // The finite field of order 3
//! let ring_operator   =   oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator::new( 3 );
//!
//! // The field of rational numbers
//! use num::rational::Ratio; // the OAT ring operator uses the built-in Rust type for rational numbers
//! let ring_operator   =   oat_rust::algebra::rings::operator_structs::ring_native::DivisionRingNative::< Ratio<i64> >::new();
//! 
//! // The field of real numbers, represented by floating points
//! // (unsafe, because arithmetic operations have numerical error)
//! let ring_operator   =   oat_rust::algebra::rings::operator_structs::ring_native::DivisionRingNative::< f64 >::new();
//! ```
//! 
//! # Create new ring types
//! 
//! OAT makes it easy to define new ring types.  To do so, simply define a struct, e.g. `MyRingOperator`,
//! and implement one or more of the following traits: [rings](`crate::algebra::rings::operator_traits::Ring`), [semirings](`crate::algebra::rings::operator_traits::Semiring`), and [division rings](``crate::algebra::rings::operator_traits::DivisionRing``).
//! 
//! 
//! 
// //! When you define a new ring, you'll actually use two types: one for the elements of the ring, 
// //! and another for 
// //! the ring *operations* (addition, multiplication, etc.).   You can use any object as a ring operator, as long as it implements the proper traits.
// //! There are separate traits for [rings](`crate::algebra::rings::operator_traits::Ring`), [semirings](`crate::algebra::rings::operator_traits::Semiring`), and [division rings](``crate::algebra::rings::operator_traits::DivisionRing``).
// //! 
// //! # Note for developers
// //!
// //! OAT differs from many other platforms in its approach to rings.  
// //! This choice stems from issues we have encountered in the past, concerning the relationship between strongly typed languages and infinite families of rings (like the family of all prime-order fields).  The
// //! main advantage of our approach is that it allows one to effectively work with infinitely many
// //! rings without defining infinitely many types.

pub mod operator_traits;
pub mod operator_structs;
