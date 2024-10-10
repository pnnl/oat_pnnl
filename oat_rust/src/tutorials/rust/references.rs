//! References and lifetimes are integral to Rust.  This file contains a few hard-to-find facts to make your experience easier
//! 
//! # Lifetime annotations describe the referenced object, not the reference itself
//! 
//! Consider the following example.  The first function delcaration does not result in an error
//! because the reference `& self.worm` points at an object with lifetime `'a`.  The second function declaration 
//! results in an error because `&'a self.worm` *points* at something with lifetime `'a`, but does not have
//! lifetime `'a` itself.
//! 
//! ```
//! // an apple containing a worm
//! pub struct AppleVal{ 
//!     worm: usize
//! }
//! 
//! impl AppleVal {  
//! 
//!     fn worm_ref< 'a >( &'a self ) -> &'a usize { & self.worm }
//! 
//!     // fn worm_ref_ref< 'a >( &'a self ) -> &'a &'a usize { && self.worm }
//! }
//! ```
//! 
//! A similar phenomenon can be found below.  Here only the third function declaration,
//! `worm_ref_ref_ref`, results in an error.
//! 
//! ```
//! // apple containing reference
//! pub struct AppleRef<'a>{ 
//!     worm: &'a usize
//! }
//! 
//! impl <'a> AppleRef < 'a > {  
//! 
//!     fn worm_ref( &'a self) -> &'a usize { self.worm }
//! 
//!     fn worm_ref_ref( &'a self ) -> &'a &'a usize { & self.worm }
//! 
//!     // fn worm_ref_ref_ref( &'a self ) -> &'a &'a &'a usize { && self.worm }    
//! }
//! ```

