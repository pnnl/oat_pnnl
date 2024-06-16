//! Conventions for comments, naming, format, and testing
//!
//! ## Comments
//! 
//! -  use !!! to tag items that developers should inspect, in the future
//!    - review all occurrences of `!!!` in the comments before making a commit
//! 
//! ## Naming conventions
//! 
//! - when naming objects, functions, etc., arrange descriptive terms in order of "frequency," low to high
//!   - example: 
//!     - yes: `vector_long, vector_short, vector_medium`
//!     - no: `long_vector`, `short_vector`, `medium_vector`
//! 
//! ## Code format 
//! 
//! - indentation
//!   - current convention is
//!     - if a close prenthesis / bracket occurs on a line different from its open bracket, indent the closed bracket once
//!     - in this case indent each line between the two brackets twice
//!     - this makes it visually clear where the parenthetical content begins/ends
//!     - example:
//!         ```
//!         let a   =   vec![
//!                                 vec![ 1, 2, 3, ]
//!                             ];
//!         let b   =   vec![
//!                                 vec![ 4, 5, 6, ]
//!                             ];
//!         ```
//! 
//! # Testing
//! 
//! - Ensuring panic
//!     - We use the following template for tests that ensure that a piece of code which is intended to panic
//!   does in fact do so.  The original code for this template came from .
//!
//!   ```
//!   // this test will pass
//!   
//!   fn auto_panic() { panic!("This function does nothing but panic.") }
//!   
//!   let result = std::panic::catch_unwind(|| { auto_panic() } );
//!   assert!(result.is_err());  
//!   
//!   // if we uncomment the following segmet of code, then the test will fail,
//!   // since the expression inside `{ .. }`, namely `false`, does not panic.
//!   // let result = std::panic::catch_unwind(|| { false } );
//!   // assert!(result.is_err()); 
//!   ```
//!     - For an alternative approach, see the module [`assert_panic`](https://docs.rs/assert-panic/latest/assert_panic/macro.assert_panic.html)
//! 
//! 
