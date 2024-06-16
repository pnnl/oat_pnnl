//! # Technical insights for working with clones and copies
//! 
//! ## derive
//! 
//! - `#[derive(Clone)]` can be used to automatically implement `Clone` on a struct.  The macro will check that all the types associated with the struct implement Clone (possibly even types that aren't associated with any field); then it will implement a field-wise clone, c.f. [this website](https://stegosaurusdormant.com/understanding-derive-clone/)
//! 
//! - the condition that every generic parameter associated with a struct must implement `Clone` can be overly restrictive.  If `#[derive(Clone)]` fails for this reason, you can always manually implement `Clone`
//! 
//! - below is an example from [this website](https://hashrust.com/blog/moves-copies-and-clones-in-rust/)
//! 
//!   ```ignore
//!   //error:the trait `Copy` may not be implemented for this type
//!   //because its nums field does not implement `Copy`
//!   #[derive(Copy, Clone)]
//!   struct Numbers {
//!       nums: Vec<i32>
//!   }
//!   ```
//! 
//! ## "manual" implementation
//! 
//! from [this website](https://hashrust.com/blog/moves-copies-and-clones-in-rust/)
//! 
//!   ```
//!   struct Point {
//!       x: i32,
//!       y: i32,
//!   }
//! 
//!   //no method in Copy because it is a marker trait
//!   impl Copy for Point {}
//! 
//!   impl Clone for Point {
//!       fn clone(&self) -> Point {
//!           *self
//!       }
//!   }
//!   ```
//! 
//! ## general principles
//! 
//! ### `Copy`
//! 
//! - corresponds to a direct "bitwise" copy of an object
//! - you can implement manually as shown in the example above, but **there are no methods** in the `Copy` trait; you have to let rust figure out everything; also your struct has to meet certain requirements
//! - **references**
//!   - quote from ryan:
//! 
//!   > and I think &T is Copy if T is Copy, so if it was in a struct, the struct should also be copy, but in this case it would just copy the reference itself, not the underlying data
//! 
//!   > which is why &mut T is not Copy regardless of whether T is Copy or not.... because you cannot have multiple mutable references to the same data
//! 
//! 
//! 
//! ### `Clone`
//!   
//! - not necessarily a "deep copy" in the general sense of the term
//! - in words of Ryan "sometimes it is a deep copy, sometimes it isn't"