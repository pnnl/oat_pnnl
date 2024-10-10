//! Sparse vectors represented by iterators.
//!
//! 
//! In OAT, any iterator that runs over [sparse vector entries](crate::algebra::vectors::entries) can be treated as a sparse vector.
//! 
//! This module contains a number of tools for working with sparse vectors, organized into two groups: [operations], and [linear combinations](linear_combinations).
//! 
// //! Recall that an *entry* is an object that stores the data of a `(key, value)` pair.
// //! Not every entry is a tuple; an entry can be 
// //! any struct that implements the [KeyValGet](crate::algebra::vectors::entries::KeyValGet) trait.
// //! 
// 
// //! 
// //!
// //! **Example**
// //! 
// //! ```
// //! let v = vec![ (1, 1.5), (2, 2.5) ]; // this is a standard Rust command; it's not special to OAT
// //! 
// //! for (index, coeff) in v {
// //!     println!("(index, coeff) = {:?}", (index, coeff));
// //! }
// //! 
// //! // Should display:
// //! // (index, coeff) = (1, 1.5)
// //! // (index, coeff) = (2, 2.5)
// //! ```
//! 
//! 
//! 
//! 
// I THINK THE FOLLOWING IS NOW DEPRECATED
// //! # Working with sparse vectors
// //!
// //! -  Simplify< HitMerge > is the preferred method of working with a linear combination of vectors where
// //!     - each vector returns entries in sorted order
// //!     - you want the sum of the entries to be returned in sorted order
// //!     - you want to create the sum now and modify it later, e.g. by adding new vectors
// //! - The `Pair< Index, Coeff >` items may outperform tuples of form `( Index, Coeff )` (at least
// //! for now) because tuples demand a certain memory structure that may impede performance, c.f. 
// //! [rewriting in memory](https://www.reddit.com/r/rust/comments/79ry4s/tuple_performance/).


// pub mod svi;
pub mod entries;
pub mod operations;
// pub mod linear_combinations;
// pub mod svi_discussion;








             


