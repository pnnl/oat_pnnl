//! Vec-of-Vec sparse matrices.
//! 
//! A "vec-of-vec" is a vector of vectors representing the nonzero entries of a sparse matrix.
//! For example, if `A` is the following matrix
//! 
//! ```text
//! 5 6
//! 7 0
//! ```
//! 
//! Then the row-major vec-of-vec representation of `A` would be 
//! 
//! ```
//! vec![
//!     vec![ (0, 5), (1, 6) ],
//!     vec![ (0, 7) ],
//! ];
//! ```
//! 
//! and the column-major vec-of-vec representation of `A` would be
//! 
//! ```
//! vec![
//!     vec![ (0, 5), (1, 7) ],
//!     vec![ (0, 6) ],
//! ];
//! ```
//! 
//! This module contains several variants of the vec-of-vec data structure for use
//! with the [matrix oracle traits](crate::algebra::matrices::matrix_oracle_traits).
//! 
//! # Design elements
//! 
//! There are three separate vec-of-vec data structures, each motivated by a different use.
//! 
//! - [Sorted entries (defualt order)](sorted) 
//!   - Raison D'Etre  
//!     
//!     - This format is used extensively throughout the documentation and unit tests in this library
//!       because it is one of the most human readable sparse vector formats available.  However, this struct
//!       has a number of type parameters; more type parameters make code less readable and examples/tests harder to write.
//!                   The [VecOfVec](VecOfVec) struct has fewer type parameters; it is easer to read, construct,
//!                   and analyze.  (The price you pay for this simplicity is flexibility/generality, but in unit tests
//!         this matters less).
//! 
//! - [Sorted entries, custom order](sorted_custom)
//!   - Raison D'Etre
//!     - *Lifetimes are required by the matrix oracle traits.*  
//!       To the best of their knowledge (please let them know if information is found to the contrary), the developers believe that a `Vec<Vec<EntryType>>` object does not have an explicit
//!       lifetime parameter.  Wrapping a `Vec<Vec<EntryType>>` in a `VecOfVec<..>`  struct
//!       provides a lifetime parameter which can be used in implementing the [matrix oracle traits](crate::algebra::matrices::matrix_oracle_traits).
//!     - *Strict order is required for safe, efficient implementation of ascending and descending oracles.*  
//!                 To facilitate the implementation of the [ViewRowAscend](crate::algebra::matrices::query::ViewRowAscend) and
//!                  [ViewRowDescend](crate::algebra::matrices::query::ViewRowDescend)
//!                  traits in an efficient manner, entries in the inner vectors should appear in
//!                  in strictly sorted order according to index 
//!                  (either ascending or descending; arbitrarily, we chose ascending).
//!                  Thus, if we want a vec-of-vec data structure that can efficiently implement matrix oracle
//!                  traits, then we need some way to ensure that the inner vectors are sorted, for *safety*.
//!                  The `VecOfVec` accomplishes this by restricting user access to the `Vec<Vec<EntryType>>` which it stores internally.
//!                  One can only generate new instances of the `VecOfVec` struct by 
//!                  using the constructor `new` function, and this function panics if it receives a vec-of-vec
//!                  where one of the inner vectors is not sorted.
//!     - *Explicit order comparators make order of entries unambiguous, and afford greater flexibility.*  
//!        Since we allow matrix indices to have essentially any type, it is necessary to be explit about 
//!                  how order is defined.  This is the reason why `VecOfVec` objects store an internal 
//!                  `OrderOperator` object, which determines when one index is strictly less than another.
//!                  Typically, an order comapartor requires zero memory.  Moreover, it's necessary to provide
//!                  an order comaparator when constructig a `VecOfVec` (so that the `new` function can
//!                  check that entries are sorted), so requiring an order comparator as an argument in the 
//!                  `new` constructor adds no additional burden to the user.
//!
//!   - *Methods and challenges of editing data stored in a `VecOfVec`.*  
//!     - Due to the requirements of object
//!         safety, the user is highly limited in their ability to modify that data stored internally by a `VecOfVec`.
//!         The most general route is to retrieve the internally stored `Vec<Vec<EntryType>>`, modify it, then 
//!         wrap it in a new `VecOfVec` struct.  Note, however, that this will incur the cost of (i) re-sorting
//!         each internal vector, or (ii) verifying that each internal vector is sorted, if it is already sorted.
//!   - *Alternatives*
//!     - If the restrictions imposed by safety requirements on the `VecOfVec` struct are overly onerous, 
//!       consider using [VecOfVecUnsorted](VecOfVecUnsorted).  This does not implement the [ViewRowAscend](crate::algebra::matrices::query::ViewRowAscend) or
//!       [ViewRowDescend](crate::algebra::matrices::query::ViewRowDescend) methods, but it is much easier to modifiy.
//! 
//! - [Unsorted entries](unsorted)
//!   - Raison D'Etre
//!     - The requirement that internal vectors be strictly sorted may place an undue borden on the
//!                  user in some cases, for example
//!       - when one wishes to make frequent updates to the matrix, without re-sorting and re-checking each internal vector each time
//!       - when defining an order comparator is difficult or onerous
//!     Nevertheless, one may still wish to implement `ViewRow` on 
//!                  Unlike a struct of type `Vec< Vec< EntryType > >`, a struct of type [VecOfVecUnsorted](VecOfVecUnsorted) has a lifetime
//!                  parameter.  We use that parameter in the implementation of the 
//!                  [matrix oracle traits](crate::algebra::matrices::matrix_oracle_traits).
//!     - **Issue** Could the same be achieved with a reference, `&'a Vec<Vec<EntryType>>`?
//!                  
//! 
//! 
//! 

pub mod sorted;
pub mod sorted_ref;
pub mod sorted_custom;
pub mod unsorted;