//! Traits and tools to work with functions (composition, evaluation, etc.)
//! 
//! The RUST language offers a range of traits to characterize functions,
//! e.g. `Fn`, `FnOnce`, and `FnMut`.  However, these traits don't provide all the
//! functionality one might wish; for example, one can't manually implement `Fn` on a struct,
//! in general.  This module provides traits to supplement some of this missing functionality,
//! as well as tools for basic operations, such as composition of functions.

pub mod evaluate;
pub mod compose;
pub mod misc_functions;