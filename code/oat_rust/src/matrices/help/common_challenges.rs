//! Common challenges, and how to address them
//!
//! # How to pass a matrix to a function without consuming it
//! 
//! Many functions take matrices as arguments.  However, once you pass a matrix to one of these
//! functions, it is often "consumed" (no longer available for use).  
//! To avoid this, it is sometimes possible to pass `&T` as an argument, rather than `T`.
//! This solution will not work, however, if `&T` does not implement the same matrix oracle traits that `T` does.
//! In this case you can use the [`OracleRef`](oat_rust::matrices::matrix_types::oracle_ref::OracleRef) object:
//! 
//! ```ignore
//! let oracle_ref = OracleRef::new( &T );
//! my_function( oracle_ref );
//! ```
//! 
//! The `OracleRef` struct is just a wrapper that contains a reference to a matrix oracle, i.e. an `&T`.  Unlike `&T`, however, an `OracleRef` is guaranteed to implement every matrix oracle trait that `T` implements.