//! Common error messages
//!
//! # Error message: Does not implement [`OracleRefInherit`](oat_rust::marices::matrix_oracle_traits::OracleRefInherit)
//! 
//! ```ignore
//! the trait bound `< YourType >: OracleRefInherit` is not satisfied the trait `OracleRefInherit` is not implemented for `YourType`rustcE0277
//! ```
//! 
//! This error message is somewhat misleading.  It typically occurs when there is a type, `YourType`, such that (i) `YourType` implements one of the oracle traits, but (ii) `& YourType` does not.
//! The error is really trying to express the idea that `& YourType` needs to implement an oracle trait, but it does not.
//! You can correct the error by implementing the necessary oracle trait; implementing the `OracleRefInherit` trait on `YourType` is only one of many
//! different ways to achieve this (implementing `OracleRefInherit` will automatically implement all oracle traits on `& YourType` which are currently impelemented for `YourType`).
//! 
//! 
