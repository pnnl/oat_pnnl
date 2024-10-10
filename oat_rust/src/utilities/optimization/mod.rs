//! Mathematical optimization
//! 
//! 
//! This module provides tools for mathematical optimization, which can be used to [tighten cycle representatives](https://www.frontiersin.org/articles/10.3389/frai.2021.681117/full).
//! 
//! 

pub mod minimize_l1;

// This conditionally includes a module which implements WEBP support.
#[cfg(feature = "gurobi")]
pub mod gurobi;