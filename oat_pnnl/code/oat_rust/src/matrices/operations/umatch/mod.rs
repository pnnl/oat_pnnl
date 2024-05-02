//! Umatch factorization
//!
//! U-match is a form of matrix factorization described in the paper [U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co)homology](https://arxiv.org/abs/2108.08831), by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.
//! This factorization is expressed as a matrix equation
//! 
//! ```ignore
//! RM = DC
//! ```
//! where `D` is an `m` x `n` matrix, `R` and `C` are upper unitriangular matrices, and `M` is a generalized matching matrix.  We call
//! 
//! - `D` the mapping array
//! - `R` the codomain COMB
//! - `C` the domain COMB
//! - `M` the matching matrix
//! 
//! U-match factorization can be used to solve several core problems in linear algebra:
//! 
//! - solving `Dx = b` for `x`
//! - computing bases for the image and (co)kernel of `D`
//! 
//! for example, `{R[:,i] : M[i,j] != 0 }` is a basis for column space of `D`, and `{ C[:,j] : M[i,j] != 0 }` is a basis for the kernel of `D`.
//! 
//! oat_rust provides functionality to compute a U-match factorization, via the function  `UmatchRowMajor::new`.
//! This produces an `UmatchRowMajor` struct that provides the user with access to `M`, `R`, `R^{-1}`, `C`, and `C^{-1}`.
//! 
//! # Memory-efficiency
//! 
//! When oat_rust computes a U-match factorization, it stores only three items in memory: the mapping array `D`, the generalized matching array `M`, and the square submatrix of `R^{-1}` indexed by row-pivot indices.
//! It can be [shown](https://arxiv.org/abs/2108.08831) that this information suffices to reconstruct any row or column of `M`, `R`, `R^{-1}`, `C`, and `C^{-1}` efficiently (row access is generally much faster than column access).
//! This makes U-match factorization highly memory efficient, in practice.
//! 



// pub mod row_major_only;
pub mod row_major;
pub mod developer;
// pub mod diagnostics;