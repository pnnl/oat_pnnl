//! Solve `Ax = b` for `x`, where `A` is a matrix and `b` is a vector.
//! 
//! We recommend
//! 
//! - [triangle] for triangular matrices
//! - [echelon] for matrices that are (partially) in echelon form
//! - [umatch](crate::algebra::matrices::operations::umatch) for all other matrices ([`Ax = b`](crate::algebra::matrices::operations::umatch::row_major::Umatch::solve_dx_equals_b) and [`xA = b`](crate::algebra::matrices::operations::umatch::row_major::Umatch::solve_xd_equals_b))
//! 
//! It's also possible to solve `Ax = b` by [inverting](crate::algebra::matrices::operations::invert) a triangular matrix `A`, but overwhelming empirical evidence suggests that this is much less efficient than [solving](triangle) directly.


pub mod triangle;
pub mod echelon;