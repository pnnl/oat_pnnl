//! General advice for computing PH with Umatch.
//! 
//! - You can compute PH by calculating a Umatch factorization of an (appropriately ordered) boundary matrix.
//!   - Recall that to creat a new Umatch factorization, you need a matrix oracle and an iterator.
//!   - The iterator has to **run over row indices in reverse order of birth**,
//!     otherwise the factorization may produce an invalid result.  This has nothing to do with factorization
//!     algorithm -- it's a structural property of the factorization itself.
//!     - Note, however, that this rule only applies to simplices of equal dimension.  Simplices of different
//!       dimensions can be placed in any order.
//!   - The computation can run much faster if you use the clearning optimization.  To make use of this
//!     optimization, you'll need to order your row index iterator first by dimension (ascending) then
//!     (for indices of equal dimension) by birth time (descending)
//!   - There is a [crate::matrices::operations::umatch::row_major::betti_numbers] function to compute betti numbers