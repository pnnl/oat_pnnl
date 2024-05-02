//! Matrix traits, objects, and operations.
//! 
//! oat_rust provides a variety of tools to work with sparse matrices.
//! The most important are the [entry lookup traits](crate::matrices::matrix_oracle_traits),
//! which provide a standardized set of commands to look up the nonzero entries
//! in rows and columns.
//! 
// //! **In oat_rust, a "sparse matrix" is considered to be any struct that implements one or more of [entry lookup traits](crate::matrices::matrix_oracle_traits).**
//! 
// //! There are many ways to store a sparse matrix in memory (triplets, compressed sparse row, compressec sparse column, etc.).
// //! However, the specific storage strategy is usually 
// //! 
// //! 
// //!
// //! Two of the most common tasks, when working with a sparse matrix `M`, are to
// //! 
// //! * look up the nonzero entries in a row of `M`
// //! * look up the nonzero entries in a column of `M`
//! 
//! 
//! # Example
//! 
//! Let's construct a matrix, then look up the entries in its rows and columns.
//! 
//! The commands we use to perform the lookup are `view_major` and `view_minor`.
//! See [here](crate::matrices::matrix_oracle_traits) for an explanation of what these terms mean.
//!    
//! ```
//! // Import a package that defines VecOfVecSimple, a type of sparse matrix that stores nonzero entries in vector-of-vectors format
//! use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
//! 
//! // Import some traits that define look-up commands for rows and columns (namely, view_major and view_minor)
//! use oat_rust::matrices::matrix_oracle_traits::{OracleMajor, OracleMinor};
//!  
//! // Create a sparse 2x2 upper triangular matrix with 1's above the diagonal.  
//! // oat_rust has many different pre-defined types for sparse matrices.  One 
//! // example is `VecOfVecSimple`.  For technical reasons, the matrix oracle 
//! // traits aren't defined on `VecOfVecSimple` directly.  Instead, the oracle 
//! // traits are implemented on references to `VecOfVecSimple`.  References in 
//! // Rust are denoted by an `&` symbol.
//! let matrix_data =    VecOfVecSimple::new(
//!                         vec![   vec![   (0, 1.),   (1, 1.)  ], 
//!                                 vec![              (1, 1.)  ] 
//!                             ], 
//!                      );
//! let matrix  =   & matrix_data; 
//!  
//! // Look up the nonzero entries in row 0.
//! // To do so, we use the `view_major` method associated with the `OracleMajor` trait.
//! let view_of_row_0    =   matrix.view_major( 0 ); 
//! itertools::assert_equal( view_of_row_0, vec![(0, 1.), (1, 1.)] );
//! 
//! // Look up the nonzero entries in column 1.
//! // To do so, we use the `view_minor` method associated with the `OracleMinor` trait.
//! let view_of_col_1    =   matrix.view_minor( 1 ); 
//! itertools::assert_equal( view_of_col_1, vec![(0, 1.), (1, 1.)] );
//! ```
//! 
// //! **Vocabulary: major and minor dimensions** In many, though not all, sparse 
// //!     matrices, it's easier to look up rows than to look up columns, or vice versa.
// //!     We call the easy dimension the *major dimension*.  When we need to look up a row or column 
// //!     of a matrix, we do so by indexing into the appropriate major or minor dimension; see 
// //!     [here](crate::matrices::matrix_oracle_traits) for more details.*
//! 
//! # Help with matrices
//! 
//! See [help] for help creating new matrices, trouble shooting, and overcoming common design challenges.
//! 
//! 


pub mod operations;
pub mod matrix_oracle_traits; 
pub mod matrix_types; 
pub mod random_constructors;
pub mod debug;
pub mod display;
pub mod help;
pub mod developer;