//! Sparse matrices represented by oracles that look up rows, columns, and entries.
//! 
//! - [About](#about)
//! - [Tools and help](#tools-and-help)
//! - [Build your own matrix](#build-your-own-matrix)
//! - [Development](#development)
//! 
//! 
//! 
//! # About
//!
//! OAT provides a variety of tools to work with sparse matrices.  Key features include
//! - flexible indexing: matrices can be indexed by arbitrary keys; for example, the boundary matrix of a simplicial complex can be indexed by simplices; this feature is powerful in a variety of homology computations, where explit enumeration of row and column indices is expensive
//! - lazy look up: rows and columns of matrices can be constructed in a lazy fashion, consistent with current state-of-the-art practices in large scale PH computation
//! - extensive unit testing
//! 
//! # Tools and help
//! 
//! To explore the complete set of matrix tools, check out the [modules](#modules) section at the bottom of this page.
//! 
//! For help, including creating new matrices, trouble shooting, and overcoming common design challenges check out [help](crate::algebra::matrices::help).
//! 
//! # Build your own matrix
//!
//! Users can work with any object *as if it were a matrix*, provided the object implements the right [query traits](crate::algebra::matrices::query).  This means, for example, that you can write your own struct to represent the boundary matrix of a simplicial complex.
//!
//! Let's illustrate with an example.  We'll define an object that represents an identity matrix, but throw in a twist:
//! - rows are indexed integers
//! - columns are indexed by strings
//! - for example, row 0 has a single nonzero entry, in column "0".  Note the difference: the row index is a number, but the column index is a string.
//! 
//! This object uses no memory at all!
//! 
//! **NB** This example uses the idea of major and minor views.  Check out the [query traits](crate::algebra::matrices::query) page for an explanation of what these words mean.
//! 
//! ```
//! use oat_rust::algebra::matrices::query::{IndicesAndCoefficients};
//! use oat_rust::algebra::matrices::query::{MatrixEntry, ViewRow, ViewCol};
//! 
//! //  ---------------------------------------------
//! //  Define the new struct, and create an instance
//! //  ---------------------------------------------
//! 
//! pub struct MyMatrix{} // this struct contains no data, and uses zero memory
//! let matrix = MyMatrix{}; // here's an instance
//! 
//! //  ----------------------------------
//! //  Implement `IndicesAndCoefficients`
//! //  ----------------------------------
//! 
//! // This boils down to declaring a specific type for RowIndex, ColIndex, and Coefficient
//! impl IndicesAndCoefficients for MyMatrix {
//!     type EntryMajor     =   ( String, f64 );    // entries of major view are tuples
//!     type EntryMinor     =   ( u64, f64 );       // entries of minor views are tuples
//!     type RowIndex       =   u64;                // major keys are integers
//!     type ColIndex       =   String;             // minor keys are strings
//!     type Coefficient    =   f64;                // coefficients are 64 bit floats
//! }
//! 
//! //  -------------------------------------------------------
//! //  Implement `MatrixEntry`, to generate individual entries
//! //  -------------------------------------------------------
//! 
//! impl MatrixEntry for MyMatrix {
//!     fn entry_major_at_minor( &self, keymaj: Self::RowIndex, keymin: Self::ColIndex, ) -> Option< Self::Coefficient > {
//!         if keymin == keymaj.to_string() {
//!             // if the row and column indices agree, then return 1.0
//!             return Some(1.0)
//!         } else {
//!             // otherwise return None to indicate that the entry is structurally zero
//!             return None
//!         }
//!     }
//! }
//! 
//! //  ----
//! //  Test
//! //  ----
//! 
//! assert_eq!( matrix.entry_major_at_minor( 0, 0.to_string() ),  Some(1.0) );
//! assert_eq!( matrix.entry_major_at_minor( 0, 1.to_string() ),  None      );
//! 
//! //  ----------------------------------------------------
//! //  Implement `ViewRow`, to generate major views
//! //  ----------------------------------------------------
//! 
//! impl ViewRow for MyMatrix {
//! 
//!     // declare that major views are lists of string-float tuples
//!     type ViewMajor            =   Vec< ( String, f64 ) >;
//!     
//!     // declare that calling `.into_iter()` on `Vec< ( String, f64 ) >` produces an object of the following type
//!     type ViewMajorIntoIter    =   std::vec::IntoIter<( String, f64 )>;
//! 
//!     // define the function that eats a row index and returns the entries of the row
//!     fn view_major( &self, index: Self::RowIndex ) -> Self::ViewMajor {
//!         // convert the row index to a column index
//!         let col_index = index.to_string();
//!         // return the nonzero entry in the specified row
//!         return vec![ ( col_index, 1.0 ) ]
//!     }
//! }
//! 
//! //  ----
//! //  Test
//! //  ----
//! 
//! assert_eq!( matrix.view_major( 0 ), vec![ ( 0.to_string(), 1.0 ) ] );
//! assert_eq!( matrix.view_major( 1 ), vec![ ( 1.to_string(), 1.0 ) ] );
//! 
//! //  --------------------------------------------------------------
//! //  Implement `ViewCol`, to generate minor views
//! //  --------------------------------------------------------------
//! 
//! impl ViewCol for MyMatrix {
//! 
//!     // declare that major views are lists of string-float tuples
//!     type ViewMinor            =   Vec< ( u64, f64 ) >;
//!
//!     // declare that calling `.into_iter()` on `Vec< ( u64, f64 ) >` should produce an object of the following type
//!     type ViewMinorIntoIter    =   std::vec::IntoIter<( u64, f64 )>;
//! 
//!     // define the function that eats a column index and returns the entries of the column
//!     fn view_minor( &self, index: Self::ColIndex ) -> Self::ViewMinor {
//!         // convert the column index to a row index
//!         let row_index = index.parse::<u64>().unwrap();
//!         // return the nonzero entry in the specified column
//!         return vec![ ( row_index, 1.0 ) ]
//!     }
//! }
//! 
//! //  ----
//! //  Test
//! //  ----
//! 
//! assert_eq!( matrix.view_minor( 0.to_string() ), vec![ ( 0, 1.0 ) ] );
//! assert_eq!( matrix.view_minor( 1.to_string() ), vec![ ( 1, 1.0 ) ] );
//! ```
//!
//! 
//! 



//!
//!
// // ! # Query
// // !
// // ! The most important are the [entry lookup traits](crate::algebra::matrices::matrix_oracle_traits),
// // ! which provide a standardized set of commands to look up the nonzero entries
// // ! in rows and columns.
// // !
// // //! **In OAT, a "sparse matrix" is considered to be any struct that implements one or more of [entry lookup traits](crate::algebra::matrices::matrix_oracle_traits).**
// // !
// // //! There are many ways to store a sparse matrix in memory (triplets, compressed sparse row, compressec sparse column, etc.).
// // //! However, the specific storage strategy is usually
// // //!
// // //!
// // //!
// // //! Two of the most common tasks, when working with a sparse matrix `M`, are to
// // //!
// // //! * look up the nonzero entries in a row of `M`
// // //! * look up the nonzero entries in a column of `M`
// // !
// // !
// // !
// // ! Let's construct a matrix, then look up the entries in its rows and columns.
// // !
// // ! The commands we use to perform the lookup are `view_major` and `view_minor`.
// // ! See [here](crate::algebra::matrices::matrix_oracle_traits) for an explanation of what these terms mean.
// // !
// // ! ```
// // ! // Import a package that defines VecOfVec, a type of sparse matrix that stores nonzero entries in vector-of-vectors format
// // ! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
// // !
// // ! // Import some traits that define look-up commands for rows and columns (namely, view_major and view_minor)
// // ! use oat_rust::algebra::matrices::query::{ViewRow, ViewCol};
// // !
// // ! // Create a sparse 2x2 upper triangular matrix with 1's above the diagonal.
// // ! // OAT has many different pre-defined types for sparse matrices.  One
// // ! // example is `VecOfVec`.  For technical reasons, the matrix oracle
// // ! // traits aren't defined on `VecOfVec` directly.  Instead, the oracle
// // ! // traits are implemented on references to `VecOfVec`.  References in
// // ! // Rust are denoted by an `&` symbol.
// // ! let matrix_data =    VecOfVec::new(
// // !                         vec![   vec![   (0, 1.),   (1, 1.)  ],
// // !                                 vec![              (1, 1.)  ]
// // !                             ],
// // !                      );
// // ! let matrix  =   & matrix_data;
// // !
// // ! // Look up the nonzero entries in row 0.
// // ! // To do so, we use the `view_major` method associated with the `ViewRow`.
// // ! let view    =   matrix.view_major( 0 );
// // ! itertools::assert_equal( view, vec![(0, 1.), (1, 1.)] );
// // !
// // ! // Look up the nonzero entries in column 1.
// // ! // To do so, we use the `view_minor` method associated with the `ViewCol`.
// // ! let view    =   matrix.view_minor( 1 );
// // ! itertools::assert_equal( view, vec![(0, 1.), (1, 1.)] );
// // ! ```
// // !
// // //! **Vocabulary: major and minor dimensions** In many, though not all, sparse
// // //!     matrices, it's easier to look up rows than to look up columns, or vice versa.
// // //!     We call the easy dimension the *major dimension*.  When we need to look up a row or column
// // //!     of a matrix, we do so by indexing into the appropriate major or minor dimension; see
// // //!     [here](crate::algebra::matrices::matrix_oracle_traits) for more details.*
//!
//!
//!
//! # Development
//!
//! Matrices are deeply integrated with many aspects of this library, so every aspect of the matrix API has design implications.  It's still evolving!  See the design documents (developer/design/matrix.rs) in the OAT git repository for notes and discussion.  If you have thoughts, we'd love to hear from you!


pub mod operations;
pub mod query;
pub mod types;
pub mod debug;
pub mod display;
pub mod help;








#[cfg(test)]
mod tests {



    // Note this useful idiom: importing names from outer (for mod tests) scope.

    use crate::algebra::matrices::query::MatrixEntry;

    #[test]
    fn doc_test_build_your_own_matrix() {

        use crate::algebra::matrices::query::IndicesAndCoefficients;
        use crate::algebra::matrices::query::{ViewRow, ViewCol};

        //  ---------------------------------------------
        //  Define the new struct, and create an instance
        //  ---------------------------------------------

        pub struct MyMatrix{} // this struct contains no data, and uses zero memory
        let matrix = MyMatrix{}; // here's an instance

        //  ----------------------------------
        //  Implement `IndicesAndCoefficients`
        //  ----------------------------------

        // This boils down to declaring a specific type for RowIndex, ColIndex, and Coefficient
        impl IndicesAndCoefficients for MyMatrix {
            type EntryMajor =   ( String, f64 );  
            type EntryMinor       =   ( u64,    f64 );                      
            type RowIndex = u64;      // major keys are integers
            type ColIndex = String;   // minor keys are strings
            type Coefficient = f64;      // coefficients are 64 bit floats
        }

        //  -------------------------------------------------------
        //  Implement `MatrixEntry`, to generate individual entries
        //  -------------------------------------------------------

        impl MatrixEntry for MyMatrix {
            fn entry_major_at_minor( &self, keymaj: Self::RowIndex, keymin: Self::ColIndex, ) -> Option< Self::Coefficient > {
                if keymin == keymaj.to_string() {
                    // if the row and column indices agree, then return 1.0
                    Some(1.0)
                } else {
                    // otherwise return None to indicate that the entry is structurally zero
                    None
                }
            }
        }

        //  ----
        //  Test
        //  ----

        assert_eq!( matrix.entry_major_at_minor( 0, 0.to_string() ),  Some(1.0) );
        assert_eq!( matrix.entry_major_at_minor( 0, 1.to_string() ),  None      );

        //  ----------------------------------------------------
        //  Implement `ViewRow`, to generate major views
        //  ----------------------------------------------------

        impl ViewRow for MyMatrix {

            // major views are lists of string-float tuples
            type ViewMajor            =   Vec< ( String, f64 ) >;
            // calling `.into_iter()` on `Vec< ( String, f64 ) >` produces an object of the following type
            type ViewMajorIntoIter    =   std::vec::IntoIter<( String, f64 )>;

            fn view_major( &self, index: Self::RowIndex ) -> Self::ViewMajor {
                // convert the row index to a column index
                let col_index = index.to_string();
                // return the nonzero entry in the specified row
                vec![ ( col_index, 1.0 ) ]
            }
        }

        //  ----
        //  Test
        //  ----

        assert_eq!( matrix.view_major( 0 ), vec![ ( 0.to_string(), 1.0 ) ] );
        assert_eq!( matrix.view_major( 1 ), vec![ ( 1.to_string(), 1.0 ) ] );

        //  --------------------------------------------------------------
        //  Implement `ViewCol`, to generate minor views
        //  --------------------------------------------------------------

        impl ViewCol for MyMatrix {

            // major views are lists of string-float tuples
            type ViewMinor            =   Vec< ( u64, f64 ) >;
            // calling `.into_iter()` on `Vec< ( u64, f64 ) >` produces an object of the following type
            type ViewMinorIntoIter    =   std::vec::IntoIter<( u64, f64 )>;

            fn view_minor( &self, index: Self::ColIndex ) -> Self::ViewMinor {
                // convert the column index to a row index
                let row_index = index.parse::<u64>().unwrap();
                // return the nonzero entry in the specified column
                vec![ ( row_index, 1.0 ) ]
            }
        }

        //  ----
        //  Test
        //  ----

        assert_eq!( matrix.view_minor( 0.to_string() ), vec![ ( 0, 1.0 ) ] );
        assert_eq!( matrix.view_minor( 1.to_string() ), vec![ ( 1, 1.0 ) ] );
    }
}
