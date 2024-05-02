


//! First steps with matrix oracles.
//! 
//! This tutorial is under construction; for an introduction to matrix oracles, see the [matrices](crate::matrices) module.




















//      ==========================================================================================
//      THE FOLLOWING DOCUMETNATION IS DEPRECATED, BUT MAY BE USEFUL RAW MATERIAL FOR THE TUTORIAL
//      ==========================================================================================
// //! 
// //! Views and major dimensions are the two basic concepts for understanding matrix 
// //! oracles in oat_rust.  A good way to 
// //! 
// //! can provide two types of information:
// //! 
// //! * **Rows and columns** What are the entries in the ith row or jth column?  
// //! * **Major dimension** With many (but not all) sparse matrices, it's easier to access
// //! one dimension than another, e.g. easier to access rows than columns.  It's important
// //! to know which is easier when you're working with the matrix.  
// //! 
// //! When you ask a matrix oracle for information about a row or column vector, it returns
// //! something called a "view."  A view is an iterator that runs over the entries of that
// //! row/column.  We call it a view because it doesn't necessarily allow you to re-write 
// //! the entries that it exposes.
// //! 
// //! **Example** If a matrix is row-major, then 
// //! 
// //! * the ith major view is the ith row of the matrix, 
// //! * the jth minor view is the jth column of the matrix.
// //! 
// //! 
// //! 
// //! # How to write your own matrix oracle
// //! 
// //! Writing your own matrix oracle is easier than you might think.
// //! 
// //! ### Example of scary source code
// //! 
// //! The first time you see source code for a matrix oracle, it can look
// //! a bit daunting.  For example, the following is an excerpt from 
// //! code that defines a vector-of-vector matrix, and implements the 
// //! `OracleMajor` trait; it looks like a real mess.
// //! 
// //!
// //! ```ignore
// //! // Define the object
// //! 
// //! pub struct VecOfVec
// //! 
// //!     < 'a, IndexCoeffPair >
// //! 
// //!     where   IndexCoeffPair:    KeyValGet,
// //!             Self:           'a
// //! 
// //! {
// //!     pub major_dimension: MajorDimension, 
// //!     pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
// //!     pub phantom: PhantomData<&'a IndexCoeffPair >
// //! }
// //! 
// //! // Implement the trait
// //! 
// //! impl < 'a, IndexCoeffPair > 
// //!     
// //!     OracleMajor
// //!     <   
// //!         'a,
// //!         usize, 
// //!         < IndexCoeffPair as KeyValGet >::Key, 
// //!         < IndexCoeffPair as KeyValGet >::Val, 
// //!     > 
// //!     
// //!     for 
// //!     
// //!     VecOfVec < 'a, IndexCoeffPair > 
// //! 
// //!     where   IndexCoeffPair:    KeyValGet + Clone + 'a,
// //!             Self: 'a
// //! {
// //!     type PairMajor = IndexCoeffPair;
// //!     type ViewMajor = Cloned <std::slice::Iter <'a, IndexCoeffPair>>; 
// //!         
// //!     fn view_major <'b: 'a>( &'b self, index: usize ) -> ViewMajor {
// //!         return self.vec_of_vec[index].iter().cloned()
// //!     } 
// //! }
// //! ```
// //!
// //! 
// //! ### Easily modify the example to do what you want
// //!
// //! Suppose we want to write a matrix oracle that represents a scalar matrix. 
// //! The only data we need to define this matrix are: (i) the scalar value, `alpha`,
// //! and (ii) the major dimension of the matrix.  Let's define a struct to house
// //! this data:
// //! 
// //! ```
// //! // Import the object that formally encodes the two symbols for major dimension (row and col)
// //! use oat_rust::matrices::matrix_oracle_traits::MajorDimension; 
// //! 
// //! // Define the struct that represents a scalar matrix with a specified major dimension.
// //! pub struct ScalarMatrixDemo
// //! {
// //!     scalar: f64,                            // the scalar must be a float
// //!     major_dimension: MajorDimension,        // row-major or col-major
// //! }
// //! ```
// //! 
// //! 
// //! The `m`th row or column of a scalar matrix is equal to `alpha` times
// //! the `m`th standard unit vector.  This vector has at most one nonzero entry, 
// //! namely `(m, alpha)`.  So an oracle representing this matrix should return an 
// //! iterator that runs over `(m, alpha)` exactly once, for any `m` (for convenience
// //! we'll assume the matrix has infinite size, so any nonnegative integer `m` is 
// //! allowed). We can write a function that returns such an iterator, for any `m`:
// //! 
// //! ```
// //! # // Import the object that formally encodes the two symbols for major dimension (row and col)
// //! # use oat_rust::matrices::matrix_oracle_traits::MajorDimension; 
// //! # 
// //! # // A struct representing a scalar matrix.
// //! # pub struct ScalarMatrixDemo
// //! # {
// //! #     scalar: f64,                              // the scalar must be a float
// //! #     major_dimension: MajorDimension,          // row-major or col-major
// //! # }
// //! 
// //! /// Given a scalar matrix `M` and an index `i`, return a view of the `i`th 
// //! /// row or column vector of `M`.
// //! fn get_vector( matrix: &ScalarMatrixDemo, index: usize ) -> Vec< (usize, f64) > 
// //! {
// //!     let alpha = matrix.scalar.clone();          // make a copy of the scalar
// //!     return vec![ (index, alpha) ]
// //! }  
// //! ```
// //! 
// //! 
// //! Now we can modify the source code in the example above (used for vec-of-vec matrices)
// //! to implement the `OracleMajor` trait for our new struct.
// //! 
// //! ```
// //! // ORIGINAL CODE
// //! // impl < 'a, IndexCoeffPair > 
// //! //     
// //! //     OracleMajor
// //! //     <   
// //! //         'a,
// //! //         usize, 
// //! //         < IndexCoeffPair as KeyValGet >::Key, 
// //! //         < IndexCoeffPair as KeyValGet >::Val, 
// //! //     > 
// //! //     
// //! //     for 
// //! //     
// //! //     VecOfVec < 'a, IndexCoeffPair > 
// //! // 
// //! //     where   IndexCoeffPair:    KeyValGet + Clone + 'a,
// //! //             Self: 'a
// //! // {
// //! //     type PairMajor = IndexCoeffPair;
// //! //     type ViewMajor = Cloned <std::slice::Iter <'a, IndexCoeffPair>>; 
// //! //         
// //! //     fn view_major <'b: 'a>( &'b self, index: usize ) -> ViewMajor {
// //! //         return self.vec_of_vec[index].iter().cloned()
// //! //     } 
// //! // }
// //! 
// //! // MODIFIED CODE
// //! 
// //! # // Import the object that formally encodes the two symbols for major dimension (row and col)
// //! # use oat_rust::matrices::matrix_oracle_traits::*;
// //! # 
// //! # // A struct representing a scalar matrix.
// //! # pub struct ScalarMatrixDemo
// //! # {
// //! #     scalar: f64,                            // the scalar must be a float
// //! #     major_dimension: MajorDimension,        // row-major or col-major
// //! # }
// //! 
// //! impl < 'a >  // delete `IndexCoeffPair`, since our scalar matrix doesn't use this
// //!     
// //!     OracleMajor
// //!     <   
// //!         'a,     // we don't have to worry about this
// //!         usize,  // our major dimension is indexed by keys of type `usize`
// //!         usize,  // our minor dimension is indexed by keys of type `usize`
// //!         f64,    // our coefficients are f64
// //!     > 
// //!     
// //!     for 
// //!     
// //!     ScalarMatrixDemo
// //! 
// //!     where   Self: 'a    // we deleted `IndexCoeffPiar` so we remove the associated type constraints
// //! {
// //!     type PairMajor = (usize, f64);              // our vector entries are represented by objects of type `(usize, f64)`
// //!     type ViewMajor = Vec< (usize, f64) >;       // our vectors are represented by objects of type `Vec< (usize, f64) >`
// //!         
// //!     // To define the `major_view` function, we essentially copy/paste the body of
// //!     // our `get_vector` into the body of the original `major_view` function.  
// //!     // Note that we replace `matrix` with `self`.
// //!     fn view_major <'b: 'a>( &'b self, index: usize ) -> ViewMajor {
// //!         let alpha = self.scalar.clone();        // make a copy of the scalar
// //!         return vec![ (index, alpha) ]  
// //!     } 
// //! }
// //! ```
// //! 
// //! That's it!  Other traits can be implemented similarly.
// //! 
// //! **Note** Most functions that take matrix oracles as inputs do not require 
// //! their inputs to implement *all* of the oracle traits -- only a *subset*.