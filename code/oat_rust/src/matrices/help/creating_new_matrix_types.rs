//! Want to create your own oracle?  Perfect!  That is what oat_rust's made for!
//! 
//! 
//! 
//! Here are some resources to help you get started.  Please feel free to reach out to the developers if you find you are getting stuck.
//! 
//! - A matrix is any struct that implements one of the [entry lookup traits](crate::matrices::matrix_oracle_traits).  In order to implement these traits, you will also ned to implement [IndicesAndCoefficients](crate::matrices::matrix_oracle_traits::IndicesAndCoefficients).
//! This means that you can create new type of matrix by defining a new struct, then implementing a lookup trait.
//! - The folder `src/matrices/matrix_types` contains several examples of objects that implement the oracle traits.
//!   - The [`scalar_matrices`](crate::matrices::matrix_types::scalar_matrices) module may be a good place to start, 
//! because the [`ScalarMatrix`](crate::matrices::matrix_types::scalar_matrices::ScalarMatrix) struct is one of the simplest oracles.
//! - The following examples may give ideas / inspiration, and highlight some common issues in implementation.
//! 
//! ## Examples
//! 
//! From time to time, you may want to define a matrix oracle whose `views` include a 
//! [generic lifetime parameter](https://doc.rust-lang.org/book/ch10-03-lifetime-syntax.html).
//! 
//! There are several ways to approach this type of situation, in practice.  The following code
//! gives several examples.  We also include some "nonexamples" -- approaches which we were
//! not able to turn into feasible solutions.  These aren't necessarily dead ends, but we hope
//! they will help illustrate some common stumbling blocks.
//! 
//! 
//! ```
//! // import crates
//! use oat_rust::matrices::matrix_oracle_traits::{ OracleMajor, IndicesAndCoefficients };          
//! 
//! 
//! //  -----------------------------------------------------------------------------------------
//! //  implement OracleMajor on a WRAPPER FOR &'a Vec< Vec< usize > > 
//! //  -----------------------------------------------------------------------------------------
//! 
//! // NB:  (1) the salient feature of this segment are that 
//! //          (a) it uses an immutable reference to a vector-of-vectors.  This is because the
//! //              `lifetime` of the reference can help us to overcome some issues concerning
//! //              lifetimes and major views.
//! //          (b) it introduces and uses a *wrapper* struct that contains a `&'a Vec< Vec< usize > >` 
//! //              This is beecause Rust's "orphan rule" prevents the implementation of a
//! //              foreign trait on a foreign type (see the following segment for details)
//! 
//! // This struct is just a wrapper around `&'a Vec< Vec< (usize,usize)`
//! struct VecOfVecReferenceWrapper< 'a > { vec_of_vec: &'a Vec< Vec< (usize,usize) > > }
//! 
//! // Here we implement a trait with no methods; the trait doesn't do anything except specify three types: KeyMin, KeyMan, and SnzVal 
//! impl < 'a > IndicesAndCoefficients for VecOfVecReferenceWrapper< 'a > { 
//!     type KeyMin = usize;  type KeyMaj = usize;  type SnzVal = usize;
//! }
//! 
//! impl < 'a >
//! 
//!         OracleMajor for 
//! 
//!         VecOfVecReferenceWrapper< 'a >
//! 
//! {
//!     type ViewMajor          = &'a [(usize,usize)];
//!     type ViewMajorIntoIter  = std::slice::Iter< 'a, (usize,usize) >;
//!     type ViewMajorEntry     = &'a (usize,usize);
//! 
//!     fn view_major( & self, index: Self::KeyMaj ) -> Self::ViewMajor {
//!         return self.vec_of_vec[index].as_slice()
//!     } 
//! }
//! 
//! // Get a major view
//! let matrix = VecOfVecReferenceWrapper{ vec_of_vec: &vec![ vec![ (1,1), (2,2) ]  ] };
//! let row = matrix.view_major( 0 );
//! itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );
//! 
//! // Get a major view from a *nested* immutable reference.
//! let _row = (&(& matrix)).view_major( 0 );      
//! 
//! 
//  THE FOLLOWING SEGMENTS NO LONGER WORK (POSSIBLY DUE TO ORPHAN RULES?)
//
//
// //! //  -----------------------------------------------------------------------------------------
// //! //  implement OracleMajor on &'a Vec< Vec< usize > > 
// //! //  -----------------------------------------------------------------------------------------
// //! 
// //! 
// //! // NB:  SO FAR AS WE ARE CURRENTLY AWARE it is not possible to implement OracleMajor on 
// //! //      &'a Vec< Vec< usize > > **except when the code resides in the same file where this
// //! //      trait is defined.**  The reason for this is Rust's orphan rule 
// //! //      (c.f. https://doc.rust-lang.org/error-index.html#E0117), which prevents the
// //! //      implementation of a foreign trait on a foreign type (where "foreign" means a trait or
// //! //      type that is defined outside the file where you want to implement the trait).
// //! //      Therfore the following code compiles in the unit tests for the file
// //! //      src/oat_rust/matrices/matrix_oracle_traits.rs, but it does not compile in doc string tests or
// //! //      unit tests defined in other files.
// //! 
// //! impl < 'a > IndicesAndCoefficients for &'a Vec < Vec < (usize,usize) > > { 
// //!     type KeyMin = usize;  type KeyMaj = usize;  type SnzVal = usize;
// //! }
// //! 
// //! impl < 'a >
// //! 
// //!         OracleMajor for 
// //! 
// //!         &'a Vec < Vec < (usize,usize) > > 
// //! 
// //! {
// //!     type ViewMajor          =   &'a [(usize,usize)];
// //!     type ViewMajorIntoIter  =   std::slice::Iter< 'a, (usize,usize) >;
// //!     type ViewMajorEntry     =   &'a (usize,usize);
// //! 
// //!     fn view_major( & self, index: Self::KeyMaj ) -> Self::ViewMajor {
// //!         return self[index].as_slice()
// //!     } 
// //! }
// //! 
// //! // Get a major view from an immutable reference.     
// //! let mut matrix = vec![ vec![ (1,1), (2,2) ]  ];
// //! let row = (& matrix).view_major( 0 );
// //! itertools::assert_equal( row, vec![ (1,1), (2,2) ].iter() );
// //! 
// //! // Get a major view from a *nested* immutable reference.
// //! let _row = (&(& matrix)).view_major( 0 );      
// //! 
// //! // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
// //! // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
// //! let _row_a = (& matrix).view_major( 0 );
// //! matrix.push( vec![ (4,4), (5,5)] );
// //! let row_b = (& matrix).view_major( 1 );   
// //! 
// //! itertools::assert_equal(row_b, vec![ (4,4), (5,5)].iter() );
// //! // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );        
// //! 
// //! 
// //! //  -----------------------------------------------------------------------------------------
// //! //  implement OracleMajor on &'a mut Vec< Vec< usize > > 
// //! //  -----------------------------------------------------------------------------------------
// //! 
// //! // NB:  SO FAR AS WE ARE CURRENTLY AWARE it is not possible to implement OracleMajor on &'a mut Vec< Vec< usize > > 
// //! //      in a manner similar to the example above.  Indeed, the following code is a facsimile of the 
// //! //      section above, with `&'a mut` replaceing `&'a`.  It does not seem to compile, due to 
// //! //      a lifetime conflict.
// //! 
// //! // impl < 'a, 'b, (usize,usize) >
// //! 
// //! //         OracleMajor
// //! //         <  usize,  &'a [(usize,usize)]  > for 
// //! 
// //! //         &'a mut Vec < Vec < (usize,usize) > > 
// //! 
// //! //         where   (usize,usize): Clone,
// //! //                 'b: 'a,
// //! // {
// //! //     fn view_major( & self, index: usize ) -> &'a [(usize,usize)] {
// //! //         return self[index].as_slice()
// //! //     } 
// //! // }
// //! 
// //! // // Get a major view from an immutable reference.     
// //! // let mut matrix = vec![ vec![ 1, 2, 3 ]  ];
// //! // let row = (& matrix).view_major( 0 );
// //! // itertools::assert_equal( row, vec![ vec![ (1,1), (2,2) ]  ].iter() );
// //! 
// //! // // Get a major view from a *nested* immutable reference.
// //! // let _row = (&(& matrix)).view_major( 0 );      
// //! 
// //! // // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
// //! // // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
// //! // let _row_a = (& matrix).view_major( 0 );
// //! // matrix.push( vec![ (4,4), (5,5) ] );
// //! // let row_b = (& matrix).view_major( 0 );   
// //! 
// //! // itertools::assert_equal(row_b.iter().cloned() , vec![ 1,2,3 ] );
// //! // // itertools::assert_equal( _row_a.iter().cloned() , vec![ 1,2,3] );           
// //! 
// //! 
//! //  -----------------------------------------------------------------------------------------
//! //  implement OracleMajor on a struct with a lifetime generic parameter
//! //  -----------------------------------------------------------------------------------------
//! 
//! //  NB: SO FAR AS WE ARE AWARE it is not straightforward to implement `OracleMajor` in a way
//! //      that leverages a generic lifetime parameter associated with a struct.  This is in 
//! //      contrast to some of the example shown below, which *do* leverage the lifetime
//! //      parameter of a struct (it is key to note that different traits are implemented, and
//! //      in the later examples the trait explicitly involves a lifetime parameter).  
//! //      
//! //  The following is a record of our "best attempt" at an implementation of `OracleMajor` that
//! //  does make use of a lifetime parameter.  The attempt failed, but others may have more
//! //  success in the future.
//! 
//! // use std::marker::PhantomData;   
//! 
//! // // Define the struct
//! // pub struct VecOfVecWithLifetime
//! // < 'a, IndexCoeffPair >
//! // {
//! //     pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
//! //     phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
//! // }
//! 
//! // // Implement the trait
//! // impl < 'a >
//! 
//! //         OracleMajor
//! //         <  usize, Cloned<Rev<std::slice::Iter<'a, (usize,usize)>>>  > for 
//! 
//! //         VecOfVecWithLifetime< 'a, (usize,usize) >
//! 
//! //         where   (usize,usize): Clone,
//! // {
//! //     fn view_major( & self, index: usize ) -> Cloned<Rev<std::slice::Iter<'a, (usize,usize)>>> {
//! //         return self.vec_of_vec[index].iter().rev().cloned()
//! //     } 
//! // }
//! 
//! 
//! //  -----------------------------------------------------------------------------------------
//! //  -----------------------------------------------------------------------------------------
//! //  Documentation of alternate oracle trait: with lifetimes
//! //  -----------------------------------------------------------------------------------------
//! //  -----------------------------------------------------------------------------------------        
//! use std::marker::PhantomData;
//! 
//! /// Entries may not appear in sorted order.
//! pub trait OracleMajorWithLifetime< 'a, MajKey, ViewMajor >
//! {
//!     /// Get a major vector.
//!     ///
//!     /// The order in which terms appear should be the same every time the
//!     /// function is called; however, the order need not be sorted.
//!     fn   view_major_with_life<'b: 'a>( &'b self, index: MajKey ) -> ViewMajor;
//! }
//! 
//! //  -----------------------------------------------------------------------------------------
//! //  implement OracleMajorWithLifetime on &'a mut Vec< Vec< (usize,usize) > >
//! //  -----------------------------------------------------------------------------------------
//! 
//! impl < 'a >
//! 
//!         OracleMajorWithLifetime
//!         <  'a, usize,  &'a [(usize,usize)]  > for 
//! 
//!         &'a mut Vec < Vec < (usize,usize) > > 
//! 
//!         where   (usize,usize): Clone,
//! {
//!     fn view_major_with_life<'b: 'a>( &'b self, index: usize ) -> &'a [(usize,usize)] {
//!         return self[index].as_slice()
//!     } 
//! }
//! 
//! // Define a mutable reference to a matrix
//! let mut matrix_core = vec![ vec![ (1,1), (2,2) ]  ];
//! let matrix_ref = &mut matrix_core;
//! 
//! // Get a major view from the mutable reference.     
//! let row = matrix_ref.view_major_with_life( 0 );
//! itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );
//! 
//! // Get a major view from a *nested* immutable reference.
//! let _row = (& matrix_ref).view_major_with_life( 0 );      
//! 
//! // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
//! // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` 
//! // below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
//! let _row_a = matrix_ref.view_major_with_life( 0 );
//! matrix_ref.push( vec![ (4,4), (5,5) ] );
//! let row_b = matrix_ref.view_major_with_life( 1 );   
//! 
//! itertools::assert_equal(row_b.iter().cloned() , vec![ (4,4), (5,5) ] );
//! // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] ); 
//! 
//! 
//! //  -----------------------------------------------------------------------------------------
//! //  implement OracleMajorWithLifetime on VecOfVecWithLifetime< 'a, (usize,usize) >
//! //  -----------------------------------------------------------------------------------------
//! 
//! // DEFINE THE STRUCT
//! pub struct VecOfVecWithLifetime
//!             < 'a, IndexCoeffPair >
//! 
//! {
//!     // pub major_dimension: MajorDimension, 
//!     pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
//!     pub phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
//! }
//! 
//! // DEFINE ITS METHODS
//! impl    < 'a, IndexCoeffPair >
//!         
//!         VecOfVecWithLifetime 
//!         < 'a, IndexCoeffPair > 
//! 
//! {
//!     // Make new (empty) VecOfVec. 
//!     pub fn new( vecvec: Vec < Vec < IndexCoeffPair > > ) -> Self  
//!     {
//!         VecOfVecWithLifetime{   
//!                     // major_dimension: major_dimension,
//!                     vec_of_vec:     vecvec,                    
//!                     phantom_kvpair: PhantomData,                   
//!                 }
//!     }
//! }
//! 
//! // IMPLEMENT OracleMajorWithLifetime
//! impl < 'a >
//! 
//!         OracleMajorWithLifetime
//!         <  'a, usize,  &'a Vec< (usize,usize) >  > for 
//! 
//!         VecOfVecWithLifetime< 'a, (usize,usize) >
//! 
//!         where   (usize,usize): Clone,
//! {
//!     fn view_major_with_life<'b: 'a>( &'b self, index: usize ) -> &'a Vec< (usize,usize) > {
//!         return &self.vec_of_vec[index]
//!     } 
//! }        
//! 
//! // Define a mutable reference to a matrix
//! let mut matrix = VecOfVecWithLifetime::new( vec![ vec![ (1,1), (2,2) ]  ] );
//! 
//! // Get a major view from the mutable reference.     
//! let row = matrix.view_major_with_life( 0 );
//! itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );
//! 
//! // Get a major view from a *nested* immutable reference.
//! let _row = (& matrix).view_major_with_life( 0 );      
//! 
//! // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
//! // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
//! let _row_a = matrix.view_major_with_life( 0 );
//! matrix.vec_of_vec.push( vec![ (4,4), (5,5) ] );
//! let row_b = matrix.view_major_with_life( 1 );   
//! 
//! itertools::assert_equal(row_b.iter().cloned() , vec![ (4,4), (5,5) ] );
//! // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] ); 
//! ```
//! 
//
// THE FOLLOWING SEGMENT OF DOCUMENTATION DATES BACK TO A TIME BEFORE WE USED ASSOCIATED TYPES FOR THE MATRIX ORACLES, AND HENCE WHEN TYPE INFERENCE WAS HARDER; IF WE EVER GO BACK TO THAT WAY OF DOING THINGS, WE MAY WANT TO REVIVE THIS EXAMPLE
//
// //! # Partial type annotations
// //! 
// //! Matrix oracles tend to generate a number of generic type parameters, especially when you start to
// //! place them in wrappers.  In this situation Rust may get confused and fail to infer the right
// //! parameters.  In these situations, you can often clear things up by specifying **some but not all*
// //! of the type paramers associated with a struct.  Here's an example:
// //! 
// //! ```
// //! use std::slice::Iter;
// //! use std::iter::Cloned;
// //! use itertools::Itertools;
// //! use oat_rust::matrices::random_constructors::random_upper_unitriangular;
// //! use oat_rust::matrices::operations::transform_vector_wise::VecWiseTransformed;
// //! 
// //! // construct an upper unitriangular matrix
// //! let matrix_size         =   3;
// //! let modulus             =   7;
// //! let array_mapping = random_upper_unitriangular( matrix_size, modulus );
// //! let array_mapping_ref = & array_mapping;
// //! 
// //! 
// //! // compute the codomain COMB of the matrix WHOSE COLUMNS APPEAR IN REVERSE ORDER
// //! let vector_transformer = |x: Cloned< Iter< '_, (usize,usize) > >| 
// //!                                                             x
// //!                                                                 .rev()
// //!                                                                 .map(   |(a,b)| 
// //!                                                                         ( modulus - 1 - a, b) 
// //!                                                                     )
// //!                                                                 .collect_vec() ;
// //! 
// //! let matrix_transformed: VecWiseTransformed< _, Cloned< Iter< '_, (usize,usize) > >, _>      =   VecWiseTransformed::new( 
// //!                                     & array_mapping_ref,
// //!                                     vector_transformer,
// //!                                 );
// //! ```