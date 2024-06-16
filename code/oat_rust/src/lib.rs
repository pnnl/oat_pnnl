
//! Sparse linear algebra in Rust (oat_rust).
//! 
//! # Quick start tutorials
//! 
//! This library provides a variety of [examples and demonstrations](crate::tutorials), including a
//! [quick start guide to writing your first program with oat_rust](crate::tutorials::oat_rust_quick_start).
//! 
//! For a quick start guide to computing homology and persistent homology, refer to the sohar_rust library.
//! 
// //! .  These include:
// //! 
// //! - [Write your first program with oat_rust](crate::tutorials::oat_rust_quick_start)
// //! - [Compute homology](crate::tutorials::compute_homology)
// //! - [Compute persistent homology](crate::tutorials::compute_ph)
//! 
//! # About
//! 
//! oat_rust is a scientific computing library for linear and homological algebra.  It provides a flexible framework to construct and analyze
//! 
//! - [coefficient rings](crate::rings) 
// //!   - integers, finite fields, rational numbers, etc.; users may also define their own rings
//! - [vectors](crate::vectors) and [vector entries](crate::entries)
//! - [matrices](crate::matrices)
//! 
//! and time and memory efficient tools for operations such as
//! 
//! - [accessing entries in rows and columns](crate::matrices::matrix_oracle_traits)
//! - [matrix and vector multiplication](crate::matrices::operations::multiply)
//! - [solving systems of equations](crate::matrices::operations::solve)
//! - [matrix factorization](crate::matrices::operations::umatch)
//! - [matrix inversion](crate::matrices::operations::invert)
//! - [**computing bases for images and kernels**](crate::matrices::operations::umatch)
//! - [transforming and combining vectors (scale, add, drop zeros, reindex, etc.)](crate::vectors::operations)
//! 
//! Key features include
//! 
//! - flexible indexing: matrices can be indexed by arbitrary keys; for example, the boundary matrix of a simplicial complex can be indexed by simplices; this feature is powerful in a variety of homology computations, where explit enumeration of row and column indices is expensive
//! - lazy look up: rows and columns of matrices can be constructed in a lazy fashion, consistent with current state-of-the-art practices in large scale PH computation
//! - extensive unit testing
//! 
//!
//! # Resources
//! 
//! **oat_rust**
//! 
// //! There are two primary resources to learn about the oat_rust library:
// //! 
// //! * **Overview** This page gives a high-level introduction to the package.
// //! 
// //! * **Documentation** A full API, including objects, functions, etc. can be found at the bottom of this page.
// //! 
//! 
//! Full documentation for oat_rust, including objects, functions, etc. can be found at the bottom of this page.
//! If you are not able to find what you need, you are welcome to reach out to the development team on our [Github page](https://github.com/ExHACT/oat_rust).
//! 
//! **RUST**
//! 
//! Rust is a low-level programming language with powerful features for scientific computing, e.g. memory safety and helpful error messages.  It has been voted [most-loved language](https://insights.stackoverflow.com/survey/2021) by the worldwide developer community since 2015.
//! 
//! * **Installation** See the Rust website for directions on [installation](https://www.rust-lang.org/learn/get-started)
//!
//! * **Visual studio code editor** If you are new to Rust, you may find the [VS Code editor](https://code.visualstudio.com/docs/languages/rust) helpful.  It has a very useful [`rust-analyzer`](https://marketplace.visualstudio.com/items?itemName=rust-lang.rust-analyzer) extension, with syntax highlighting, help for debugging, code completion, and more.
//! 
//! * **Debugging** Rust is very good about providing helpful error messages, as a rule.
//!   It is often possible to find help just by copying these messages into a web search.
//!   The oat_rust developers have also collected a short list of [debugging tips and tricks](crate::developer::rust_debugging), which you may find useful.
//! 
// //! # User guide
// //! 
// //! 
// //! oat_rust revolves around three objects: rings, vectors, and matrices.  However, the library does not 
// //! define a ring object, vector object, and matrix object.  Instead, it uses the Rust notion of a 
// //! [trait](https://doc.rust-lang.org/book/ch10-02-traits.html).  You can write
// //! your own object, and as long as it implements the appropriate traits, oat_rust will be 
// //! able to use it as a ring, vector, or matrix.  Here are key examples of how this works in practice:
// //! 
// //! *  **Ring** traits
// //! 
// //!     Most programming lanuages have "types" for working with certain rings.  For example, there 
// //!     might be a `float` type, an `integer` type, and a `bool` type.  oat_rust makes it easy
// //!     to define your *own* ring type.  To do so, you'll actually use two types: one for the elements of the ring, 
// //!     and another for 
// //!     the ring *operations* (addition, multiplication, etc.).  
// //! 
// //!     ```
// //!     // Let's try arithmetic with the two-element (Boolean) field.  You can code your own 
// //!     // implementation of this field, but for brevity we'll use oat_rust's pre-built 
// //!     // implementation.  
// //! 
// //!     // First get some elements. 
// //!     let x = true; // the type of elements is `bool`
// //!     let y = false;
// //! 
// //!     // Then define the ring operator.
// //!     use oat_rust::rings::operator_structs::field_prime_order::BooleanFieldOperator; // this line imports the contstructor that will build the ring operator
// //!     use oat_rust::rings::operator_traits::{Semiring, Ring, DivisionRing}; // this line imports the traits for rings, semirings, and division rings
// //!     let ring_operator = BooleanFieldOperator{}; // the type of the ring operator is `BooleanFieldOperator`
// //! 
// //!     // Finally, use the ring operator to perform basic ring operations:
// //!     assert_eq!( true,   ring_operator.add( x, y ) );
// //!     assert_eq!( true,   ring_operator.subtract( x, y ) );
// //!     assert_eq!( false,  ring_operator.multiply( x, y ) );
// //!     assert_eq!( false,  ring_operator.divide( y, x ) );
// //!     ```
// //!     
// //!     You can use any object as a ring operator, as long as it implements the proper traits.
// //!     There are separate traits for [rings](`crate::rings::operator_traits::Ring`), [semirings](`crate::rings::operator_traits::Semiring`), and [division rings](``crate::rings::operator_traits::DivisionRing``).
// //! 
// //! * **Vector entry**  traits
// //! 
// //!     An entry in a vector `v`, is a pair `(i, a)` such that `v[i] = a`.  
// //!     There are many ways to store a vector entry in computer memory, e.g. in a tuple, 
// //!     a list, a dictionary, etc.  There are some operations we'd like to perform on an entry,
// //!     no matter the data structure stores it:
// //! 
// //!     * [KeyValGet](crate::entries::KeyValGet) 
// //!         allows one to determine the value of  `i` or `a`.  
// //! 
// //! 
// //!     * [KeyValSet](crate::entries::KeyValSet)
// //!         allows one to change the value of  `i` or  `a`.  
// //! 
// //!     Here are some examples:
// //!     
// //!     ```
// //!     // Import the KeyValGet and KeyValSet traits, so that we can use them.
// //!     use oat_rust::entries::{KeyValGet, KeyValSet}; 
// //!
// //!     // Define a vector entry.
// //!     // Every tuple of length 2 (that meets certain requirements) implements the KeyValGet and KeyValSet traits, automatically.
// //!     let mut vector_entry = (1, 0.);
// //!
// //!     // Use the methods associated with these traits.
// //!     // The associated methods are `.key()`, `.val()`, `.set_key`, and `.set_val`
// //!     assert_eq!( vector_entry.key(), 1  ); // the .key() method retreives the index
// //!     assert_eq!( vector_entry.val(), 0. ); // the .val() method retreives the coefficient
// //!
// //!     vector_entry.set_key( 2 );            // the .set_key() method sets the index
// //!     assert_eq!( vector_entry.key(), 2  ); 
// //!     ```
// //! 
// //! * **Vector** traits
// //! 
// //!     oat_rust requires no traits at all to represent sparse vectors: 
// //!     any [iterator](https://doc.rust-lang.org/book/ch13-02-iterators.html) that runs over sparse vector entries 
// //!     (an entry is an object that implements the [KeyValGet](crate::entries::KeyValGet) trait) 
// //!     can be used to represent a  sparse vector.
// //! 
// //! 
// //! * **Matrix**  traits
// //!       
// //! 
// //!     In oat_rust, the most common task that involves a sparse matrix is to look up the nonzero entries in one of its
// //!     rows or columns.  The so-called "oracle" traits provide a simple framework for performing this type of look-up operation. 
// //!     Details about the oracle traits can be found in the [matrix_oracle](matrices::matrix_oracle_traits) module.
// //! 
// //!    *[Vocabulary: major and minor dimensions] In many, though not all, sparse 
// //!         matrices, it's easier to look up rows than to look up columns, or vice versa.
// //!         We call the easy dimension the *major dimension*.  When we need to look up a row or column 
// //!         of a matrix, we do so by indexing into the appropriate major or minor dimension; see 
// //!         [here](crate::matrices::matrix_oracle_traits) for more details.*
// // //!         The [WhichMajor](matrices::matrix_oracle_traits::WhichMajor) traight indicates whether
// // //!         should think of the "easy" dimension as rows or as columns. 
// //!     
// //! 
// //!     The most commonly used oracle traits are:
// //! 
// //!     [OracleMajor](matrices::matrix_oracle_traits::OracleMajor): returns the entries in a row (if the matrix is row-major) or 
// //!     column (if the matrix is column-major).  Entries may not appear in sorted order. <br />
// //!     [OracleMajorAscend](matrices::matrix_oracle_traits::OracleMajorAscend): returns entries in ascending order of index <br />
// //!     [OracleMajorDescend](matrices::matrix_oracle_traits::OracleMajorDescend): returns entries in descending order of index <br />
// //!     [OracleMinor](matrices::matrix_oracle_traits::OracleMinor): returns the entries in a row (if the matrix is column-major) or 
// //!     column (if the matrix is row-major).  Entries may not appear in sorted order. <br />
// //!     [OracleMinorAscend](matrices::matrix_oracle_traits::OracleMinorAscend): returns entries in ascending order of index, <br />
// //!     [OracleMinorDescend](matrices::matrix_oracle_traits::OracleMinorDescend): returns entries in descending order of index, <br />
// //!        
// //!        
// //!     ```
// //!     // Import some packages to work with sparse vec-of-vec matrices, and other traits.
// //!     use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
// //!     use oat_rust::matrices::matrix_oracle_traits::{MajorDimension, OracleMajor};
// //!     use std::iter::FromIterator;
// //! 
// //!     // Create a vector-of-vectors sparse matrix.  
// //!     // In particular, we will construct a 2x2 upper triangular matrix with 1's above the diagonal.  
// //!     // The matrix is row-major; for vec-of-vec matrices, that means that each vector represents a row.
// //!     let matrix_data =    VecOfVecSimple::new(
// //!                         vec![ vec![(0, 1.), (1, 1.)], vec![(1, 1.)] ], // the vector of vectors,
// //!                     );
// //!     let matrix  =   & matrix_data;
// //! 
// //!     // Access a row.
// //!     // Since this matrix is row-major, we use the OracleMajor trait.  This trait accesses
// //!     // vectors along the major dimension via the command `view_major`
// //!     let row_iterator    =   matrix.view_major( 0 ); // access the 0th row.  the result is an iterable, i.e. a struct that can be easily transformed into an iterator.
// //!     let row_vector      =   Vec::from_iter( row_iterator ); // collect the elements of the iterator into a Rust vector
// //!     assert_eq!( row_vector, vec![(0, 1.), (1, 1.)] );
// //! 
// //!     ```
// 
// //! 
// //! 
// //! 
// //! ## Operations on sparse vector iterators
// //! 
// //! Let `iter_a`, `iter_b`, and `iter_c` be sparse vectors, i.e. iterators that run over 
// //! sparse matrix entries.  For example, we could define `iter_a`, `iter_b`, `iter_c` as follows
// //! 
// //! ```
// //! // First define the entries
// //! // Note that `vec!` creates a standard rust vector, which is different 
// //! // from the sort of vector we care about)
// //! let entries_a   =   vec![ (1, 1.), (4, 4.) ];
// //! let entries_b   =   vec![ (2, 2.), (3, 3.) ];
// //! let entries_c   =   vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
// //! 
// //! // Now define the sparse vector iterators
// //! // Note that `iter()` creates an interator, and `cloned()` reformats 
// //! // the entries of each iterator.
// //! let iter_a      =   entries_a.iter().cloned(); 
// //! let iter_b      =   entries_b.iter().cloned();
// //! let iter_c      =   entries_c.iter().cloned();
// //! ```
// //! 
// //! Let's also define the operator of the coefficient ring we want to work with.
// //!     
// //! ```
// //! // Load the module that allows us to define our coefficient ring.
// //! use oat_rust::rings::operator_structs::ring_native::*;
// 
// //! // Define the operator of a coefficient ring (in this case, floating point real numbers)
// //! let ring_operator = DivisionRingNative::<f64>::new();   
// //! ```
// //!     
// //! * **Scale, drop zeros, gather, simplify**
// //! 
// //!     We can scale, drop zero entries, and gather terms as follows
// //!     
// //!     ```
// //!     use oat_rust::vectors::operations::*;
// //!     use oat_rust::rings::operator_structs::ring_native::*;
// //! 
// //!     # // Define the vector
// //!     # let entries_a   =   vec![ (1, 1.), (4, 4.) ];
// //!     # let entries_c   =   vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
// //!     # let iter_a      =   entries_a.iter().cloned();
// //!     # let iter_c      =   entries_c.iter().cloned();
// //!     #
// //!     # // Define the operator of the coefficient ring (in this case, floating point real numbers)
// //!     # let ring_operator = DivisionRingNative::<f64>::new();        
// //!       
// //!     // SCALE A VECTOR BY 2.
// //!     // Example: convert [ (1, 1.), (4, 4.) ] into [ (1, 2.), (4, 8.) ]
// //!     let scaled : Vec<_> = iter_a
// //!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                             .scale( 2., ring_operator.clone() )
// //!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     assert_eq!( scaled, vec![ (1, 2.), (4, 8.) ]);
// //!       
// //!     // DROP ZERO ENTRIES
// //!     // Example: convert [ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ] into [ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]
// //!     let dropped : Vec<_> = iter_c
// //!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                             .drop_zeros( ring_operator.clone() )
// //!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     
// //!     assert_eq!( dropped, vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]);
// //!       
// //!     // GATHER CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX
// //!     // The resulting vector has no repeating consecutive indices; each index gets 
// //!     // the sum of the corresponding coefficients.
// //!     // Example: convert [(1,1.), (1,0.5), (2,0.), (1,0.)] into [(1,1.5), (2,0.), (1,0.)]
// //!     let gathered : Vec<_> = iter_c
// //!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                             .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
// //!                             .gather( ring_operator.clone() )
// //!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     assert_eq!( gathered, vec![ (1, 1.), (2, 2.), (3, 6.), (4, 0.) ]);   
// //! 
// //!     // GATHER CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX, AND DROP RESULTING ZEROS
// //!     let simplified : Vec<_> = iter_c
// //!                                 .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                                 .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
// //!                                 .simplify( ring_operator.clone() )
// //!                                 .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     assert_eq!( simplified, vec![ (1, 1.), (2, 2.), (3, 6.), ]);  
// //!     ```
// //! * **Combine iterators in sorted order** (basic)
// //! 
// //!   We can combine two iterators, `A` and `B`, into a single iterator `C` using the 
// //!     [merge](https://docs.rs/itertools/0.7.2/itertools/fn.merge.html) function from [itertools](https://docs.rs/itertools/latest/itertools/). 
// //!     The resulting iterator, `C`, is a 
// //!     [Merge struct](https://docs.rs/itertools/0.7.8/itertools/structs/struct.Merge.html).
// //!     Iterator `C` will iterate over all the entries in `A` and `B`.
// //!     If the items of `A` and `B` appear in sorted order, then the items of `C` will also 
// //!     appear in sorted order.
// //!     ```
// //!     # let entries_a   =   vec![ (1, 1.), (4, 4.) ];
// //!     # let entries_b   =   vec![ (2, 2.), (3, 3.) ];
// //!     # let iter_a      =   entries_a.iter().cloned(); 
// //!     # let iter_b      =   entries_b.iter().cloned();
// //! 
// //!     use itertools::merge;
// //!     use std::iter::FromIterator;
// //!     
// //!     // Merge [ (1, 1.), (4, 4.) ] and [ (2, 2.), (3, 3.) ].
// //!     // The entries in these vectors are in sorted order, so the resulting iterator will 
// //!     // also produce items in sorted order.
// //!     let iter_merged   =   merge(iter_a, iter_b);
// //!     let entries_mrgd  =   Vec::from_iter(iter_merged);
// //!     assert_eq!( entries_mrgd, vec![ (1, 1.), (2, 2.), (3, 3.), (4, 4.) ])
// //!     ```
// //! 
// //!     We can merge any `k` iterators of the same kind using the 
// //!     [merge](https://docs.rs/itertools/0.7.2/itertools/fn.merge.html) 
// //!     function from 
// //!     [itertools](https://docs.rs/itertools/latest/itertools/).  
// //! 
// //! * **Combine iterators in sorted order** (advanced)
// //!  
// //!     For advanced usage (eg matrix reduction), we also provide a
// //!     customized merge process in the [hit_merge](utilities::iterators::hit_merge) module.
// //! 
// //! * **Add**
// //! 
// //!     We can add the vectors represented by `iter_a` and `iter_b` by 
// //!     first combining (e.g., with the `merge` function discussed above), 
// //!     then applying the `gather` method.
// //! 
// //! * **Subtract**
// //! 
// //!     We can subtract  `iter_a` from `iter_b` by first scaling `iter_a` by `-1`, then adding.
// //!̦
// //! 


pub mod rings;
pub mod vectors;
pub mod matrices;
pub mod utilities;
pub mod entries;
pub mod developer;
pub mod tutorials;
