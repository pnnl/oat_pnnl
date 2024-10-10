
//! # Open Applied Topology
//! 
//! Open Applied Topology (OAT) is software for fast, user-friendly algebra and topology.
//! 
//! - [Welcome!](#welcome)
//! - [Community](#community)
//! - [Values](#values)
//! - [Mission](#features)
//! - [Get started](#get-started)
//! 
//! # Welcome!
//! 
//! Welcome!  This is the <span style="color: SteelBlue;">[Rust component](./link/to/contributor/doc/explaining/name/conventions)</span> of Open Applied Topology.  It's registered as `oat_rust` on Crates.io.  OAT provides powerful tools for applied topology, including
//! 
// //! - Persistent homology
// //! - Simplicial complexes
// //! - Homological algebra
// //! - Cycle optimization
// //! - Interactive 2d and 3d visualization
//! - [Algebra](crate::algebra)
//!   - [Rings](crate::algebra::rings)
//!   - [Vectors](crate::algebra::vectors)
//!   - [Matrices](crate::algebra::matrices)
//!   - [Chain complexes](crate::algebra::chains)
//!   - [Barcodes and persistence diagrams](crate::algebra::chains::barcode)
//!   - [(Persistent) (relative) (co)homology](crate::algebra::chains::jordan)
//!   - [Optimal cycle representatives](crate::algebra::chains::factored)
//! - [Topology](crate::topology)
//!   - [Simplicial complexes](crate::topology::simplicial)
//!   - [Boundary matrices](crate::topology::simplicial)
//!   - [Point clouds](crate::topology::point_cloud)
//! - [Utilities](crate::utilities)
//!   - [Iterators](crate::utilities::iterators)
//!   - [Optimization](crate::utilities::optimization)
//!   - [Order](crate::utilities::order)
//!   - [Sets](crate::utilities::sets)
//!   - [Indexing](crate::utilities::indexing)
//!   - and more!
//! - Interactive 2d and 3d visualization
//!   - See our sister library, OAT-Python
//! 
//! 
//! 
//! # Community
//! 
//! OAT is by and for the open source community.  <span style="color: SteelBlue;">Reach out to the developers</span> if you
//! - Need help getting started
//! - Wish for a missing feature
//! - Want to try coding
//! 
//! A collaboration of 20 research centers at colleges, universities, private, and public organizations support OAT's
//! development. The founding developers are Princton University, Macalester College, and the University of Delaware
//! The National Science Foundation granted seed funding for OAT in
//! 2019 under the [ExHACT]((https://www.nsf.gov/awardsearch/showAward?AWD_ID=1854748&HistoricalAwards=false))
//! project, and [Pacific Northwest National Laboratory](https://www.pnnl.gov/) (PNNL) provides continuing financial
//! support.  PNNL now coordinates development.  See <span style="color: SteelBlue;">[here](./ATTRIBUTIONS.md)</span>
//! for further details.
//! 
//! # <span style="color: SteelBlue;">[Values](./CODE_OF_CONDUCT.md)</span>
//! 
//! Our <span style="color: SteelBlue;">[shared values](./CODE_OF_CONDUCT.md)</span> are
//! 
//! - Inclusion
//! - Respect, and 
//! - A shared passion to expand human knowledge, through algebraic topology
//! 
//! 
//! # Mission
// //! 
// // //! OAT offers powerful features for beginners through advanced users.
// //!  
// //! <details>
// //! <summary>Performance</summary>
// //! <br>
// //! OAT is a first-class engine for cutting-edge applications.  It is specially suited to large, sparse // matrices.
// //!     The core library is written in Rust, a low-level systems programming language designed for safety and // performance.
// //!     High-level wrappers are available in Python. 
// //! </details>
// //! 
// //! <details>
// //! <summary>Reliability</summary>
// //! <br>
// //!   - Verification and test coverage: the OAT library is extensively tested.  In addition, the modular design // of the library makes it easy for users to generate their own certificates of correctness.
// //!   - Safety: OAT inherits strong safety guarantees from the features of the Rust compiler, especially in the // open source development process
// //! </details>
// //! 
// //! <details>
// //! <summary>Transparency</summary>
// //! <br>
// //!   - Documentation:  emphasizes clarity and accessibility for users with all backgrounds. OAT docs provide // explicit descriptions of both code *and* the underlying mathematical concepts. 
// //!   - [Tutorials](crate::tutorials) offer examples and helpful tips for beginners through advanced Rust users.
// //!   - Indexing: is one of the most pervasive challenges to writing transparent, interpretable code in // computational topology.  OAT's matrices and vectors can be indexed by simplices, cubes, and other user-defined // data structures, in addition to integers.
// //! </details>
// //! 
// //! <details>
// //! <summary>Modularity</summary>
// //! <br>
// //! 
// // //! Creative recombination of building blocks is among the most important ways we innovate.  
// //! 
// //!   OAT breaks problems into the same basic building blocks that topologists use when writing on a chalk // board. Users can mix and match those blocks with confidence, with a simple, streamlined interface.  They can // even create new components that work seemlessly with the rest of the library, including coefficient rings, // sparse matrix data structures, and customized filtrations on simplicial complexes.
// //! </details>
//! 
//! 
//! 
//! **Performance**
//!     
//! OAT is a first-class solver for cutting-edge applications.  It is ideally suited to large, sparse data sets.
//!     The core library is written in Rust, a low-level systems programming language designed for safety and performance.
//!     High-level wrappers are available in Python. 
//! 
//! \
//! 
//! **Reliability**
//! 
//! OAT has more unit tests than type definitions and function definitions, combined.
//!   Its modular design enables end users to write their own checks for correctness, with ease.
//!   The library inherits strong safety guarantees from the the Rust compiler.
//! 
//! 
//! **Transparency**
//! 
//! OAT documentation emphasizes clarity and accessibility for users with all backgrounds.  It includes more than 180 working examples, and describes both code and underlying mathematical concepts in detail.
//! [Tutorials](crate::tutorials) illustrate how to combine multiple tools into larger applications.
//! The platform's modular design breaks large solvers into basic components, which it exposes to the user for inspection.  In addition, the library provides powerful methods to inspect and analyze objects, consistent with the way humans naturally think about problems; for example, you can look up rows and columns of boundary matrices using *cubes*, *simplices*, or *cells* as keys.
//!   
//! 
//! **Modularity**
//! 
// //! Creative recombination of building blocks is among the most important ways we innovate.  
//! 
//!   OAT reduces complex problems to the same basic building blocks that topologists use when writing on a chalk board. Users can mix and match those blocks with confidence, using a simple, streamlined interface.  They can even create new components that work seemlessly with the rest of the library, including coefficient rings, sparse matrix data structures, and customized filtrations on simplicial complexes.
//! 
//! 
//! 
//! 
//! # Get Started
//! 
// //! #### Contents
// //! 
// // ! OAT includes, *but is not limited to*
// // ! 
// //! - [Algebra](crate::algebra)
// //!   - [Rings](crate::algebra::rings)
// //!   - [Vectors](crate::algebra::vectors)
// //!   - [Matrices](crate::algebra::matrices)
// //!   - [Chain complexes](crate::algebra::chains)
// //!   - [Barcodes and persistence diagrams](crate::algebra::chains::barcode)
// //!   - [(Persistent) (relative) (co)homology](crate::algebra::chains::jordan)
// //!   - [Optimal cycle representatives](crate::algebra::chains::factored)
// //! - [Topology](crate::topology)
// //!   - [Simplicial complexes](crate::topology::simplicial)
// //!   - [Boundary matrices](crate::topology::simplicial)
// //!   - [Point clouds](crate::topology::point_cloud)
// //! - [Utilities](crate::utilities)
// //!   - [Iterators](crate::utilities::iterators)
// //!   - [Optimization](crate::utilities::optimization)
// //!   - [Order](crate::utilities::order)
// //!   - [Sets](crate::utilities::sets)
// //!   - [Indexing](crate::utilities::indexing)
// //!   - and more!
//! 
//! **Try the tutorials**
//! 
//! - Check out the <span style="color: SteelBlue;">Python jupyter notebook tutorials</span>
//! - Try the [Rust tutorials](crate::tutorials), including a
//! [quick start guide to writing your first program with OAT](crate::tutorials::oat_quick_start).
//! 
//! **Find answers in the documentation**
//! 
//! Rust documentation places lists (of objects, functions, etc.) in two places: either the bottom of a page, or in the
//! menu bar on the left. You can also use the search bar at the top to pull up a list of related terms.  The question mark button to the right of the bar gives other helpful tips (for example, you can search for functions based on their type signature).
//! 
//! **Get started with Rust**
//! 
//! See this  <span style="color: SteelBlue;">orientation</span> for help
//! 
//! - installing
//! - debugging
//! - avoiding the most time-consuming parts of the coding process
//! 
//! 
//! 
//!
//! 
// //! # User guide
// //! 
// //! 
// //! OAT revolves around three objects: rings, vectors, and matrices.  However, the library does not 
// //! define a ring object, vector object, and matrix object.  Instead, it uses the Rust notion of a 
// //! [trait](https://doc.rust-lang.org/book/ch10-02-traits.html).  You can write
// //! your own object, and as long as it implements the appropriate traits, OAT will be 
// //! able to use it as a ring, vector, or matrix.  Here are key examples of how this works in practice:
// //! 
// //! *  **Ring** traits
// //! 
// //!     Most programming lanuages have "types" for working with certain rings.  For example, there 
// //!     might be a `float` type, an `integer` type, and a `bool` type.  OAT makes it easy
// //!     to define your *own* ring type.  To do so, you'll actually use two types: one for the elements of the ring, 
// //!     and another for 
// //!     the ring *operations* (addition, multiplication, etc.).  
// //! 
// //!     ```
// //!     // Let's try arithmetic with the two-element (Boolean) field.  You can code your own 
// //!     // implementation of this field, but for brevity we'll use OAT's pre-built 
// //!     // implementation.  
// //! 
// //!     // First get some elements. 
// //!     let x = true; // the type of elements is `bool`
// //!     let y = false;
// //! 
// //!     // Then define the ring operator.
// //!     use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator; // this line imports the contstructor that will build the ring operator
// //!     use oat_rust::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing}; // this line imports the traits for rings, semirings, and division rings
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
// //!     There are separate traits for [rings](`crate::algebra::rings::operator_traits::Ring`), [semirings](`crate::algebra::rings::operator_traits::Semiring`), and [division rings](``crate::algebra::rings::operator_traits::DivisionRing``).
// //! 
// //! * **Vector entry**  traits
// //! 
// //!     An entry in a vector `v`, is a pair `(i, a)` such that `v[i] = a`.  
// //!     There are many ways to store a vector entry in computer memory, e.g. in a tuple, 
// //!     a list, a dictionary, etc.  There are some operations we'd like to perform on an entry,
// //!     no matter the data structure stores it:
// //! 
// //!     * [KeyValGet](crate::algebra::vectors::entries::KeyValGet) 
// //!         allows one to determine the value of  `i` or `a`.  
// //! 
// //! 
// //!     * [KeyValSet](crate::algebra::vectors::entries::KeyValSet)
// //!         allows one to change the value of  `i` or  `a`.  
// //! 
// //!     Here are some examples:
// //!     
// //!     ```
// //!     // Import the KeyValGet and KeyValSet traits, so that we can use them.
// //!     use oat_rust::algebra::vectors::entries::{KeyValGet, KeyValSet}; 
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
// //!     OAT requires no traits at all to represent sparse vectors: 
// //!     any [iterator](https://doc.rust-lang.org/book/ch13-02-iterators.html) that runs over sparse vector entries 
// //!     (an entry is an object that implements the [KeyValGet](crate::algebra::vectors::entries::KeyValGet) trait) 
// //!     can be used to represent a  sparse vector.
// //! 
// //! 
// //! * **Matrix**  traits
// //!       
// //! 
// //!     In OAT, the most common task that involves a sparse matrix is to look up the nonzero entries in one of its
// //!     rows or columns.  The so-called "oracle" traits provide a simple framework for performing this type of look-up operation. 
// //!     Details about the oracle traits can be found in the [matrix_oracle](algebra::matrices::matrix_oracle_traits) module.
// //! 
// //!    *[Vocabulary: major and minor dimensions] In many, though not all, sparse 
// //!         matrices, it's easier to look up rows than to look up columns, or vice versa.
// //!         We call the easy dimension the *major dimension*.  When we need to look up a row or column 
// //!         of a matrix, we do so by indexing into the appropriate major or minor dimension; see 
// //!         [here](crate::algebra::matrices::matrix_oracle_traits) for more details.*
// // //!         The [WhichMajor](algebra::matrices::query::WhichMajor) traight indicates whether
// // //!         should think of the "easy" dimension as rows or as columns. 
// //!     
// //! 
// //!     The most commonly used oracle traits are:
// //! 
// //!     [ViewRow](algebra::matrices::query::ViewRow): returns the entries in a row (if the matrix is row-major) or 
// //!     column (if the matrix is column-major).  Entries may not appear in sorted order. <br />
// //!     [ViewRowAscend](algebra::matrices::query::ViewRowAscend): returns entries in ascending order of index <br />
// //!     [ViewRowDescend](algebra::matrices::query::ViewRowDescend): returns entries in descending order of index <br />
// //!     [ViewCol](algebra::matrices::query::ViewCol): returns the entries in a row (if the matrix is column-major) or 
// //!     column (if the matrix is row-major).  Entries may not appear in sorted order. <br />
// //!     [ViewColAscend](algebra::matrices::query::ViewColAscend): returns entries in ascending order of index, <br />
// //!     [ViewColDescend](algebra::matrices::query::ViewColDescend): returns entries in descending order of index, <br />
// //!        
// //!        
// //!     ```
// //!     // Import some packages to work with sparse vec-of-vec matrices, and other traits.
// //!     use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
// //!     use oat_rust::algebra::matrices::query::{MajorDimension, ViewRow};
// //!     use std::iter::FromIterator;
// //! 
// //!     // Create a vector-of-vectors sparse matrix.  
// //!     // In particular, we will construct a 2x2 upper triangular matrix with 1's above the diagonal.  
// //!     // The matrix is row-major; for vec-of-vec matrices, that means that each vector represents a row.
// //!     let matrix_data =    VecOfVec::new(
// //!                         vec![ vec![(0, 1.), (1, 1.)], vec![(1, 1.)] ], // the vector of vectors,
// //!                     );
// //!     let matrix  =   & matrix_data;
// //! 
// //!     // Access a row.
// //!     // Since this matrix is row-major, we use the ViewRow trait.  This trait accesses
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
// //! // Note that `vec!` creates a rust Vec, which is a general data structure for storing lists
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
// //! use oat_rust::algebra::rings::operator_structs::ring_native::*;
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
// //!     use oat_rust::algebra::vectors::operations::*;
// //!     use oat_rust::algebra::rings::operator_structs::ring_native::*;
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


pub mod algebra;
pub mod topology;
pub mod utilities;
// pub mod developer;
pub mod tutorials;
