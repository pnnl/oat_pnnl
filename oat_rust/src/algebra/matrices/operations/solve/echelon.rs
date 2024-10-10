//! Solve `Ax = b` for `x`, where `A` is a matrix partially in echelon form.
//! 
//! - [Background](#background)
//! - [Solve xA = b](#solvexab), with several [alternate formulations](#alternate-formulations)
//! - [Solve Ay = c](#solveayc)
//! - [Quotients, remainders, and solutions](#qrsolutions)
//! - [Panics](#panics) will (eventually) occur when `b` or `c` aren't properly sorted.
//! 
//! # Background
//! 
//! Suppose we have an invertible upper triangular matrix `A` and a vector `c`.
//! 
//! We can solve `xA = c` for `c` using [back substitution](https://algowiki-project.org/en/Backward_substitution).
//! This process adds linear multiples of the rows of `A` to `c` iteratively, to eliminate leading entries, until
//! the vector has been reduced to zero. 
//! 
//! If we permute the rows of `A` we can still solve `xA = c` via back-substitution.  To do this efficiently, we
//! just need an efficient way to look up the row with leading entry `i`, for each column index `i = 0, 1, 2, ..`.
//! We refer to the map sending column indices to the corresponding row indices the **leading entry bijection**.
//! 
//! If we delete some rows from the permuted matrix `A` then it may no longer be possible to solve `xA = c` for x.
//! However, we can still apply the back substitution algorithm using just the rows that haven't been deleted --
//! stopping when we encounter a leading entry that can't be eliminated with the leading entry of a row from `A`.
//! This algorithm will produce a vector `q`.  We call
//! 
//! - `q` the **quotient**, and
//! - `c - Aq` the **remainder**
//! 
//! Since we've deleted some rows, there's no longer a nice bijection from column indices to row indices.  But there's
//! still a bijection from a subset of column indices to the row indices that haven't been deleted.  We call this the
//! 
//! - **partially defined leading entry bijection**
//! 
//! 
// //! We start with a matrix `A` and a vector `b`, were `b` is an iterator that runs over the
// //! structural nonzero entries of a sparse vector.  
// //! 
// //! Here is a conceptual example of what is meant by "echelon solve".  Suppose we have
// //! a matrix `A` and vector `b` with the sparse structure show below (`*` indicates a nonzero entry).
// // //! It can be shown that `b` lies in the row space of `A` iff `b[1:3]` lies in the
// // //! row space of the submatrix `A[2:4, 1:3]`.  
// //! If there exists a 
// //! vector `x` such that `xA = b`, then we can find such an `x` via "back substituion": first we add a scalar
// //! multiple of `A[2,:]` to clear the first entry in `b`, resulting in a vector `c`.
// //! Then we add a scalar multiple of `A[3,:]` to clear the first nonzero entry in `c` (which occurs in column 2).
// //! This process continues until we have eliminated all nonzero entries.  The
// //! solution `x` can then be constructed from the collection of scalar multipliers
// //! used to perform the elimination.
// //! 
// //! 
// //! ```text
// //!              1  2  3  4  5  6  7
// //!           -------------------------
// //!         1 |                        |
// //!         2 |  *  *  *  *  *  *  *   |
// //!   A  =  3 |     *  *  *  *  *  *   |
// //!         4 |        *  *  *  *  *   |
// //!         5 |                        |
// //!           -------------------------
// //! 
// //!              1  2  3  4  5  6  7
// //!            -------------------------
// //!   b  =    |  *  *  *  *  *  *  *   |
// //!            -------------------------
// //! ```
// //! 
// //! We can use an analogous method whenever we have the following data
// //! 
// //! - a matrix `A`
// //! - a vector `b`
// //! - a set of row indices `I` and a set of column indices `J` such that `A[I,J]` is upper triangular and invertible.
// //! 
// //! we require, in addition, that 
// //! 
// //! - if `[i,j]` is a diagonal entry of `A[I,J]`, then all entries to the left of `[i,j]` in `A` vanish.
// //!   - this applies only in cases where we we add major views of `A` to `b`
// //! - if `[i,j]` is a diagonal entry of `A[I,J]`, then all entries below `[i,j]`  in `A` vanish.
// //!   - this applies only in cases where we we add minor views of `A` to `b`
// //! 
// //! In practice, when we use OAT to solve `xA = b`, we do not provide `I` and `J` explicitly;
// //! rather we note that there are (mutually inverse) bijections `f: I -> J` and `g: J -> I` such that
// //! `[i, f(i)]` and `[g(j), j]` lie on the diagonal of `A[I,J]`.  We pass these bijections to the solver, instead.
//! 
//! 
//! # <a name="solvexab">Solve xA = b for y</a>
//! 
//! OAT's echelon solvers will compute a quotient and remainder, given the following inputs
//! 
//! - a matrix `A` where no two rows share a leading nonzero entry in the same column
//! - a partially defined leading entry bijection
//! - an order operator (used to specify the order of row entries)
//! 
//! 
//!  
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified};        
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanFieldOperator::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                         vec![ 
//!                             vec![                                     ], 
//!                             vec![                                     ],                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );    
//! 
//! // partial bijection from the column indices of the diagonal elements to their row indices
//! let partial_bijection = |x: usize| { if x < 3 { Some(x+2) } else { None } };
//! 
//! // the problem vector b
//! let b = vec![ (0,true), (1,true) ];
//! 
//! // solver; attempts to solve xA = b for x
//! let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
//!                                                                 b.clone(), // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* minor keys to matched major keys
//!                                                                 BooleanFieldOperator::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! // multiply the solution with A, and check that the product equals b                               
//! let product = division.solution().unwrap().multiply_matrix_major_ascend( &matrix, ring_operator, OrderOperatorByKey::new() );
//! assert_eq!( product.collect_vec(), b );   // check the solution, i.e. check that xA = b  
//! ```
//! 
//! # Alternate formulations
//! 
//! There are many ways to specify the partially defined leading entry bijection.
//! You can use any object that implements the right version of the [EvaluateFunction] trait.  For example,
//! you can replace `let partial_bijection = |x: usize| { if x < 3 { Some(x+2) } else { None } };` in the preceding
//! example with a `Vec`:
//! 
//! ```
//! # use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! # use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! # use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified};        
//! # use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! # use oat_rust::algebra::vectors::operations::VectorOperations;
//! # use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! # use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! # use std::collections::HashMap;
//! # use itertools::Itertools;
//! # use assert_panic::assert_panic;
//! # 
//! # // a matrix A with an invertible upper triangular submatrix
//! # let matrix  =   VecOfVec::new(
//! #                         vec![ 
//! #                             vec![                                     ], 
//! #                             vec![                                     ],                                     
//! #                             vec![ (0,true),  (1,true), (2,true)       ], 
//! #                             vec![            (1,true), (2,true)       ],
//! #                             vec![                      (2,true)       ],                                       
//! #                         ],
//! #                     );  
//! # 
//! # // the problem vector b
//! # let b = vec![ (0,true), (1,true) ];
//! # 
//! # // the partial bijection
//! let partial_bijection = vec![2,3,4]; // or as a vector
//! # 
//! # // create a solver to solve xA = b for x
//! # let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
//! #                                                                 b.clone(), // entries must appear in strictly ASCENDING order
//! #                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//! #                                                                 partial_bijection, // maps *matched* minor keys to matched major keys
//! #                                                                 BooleanFieldOperator::new(), // defines the ring operations
//! #                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//! #                                                             );
//! #                                                         
//! # // check the solution, i.e. check that xA = b
//! # let product = division.solution().unwrap().multiply_matrix_major_ascend( &matrix, BooleanFieldOperator::new(), OrderOperatorByKey::new() );
//! # assert_eq!( product.collect_vec(), b.clone() );  
//! ```
//! 
//! or a `HashMap`: 
//! 
//! ```
//! # use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! # use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! # use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified};        
//! # use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! # use oat_rust::algebra::vectors::operations::VectorOperations;
//! # use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! # use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! # use std::collections::HashMap;
//! # use itertools::Itertools;
//! # use assert_panic::assert_panic;
//! # 
//! # // a matrix A with an invertible upper triangular submatrix
//! # let matrix  =   VecOfVec::new(
//! #                         vec![ 
//! #                             vec![                                     ], 
//! #                             vec![                                     ],                                     
//! #                             vec![ (0,true),  (1,true), (2,true)       ], 
//! #                             vec![            (1,true), (2,true)       ],
//! #                             vec![                      (2,true)       ],                                       
//! #                         ],
//! #                     );  
//! # 
//! # // the problem vector b
//! # let b = vec![ (0,true), (1,true) ];
//! # 
//! # // the partial bijection
//! let partial_bijection = HashMap::from( [(0,2),(1,3),(2,4)] ); // or as a hashmap
//! # 
//! # // create a solver to solve xA = b for x
//! # let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
//! #                                                                 b.clone(), // entries must appear in strictly ASCENDING order
//! #                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//! #                                                                 partial_bijection, // maps *matched* minor keys to matched major keys
//! #                                                                 BooleanFieldOperator::new(), // defines the ring operations
//! #                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//! #                                                             );
//! #                                                         
//! # // check the solution, i.e. check that xA = b
//! # let product = division.quotient().multiply_matrix_major_ascend( &matrix, BooleanFieldOperator::new(), OrderOperatorByKey::new() );
//! # assert_eq!( product.collect_vec(), b.clone() );  
//! ```
//! 
//! 
//! 
//! # <a name="solveayc">Solve Ay = c for y</a>
//! 
//! When solving for columns we have to eliminate **trailing** entries, rather than leading entries, in order to take advantage
//! of the upper triangular structure.  Therefore the entries of `c` have to appear in **strictly descending order**.
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified};        
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use std::collections::HashMap;
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! let matrix  =   VecOfVec::new(
//!                         vec![ 
//!                             vec![                                     ], 
//!                             vec![                                     ],                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );  
//! 
//! let c = vec![ (3,true), (2,true) ];   // define a sparse vector c                                                                                                                          
//! let division =  EchelonSolverMinorDescendWithMinorKeys::solve(
//!                                                                 c.clone(), // entries must appear in strictly DESCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 HashMap::from( [(2,0),(3,1),(4,2)] ), // or as a hashmap, // maps *matched* major keys to matched minor keys
//!                                                                 BooleanFieldOperator::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
//!                                                             );
//!                                                         
//! let product = division.quotient().multiply_matrix_minor_descend( &matrix, BooleanFieldOperator::new(), OrderOperatorByKey::new() );
//! 
//! // solution satisfies Ay = c      
//! assert_eq!( 
//!     product.collect_vec(), 
//!     c.clone() 
//! );   
//! ```
//! 
//! # <a name="qrsolutions">Quotients, remainders, and solutions</a>
//! 
//! The solvers return a [DivisionWithRemainder] object, that can be used to extract quotients, remainders, and solutions.
//! 
//! **Fact** a quotient is a solution iff the corresponding remainder is zero.
//! 
//! Calling `.solution()` on a [DivisionWithRemainder] object will therefore return `Some(quotient)` if the remainer is zero, and `None`
//! otherwise.
//! 
//! 
//! #### Case 1: There exists a solution
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified, };        
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use std::collections::HashMap;
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanFieldOperator::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                         vec![ 
//!                             vec![                                     ], 
//!                             vec![                                     ],                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );  
//! 
//! // create a solver to solve xA = b for x
//! let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
//!                                                                 vec![ (0,true), (1,true) ], // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* minor keys to matched major keys
//!                                                                 BooleanFieldOperator::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! // the remainder is zero
//! assert_eq!( 
//!     division.clone().remainder().collect_vec(),
//!     vec![ ]
//! );
//!
//! // therefore the solution equals the quotient
//! assert_eq!( 
//!     division.clone().solution(), 
//!     Some( division.clone().quotient().collect_vec() ) 
//! ); 
//! 
//! // we can also split the division object into a separate quotient and remainder
//! let qr = division.clone().quotient_remainder();
//! assert_eq!( 
//!     qr.remainder.clone().collect_vec(), 
//!     vec![] 
//! );
//! assert_eq!( 
//!     division.clone().solution(), 
//!     Some( qr.quotient.clone() ),
//! );
//! ```
//! 
//! 
//! #### Case 2: There exists no solution
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified};        
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use std::collections::HashMap;
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanFieldOperator::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                         vec![ 
//!                             vec![                                     ], 
//!                             vec![                                     ],                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );  
//! 
//! // create a solver to solve xA = b for x
//! let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
//!                                                                 vec![ (0,true), (1,true), (2,true), (3,true) ], // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* minor keys to matched major keys
//!                                                                 BooleanFieldOperator::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! // the remainder is nonzero
//! assert_eq!( 
//!     division.clone().remainder().collect_vec(),
//!     vec![ (3,true) ]
//! );
//!
//! // therefore the solution equals None
//! assert_eq!( 
//!     division.clone().solution(), 
//!     None
//! );                                     
//! ```
//! 
//! 
//! # Panics
//! 
//! 
//! The solver will panic if we pass an unsorted input, and pull all entries from the quotient and remainder
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified};        
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanFieldOperator::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                         vec![ 
//!                             vec![                                     ], 
//!                             vec![                                     ],                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );  
//! 
//! let d = vec![ (1,true), (0,true), (3,true), (2,true), ];
//! let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
//!                                                                 d.clone(), // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* minor keys to matched major keys
//!                                                                 BooleanFieldOperator::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );                                                                
//! assert_panic!( for _ in division.solver() {} );   // the elimination procedure clears the first three columns, but cannot clearn the fourth
//! ```
//! 
//! The solver *may possibly* compute a quotient without throwing an error,
//! if it finds an entry it can't eliminate before hitting the first consecutive out-of-order pair
//! 
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
//! use oat_rust::algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified};        
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanFieldOperator::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                         vec![ 
//!                             vec![                                     ], 
//!                             vec![                                     ],                                     
//!                             vec![ (0,true),  (1,true), (2,true)       ], 
//!                             vec![            (1,true), (2,true)       ],
//!                             vec![                      (2,true)       ],                                       
//!                         ],
//!                     );  
//! 
//! 
//! let d = vec![ (0,true), (1,true), (2,true), (3,true), (4,true), (2,true), ];
//! let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
//!                                                                 d.clone(), // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* minor keys to matched major keys
//!                                                                 BooleanFieldOperator::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! let mut qr = division.quotient_remainder();
//! 
//! // no panic when computing the quotient
//! assert_eq!( 
//!     qr.quotient.clone(), 
//!     vec![ (2,true) ] 
//! ); 
//! 
//! // no panic when computing the first element of the remainder 
//! assert_eq!( 
//!     qr.remainder.next(), 
//!     Some( (3,true) ) 
//! ); 
//! 
//! // panics once we exhaust the problem vector
//! assert_panic!( 
//!     for _ in qr.remainder {} 
//! );   
//! ```



//  ===========================================================================================
//  ===========================================================================================
//  ===========================================================================================
//  ===========================================================================================



use derive_getters::Dissolve;
use derive_new::new;


use crate::algebra::matrices::types::transpose::AntiTranspose;
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing, };
use crate::algebra::vectors::entries::{KeyValSet, KeyValGet};
use crate::algebra::vectors::operations::{VectorOperations, LinearCombinationSimplified};

use crate::utilities::iterators::merge::hit::{hit_bulk_insert, hit_merge_by_predicate};
use crate::utilities::order::{JudgePartialOrder, ReverseOrder};
use crate::utilities::functions::evaluate::EvaluateFunction;
use crate::utilities::iterators::general::{IterTwoType, HeadTail, RequireStrictAscentWithPanic, TransformIter};
// use debugit::DebugIt as D;


/// Standardizes the method of extracting a remainder from a solver
/// 
/// See the documentaiton for [DivisionWithRemainder], for an explanation.
pub trait IntoRemainder{
    type Remainder;
    fn into_remainder( self ) -> Self::Remainder;
}

/// Wrapper for the quotient and remainder of a division b/A
/// 
/// See the documentaiton for [DivisionWithRemainder], for an explanation of terms.
#[derive(Clone,Copy,Debug)]
pub struct QuotientRemainder< Q, R >{
    pub quotient: Q,
    pub remainder: R,
}

// /// Splits a solver into a quotient and a solution
// /// 
// /// See [QuotientRemainder] for the full interpretation.
// pub trait IntoQuotientRemainder< Remainder >
//     where 
//         Self:       Iterator + IntoRemainder,
// {
//     /// Splits the solution into a quotient and a remainder
//     /// 
//     /// Operates by draining the iterator 
//     /// 
//     /// See [QuotientRemainder] for the full interpretation.
//     fn into_quotient_remainder( self ) -> QuotientRemainder< Vec< Self::Item >, Remainder > {
//         // let mut q = Vec::new();
//         // while let Some(i) = self.next() { q.push(i) }
//         // q.shrink_to_fit();
//         // return (q, self)
//     }
// }


/// Safe wrapper for iterators that solve linear systems, with quotient and remainder parts
/// 
/// Echelon solvers use an incremental approach that involves iterative
/// elimination of leading entries of a vector `r` by adding scalar multiples of rows or columns from a matrix `A`.  The solver terminates when either
/// 
/// - no nonzero entries remain, OR
/// - the solver cannot eleiminate a leading entry
/// 
/// We call the linear coefficients used to combine rows/columns of `A` the **quotient**, and we call whatever remains after the elimination the **remainder**
/// 
/// A solution exists iff the remainder equals zero; in this case the quotient is a solution.
/// 
/// This struct provides methods for extracting quotients, remainders, and solutions from a solver.

#[derive(Copy,Debug,new)]
pub struct DivisionWithRemainder< I >
    where I:    Iterator + IntoRemainder,
{
    solver:  I
}

impl < I > DivisionWithRemainder< I >
    where 
        I:              Iterator + IntoRemainder,
        I::Remainder:   Iterator,

{
    /// Splits apart the quotient and remainder of an echelon solver
    /// 
    /// The quotient is a bona fide solution to `Ax = b` or `xA = b` iff the remainder is empty.
    /// 
    /// See [DivisionWithRemainder] for details.
    pub fn quotient_remainder( mut self ) -> QuotientRemainder< Vec<I::Item>, I::Remainder > {
        let mut quotient = Vec::new();
        for i in self.solver.by_ref() { quotient.push(i) }
        quotient.shrink_to_fit();
        let remainder = self.solver.into_remainder();
        QuotientRemainder{ quotient, remainder }       
    }

    /// Returns the remainder
    /// 
    /// See [DivisionWithRemainder] for details.
    pub fn remainder( mut self ) -> I::Remainder {
        for _ in self.solver.by_ref() {} // drain the quotient
        self.solver.into_remainder() // what remains is the remainder
    }    

    /// Returns the quotient, which is a solution iff the remainder is zero
    /// 
    /// See [DivisionWithRemainder] for furhter details.
    /// 
    /// This is equivalent to calling `self.solver()`, because the solver iterates over elements of the quotient.
    pub fn quotient( self ) -> I { self.solver }

    /// Returns the solver, which iterates over elements of the quotient
    /// 
    /// See [DivisionWithRemainder] for furhter details.
    pub fn solver( self ) -> I { self.solver }    

    /// Returns a solution formated as `Vec<T>`, if a solution exists
    /// 
    /// Concretely, returns the quotient of `b/A` if the remainder is zero, and returns `None` otherwise.
    /// 
    /// Alternate methods
    /// 
    /// - [DivisionWithRemainder::quotient] provides an iterator for the quotient; this can be faster and more memory efficient, however **the quotient is only a bona fide solution when the remainder is zero**
    /// - [DivisionWithRemainder::quotient_remainder], provides the same information, plus information about the remainder.  However, take care to remember that the quotient returned by that method is only a bona fide solution if the remainder is empty.
    /// 
    /// See [DivisionWithRemainder] for furhter details.
    pub fn solution( self ) -> Option< Vec< I::Item > > {   
        let mut qr = self.quotient_remainder();
        match qr.remainder.next().is_none() {
            true    =>  Some( qr.quotient ),
            false   =>  None,
        }
    }    
}


impl < I >

    Clone for 
    
    DivisionWithRemainder< I >

    where I:    Clone + Iterator + IntoRemainder,
{
    fn clone(&self) -> Self {
        DivisionWithRemainder { solver: self.solver.clone() }
    }
}


//  ======================================================================================
//  ======================================================================================
//  ======================================================================================
//  ASCENDING MAJOR SOLVE
//  ======================================================================================
//  ======================================================================================
//  ======================================================================================



//  ======================================================================================
//  ECHELON SOLVE ASCEND WITH MINOR KEYS
//  ======================================================================================




// ECHELON MAJOR ASCEND SOLVE W/ MINOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------

/// Iterates over the entries of the solution `x` to a matrix equation `xA = b`, where
/// (i) each entry of `x` representing a pair of form `(major_key, coefficient)` is replaced by an entry representing `(match(major_key), coefficient)`, where
/// `match(major_key)` is the matched minor key, and
/// (ii) entries appear in ascending order, according to minor index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the minor keys of `A`, and the entries of `b` appear in ascending order of index
/// - There is a partial bijection `match: keymin -> keymaj` from the minor keys of `A` to the major keys of `A`; concretely, this means a bijection from a subset of the minor keys to a subset of the major keys
/// - Whenever `keymaj = match( keymin )`, the leading entry of `A.view_major_ascend( keymaj )` has index `keymin`
/// - A solution to `xA = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymin, -alpha) )`.
#[derive(Dissolve)]
pub struct EchelonSolverMajorAscendWithMinorKeys< 
                    ProblemVector,
                    MatchingFromKeyMinToKeyMaj,
                    Matrix,                   
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone + PartialEq,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMajorAscend:        IntoIterator,        
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,                   

{   
    ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
    matrix:                             Matrix, // the `A` in `Ax = b`
    entries_to_eliminate:               HeadTail<
                                                LinearCombinationSimplified<   
                                                        IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                                                < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                                                RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >,
                                                                // Matrix::EntryMajor,
                                                            >,
                                                        Matrix::ColIndex, Matrix::Coefficient, RingOperator, OrderOperator 
                                                    >
                                            >,
    matching_from_keymin_to_keymaj:     MatchingFromKeyMinToKeyMaj,
    // phantom_lifetime:                   PhantomData< Matrix::RowIndex > 
} 




impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        EchelonSolverMajorAscendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone + PartialEq,
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                         Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,  

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolverMajorAscendWithMinorKeys`].
    pub fn solve(
                b:                              ProblemVector,                
                a:                              Matrix,
                matching_from_keymin_to_keymaj: MatchingFromKeyMinToKeyMaj,
                ring_operator:                  RingOperator,
                order_operator:                 OrderOperator,
            )
            ->
            DivisionWithRemainder<
                    EchelonSolverMajorAscendWithMinorKeys<
                            ProblemVector,
                            MatchingFromKeyMinToKeyMaj,
                            Matrix,                    
                            RingOperator,
                            OrderOperator,
                        >
                >
    {
        // we will eliminate entries by adding other vectors; addition is performed by merging in new iterators;
        // because the vectors we merge in may have a different type than the problem vector b, we need to do a
        // "type union" by wrapping `b` inside an `IterTwoType`
        let entries_to_eliminate
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    IterTwoType::Iter2( b.into_iter().require_strict_ascent_with_panic( order_operator.clone() ) )   // wrap b in an IterTwoType enum so that it can be easily combined with  
                                        .scale( ring_operator.minus_one(), ring_operator.clone() )
                                ),
                            order_operator,
                        )
                        .simplify( ring_operator.clone() );
        // load the entries in the tail of a HeadTail
        let entries_to_eliminate = HeadTail { head: None, tail: entries_to_eliminate };
        DivisionWithRemainder::new(
            EchelonSolverMajorAscendWithMinorKeys{
                    ring_operator,
                    matching_from_keymin_to_keymaj,
                    matrix:                             a,
                    entries_to_eliminate, 
                    // phantom_lifetime:                   PhantomData,
                }
        )
    }
  
}



impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        IntoRemainder for

        EchelonSolverMajorAscendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone + PartialEq,
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                         Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,  


{
    type Remainder = 
            HeadTail<
                    LinearCombinationSimplified<   
                            IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                    < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                    RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >, // the vector we have tried to eliminate
                                >,
                            Matrix::ColIndex, Matrix::Coefficient, RingOperator, OrderOperator 
                        >        
                >;


        /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
        /// in the process.
        /// 
        /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
        /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
        /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
        /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
        /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder {
        self.entries_to_eliminate
    }

}




impl    <
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        EchelonSolverMajorAscendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:               Clone + PartialEq,
        Matrix::ViewMajorAscend:        IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:            Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,  

{
    type Item = Matrix::EntryMajor;

    fn next( &mut self ) -> Option<Self::Item> { 

        // THIS ITERATOR PROCEEDS BY ELIMINATING ENTRIES IN `self.entries_to_eliminate`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_eliminate` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( mut entry_to_eliminate )       =   self.entries_to_eliminate.next() {

            // println!("EchelonSolverMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: STARTING");

            // KEY TO LOOK UP THE VECTOR THAT WILL EXECUTE THE ELIMINATION
            // let keymaj_of_eliminating_viewmaj = self.matching_from_keymin_to_keymaj.evaluate_function( entry_to_eliminate.key() );

            match self.matching_from_keymin_to_keymaj.evaluate_function( entry_to_eliminate.key() ) {
                None =>  {
                    // println!("!!!!!!!!!! EXITING BRANCH WITH REMAINDER");
                    // put back the entry
                    self.entries_to_eliminate.head = Some( entry_to_eliminate );
                    None
                } Some( keymaj_of_eliminating_viewmaj ) => {
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: WAYPOINT 1");            

                    // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
                    let mut seed_of_eliminating_iterator    =   self.matrix.view_major_ascend( keymaj_of_eliminating_viewmaj ).into_iter();
                    
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: WAYPOINT 2");            

                    // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
                    let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

                    // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
                    let scale_factor     =    self.ring_operator.negate(
                                                            self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                        );
                    
                    // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
                    // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate
                    let eliminating_iterator       
                            =   IterTwoType::Iter1( seed_of_eliminating_iterator ) 
                                    .scale( scale_factor.clone(), self.ring_operator.clone() );

                    // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate`
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): hit_bulk_insert: STARTING");
                    hit_bulk_insert( &mut self.entries_to_eliminate.tail.unsimplified, vec![eliminating_iterator] ); 
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): hit_bulk_insert: ENDING");

                    // NOW RE-USE `entry_to_eliminate` AS A CONTAINER FOR THE NEXT ENTRY IN THE INVERSE MATRIX
                    entry_to_eliminate.set_val( scale_factor );

                    Some( entry_to_eliminate ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function
                }
            }

        } else {           
            None
        }
    }
}  



impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Clone for    

        EchelonSolverMajorAscendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        ProblemVector::IntoIter:        Clone,
        MatchingFromKeyMinToKeyMaj:     Clone + EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         Clone + ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone + PartialEq,
        Matrix::ViewMajorAscend:                IntoIterator,
        Matrix::ViewMajorAscendIntoIter: Clone,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                         Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,              
{
    fn clone(&self) -> Self {
        EchelonSolverMajorAscendWithMinorKeys{ 
                ring_operator: self.ring_operator.clone(),
                matrix: self.matrix.clone(),
                entries_to_eliminate: self.entries_to_eliminate.clone(),
                matching_from_keymin_to_keymaj: self.matching_from_keymin_to_keymaj.clone(),
            }
    }
}



















//  ======================================================================================
//  ECHELON MAJOR ASCEND SOLVE ASCEND WITH MAJOR KEYS
//  ======================================================================================




// ECHELON SOLVE W/ MAJOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------

/// Iterates over the entries of the solution `x` to a matrix equation `xA = b`, where entries appear in ascending order, according to (major) index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the minor keys of `A`, and the entries of `b` appear in ascending order of index
/// - There is a partial bijection `match: keymin -> keymaj` from the minor keys of `A` to the major keys of `A`; concretely, this means a bijection from a subset of the minor keys to a subset of the major keys
/// - Whenever `keymaj = match( keymin )`, the leading entry of `A.view_major_ascend( keymaj )` has index `keymin`
/// - A solution to `xA = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymaj, -alpha) )`.
#[derive(Dissolve)]
pub struct EchelonSolverMajorAscendWithMajorKeys<
                    ProblemVector,
                    MatchingFromKeyMinToKeyMaj,
                    Matrix,                    
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                 Clone + PartialEq,
        Matrix::Coefficient:                 Clone,   
        Matrix::ViewMajorAscend:        IntoIterator,        
        Matrix::EntryMajor:         Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                Clone + JudgePartialOrder <  Matrix::EntryMajor >,                   

{   
    ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
    matrix:                             Matrix, // the `A` in `Ax = b`
    entries_to_eliminate:               HeadTail<
                                                LinearCombinationSimplified<   
                                                        IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                                                < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                                                RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >,
                                                                // Matrix::EntryMajor,
                                                            >,
                                                        Matrix::ColIndex, Matrix::Coefficient, RingOperator, OrderOperator 
                                                    >
                                            >,
    matching_from_keymin_to_keymaj:     MatchingFromKeyMinToKeyMaj,
} 


impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        EchelonSolverMajorAscendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:               Clone + PartialEq,
        Matrix::RowIndex:               Clone,
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                         Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,  

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolverMajorAscendWithMajorKeys`].
    pub fn solve(
                b:                              ProblemVector,                
                a:                              Matrix,
                matching_from_keymin_to_keymaj: MatchingFromKeyMinToKeyMaj,
                ring_operator:                  RingOperator,
                order_operator:               OrderOperator,
            )
            ->
            DivisionWithRemainder<
                    EchelonSolverMajorAscendWithMajorKeys<
                            ProblemVector,
                            MatchingFromKeyMinToKeyMaj,
                            Matrix,                    
                            RingOperator,
                            OrderOperator,
                        >  
                >
    {
        // we will eliminate entries by adding other vectors; addition is performed by merging in new iterators;
        // because the vectors we merge in may have a different type than the problem vector b, we need to do a
        // "type union" by wrapping `b` inside an `IterTwoType`
        let entries_to_eliminate
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    IterTwoType::Iter2( b.into_iter().require_strict_ascent_with_panic(order_operator.clone()) )   // wrap b in an IterTwoType enum so that it can be easily combined with  
                                        .scale( ring_operator.minus_one(), ring_operator.clone() ) 
                                ),
                            order_operator,
                        )
                        .simplify( ring_operator.clone() );
        // load the entries in the tail of a HeadTail
        let entries_to_eliminate = HeadTail { head: None, tail: entries_to_eliminate }; 
        DivisionWithRemainder::new(
            EchelonSolverMajorAscendWithMajorKeys{
                ring_operator,
                matching_from_keymin_to_keymaj,
                matrix:                             a,
                entries_to_eliminate, 
            } 
        )
    }
}





impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        IntoRemainder for

        EchelonSolverMajorAscendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone + PartialEq,
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                         Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,  

{

    type Remainder = 
        HeadTail<
            LinearCombinationSimplified<   
                    IterTwoType< // this enum allows us to treat two different iterator types as a single type
                            < Matrix::ViewMajorAscend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                            RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >, // the vector we have tried to eliminate
                        >,
                    Matrix::ColIndex, Matrix::Coefficient, RingOperator, OrderOperator 
                >
            >;


    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder
    {
        self.entries_to_eliminate
    }   
}








impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        EchelonSolverMajorAscendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromKeyMinToKeyMaj:     EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone + PartialEq,
        Matrix::RowIndex:                         Clone, // this is necessary to avoid a "move" error, since the following function uses a major key once to look up a major view (which consums the keymaj) and once to return the value of the keymaj to the user
        Matrix::ViewMajorAscend:                IntoIterator,
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                         Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:             Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,  

{
    type Item = (Matrix::RowIndex, Matrix::Coefficient);

    fn next( &mut self ) -> Option<Self::Item> { 

    // unsafe {
    //     println!("about to invoke EchelonSolverMajorAscendWithMajorKeys.nex(); problem vector =  {:?}", D(&self.entries_to_eliminate) );
    //     std::thread::sleep(std::time::Duration::from_secs(1));

    // }
    

    // THIS ITERATOR PROCEED BY ELIMINATING ENTRIES IN `self.entries_to_eliminate`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_eliminate` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( entry_to_eliminate )       =   self.entries_to_eliminate.next() {

            // unsafe {
            //     println!("entry to eliminate =  {:?}", D(&entry_to_eliminate) );
            // }
                            

            // INDEX OF THE VECTOR THAT WILL EXECUTE THE ELIMINATION
            // let keymaj_of_eliminating_viewmaj = self.matching_from_keymin_to_keymaj.evaluate_function( entry_to_eliminate.key() );

            match self.matching_from_keymin_to_keymaj.evaluate_function( entry_to_eliminate.key() ) {
                None => {
                    self.entries_to_eliminate.head = Some( entry_to_eliminate ); // put back the entry; we cannot clear it
                    None
                }, Some( keymaj_of_eliminating_viewmaj ) => {
                    // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
                    let mut seed_of_eliminating_iterator    =   self.matrix.view_major_ascend( keymaj_of_eliminating_viewmaj.clone() ).into_iter();
                    
                    // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
                    let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

                    // unsafe {
                    //     println!("keymaj_of_eliminating_viewmaj =  {:?}", D(&keymaj_of_eliminating_viewmaj) );                
                    //     println!("eliminating entry =  {:?}", D(&eliminating_entry) );
                    //     let seed_of_eliminating_iterator_2    =   self.matrix.view_major_ascend( keymaj_of_eliminating_viewmaj.clone() ).into_iter().collect_vec();
                    //     println!("eliminating iterator entries =  {:?}", D(&seed_of_eliminating_iterator_2) );
                    // }            

                    // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
                    let scale_factor     =    self.ring_operator.negate(
                                                            self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                        );
                    
                    // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
                    // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate
                    let eliminating_iterator       
                            =   IterTwoType::Iter1( seed_of_eliminating_iterator ) 
                                    .scale( scale_factor.clone(), self.ring_operator.clone() );


                    // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate`
                    hit_bulk_insert( &mut self.entries_to_eliminate.tail.unsimplified, vec![eliminating_iterator] ); 

                    Some( (keymaj_of_eliminating_viewmaj, scale_factor) ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function
                }
            }

        } else {
            None
        }
    }
}  






impl    < 
            ProblemVector,
            MatchingFromKeyMinToKeyMaj,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Clone for    

        EchelonSolverMajorAscendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMinToKeyMaj,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                      IntoIterator< Item = Matrix::EntryMajor >, // the iterator runs over entries of the same type as the matrix        
        ProblemVector::IntoIter:            Clone,
        MatchingFromKeyMinToKeyMaj:         Clone + EvaluateFunction< Matrix::ColIndex, Option< Matrix::RowIndex > >,        
        Matrix:                             Clone + ViewRowAscend + IndicesAndCoefficients,
        Matrix::ColIndex:                   Clone + PartialEq,
        Matrix::RowIndex:                   Clone, // this is necessary to avoid a "move" error, since the following function uses a major key once to look up a major view (which consums the keymaj) and once to return the value of the keymaj to the user
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::ViewMajorAscendIntoIter:    Clone,
        RingOperator:                       Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        Matrix::Coefficient:                Clone,
        OrderOperator:                      Clone + JudgePartialOrder <  Matrix::EntryMajor >,               
        Matrix::EntryMajor:                 Clone + KeyValGet < Matrix::ColIndex, Matrix::Coefficient > + KeyValSet < Matrix::ColIndex, Matrix::Coefficient >,              
{
    fn clone(&self) -> Self {
        EchelonSolverMajorAscendWithMajorKeys{ 
                ring_operator: self.ring_operator.clone(),
                matrix: self.matrix.clone(),
                entries_to_eliminate: self.entries_to_eliminate.clone(),
                matching_from_keymin_to_keymaj: self.matching_from_keymin_to_keymaj.clone(),
            }
    }
}










//  ======================================================================================
//  ======================================================================================
//  ======================================================================================
//  DESCENDING MINOR SOLVE
//  ======================================================================================
//  ======================================================================================
//  ======================================================================================







//  ======================================================================================
//  ECHELON MINOR DESCEND SOLVE WITH MAJOR KEYS
//  ======================================================================================




// ECHELON SOLVE W/ MAJOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------


/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, where
/// (i) each entry of `x`, which represents a pair of form `(minor_key, coefficient)`, is replaced by an entry representing `(match(minor_key), coefficient)`, where
/// `match(minor_key)` is the matched major key, and
/// (ii) entries appear in descending order, according to major index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the major keys of `A`, and the entries of `b` appear in descending order of index
/// - There is a partial bijection `match: keymaj -> keymin` from the major keys of `A` to the minor keys of `A`; concretely, this means a bijection from a subset of the major keys to a subset of the minor keys
/// - Whenever `keymin = match( keymaj )`, the leading entry of `A.view_minor_descend( keymin )` has index `keymaj`
/// - A solution to `Ax = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymaj, -alpha) )`.
/// 
/// Under the hood, `EchelonSolverMinorDescendWithMajorKeys` is simply a wrapper for a struct `EchelonSolverMajorAscendWithMinorKeys` which
/// contains a reference to the antitranspose of the matrix; see source code for full details.
#[derive(Dissolve)]
pub struct EchelonSolverMinorDescendWithMajorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:               Clone + PartialEq,
        Matrix::Coefficient:                 Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,                   

{   
    antitransposed_solver:      EchelonSolverMajorAscendWithMinorKeys< // the minor keys of the antitransposed matrix are the major keys of the un-antitransposed matrix
                                        ProblemVector,
                                        MatchingFromKeyMajToKeyMin,
                                        AntiTranspose< Matrix >,
                                        RingOperator,
                                        ReverseOrder< OrderOperator >,
                                    >,
} 

impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperatorUnreversed,
        >  

        EchelonSolverMinorDescendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderOperatorUnreversed,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:                 Clone + PartialEq,
        Matrix::Coefficient:                 Clone,   
        Matrix::ViewMinorDescend:       IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperatorUnreversed:        Clone + JudgePartialOrder <  Matrix::EntryMinor >,   

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolverMajorAscendWithMinorKeys`].
    pub fn solve(
                b:                              ProblemVector,                
                a:                              Matrix,
                matching_from_keymaj_to_keymin: MatchingFromKeyMajToKeyMin,
                ring_operator:                  RingOperator,
                order_operator_unreversed:    OrderOperatorUnreversed,
            )
            ->
            DivisionWithRemainder<
                    EchelonSolverMinorDescendWithMajorKeys< 
                            ProblemVector,
                            MatchingFromKeyMajToKeyMin,
                            Matrix,                     
                            RingOperator,
                            OrderOperatorUnreversed,
                        >             
                > 
    {
        DivisionWithRemainder::new(
            EchelonSolverMinorDescendWithMajorKeys{
                    // Note that we use EchelonSolverMajorAscendWithMinorKeys, since the minor keys of the antitranspose are the major keys of the original matrix
                    antitransposed_solver:   EchelonSolverMajorAscendWithMinorKeys::solve(
                                                    b,
                                                    AntiTranspose::new( a ),
                                                    matching_from_keymaj_to_keymin,
                                                    ring_operator,
                                                    ReverseOrder::new( order_operator_unreversed ),
                                                ).quotient()
                }
        )
    }

}








impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperatorUnreversed,
        >  

        IntoRemainder for

        EchelonSolverMinorDescendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderOperatorUnreversed,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:                 Clone + PartialEq,
        Matrix::Coefficient:                 Clone,   
        Matrix::ViewMinorDescend:       IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperatorUnreversed:        Clone + JudgePartialOrder <  Matrix::EntryMinor >,   

{

    type Remainder = 
        HeadTail<
                LinearCombinationSimplified<   
                        IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                < Matrix::ViewMinorDescend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                RequireStrictAscentWithPanic< ProblemVector::IntoIter, ReverseOrder< OrderOperatorUnreversed > >, // the vector we have tried to eliminate
                            >,
                        Matrix::RowIndex, Matrix::Coefficient, RingOperator, ReverseOrder< OrderOperatorUnreversed >
                    >
            >;

    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder                 
    {
        self.antitransposed_solver.into_remainder()
    }    
}









impl    <
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        EchelonSolverMinorDescendWithMajorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:                         Clone + PartialEq,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,   

{
    type Item = Matrix::EntryMinor;

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitransposed_solver.next()
    }
}  






impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Clone for    

EchelonSolverMinorDescendWithMajorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        ProblemVector::IntoIter:        Clone,
        MatchingFromKeyMajToKeyMin:     Clone + EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         Clone + ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:                         Clone + PartialEq,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMinorDescend:       IntoIterator,        
        Matrix::ViewMinorDescendIntoIter:   Clone,
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,                    
{
    fn clone(&self) -> Self {
        EchelonSolverMinorDescendWithMajorKeys{ 
                antitransposed_solver: self.antitransposed_solver.clone()
            }
    }
}
















//  ======================================================================================
//  ECHELON MINOR DESCEND SOLVE WITH MINOR KEYS
//  ======================================================================================




// ECHELON SOLVE W/ MINOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------


/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, where entries appear in descending order, according to (minor) index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the major keys of `A`, and the entries of `b` appear in descending order of index
/// - There is a partial bijection `match: keymaj -> keymin` from the major keys of `A` to the minor keys of `A`; concretely, this means a bijection from a subset of the major keys to a subset of the minor keys
/// - Whenever `keymin = match( keymaj )`, the leading entry of `A.view_minor_descend( keymin )` has index `keymaj`
/// - A solution to `Ax = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `keymin` be the index of this entry.  (iii) if `keymin` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `keymaj` be the major key matched to `keymin`.  Add a scalar multiple of `A.view_major_ascend( keymaj )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (keymaj, -alpha) )`.
/// 
/// Under the hood, `EchelonSolverMinorDescendWithMajorKeys` is simply a wrapper for a struct `EchelonSolverMajorAscendWithMinorKeys` which
/// contains a reference to the antitranspose of the matrix; see source code for full details.
pub struct EchelonSolverMinorDescendWithMinorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:                         Clone + PartialEq,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,                   

{   
    antitransposed_solver:      EchelonSolverMajorAscendWithMajorKeys< // the major keys of the antitransposed matrix are the minor keys of the un-antitransposed matrix
                                        ProblemVector,
                                        MatchingFromKeyMajToKeyMin,
                                        AntiTranspose< Matrix >,                 
                                        RingOperator,
                                        ReverseOrder< OrderOperator >,
                                    >,
} 

impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        EchelonSolverMinorDescendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:                         Clone + PartialEq,
        Matrix::ColIndex:               Clone,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,   

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolverMajorAscendWithMinorKeys`].
    pub fn solve(
                b:                              ProblemVector,                
                a:                              Matrix,
                matching_from_keymaj_to_keymin: MatchingFromKeyMajToKeyMin,
                ring_operator:                  RingOperator,
                order_operator:               OrderOperator,
            )
            ->
            DivisionWithRemainder<
                    EchelonSolverMinorDescendWithMinorKeys< 
                            ProblemVector,
                            MatchingFromKeyMajToKeyMin,
                            Matrix,                     
                            RingOperator,
                            OrderOperator,
                        >             
                > 
    {
        DivisionWithRemainder::new(
            EchelonSolverMinorDescendWithMinorKeys{
                // Note that we use EchelonSolverMajorAscendWithMajorKeys, since the major keys of the antitranspose are the minor keys of the original matrix
                antitransposed_solver:   EchelonSolverMajorAscendWithMajorKeys::solve(
                                                b,
                                                AntiTranspose::new( a ),
                                                matching_from_keymaj_to_keymin,
                                                ring_operator,
                                                ReverseOrder::new( order_operator ),
                                            ).quotient()
            }
        )
    }
}





impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        IntoRemainder for

        EchelonSolverMinorDescendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::RowIndex:                         Clone + PartialEq,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,   

{

    type Remainder = 
            HeadTail<
                    LinearCombinationSimplified<   
                            IterTwoType< // this enum allows us to treat two different iterator types as a single type
                                    < Matrix::ViewMinorDescend as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                    RequireStrictAscentWithPanic< ProblemVector::IntoIter, ReverseOrder< OrderOperator > >, // the vector we have tried to eliminate
                                >,
                            Matrix::RowIndex, Matrix::Coefficient, RingOperator, ReverseOrder< OrderOperator >
                        >
                >;

    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder     
    {
        self.antitransposed_solver.into_remainder()
    }    
}





impl    <
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        EchelonSolverMinorDescendWithMinorKeys<
                ProblemVector,
                MatchingFromKeyMajToKeyMin,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        MatchingFromKeyMajToKeyMin:     EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         ViewColDescend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone,
        Matrix::RowIndex:                         Clone + PartialEq,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMinorDescend:               IntoIterator,        
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,   

{
    type Item = (Matrix::ColIndex, Matrix::Coefficient);

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitransposed_solver.next()
    }
}  





impl    < 
            ProblemVector,
            MatchingFromKeyMajToKeyMin,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Clone for    

EchelonSolverMinorDescendWithMinorKeys<
                    ProblemVector,
                    MatchingFromKeyMajToKeyMin,
                    Matrix,                     
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                  IntoIterator< Item = Matrix::EntryMinor >, // the iterator runs over entries of the same type as the matrix
        ProblemVector::IntoIter:        Clone,
        MatchingFromKeyMajToKeyMin:     Clone + EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColIndex > >,
        Matrix:                         Clone + ViewColDescend + IndicesAndCoefficients,
        Matrix::ColIndex:                         Clone,
        Matrix::RowIndex:                         Clone + PartialEq,
        Matrix::Coefficient:                         Clone,   
        Matrix::ViewMinorDescend:       IntoIterator,        
        Matrix::ViewMinorDescendIntoIter:   Clone,
        Matrix::EntryMinor:             Clone + KeyValGet < Matrix::RowIndex, Matrix::Coefficient > + KeyValSet < Matrix::RowIndex, Matrix::Coefficient >,      
        RingOperator:                   Clone + Semiring< Matrix::Coefficient > + Ring< Matrix::Coefficient > + DivisionRing< Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::EntryMinor >,                     
{
    fn clone(&self) -> Self {
        EchelonSolverMinorDescendWithMinorKeys{ 
                antitransposed_solver: self.antitransposed_solver.clone()
            }
    }
}












//  =========================================================================
//  =========================================================================
//  =========================================================================
//  IMPLEMENTATIONS OF CLONE
//  =========================================================================
//  =========================================================================
//  =========================================================================












































//  TESTS
//  ===========================================================================


#[cfg(test)]
mod doctring_tests {



    #[test]
    fn doc_test_module() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys, EchelonSolverMinorDescendWithMinorKeys};
                
        use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
        use crate::utilities::order::{OrderOperatorByKey};
        use itertools::Itertools;
        use assert_panic::assert_panic;

        // define the ring operator
        let ring_operator = BooleanFieldOperator::new();

        // build a matrix A with an invertible upper triangular submatrix
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![                                     ], 
                                    vec![                                     ],                                     
                                    vec![ (0,true),  (1,true), (2,true)       ], 
                                    vec![            (1,true), (2,true)       ],
                                    vec![                      (2,true)       ],                                       
                                ],
                            );    

        // ======================================================================================                        
        // SOLVE xA = b FOR x
        // ======================================================================================
                        
        let b = vec![ (0,true), (1,true) ];   // define a sparse vector b           

        // define the key-mapping with a wrapper around a closure operator (|x| -> expression(x))      
        // --------------------------------------------------------------------------------------
                        
        let partial_bijection = |x: usize| { if x < 3 { Some(x+2) } else { None } };
        let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
                                                                        b.clone(), // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* minor keys to matched major keys
                                                                        BooleanFieldOperator::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );
                                                                
        let product = division.quotient().multiply_matrix_major_ascend( &matrix, ring_operator, OrderOperatorByKey::new() );
        assert_eq!( product.collect_vec(), b );   // check the solution, i.e. check that xA = b  


        // define the key-mapping with vector
        // --------------------------------------------------------------------------------------   
                        
        // create a solver to solve xA = b for x
        let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
                                                                        b.clone(), // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* minor keys to matched major keys
                                                                        BooleanFieldOperator::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );
                                                                
        // check the solution, i.e. check that xA = b
        let product = division.clone().quotient().multiply_matrix_major_ascend( &matrix, ring_operator, OrderOperatorByKey::new() );
        assert_eq!( product.collect_vec(), b );  

        // different look-up methods
        // --------------------------------------------------------------------------------------             

        assert_eq!( division.clone().solution(), Some( division.clone().quotient().collect_vec() ) );        

        // ======================================================================================
        // SOLVE Ay = c FOR y
        // ======================================================================================

        let c = vec![ (3,true), (2,true) ];   // define a sparse vector c  

        // define the key-mapping with a wrapper around a closure operator (|x| -> expression(x))      
        // --------------------------------------------------------------------------------------        
                                                                                                                                
        let partial_bijection = |x: usize| { if (2..=4).contains(&x) { Some(x-2) } else { None } };
        let division =  EchelonSolverMinorDescendWithMinorKeys::solve(
                                                                        c.clone(), // entries must appear in strictly DESCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* major keys to matched minor keys
                                                                        BooleanFieldOperator::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
                                                                    );
                                                                
        let product = division.quotient().multiply_matrix_minor_descend( &matrix, ring_operator, OrderOperatorByKey::new());
        assert_eq!( product.collect_vec(), c );   // check the solution, i.e. check that Ay = c      

        // ======================================================================================                        
        // GET A REMAINDER
        // ======================================================================================


        // quotienting out columns
        // --------------------------------------------------------------------------------------   

        let d = vec![ (3,true), (2,true), (1,true), (0,true)  ];

        let partial_bijection = |x: usize| { if (2..=4).contains(&x) { Some(x-2) } else { None } };
        let division =  EchelonSolverMinorDescendWithMinorKeys::solve(
            d, // entries must appear in strictly DESCENDING order
            & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
            EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* major keys to matched minor keys
            BooleanFieldOperator::new(), // defines the ring operations
            OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
        );    
        let remainder = division.remainder().collect_vec();                                                        
        assert!( remainder.eq( &vec![ (1,true), (0,true) ] ) );   // the elimination procedure clears all but the first two entries

        // quotienting out rows
        // --------------------------------------------------------------------------------------   

        let d = vec![ (0,true), (1,true), (2,true), (3,true), ];
        let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
                                                                        d, // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* minor keys to matched major keys
                                                                        BooleanFieldOperator::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );                                                                
        assert!( division.clone().remainder().eq( vec![ (3,true) ] ) );   // the elimination procedure clears the first three columns, but cannot clearn the fourth

        // alternate look-up method
        // --------------------------------------------------------------------------------------   
                                                          
        assert!( division.clone().quotient_remainder().remainder.eq( vec![ (3,true) ] ) );   // the elimination procedure clears all but the first two entries        
         

        // ======================================================================================                        
        // PANIC FOR UNSORTED INPUTS
        // ======================================================================================

        // panics if we pass an unsorted input, and pull all entries from the quotient and remainder
        // --------------------------------------------------------------------------------------   

        let d = vec![ (1,true), (0,true), (3,true), (2,true), ];
        let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
                                                                        d, // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* minor keys to matched major keys
                                                                        BooleanFieldOperator::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );                                                                
        assert_panic!( for _ in division.solver() {} );   // the elimination procedure clears the first three columns, but cannot clearn the fourth


        // the solver *may possibly* compute a quotient without throwing an error, if it finds an entry it can't eliminate before hitting the first consecutive out-of-order pair
        // --------------------------------------------------------------------------------------   

        let d = vec![ (0,true), (1,true), (2,true), (3,true), (4,true), (2,true), ];
        let division =  EchelonSolverMajorAscendWithMajorKeys::solve(
                                                                        d, // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* minor keys to matched major keys
                                                                        BooleanFieldOperator::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );
        
        let mut qr = division.quotient_remainder();
        assert!( qr.quotient.clone().eq( &vec![ (2,true) ] ) ); // no panic when computing the quotient
        assert_eq!( qr.remainder.next(), Some( (3,true) ) ); // no panic when computing the first element of the remainder 
        assert_panic!( for _ in qr.remainder {} );   // panics once we exhaust the problem vector

    }

}










//  TESTS
//  ===========================================================================


#[cfg(test)]
mod tests {

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    
    use itertools::Itertools;
    use crate::{algebra::matrices::types::vec_of_vec::sorted::VecOfVec, algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator};


    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Verifies that the output of a solver actually solves the problem Ax = b
    /// **IT IS IMPORTANT THAT THE MATRIX BE IN PARTIAL ECHECLON FORM**
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_echelon_solve_on_invertible_uppertriangular_matrix< 'a, MatchingFromKeyMinToKeyMaj, MatchingFromKeyMajToKeyMin, Coefficient, RingOperator >( 
                matrix:                             & VecOfVec< usize, Coefficient >, 
                matching_from_keymin_to_keymaj:     MatchingFromKeyMinToKeyMaj,
                matching_from_keymaj_to_keymin:     MatchingFromKeyMajToKeyMin,
                ring_operator:                      RingOperator, 
                matrix_size:                        usize 
            ) 
        where   Coefficient:                         Clone + PartialEq + Ord + std::fmt::Debug,
                RingOperator:                   Semiring< Coefficient > + Ring< Coefficient > + DivisionRing< Coefficient > + Clone,
                MatchingFromKeyMinToKeyMaj:     Clone + EvaluateFunction< usize, Option< usize > >,
                MatchingFromKeyMajToKeyMin:     Clone + EvaluateFunction< usize, Option< usize > >,                
    {
        use crate::{algebra::matrices::operations::multiply::{vector_matrix_multiply_major_ascend_simplified, vector_matrix_multiply_minor_descend_simplified}, utilities::order::OrderOperatorAuto};

        


        // generate random vectors of 0's and 1's; try to get one vector each with 0, 1, 2, 3, and `matrix_size` nonzero entries in total
        for num_nonzeros in [0, 1, 2, 3, matrix_size].iter().map(|x| std::cmp::min(x, &matrix_size) ) {
            for vector_support in (0 .. matrix_size).combinations( *num_nonzeros ) {
                let mut vec = vector_support
                                .iter()
                                .cloned()
                                .map(|x| ( x, RingOperator::one() ))
                                .collect_vec();
                vec.sort();


                //  view-major solutions (corresponding to xA = b, where A is row-major)
                //  --------------------------------------------------------------------                

                // compute a view-major solution with minor keys
                let solution_major_with_minor_keys = EchelonSolverMajorAscendWithMinorKeys::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymin_to_keymaj.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with minor keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_major_ascend_simplified( 
                            solution_major_with_minor_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );

                // compute a view-major solution with major keys
                let solution_major_with_major_keys = EchelonSolverMajorAscendWithMajorKeys::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymin_to_keymaj.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with major keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_major_ascend_simplified( 
                            solution_major_with_major_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );    


                //  view-minor solutions (corresponding to Ax = b, where A is row-major)
                //  --------------------------------------------------------------------

                // REVERSE THE ORDER OF ENTRIES, SINCE WE ARE NOW DEALING WITH DESCENDING VIEWS
                vec.reverse(); 

                // compute a view-minor solution with minor keys
                let solution_minor_with_minor_keys = EchelonSolverMinorDescendWithMinorKeys::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymaj_to_keymin.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with minor keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_minor_descend_simplified( 
                            solution_minor_with_minor_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );

                // compute a view-minor solution with major keys
                let solution_major_with_major_keys = EchelonSolverMinorDescendWithMajorKeys::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_keymaj_to_keymin.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with major keys
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    vector_matrix_multiply_minor_descend_simplified( 
                            solution_major_with_major_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );                                
                

            }
        }                                                               
    }


    
    #[test]
    fn test_echelon_solve_on_specific_matrices() {
        use num::rational::Ratio;        
        // use crate::algebra::matrices::query::MajorDimension;
        use crate::algebra::rings::operator_structs::ring_native::DivisionRingNative;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;

        // Define the integer type for the rational numbers we'll use
        type IntegerType = isize;
        let modulus = 1049;

        // Define shorthand functions for creating new rational numbers
        let r = Ratio::new;    
        let q = Ratio::from_integer;   

        // Define matching function from Matrix::ColIndex (=usize) to Matrix::RowIndex (=usize)
        let matching_keymin_to_keymaj_closure = | x: usize | -> Option< usize > { Some( x ) };
        let matching_keymin_to_keymaj = EvaluateFunctionFnMutWrapper::new( matching_keymin_to_keymaj_closure );
        
        // Define the ring operators
        let ring_operator_q  =   DivisionRingNative::< Ratio<IntegerType> >::new();  
        let ring_operator_p  =   PrimeOrderFieldOperator::new(modulus);                  

        // Test individual matrices        
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)) ]
                                ],
        );
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q, 1 );       

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(1)), ], 
                                    vec![            (1,q(1)), ],
                                ],
        );
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q, 2 );       

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
        );     
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q, 3 );        
        
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
        );        
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q, 3 );           

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
        );     
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj.clone(), ring_operator_q, 4 );         
        
        // MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        use rand::Rng;        // we use this module to generate random elements
        let matrix_size =   20;

        let mut rng = rand::thread_rng(); // this generates random integers
        let mut vec_of_vec = vec![];
        for majkey in 0 .. matrix_size {
            let coefficient_leading         =   rng.gen_range( 1 .. modulus );
            let mut new_vec     =   vec![ (majkey, coefficient_leading) ]; // start with a nonzero diagonal element            
            for q in majkey+1 .. matrix_size { // fill out the rest of this row of the matrix
                let coefficient   = rng.gen_range( 0 .. modulus );
                let flag = rng.gen_range(0usize .. 3usize);
                if      flag == 0 { new_vec.push( ( q, 0 )           ) }
                else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
                else              { continue }
            }
            vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        }

        let matrix  =   VecOfVec::new(vec_of_vec); // formally wrap the matrix in a VecOfVec struct
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_keymin_to_keymaj.clone(), matching_keymin_to_keymaj, ring_operator_p, matrix_size );                 

    }

}
