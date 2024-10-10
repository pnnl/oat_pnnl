//! Computing persistent homology
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
//!
//! # Example
//! 
//! Under construction

// //! An example persistent homology calculation.
// //! 
// //! Here we compute a basis for persistent homology in dimension 1,
// //! for the simplicial complex on vertex set `{0,1,2,3}`
// //! whose maximal faces are `{0,1,2}, {0,3}, {1,3}, {2,3}`.
// //! 
// //! 
// //! ## Mathematical overview
// //! 
// //! We compute PH using a matrix facotrization technique called U-match factorization,
// //! which is closely related to other factorization methods in PH, such as `R=DV` and `RU = D`.  We apply
// //! this method to the boundary matrix `D`, whose `k`th row and column correspond to the `k`th simplex
// //! added to the filtration.  This matrix has a row/column for every simplex of every dimension.
// //! 
// //! A U-match factorization of `D` is an equation `RM = DC`, where `R` and `C` are square upper-unitriangular
// //! matrices and `M` is a generalized matching matrix
// //! (see the [U-match](oat_rust::matrices::operations::umatch) module for more details).   To compute the barcode in dimension `d` (and, optionally, a basis of cycle representatives)
// //! you need only two facts: 
// //! 
// //! * **Finite bars** The barcode has one interval of form `[birth(s), birth(t))` for every `d`-simplex `s` 
// //! and `d+1`-simplex `t` such that `M[s,t]` is nonzero.  Here `birth(s)` denotes the time when
// //! `s` enters the filtration.
// //! The corresponding basis vector is column `s` of the matrix `R`.  We call `s` and `t` birth
// //! and death simplices.  Usually we exclude empty intervals.
// //! 
// //! * **Infinite bars**  The barcode has one interval of form `[birth(s), inf)` for every `d`-simplex `s`
// //! such that `M[s,:]=0` and `M[:,s]=0`.  The corresponding basis vector is column `s` of `C`.
// //! 
// //! 
// //! 
// //! ## Compute the barcode
// //! 
// //! In this example, the complex is filtered simplex-wise (i.e. one simplex is added at a time),
// //! and simplices are added in sorted lexicographic order (by convention,
// //! simplices of lower dimension appear before simplices of higher dimension).
// //! The coefficient field is the finite field of order 3.
// //! 
// //! 
// //! To compute the barcode, enumerate birth/death pairs (see below),
// //! and calculate the corresponding intervals as per the mathematical overview.
// //! 
// //! ## Compute birth/death pairs and cycle representatives
// //! 
// //! The following code computes birth/death pairs anc cycle representatives.
// //! 
// //! To run this example on your desktop computer, first check out the
// //! [quick start tutorial in OAT]() for instructions
// //! on installing Rust and running a program.  As part of this process,
// //! you'll create a new folder that contains a file called `main.rs`.  Inside
// //! `main.rs` is some text that reads `fn main{ .. }`.  Delete everything
// //! between `{` and `}`, and paste in the following:
// //! 
// //!
// //! ```
// //! use oat_rust::topology::simplicial::from::relation::BoundaryMatrixDowker;      
// //! use oat_rust::algebra::matrices::query::ViewColDescend;
// //! use oat_rust::algebra::matrices::operations::umatch::row_major::{Umatch::factor_with_clearing};
// //! use oat_rust::utilities::order::OrderOperatorAutoLt;        
// //! use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
// //! 
// //! // Define the ring operator for the finite field of order 3.
// //! // You can use this object to perform arithmetic operations, e.g., 
// //! // to add 1 and 1, you can run `ring_operator.add(1,1)`.
// //! let ring_operator          =   PrimeOrderFieldOperator::new(3);
// //! 
// //! // We will build a dowker complex.
// //! // A dowker complex is defined by a vertex set V and a family S
// //! // of subsets of V.  A subset of V forms a simplex iff it is 
// //! // a subset of some element of S.  We refer to the elements 
// //! // of S as "dowker simplices".
// //! 
// //! // Define the maximum homology dimension we're interested in.
// //! let maxdim    =   2;
// //! 
// //! // Here is one storage format for the set S.
// //! let dowker_simplices =  
// //!         vec![ 
// //!                 vec![0,1,2],
// //!                 vec![0,3],                      
// //!                 vec![1,3], 
// //!                 vec![2,3]                                           
// //!             ];      
// //! 
// //! // Build the boundary matrix.
// //! // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
// //! let boundary_matrix = BoundaryMatrixDowker::from_vectors( dowker_simplices, ring_operator.clone() );
// //! 
// //! // This iterates over simplices in descending filtration order.
// //! // When computing PH in dimension d, we only need to iterate
// //! // over simplices of dimension d and below.
// //! let iter_keymaj = boundary_matrix.row_indices_in_descending_order(maxdim);
// //! 
// //! // Compute a umatch factorization of the boundary matrix.
// //! let umatch = Umatch::factor_with_clearing(
// //!         boundary_matrix, 
// //!         iter_keymaj.clone(), 
// //!         ring_operator.clone(), 
// //!         OrderOperatorAutoLt::new(), 
// //!         OrderOperatorAutoLt::new(),
// //!     );
// //! 
// //! // Get the matrix `C` from the U-match factorization.
// //! // This matrix is also called the `domain COMB`.
// //! let comb_domain = umatch.comb_domain();
// //! 
// //! // Get the matrix `R` from the U-match factorization.
// //! // This matrix is also called the `codomain COMB`.
// //! let comb_codomain = umatch.comb_domain();
// //! 
// //! // Get the matrix `M` from the U-match factorization.  
// //! // This matrix is also called the matching array.
// //! let matching = umatch.matching_ref();
// //! 
// //! // Print the birth/death pairs and corresopnding basis vectors / cycle representatives
// //! for simplex in iter_keymaj {
// //! 
// //!     // This `if` statement handles the case where `simplex` belongs to a
// //!     // birth-death pair as described in condition (I) of the mathematical overview,
// //!     // above.
// //!     if let Some( cofacet ) = matching.keymaj_to_keymin( &simplex ) {
// //!         let birth_death_pair = (&simplex, &cofacet);
// //!         let cycle_rep: Vec<_> = comb_codomain.view_minor_descend( simplex.clone() ).collect();
// //!         println!(""); // new liine
// //!         println!("there is a birth-death pair {:?}", birth_death_pair );
// //!         println!("the corresponding cycle representative is {:?}", cycle_rep);
// //!     }
// //! 
// //!     // This handles condition (II) of the mathematical overview, above.
// //!     else if ! matching.contains_keymin( & simplex ) {
// //!         let cycle_rep: Vec<_> = comb_domain.view_minor_descend( simplex.clone() ).collect();
// //!         println!(""); // new liine
// //!         println!("there is an essential birth simplex {:?}", & simplex );
// //!         println!("the corresponding cycle representative is {:?}", cycle_rep);
// //!     }
// //! }
// //! 
// //! // This should print the following:
// //! // 
// //! // ```bash
// //! // $ there is a birth-death pair ([1, 2], [0, 1, 2])
// //! // $ the corresponding cycle representative is [([1, 2], 1), ([0, 2], 2), ([0, 1], 1)]
// //! // $ 
// //! // $ there is an essential birth simplex [1, 3]
// //! // $ the corresponding cycle representative is [([1, 3], 1), ([0, 3], 2), ([0, 1], 1)]
// //! // $ 
// //! // $ there is an essential birth simplex [2, 3]
// //! // $ the corresponding cycle representative is [([2, 3], 1), ([0, 3], 2), ([0, 2], 1)]
// //! // 
// //! // To extract the barcode, simply enumerate birth/death pairs, and calculate
// //! // the corresponding intervals.







// #[cfg(test)]
// mod doc_test_drafts {
//     use itertools::Itertools;
//     use oat_rust::utilities::sequences_and_ordinals::SortedVec;


//     #[test]
//     fn compute_homology_var() {

//         use crate::topology::simplicial::from::relation::BoundaryMatrixDowker;      
//         use crate::algebra::matrices::query::ViewColDescend;
//         use crate::algebra::matrices::operations::umatch::row_major::{Umatch::factor_with_clearing};
//         use oat_rust::utilities::order::OrderOperatorAutoLt;        
//         use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;

//         // Define the ring operator for the finite field of order 3.
//         // You can use this object to perform arithmetic operations, e.g., 
//         // to add 1 and 1, you can run `ring_operator.add(1,1)`.
//         let ring_operator          =   PrimeOrderFieldOperator::new(3);

//         // We will build a dowker complex.
//         // A dowker complex is defined by a vertex set V and a family S
//         // of subsets of V.  A subset of V forms a simplex iff it is 
//         // a subset of some element of S.  We refer to the elements 
//         // of S as "dowker simplices".

//         // Define the maximum homology dimension we're interested in.
//         let maxdim    =   2;

//         // Here is one storage format for the set S.
//         let dowker_simplices =  
//                 vec![ 
//                         vec![0,1,2],
//                         vec![0,3],                      
//                         vec![1,3], 
//                         vec![2,3]                                           
//                     ];      

//         // Build the boundary matrix.
//         // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
//         let boundary_matrix = BoundaryMatrixDowker::from_vectors( dowker_simplices, ring_operator.clone() );

//         // This iterates over simplices in descending filtration order.
//         // When computing PH in dimension d, we only need to iterate
//         // over simplices of dimension d and below.
//         let iter_keymaj = boundary_matrix.row_indices_in_descending_order(maxdim);

//         // Compute a umatch factorization of the boundary matrix.
//         let umatch = Umatch::factor_with_clearing(
//                 boundary_matrix, 
//                 iter_keymaj.clone(), 
//                 ring_operator.clone(), 
//                 OrderOperatorAutoLt::new(), 
//                 OrderOperatorAutoLt::new(),
//             );
        
//         // Get the matrix `C` from the U-match factorization.
//         // This matrix is also called the `domain COMB`.
//         let comb_domain = umatch.comb_domain();
        
//         // Get the matrix `R` from the U-match factorization.
//         // This matrix is also called the `codomain COMB`.
//         let comb_codomain = umatch.comb_domain();
        
//         // Get the matrix `M` from the U-match factorization.  
//         // This matrix is also called the matching array.
//         let matching = umatch.matching_ref();
        
//         // Print the birth/death pairs and corresopnding basis vectors / cycle representatives
//         for simplex in iter_keymaj {
        
//             // This `if` statement handles the case where `simplex` belongs to a
//             // birth-death pair as described in condition (I) of the mathematical overview,
//             // above.
//             if let Some( cofacet ) = matching.keymaj_to_keymin( &simplex ) {
//                 let birth_death_pair = (&simplex, &cofacet);
//                 let cycle_rep: Vec<_> = comb_codomain.view_minor_descend( simplex.clone() ).collect();
//                 println!(""); // new liine
//                 println!("there is a birth-death pair {:?}", birth_death_pair );
//                 println!("the corresponding cycle representative is {:?}", cycle_rep);
//             }
        
//             // This handles condition (II) of the mathematical overview, above.
//             else if ! matching.contains_keymin( & simplex ) {
//                 let cycle_rep: Vec<_> = comb_domain.view_minor_descend( simplex.clone() ).collect();
//                 println!(""); // new liine
//                 println!("there is an essential birth simplex {:?}", & simplex );
//                 println!("the corresponding cycle representative is {:?}", cycle_rep);
//             }
//         }

//         // This should print the following:
//         // 
//         // ```bash
//         // $ there is a birth-death pair ([1, 2], [0, 1, 2])
//         // $ the corresponding cycle representative is [([1, 2], 1), ([0, 2], 2), ([0, 1], 1)]
//         // $ 
//         // $ there is an essential birth simplex [1, 3]
//         // $ the corresponding cycle representative is [([1, 3], 1), ([0, 3], 2), ([0, 1], 1)]
//         // $ 
//         // $ there is an essential birth simplex [2, 3]
//         // $ the corresponding cycle representative is [([2, 3], 1), ([0, 3], 2), ([0, 2], 1)]
//         // 
//         // To extract the barcode, simply enumerate birth/death pairs, and calculate
//         // the corresponding intervals.
//     }





//     #[test]
//     fn compute_homology_projective() {

//         use std::collections::HashSet;
//         use std::iter::FromIterator;
//         use crate::topology::simplicial::from::relation::BoundaryMatrixDowker;
//         use crate::topology::simplicial::simplices::vector::{subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend};        
//         use crate::algebra::matrices::query::ViewColDescend;
//         use crate::algebra::matrices::operations::umatch::row_major::{Umatch::factor_with_clearing};
//         use oat_rust::utilities::order::OrderOperatorAutoLt;        
//         use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
//         use oat_rust::rings::operator_structs::ring_native::DivisionRingNative;

//         // Define the ring operator for the finite field of order 3.
//         // You can use this object to perform arithmetic operations, e.g., 
//         // to add 1 and 1, you can run `ring_operator.add(1,1)`.
//         // let ring_operator          =   DivisionRingNative::<Ratio<i64>>::new();
//         let ring_operator = PrimeOrderFieldOperator::new(2);

//         // We will build a dowker complex.
//         // A dowker complex is defined by a vertex set V and a family S
//         // of subsets of V.  A subset of V forms a simplex iff it is 
//         // a subset of some element of S.  We refer to the elements 
//         // of S as "dowker simplices".

//         // Define the maximum homology dimension of interest
//         let maxdim = 1;

//         // Define the Dowker simplices
//         let dowker_simplices =  
//                 vec![ 
//                         vec![0, 1, 2], vec![0, 3, 4], vec![1, 3, 5], vec![2, 4, 5], vec![0, 2, 3], vec![2, 3, 5], vec![1, 2, 4], vec![0, 1, 5], vec![1, 3, 4], vec![0, 4, 5]                                       
//                     ];      

//         // Build the boundary matrix.
//         // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
//         let boundary_matrix = BoundaryMatrixDowker::from_vectors( dowker_simplices, ring_operator.clone() );

//         // This iterates over simplices in descending order of 
//         // dimension (first) and descending lexicographic order (second).
//         // When computing homology of dimension d, we only need to iterate
//         // over simplices of dimension d and below.
//         let iter_keymaj = boundary_matrix.row_indices_in_descending_order(maxdim);

//         oat_rust::matrices::display::print_indexed_major_views( &boundary_matrix, iter_keymaj.clone() );

//         println!("now descend");       
//         for keymaj in iter_keymaj.clone() { println!("{:?}", keymaj)  };

//         println!("now ascend");
//         for keymaj in subsimplices_dim_0_thru_d_iter_ascend( & dowker_simplices_vec_format, 1 ) { println!("{:?}", keymaj)  };        


//         // Compute a umatch factorization of the boundary matrix.
//         // For details on what this factorization entails, see the paper 
//         // "U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co) homology"
//         // by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.  You can also check out the
//         // OAT documentation for `umatch`.
//         let umatch = Umatch::factor_with_clearing(
//                 boundary_matrix, 
//                 iter_keymaj, 
//                 ring_operator.clone(), 
//                 OrderOperatorAutoLt::new(), 
//                 OrderOperatorAutoLt::new(),
//             );
        
//         // Get the domain COMB (cf the  paper on umatch factorization)
//         let comb_domain = umatch.comb_domain();
        
//         // Get the matching array (cf the  paper on umatch factorization)
//         let matching = umatch.matching_ref();
        
//         // The set {columns of the domain comb that are not matched upward or downward} 
//         // forms a basis for homology
//         let dim = 1;
//         let mut betti = 0;
        
//         // Print the basis vectors
//         println!(""); // an empty line, for spacing        
//         println!("Each of the following lines represents a basis vector for homology in dimension {:?}", dim);
//         for simplex in subsimplices_dim_d_iter_ascend( &dowker_simplices_vec_format, dim ).unwrap() {
//             if matching.contains_keymaj( &simplex ) { continue }
//             if matching.contains_keymin( &simplex ) { continue }   
//             let basis_vec = comb_domain.view_minor_descend( simplex );
//             let basis_vec: Vec<_> = basis_vec.collect();
//             println!("basis vector {:?}: {:?}", betti, basis_vec);
//             betti += 1;
//         }

//         // Print the betti number
//         println!(""); // an empty line, for spacing
//         println!("The betti number in dimension {:?} is {:?}.", dim, betti);
//         println!(""); // an empty line, for spacing        

//         // Print the matching matrix scalars
//         println!("The nonzero matching entries are {:?}.", umatch.matching_ref().vec_snzval_ref() );

//         // Make sure that all changes are saved to `main.rs`, then run the
//         // program as described in the [quick start tutorial in OAT]().
//         // This should print the following:
//         // 
//         // Each of the following lines represents a basis vector for homology in dimension 1
//         // basis vector 0: [([1, 3], 1), ([0, 3], 2), ([0, 1], 2)]
//         // basis vector 1: [([2, 3], 1), ([0, 3], 2), ([0, 2], 2)]
//         // 
//         // The betti number in dimension 1 is 2
//     }
// }






