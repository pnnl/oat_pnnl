
//! An example homology calculation.
//! 
//! Here we compute a basis for homology in dimension 1,
//! for the simplicial complex on vertex set `{0,1,2,3}`
//! whose maximal faces are `{0,1,2}, {0,3}, {1,3}, {2,3}`.
//! The coefficient field is the finite field of order 3.
//! 
//! **To run this example on your desktop computer** first checkout the
//! [quick start tutorial in oat_rust]() for instructions
//! on installing Rust and running a program.  As part of this process,
//! you'll create a new folder that contains a file called `main.rs`.  Inside
//! `main.rs` is some text that reads `fn main{ .. }`.  Delete everything
//! between `{` and `}`, and paste in the following:
//! 
//!
//! ```
//! use std::collections::HashSet;
//! use std::iter::FromIterator;
//! use oat_rust::matrices::matrix_oracle_traits::OracleMinorDescend;
//! use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor_with_clearing};
//! use sohar_rust::simplicial::from::relation::DowkerComplexBoundaryMatrixRowMajor;
//! use sohar_rust::simplicial::simplices::{subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend};
//! use oat_rust::utilities::partial_order::OrderComparatorAutoLt;        
//! use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
//! 
//! // Define the ring operator for the finite field of order 3.
//! // You can use this object to perform arithmetic operations, e.g., 
//! // to add 1 and 1, you can run `ring_operator.add(1,1)`.
//! let ring_operator          =   PrimeOrderFieldOperator::new(3);
//! 
//! // We will build a dowker complex.
//! // A dowker complex is defined by a vertex set V and a family S
//! // of subsets of V.  A subset of V forms a simplex iff it is 
//! // a subset of some element of S.  We refer to the elements 
//! // of S as "dowker simplices".
//! 
//! // Each dowker simplex is represented by a vector, and we store the
//! // list of all such simplices inside a larger vector.
//! let dowker_simplices_vec_format =   vec![    
//!                                         vec![0,1,2], 
//!                                         vec![0,3], 
//!                                         vec![1,3], 
//!                                         vec![2,3]  
//!                                     ];                                
//! 
//! // For certain calculations it's better to store a simplex as a *set*,
//! // rather than a vector.  This object stores the same collection of 
//! // simplices, represented by sets.
//! let dowker_simplices_set_format: Vec<_>  =   dowker_simplices_vec_format
//!                                                         .iter()
//!                                                         .cloned()
//!                                                         .map( |x| HashSet::from_iter( x ) )
//!                                                         .collect();
//! 
//! // Build the boundary matrix.
//! // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
//! let boundary_matrix = DowkerComplexBoundaryMatrixRowMajor::new( dowker_simplices_set_format, ring_operator.clone() );
//! 
//! // This iterates over simplices in descending order of 
//! // dimension (first) and descending lexicographic order (second).
//! // When computing homology of dimension d, we only need to iterate
//! // over simplices of dimension d and below.
//! let iter_keymaj = subsimplices_dim_0_thru_d_iter_descend( &dowker_simplices_vec_format, 1 );
//! 
//! // Compute a umatch factorization of the boundary matrix.
//! // For details on what this factorization entails, see the paper 
//! // "U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co) homology"
//! // by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.  You can also check out the
//! // oat_rust documentation for `UmatchRowMajor`.
//! let umatch = new_umatchrowmajor_with_clearing(
//!         boundary_matrix, 
//!         iter_keymaj, 
//!         ring_operator.clone(), 
//!         OrderComparatorAutoLt::new(), 
//!         OrderComparatorAutoLt::new(),
//!     );
//! 
//! // Get the domain COMB (cf the  paper on umatch factorization)
//! let array_comb_domain = umatch.array_comb_domain();
//! 
//! // Get the matching array (cf the  paper on umatch factorization)
//! let array_matching = umatch.array_matching_ref();
//! 
//! // The set {columns of the domain comb that are not matched upward or downward} 
//! // forms a basis for homology
//! let dim = 1;
//! let mut betti = 0;
//! 
//! // Print the basis vectors
//! println!(""); // an empty line, for spacing        
//! println!("Each of the following lines represents a basis vector for homology in dimension {:?}", dim);
//! for simplex in subsimplices_dim_d_iter_ascend( &dowker_simplices_vec_format, dim ) {
//!     if array_matching.contains_keymaj( &simplex ) { continue }
//!     if array_matching.contains_keymin( &simplex ) { continue }        
//!     let basis_vec = array_comb_domain.view_minor_descend( simplex );
//!     let basis_vec: Vec<_> = basis_vec.collect();
//!     println!("basis vector {:?}: {:?}", betti, basis_vec);
//!     betti += 1;
//! }
//! 
//! // Print the betti number
//! println!(""); // an empty line, for spacing
//! println!("The betti number in dimension {:?} is {:?}.", dim, betti);
//! println!(""); // an empty line, for spacing    
//! ```    
//! 
//! Make sure that all changes are saved to `main.rs`, then run the
//! program as described in the [quick start tutorial in oat_rust]().
//! This should print the following:
//! 
//! ```bash
//! $ Each of the following lines represents a basis vector for homology in dimension 1
//! $ basis vector 0: [([1, 3], 1), ([0, 3], 2), ([0, 1], 2)]
//! $ basis vector 1: [([2, 3], 1), ([0, 3], 2), ([0, 2], 2)]
//! $ 
//! $ The betti number in dimension 1 is 2
//! ```
//! 
//! 
//! # Change the coefficient ring
//! 
//! oat_rust has a number of different [predefined coefficient rings](oat_rust::rings), which you can substitute into
//! the example above in order to calculate homology with different coefficients.  Simply replace the
//! line 
//! ```ignore
//! let ring_operator   =   PrimeOrderFieldOperator::new(3);
//! ``` 
//! with one of the `let ring_operator = ...` lines listed under *Predefined Rings*, [here](oat_rust::rings).
//! 
//! 
//! 
//! 
// //! 
// //! In the example above, we defined the ring operator in two steps: first, importing the definition of the operator
// //! using a `use` statement, then building the operator itself using a `new` function:
// //! 
// //! ```
// //! // import the definition
// //! use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator; 
// //! // build the operator
// //! let ring_operator   =   PrimeOrderFieldOperator::new(3); 
// //! ```
// //! 
// //! If we preffered, we could have combined these steps:
// //! 
// //! ```
// //! let ring_operator   =   oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator::new(3);
// //! ```
// //! 
// //! oat_rust has a number of different [predefined coefficient rings](oat_rust::rings) (see the section titled *Predefined rings*), which you can substitute into
// //! the example above, to calculate homology with different coefficients.



#[cfg(test)]
mod doc_test_drafts {
    use itertools::Itertools;
    use num::rational::Ratio;

    use crate::simplicial::simplices::{subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend, subsimplices_dim_0_thru_d_iter_ascend};        


    #[test]
    fn compute_homology_var() {

        use std::collections::HashSet;
        use std::iter::FromIterator;
        use crate::simplicial::from::relation::DowkerComplexBoundaryMatrixRowMajor;
        use crate::simplicial::simplices::{subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend};        
        use oat_rust::matrices::matrix_oracle_traits::OracleMinorDescend;
        use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor_with_clearing};
        use oat_rust::utilities::partial_order::OrderComparatorAutoLt;        
        use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;

        // Define the ring operator for the finite field of order 3.
        // You can use this object to perform arithmetic operations, e.g., 
        // to add 1 and 1, you can run `ring_operator.add(1,1)`.
        let ring_operator          =   PrimeOrderFieldOperator::new(3);

        // We will build a dowker complex.
        // A dowker complex is defined by a vertex set V and a family S
        // of subsets of V.  A subset of V forms a simplex iff it is 
        // a subset of some element of S.  We refer to the elements 
        // of S as "dowker simplices".

        // Each dowker simplex is represented by a vector, and we store the
        // list of all such simplices inside a larger vector.
        let dowker_simplices_vec_format =  vec![    
                                                                    vec![0,1,2], 
                                                                    vec![0,3], 
                                                                    vec![1,3], 
                                                                    vec![2,3]  
                                                                ];                                

        // For certain calculations it's better to store a simplex as a *set*,
        // rather than a vector.  This object stores the same collection of 
        // simplices, represented by sets.
        let dowker_simplices_set_format: Vec<_>  =   dowker_simplices_vec_format
                                                                .iter()
                                                                .cloned()
                                                                .map( |x| HashSet::from_iter( x ) )
                                                                .collect();

        // Build the boundary matrix.
        // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
        let boundary_matrix = DowkerComplexBoundaryMatrixRowMajor::new( dowker_simplices_set_format, ring_operator.clone() );

        // This iterates over simplices in descending order of 
        // dimension (first) and descending lexicographic order (second).
        // When computing homology of dimension d, we only need to iterate
        // over simplices of dimension d and below.
        let iter_keymaj = subsimplices_dim_0_thru_d_iter_descend( &dowker_simplices_vec_format, 1 );

        // Compute a umatch factorization of the boundary matrix.
        // For details on what this factorization entails, see the paper 
        // "U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co) homology"
        // by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.  You can also check out the
        // oat_rust documentation for `UmatchRowMajor`.
        let umatch = new_umatchrowmajor_with_clearing(
                boundary_matrix, 
                iter_keymaj, 
                ring_operator.clone(), 
                OrderComparatorAutoLt::new(), 
                OrderComparatorAutoLt::new(),
            );
        
        // Get the domain COMB (cf the  paper on umatch factorization)
        let array_comb_domain = umatch.array_comb_domain();
        
        // Get the matching array (cf the  paper on umatch factorization)
        let array_matching = umatch.array_matching_ref();
        
        // The set {columns of the domain comb that are not matched upward or downward} 
        // forms a basis for homology
        let dim = 1;
        let mut betti = 0;
        
        // Print the basis vectors
        println!(""); // an empty line, for spacing        
        println!("Each of the following lines represents a basis vector for homology in dimension {:?}", dim);
        for simplex in subsimplices_dim_d_iter_ascend( &dowker_simplices_vec_format, dim ) {
            if array_matching.contains_keymaj( &simplex ) { continue }
            if array_matching.contains_keymin( &simplex ) { continue }        
            let basis_vec = array_comb_domain.view_minor_descend( simplex );
            let basis_vec: Vec<_> = basis_vec.collect();
            println!("basis vector {:?}: {:?}", betti, basis_vec);
            betti += 1;
        }

        // Print the betti number
        println!(""); // an empty line, for spacing
        println!("The betti number in dimension {:?} is {:?}.", dim, betti);
        println!(""); // an empty line, for spacing        

        // This should print the following:
        // 
        // Each of the following lines represents a basis vector for homology in dimension 1
        // basis vector 0: [([1, 3], 1), ([0, 3], 2), ([0, 1], 2)]
        // basis vector 1: [([2, 3], 1), ([0, 3], 2), ([0, 2], 2)]
        // 
        // The betti number in dimension 1 is 2
    }





    #[test]
    fn compute_homology_projective() {

        use std::collections::HashSet;
        use std::iter::FromIterator;
        use crate::simplicial::from::relation::DowkerComplexBoundaryMatrixRowMajor;
        use crate::simplicial::simplices::{subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend};        
        use oat_rust::matrices::matrix_oracle_traits::OracleMinorDescend;
        use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor_with_clearing};
        use oat_rust::utilities::partial_order::OrderComparatorAutoLt;        
        use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;

        // Define the ring operator for the finite field of order 3.
        // You can use this object to perform arithmetic operations, e.g., 
        // to add 1 and 1, you can run `ring_operator.add(1,1)`.
        // let ring_operator          =   DivisionRingNative::<Ratio<i64>>::new();
        let ring_operator = PrimeOrderFieldOperator::new(2);

        // We will build a dowker complex.
        // A dowker complex is defined by a vertex set V and a family S
        // of subsets of V.  A subset of V forms a simplex iff it is 
        // a subset of some element of S.  We refer to the elements 
        // of S as "dowker simplices".

        // Here is one storage format for the set S.
        let dowker_simplices_vec_format =  
                vec![ 
                        vec![0, 1, 2], vec![0, 3, 4], vec![1, 3, 5], vec![2, 4, 5], vec![0, 2, 3], vec![2, 3, 5], vec![1, 2, 4], vec![0, 1, 5], vec![1, 3, 4], vec![0, 4, 5]                                       
                    ];      

        // Here is another format; it replaces inner vectors with sets; each format has its own use.
        let dowker_simplices_set_format  =  dowker_simplices_vec_format
                                                    .iter()
                                                    .cloned()
                                                    .map( |x| HashSet::from_iter( x ) )
                                                    .collect_vec();

        // Build the boundary matrix.
        // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
        let boundary_matrix = DowkerComplexBoundaryMatrixRowMajor::new( dowker_simplices_set_format, ring_operator.clone() );

        // This iterates over simplices in descending order of 
        // dimension (first) and descending lexicographic order (second).
        // When computing homology of dimension d, we only need to iterate
        // over simplices of dimension d and below.
        let iter_keymaj = subsimplices_dim_0_thru_d_iter_descend( &dowker_simplices_vec_format, 1 );



        oat_rust::matrices::display::print_indexed_major_views( &boundary_matrix, iter_keymaj.clone() );

        println!("now descend");       
        for keymaj in iter_keymaj.clone() { println!("{:?}", keymaj)  };

        println!("now ascend");
        for keymaj in subsimplices_dim_0_thru_d_iter_ascend( & dowker_simplices_vec_format, 1 ) { println!("{:?}", keymaj)  };        


        // Compute a umatch factorization of the boundary matrix.
        // For details on what this factorization entails, see the paper 
        // "U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co) homology"
        // by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.  You can also check out the
        // oat_rust documentation for `umatch`.
        let umatch = new_umatchrowmajor_with_clearing(
                boundary_matrix, 
                iter_keymaj, 
                ring_operator.clone(), 
                OrderComparatorAutoLt::new(), 
                OrderComparatorAutoLt::new(),
            );
        
        // Get the domain COMB (cf the  paper on umatch factorization)
        let array_comb_domain = umatch.array_comb_domain();
        
        // Get the matching array (cf the  paper on umatch factorization)
        let array_matching = umatch.array_matching_ref();
        
        // The set {columns of the domain comb that are not matched upward or downward} 
        // forms a basis for homology
        let dim = 1;
        let mut betti = 0;
        
        // Print the basis vectors
        println!(""); // an empty line, for spacing        
        println!("Each of the following lines represents a basis vector for homology in dimension {:?}", dim);
        for simplex in subsimplices_dim_d_iter_ascend( &dowker_simplices_vec_format, dim ) {
            if array_matching.contains_keymaj( &simplex ) { continue }
            if array_matching.contains_keymin( &simplex ) { continue }   
            let basis_vec = array_comb_domain.view_minor_descend( simplex );
            let basis_vec: Vec<_> = basis_vec.collect();
            println!("basis vector {:?}: {:?}", betti, basis_vec);
            betti += 1;
        }

        // Print the betti number
        println!(""); // an empty line, for spacing
        println!("The betti number in dimension {:?} is {:?}.", dim, betti);
        println!(""); // an empty line, for spacing        

        // Print the matching matrix scalars
        println!("The nonzero matching entries are {:?}.", umatch.array_matching_ref().vec_snzval_ref() );

        // Make sure that all changes are saved to `main.rs`, then run the
        // program as described in the [quick start tutorial in oat_rust]().
        // This should print the following:
        // 
        // Each of the following lines represents a basis vector for homology in dimension 1
        // basis vector 0: [([1, 3], 1), ([0, 3], 2), ([0, 1], 2)]
        // basis vector 1: [([2, 3], 1), ([0, 3], 2), ([0, 2], 2)]
        // 
        // The betti number in dimension 1 is 2
    }
}