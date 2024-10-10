//! U-match factorization for row-major matrices
//! 
//! # Quick start
//!
//! To factor a matrix, call [Umatch::factor].  This will produce a [Umatch], which you can use to obtain copies of the matching matrix, the COMB's, and more.
//! 
//! 
//! 
//! # Example
//! 
//! ```
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
//! use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
//! use oat_rust::algebra::matrices::operations::{umatch::row_major::Umatch, multiply::ProductMatrix};
//! use oat_rust::algebra::matrices::debug::verify_that_product_is_identity;
//! use oat_rust::algebra::matrices::query::{ViewRowAscend, ViewColDescend};
//! use oat_rust::algebra::matrices::display::print_indexed_minor_views;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::order::{OrderOperatorAuto, ReverseOrder};
//! use itertools::Itertools;
//! 
//! // DEFINE INPUTS
//! // ===============================
//! 
//! // define the ring operator and order operator
//! let modulus               =   5;
//! let ring_operator         =   PrimeOrderFieldOperator::new( modulus );        
//! let order_operator        =   OrderOperatorAuto;
//! 
//! // define the matrix we wish to factor
//! let mapping_data          =   VecOfVec::new( 
//!                                   vec![   
//!                                       vec![(0,1), (1,1), (2,1)],  // row 0
//!                                       vec![                   ],  // row 1
//!                                       vec![              (2,1)],  // row 2
//!                                   ] 
//!                               );
//! let mapping               =   & mapping_data;
//!                                 
//! // COMPUTE U-MATCH
//! // ===============================
//!                                 
//! let umatch
//!     =   Umatch::factor(
//!             mapping,  // the matrix we wish to factor
//!             (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
//!             ring_operator.clone(), // the operator for the coefficient ring
//!             order_operator.clone(), // order comparator for the entries in each row
//!             order_operator.clone(), // order comparator for the entries in each column
//!         );
//!     
//!     
//! // INSPECT FACTORIZATION
//! // ===============================
//!     
//! // extract R, R^{-1}, C, C^{-1}, and M
//! let r           =   umatch.comb_codomain();        // the codomain COMB
//! let rinv        =   umatch.comb_codomain_inv();    // inverse of the the codomain COMB
//! let c           =   umatch.comb_domain();          // the domain COMB
//! let cinv        =   umatch.comb_domain_inv();      // inverse of the domain COMB
//! let m           =   umatch.matching_ref();         // the generalized matching matrix
//!     
//!     
//! println!("\nMinor views of the codomain COMB");   print_indexed_minor_views( &r, 0..3 ); 
//! println!("\nMinor views of the   domain COMB");   print_indexed_minor_views( &c, 0..3 ); 
//! println!("\nMinor views of the matching matrix"); print_indexed_minor_views( &m, 0..3 ); 
//!     
//! // this will print the following:
//! //
//! // Minor views of the codomain COMB
//! // minor_view 0: [(0, 1)]
//! // minor_view 1: [(1, 1)]
//! // 
//! // Minor views of the   domain COMB
//! // minor_view 0: [(0, 1)]
//! // minor_view 1: [(1, 1), (0, 3)]
//! // minor_view 2: [(2, 1), (0, 2)]
//! // 
//! // Minor views of the matching matrix
//! // minor_view 0: [(0, 1)]
//! // minor_view 1: []
//! // minor_view 2: [(1, 1)]
//! 
//! // SOLVE Ax = b FOR x
//! // ===============================
//! 
//! let b   =   [ (2,1), (0,1) ]; // note we list entries in reverse order
//! let x   =   umatch.solve_dx_equals_b( b.clone() ).unwrap().collect_vec();
//! let dx  =   umatch.multiply_dv(x);
//! assert!( dx.eq( b ) );
//!     
//!     
//! // VERIFY THE CALCULATION
//! // ===============================
//!     
//! // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
//! verify_that_product_is_identity( &c, &cinv, 0..3, ring_operator.clone(), OrderOperatorAuto );
//!     
//! // check that the product of the codomain COMB with its inverse is identity: R * R^{-1} = I
//! verify_that_product_is_identity( &r, &rinv, 0..3, ring_operator.clone(), OrderOperatorAuto );
//!     
//! // check the factorization: R^{-1} * D * C = M
//! let rinv_d   = ProductMatrix::new( &rinv,   &mapping, ring_operator.clone(), OrderOperatorAuto );      
//! let rinv_d_c = ProductMatrix::new( &rinv_d, &c,       ring_operator.clone(), OrderOperatorAuto );                
//! for keymaj in 0 .. 3 { 
//!     assert_eq!(
//!         rinv_d_c.view_major_ascend( keymaj ).collect_vec(),
//!         m.view_major_ascend( keymaj ).collect_vec()
//!     ) 
//! }   
//! ```
//! 
//! 
//! # About
//! 
//! 
//! 
//! To learn more about U-match factorization, check out the [umatch](crate::algebra::matrices::operations::umatch) module, or see the paper [U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co)homology](https://arxiv.org/abs/2108.08831), by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.
//! All COMB's are calculated as described in this paper.

pub mod comb;
pub mod construction;


//  ==================================



use itertools::Itertools;



use comb::*;
use construction::*;

use crate::algebra::matrices::operations::solve::triangle::{TriangularSolverMinorDescend, TriangularSolverMajorAscend};
use crate::algebra::matrices::types::bimajor::{MatrixBimajor, MatrixBimajorData};
use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
use crate::algebra::vectors::entries::{KeyValSet, KeyValNew};
use crate::algebra::matrices::operations::transform_entry_wise::ReindexSquareMatrix;
use crate::algebra::matrices::types::matching::{GeneralizedMatchingArrayWithMajorOrdinals};
use crate::algebra::matrices::types::prepend_viewmaj::PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend;

use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients, ViewsMinorDescend};
use crate::utilities::functions::evaluate::{ EvaluateFunctionFnMutWrapper };

use crate::utilities::iterators::general::{FilterOutMembers};
use crate::utilities::iterators::merge::hit::{HitMerge};
use crate::algebra::vectors::entries::{KeyValGet};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};

use crate::utilities::order::{JudgePartialOrder, OrderOperatorByKeyCutsom, ReverseOrder};


use crate::algebra::vectors::operations::{Scale, VectorOperations, Simplify};

use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Debug;
use std::iter::{Cloned, Peekable};



use crate::algebra::matrices::operations::solve::echelon::{EchelonSolverMajorAscendWithMajorKeys};

use crate::algebra::matrices::operations::transform_vector_wise::{OnlyKeyMinInsideCollection, OnlyKeyMinOutsideCollection, OnlyKeyMajInsideCollection, OnlyKeyMajOutsideCollection};










/// Provides method to identify elements on the Pareto frontier of a matrix
pub trait ParetoShortCircuit<T> {
    /// Returns `Some(T)` when `self` is a major view, and the first element of `self` lies on the Pareto frontier.
    fn pareto_short_circuit(& self) -> Option< T >;
}





//  =========================================================================================================
//  U-MATCH OBJECT
//  =========================================================================================================





/// A (compressed representation of) a U-match factorization; all matrices concerned are represented in row-major format.
/// 
/// Internally, this object stores a copy of the matrix `R_{\rho \rho}` defined in Hang et al., 
/// "Umatch factorization: ..." 2021 (paper link)[https://arxiv.org/abs/2108.08831] **with its diagonal elements deleted**.  
/// More precisely, it contains a `VecOfVec< usize, Mapping::Coefficient >`, whose `k`th row contains the off-diagonal elements of
/// `R_{\rho \rho}`; where we replace each index `\rho_i` with the corresponding integer `i`, so that elements of the 
/// `VecOfVec` are tuples of form `(usize, Mapping::Coefficient)`.
/// 
/// 
/// # Design notes
/// 
/// **Why keep `Mapping::ViewMajorAscend` as a generic type parameter?**  
/// 
/// **Why store the "off diagonal seed" of the codomain COMB as a `VecOfVec` instead of a `VecOfVec`?**  Because we want to wrap this
/// struct in a `PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend` struct; this wrapper suppies the missing diagonal entries of,
/// the seed, thus forming the complete seed of the codomain COMB, which is a primitive and fundamental object for this 
/// computational library.  If we wanted to create an oracle for the seed that exposed entries as *references* to tuples,
/// then we would have to store the diagonal entries in memory.  However,  [`numerical experiments`](https://arxiv.org/pdf/2108.08831.pdf) 
/// have shown that the number of diagonal entries often eclipses the number of off-diagonal elements.  Thus storing 
/// the diagonal entries might incur a nontrivial memory cost.
/// 
/// **Why store the "off diagonal seed" of the codomain COMB as a `VecOfVec` instead of a `Vec< Vec< (usize, Mapping::Coefficient) > >`?**
/// Because `VecOfVec` comes with certain guarantees about order of entries, which allows us to safely implement 
/// the `ViewRowAscend` and `ViewRowDescend` traits.
/// 
/// **Remark** One can always obtain a `VecOfVecFromBorrow` from a `VecOfVec` via the
/// [`from_VecOfVec`](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVecFromBorrow::from_VecOfVec) method.
#[derive(Clone, Debug)]
pub struct Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
    where   
        Mapping:                    ViewRowAscend + IndicesAndCoefficients,
        Mapping::ColIndex:          Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:          Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::ViewMajorAscend:   IntoIterator,
        Mapping::EntryMajor:        KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
{
    mapping:                      Mapping,
    matching:                     GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
    comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:   
                                MatrixBimajorData<
                                    VecOfVec< usize, Mapping::Coefficient >,
                                    VecOfVec< usize, Mapping::Coefficient >,
                                >,
    ring_operator:              RingOperator,
    order_operator_major:       OrderOperatorRowEntries,
    order_operator_minor:       OrderOperatorColEntries,    
}


//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- CONSTRUCTORS
//  ---------------------------------------------------------------------------------------------------------



impl < Mapping, RingOperator, OrderOperatorRowIndex, OrderOperatorColIndex, >  

    Umatch 
        < 
            Mapping, 
            RingOperator, 
            OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex, >,  // recall that entries in the major views have minor keys
            OrderOperatorByKeyCutsom< Mapping::RowIndex, Mapping::Coefficient, Mapping::EntryMinor, OrderOperatorRowIndex, >, // recall that entries in the minor views have major keys
        >  
    
    where   
        Mapping::ColIndex:              Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:              Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        Mapping:                        ViewRowAscend + IndicesAndCoefficients,
        Mapping::ViewMajorAscend:       IntoIterator,
        Mapping::EntryMajor:            KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
        // OrderOperatorByKey<usize, Mapping::Coefficient, (usize, Mapping::Coefficient)>: JudgePartialOrder< (usize, Mapping::Coefficient)>
{

    /// Generate a new U-match factorization
    /// 
    /// # Arguments
    /// 
    /// - `mapping`: matrix you want to factor (or "mapping matrix")
    /// - `iter_keymaj`: an iterator that runs over the row indices (= major indices) of the matrix, in *strictly descending order*
    /// - `ring_operator`: operator for the coefficient ring
    /// - `order_operator_keymaj`: order operator for the row indices (= major keys)
    /// - `order_operator_keymin`: order operator for the column indices (= minor keys)
    /// 
    /// # Special features of this function
    /// 
    /// Although the Umatch struct can take arbitrary orders the entries in its major and minor views, this
    /// constructor requires the user to provide an order on minor and major keys.  This is a safety measure, to prevent
    /// common forms of errors.  However, if you really want a custom order on entries, there are other ways to construct
    /// a U-match factorization, besides this constructor.
    pub fn factor
            < IterRowIndex > ( 
                mapping:                    Mapping, 
                iter_keymaj:                IterRowIndex,
                ring_operator:              RingOperator,
                order_operator_keymaj:      OrderOperatorRowIndex,                   
                order_operator_keymin:      OrderOperatorColIndex,             
            ) 
        -> 
        Self

    where   Mapping:                    ViewRowAscend + IndicesAndCoefficients, 
            IterRowIndex:               Iterator < Item = Mapping::RowIndex >,
            Mapping::ColIndex:          Clone + Hash + std::cmp::Eq + Debug, 
            Mapping::RowIndex:          Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            Mapping::Coefficient:      Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            Mapping::ViewMajorAscend:   Clone, // !!! remove clone
            Mapping::ViewMajorAscendIntoIter: Clone,            
            Mapping::EntryMajor:        KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
            OrderOperatorColIndex:        Clone + JudgePartialOrder <  Mapping::ColIndex >, // !!! remove clone
            OrderOperatorRowIndex:        Clone + JudgePartialOrder <  Mapping::RowIndex >,
            
    {
        
        let ( comb_codomain_inv_off_diag_pivot_block, matching ) : ( VecOfVec<usize, Mapping::Coefficient>, GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient > )
            = get_codomain_comb_inv_off_diag_pivot_block( 
                    & mapping, 
                    iter_keymaj, 
                    ring_operator.clone(),
                    order_operator_keymin.clone(),
                    // order_operator_keymaj.clone(),
                );

        let comb_codomain_inv_off_diag_pivot_block
            =   MatrixBimajorData { 
                    matrix_minor_data: comb_codomain_inv_off_diag_pivot_block.transpose_deep( matching.num_pairs() ).unwrap(),
                    matrix_major_data: comb_codomain_inv_off_diag_pivot_block, 
                    // matrix_minor_data: transpose_deep, // the number of rows we specify is gauranteed to be correct; it makes the matrix square
                };
        
        Umatch{ 
                mapping, 
                matching, 
                comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     comb_codomain_inv_off_diag_pivot_block,   
                ring_operator,
                order_operator_major:      OrderOperatorByKeyCutsom
                                                                ::<Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex, >
                                                                ::new( order_operator_keymin ),
                order_operator_minor:     OrderOperatorByKeyCutsom
                                                                ::<Mapping::RowIndex, Mapping::Coefficient, Mapping::EntryMinor, OrderOperatorRowIndex, >
                                                                ::new( order_operator_keymaj ),                                                                
                // order_operator_minor:     OrderOperatorByKeyCutsom::new( order_operator_keymaj ),                
                // phantom_viewmajorascend:        PhantomData,
            }
        
    }



    /// Same as [Umatch::factor], but applies the clearning optimization.
    /// 
    /// Concretely, this means that when iterating over row indices, the solver skips indices that have
    /// already been identified as nonzero columns in the matching array.
    pub fn factor_with_clearing
            < IterRowIndex, KeyBoth, > ( 
                mapping:              Mapping, 
                iter_keymaj:                IterRowIndex,
                ring_operator:              RingOperator,
                order_operator_keymin:      OrderOperatorColIndex,
                order_operator_keymaj:      OrderOperatorRowIndex,                
            ) 
        -> 
        Umatch< 
                Mapping, 
                RingOperator, 
                OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex >, 
                OrderOperatorByKeyCutsom< Mapping::RowIndex, Mapping::Coefficient, Mapping::EntryMinor, OrderOperatorRowIndex >, 
            > 

    where   Mapping:                           ViewRowAscend + IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >, 
            IterRowIndex:                             Iterator < Item = Mapping::RowIndex >,
            Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq + Debug, 
            Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            Mapping::Coefficient:                   Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:                           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            Mapping::ViewMajorAscend:          ParetoShortCircuit<Mapping::EntryMajor>, // !!! remove clone
            Mapping::EntryMajor:     KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
            OrderOperatorColIndex:                    Clone + JudgePartialOrder <  Mapping::ColIndex >, // !!! remove clone
            OrderOperatorRowIndex:                    Clone + JudgePartialOrder <  Mapping::RowIndex >,
            // HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex>>: Clone // !!!! remove this        
            
    {
        
        let ( comb_codomain_inv_off_diag_pivot_block, matching ) : ( VecOfVec<usize, Mapping::Coefficient>, GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient > )
            = get_codomain_comb_inv_off_diag_pivot_block_with_clearing( 
                    & mapping, 
                    iter_keymaj, 
                    ring_operator.clone(),
                    order_operator_keymin.clone(),
                    // order_operator_keymaj.clone(),
                );

        let comb_codomain_inv_off_diag_pivot_block
            =   MatrixBimajorData { 
                    matrix_minor_data: comb_codomain_inv_off_diag_pivot_block.transpose_deep( matching.num_pairs() ).unwrap(), // the number of rows we specify is gauranteed to be correct; it makes the matrix square                
                    matrix_major_data: comb_codomain_inv_off_diag_pivot_block, 
                };                
        
        Umatch{ 
                mapping, 
                matching, 
                comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     comb_codomain_inv_off_diag_pivot_block,   
                ring_operator,
                order_operator_major:      OrderOperatorByKeyCutsom
                                                                ::<Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex, >
                                                                ::new( order_operator_keymin ),
                order_operator_minor:     OrderOperatorByKeyCutsom
                                                                ::<Mapping::RowIndex, Mapping::Coefficient, Mapping::EntryMinor, OrderOperatorRowIndex, >
                                                                ::new( order_operator_keymaj ),                                                                
                // order_operator_minor:     OrderOperatorByKeyCutsom::new( order_operator_keymaj ),                
                // phantom_viewmajorascend:        PhantomData,
            }
        
    }


}





//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- GENERAL IMPLEMENTATIONS (THERE ARE SPECIFIC IMPLEMENTATIONS FOR KEYMIN=KEMAJ BELOW)
//  ---------------------------------------------------------------------------------------------------------

     
impl < Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >  

    Umatch 
    < Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >  
    
    where   
        Mapping::ColIndex:                  Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                  Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        Mapping:                            ViewRowAscend + IndicesAndCoefficients,
        Mapping::ViewMajorAscend:           IntoIterator,
        Mapping::EntryMajor:                KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
        // OrderOperatorByKey<usize, Mapping::Coefficient, (usize, Mapping::Coefficient)>: JudgePartialOrder< (usize, Mapping::Coefficient)>
{
  
// }


    /// Generate a new U-match factorization (this constructor usually throws an error due to a technial issue re: type inference, so it is usually more convenient to use the function `new_umatchrowmajor`).
    pub fn new_custom_order
            < IterRowIndex, OrderOperatorColIndex, OrderOperatorRowIndex, EntryMinor > ( 
                mapping:                  Mapping, 
                iter_keymaj:              IterRowIndex,
                ring_operator:            RingOperator,
                order_operator_keymin:    OrderOperatorColIndex,
                order_operator_keymaj:    OrderOperatorRowIndex,                
            ) 
        -> 
        Umatch< 
                Mapping, 
                RingOperator, 
                OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex >, 
                OrderOperatorByKeyCutsom< Mapping::RowIndex, Mapping::Coefficient, EntryMinor, OrderOperatorRowIndex >, 
            > 

    where   Mapping:           ViewRowAscend + IndicesAndCoefficients,
            IterRowIndex:             Iterator < Item = Mapping::RowIndex >,
            Mapping::ColIndex:                 Clone + Hash + std::cmp::Eq + Debug, 
            Mapping::RowIndex:                 Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
            Mapping::Coefficient:                 Clone + Debug,            // !!! remove Debug eventually       
            RingOperator:           Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            Mapping::ViewMajorAscend:        Clone, // !!! remove clone
            Mapping::EntryMajor:  KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
            OrderOperatorColIndex:   Clone + JudgePartialOrder <  Mapping::ColIndex >, // !!! remove clone
            OrderOperatorRowIndex:   Clone + JudgePartialOrder <  Mapping::RowIndex >,
            HitMerge<Peekable<Scale<<Mapping::ViewMajorAscend as IntoIterator>::IntoIter, Mapping::ColIndex, RingOperator, Mapping::Coefficient>>, OrderOperatorByKeyCutsom< Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex>>: Clone // !!!! remove this        
            
    {
        
        let ( comb_codomain_inv_off_diag_pivot_block, matching ) : ( VecOfVec<usize, Mapping::Coefficient>, GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient > )
            = get_codomain_comb_inv_off_diag_pivot_block( 
                    & mapping, 
                    iter_keymaj, 
                    ring_operator.clone(),
                    order_operator_keymin.clone(),
                    // order_operator_keymaj.clone(),
                );

        let comb_codomain_inv_off_diag_pivot_block
            =   MatrixBimajorData { 
                    matrix_minor_data: comb_codomain_inv_off_diag_pivot_block.transpose_deep( matching.num_pairs() ).unwrap(), // the number of rows we specify is gauranteed to be correct; it makes the matrix square                
                    matrix_major_data: comb_codomain_inv_off_diag_pivot_block, 
                };                
        
        Umatch{ 
                mapping, 
                matching, 
                comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     comb_codomain_inv_off_diag_pivot_block,   
                ring_operator,
                order_operator_major:      OrderOperatorByKeyCutsom
                                                                ::<Mapping::ColIndex, Mapping::Coefficient, Mapping::EntryMajor, OrderOperatorColIndex, >
                                                                ::new( order_operator_keymin ),
                order_operator_minor:     OrderOperatorByKeyCutsom
                                                                ::<Mapping::RowIndex, Mapping::Coefficient, EntryMinor, OrderOperatorRowIndex, >
                                                                ::new( order_operator_keymaj ),                                                                
                // order_operator_minor:     OrderOperatorByKeyCutsom::new( order_operator_keymaj ),                
                // phantom_viewmajorascend:        PhantomData,
            }
        
    }


// //  =========================================================================================================
// //  U-MATCH REF OBJECT
// //  =========================================================================================================



// #[derive(Debug)]
// pub struct Umatch< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
//     where   
//         Mapping::ColIndex:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         Mapping::RowIndex:     Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         Mapping::ViewMajorAscend:            IntoIterator,
//         Mapping::EntryMajor:      KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
// {
//     mapping:                                          &'a Mapping,
//     matching:                                         &'b GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
//     comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     &'b VecOfVec< usize, Mapping::Coefficient >,
//     ring_operator:                                          RingOperator,
//     order_operator_major:                                 OrderOperatorRowEntries,
//     order_operator_minor:                                 OrderOperatorColEntries,    
//     phantom_viewmajorascend:                                PhantomData< Mapping::ViewMajorAscend >, // required b/c otherwise the compiler complains that the type parameter `Mapping::ViewMajorAscend` is unused
// }

// //  Implement
// //  ---------------------------------------------------------------------------------------------------------

// impl < 'a, 'b, Mapping, Mapping::ViewMajorAscend, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  

//     UmatchRowMajorWithRefs
//     < 'a, 'b, Mapping, Mapping::ViewMajorAscend, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  
    
//     where   
//         Mapping::ColIndex:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
//         Mapping::RowIndex:                 Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
//         Mapping:           ViewRowAscend + IndicesAndCoefficients,
//         Mapping::ViewMajorAscend:        IntoIterator,
//         Mapping::EntryMajor:  KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,


// {
//     pub fn new
//         ( umatch: &'b Umatch < 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  ) 
//         ->
//         UmatchRowMajorWithRefs
//              <'a, 'b, Mapping, Mapping::ViewMajorAscend, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
//         where
//             RingOperator:       Clone,
//             OrderOperatorRowEntries: Clone,
//             OrderOperatorColEntries: Clone,
//     {
//         UmatchRowMajorWithRefs{ 
//                 mapping:                                            umatch.mapping, 
//                 matching:                                         & umatch.matching,
//                 comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal:     & umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal,
//                 ring_operator:                                            umatch.ring_operator.clone(),
//                 order_operator_major:                                   umatch.order_operator_major.clone(),
//                 order_operator_minor:                                   umatch.order_operator_minor.clone(),                
//                 phantom_viewmajorascend:                                  PhantomData,
//             }
//     }

    /// Returns a copy of the ring operator
    pub fn ring_operator( &self ) -> RingOperator 
        where
            RingOperator:   Clone,
    { self.ring_operator.clone() }

    /// Returns a copy of the order comparator for index-value pairs whose index is a `Mapping::ColIndex`.
    pub fn order_operator_major( &self ) -> OrderOperatorRowEntries 
        where
            OrderOperatorRowEntries:   Clone,
    { (self.order_operator_major).clone() } 
    
    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `Mapping::ColIndex`.
    pub fn order_operator_major_reverse( &self ) -> ReverseOrder< OrderOperatorRowEntries >
        where
            OrderOperatorRowEntries:   Clone,
    { ReverseOrder::new( (self.order_operator_major).clone() ) }    
    
    
    // /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `Mapping::ColIndex`, wrapped in a struct that implements additional ordering traits.
    // pub fn order_operator_major_reverse_extended_functionality( &self ) -> InferTotalOrderFromJudgePartialOrder< ReverseOrder< OrderOperatorRowEntries > >
    //     where
    //         OrderOperatorRowEntries:   Clone,
    // { 
    //     InferTotalOrderFromJudgePartialOrder::new( // we have to construct an order comparator that spits out values of type Ord rather than type bool
    //         self.order_operator_major_reverse() 
    //     )         
    // }        

    /// Returns a copy of the order comparator for index-value pairs whose index is a `Mapping::RowIndex`.
    pub fn order_operator_minor( &self ) -> OrderOperatorColEntries 
        where
            OrderOperatorColEntries:   Clone,
    { self.order_operator_minor.clone() }  
    
    /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `Mapping::RowIndex`.
    pub fn order_operator_minor_reverse( &self ) -> ReverseOrder< OrderOperatorColEntries >
        where
            OrderOperatorColEntries:   Clone,
    { ReverseOrder::new( self.order_operator_minor.clone() ) }   


    // /// Returns a copy of the **inverted** order comparator for index-value pairs whose index is a `Mapping::RowIndex`, wrapped in a struct that implements additional ordering traits.
    // pub fn order_operator_minor_reverse_extended_functionality( &self ) -> InferTotalOrderFromJudgePartialOrder< ReverseOrder< OrderOperatorColEntries > >
    //     where
    //         OrderOperatorColEntries:   Clone,
    // { 
    //     InferTotalOrderFromJudgePartialOrder::new( // we have to construct an order comparator that spits out values of type Ord rather than type bool
    //         self.order_operator_minor_reverse() 
    //     )         
    // }       
    
    
   

    /// Returns the (row-major) codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn comb_codomain( &self ) -> CombCodomain< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  {
        CombCodomain{ umatch: self }
    }

    /// Returns the (row-major) inverse of the codomain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* rows of the mapping array.
    pub fn comb_codomain_inv( &self ) -> CombCodomainInv< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  {
        CombCodomainInv{ umatch: self }
    }  
    
    /// Returns the (row-major) domain COMB, indexed by `Mapping::ColIndex`
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn comb_domain( &self ) -> CombDomain< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  {
        CombDomain{ umatch: self }
    }

    /// Returns the (row-major) inverse of the domain COMB
    /// 
    /// # Design notes
    /// 
    /// This matrix cannot be indexed by integers without extra work, since we only assign ordinals to *matched* columns of the mapping array.
    pub fn comb_domain_inv( &self ) -> CombDomainInv< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >  {
        CombDomainInv{ umatch: self }
    }      

    /// Returns a reference to the matching array of the internally stored  U-match factorization.
    pub fn matching_ref( &self ) -> & GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient > { & self.matching }

    /// Returns a reference to the mapping array of the internally stored U-match factorization.
    pub fn mapping_ref( &self ) -> & Mapping { & self.mapping }    


    // pub fn packet_comb_codomain< 'a >( &'a self ) -> 
    //     MatrixAlgebraPacket< CombCodomain< 'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >
    //     where
    //         RingOperator:               Clone,
    //         OrderOperatorRowEntries:    Clone,
    //         OrderOperatorColEntries:    Clone,            
    // {
    //     MatrixAlgebraPacket{ matrix: self.comb_codomain(), ring: self.ring_operator(), row_entry_order: self.order_operator_major(), col_entry_order: self.order_operator_minor() }
    // }

    /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the (row-major) codomain COMB
    pub fn comb_codomain_packet( &self ) -> 
        MatrixAlgebraPacket< CombCodomain< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, RingOperator, OrderOperatorColEntries, OrderOperatorColEntries >
        where
            RingOperator:               Clone,
            OrderOperatorRowEntries:    Clone,
            OrderOperatorColEntries:    Clone,            
    {
        MatrixAlgebraPacket{ matrix: self.comb_codomain(), ring: self.ring_operator(), row_entry_order: self.order_operator_minor(), col_entry_order: self.order_operator_minor() }
    }

    /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the (row-major) inverse of the domain COMB
    pub fn comb_codomain_inv_packet( &self ) -> 
        MatrixAlgebraPacket< CombCodomainInv< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, RingOperator, OrderOperatorColEntries, OrderOperatorColEntries >
        where
            RingOperator:               Clone,
            OrderOperatorRowEntries:    Clone,
            OrderOperatorColEntries:    Clone,            
    {
        MatrixAlgebraPacket{ matrix: self.comb_codomain_inv(), ring: self.ring_operator(), row_entry_order: self.order_operator_minor(), col_entry_order: self.order_operator_minor() }
    }

    /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the (row-major) domain COMB
    pub fn comb_domain_packet( &self ) -> 
        MatrixAlgebraPacket< 
                CombDomain< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
                RingOperator, 
                OrderOperatorRowEntries, 
                OrderOperatorRowEntries
            >
        where
            RingOperator:               Clone,
            OrderOperatorRowEntries:    Clone,         
    {
        MatrixAlgebraPacket{ matrix: self.comb_domain(), ring: self.ring_operator(), row_entry_order: self.order_operator_major(), col_entry_order: self.order_operator_major() }
    }

    /// Returns a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket) for the inverse of the (row-major) domain COMB
    pub fn comb_domain_inv_packet( &self ) -> 
        MatrixAlgebraPacket< 
                CombDomainInv< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
                RingOperator, 
                OrderOperatorRowEntries, 
                OrderOperatorRowEntries,
            >
        where
            RingOperator:               Clone,
            OrderOperatorRowEntries:    Clone,          
    {
        MatrixAlgebraPacket{ matrix: self.comb_domain_inv(), ring: self.ring_operator(), row_entry_order: self.order_operator_major(), col_entry_order: self.order_operator_major() }
    }       

    /// Returns a reference to the matching array of the internally stored  U-match factorization, wrapped in a convenient convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket)
    pub fn matching_ref_packet( &self ) 
        -> MatrixAlgebraPacket< 
            & GeneralizedMatchingArrayWithMajorOrdinals< Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >,
            RingOperator, 
            OrderOperatorRowEntries, 
            OrderOperatorColEntries 
        >
        where
            RingOperator:               Clone,
            OrderOperatorRowEntries:    Clone,  
            OrderOperatorColEntries:    Clone,                      
    {
        MatrixAlgebraPacket{ matrix: self.matching_ref(), ring: self.ring_operator(), row_entry_order: self.order_operator_major(), col_entry_order: self.order_operator_minor() }
    } 

    /// Returns a reference to the mapping array of the internally stored U-match factorization, wrapped in a convenient convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket)
    pub fn mapping_ref_packet( &self ) 
        ->  MatrixAlgebraPacket< 
            & Mapping,
            RingOperator, 
            OrderOperatorRowEntries, 
            OrderOperatorColEntries 
        >
        where
            RingOperator:               Clone,
            OrderOperatorRowEntries:    Clone,  
            OrderOperatorColEntries:    Clone,                              
    {
        MatrixAlgebraPacket{ matrix: self.mapping_ref(), ring: self.ring_operator(), row_entry_order: self.order_operator_major(), col_entry_order: self.order_operator_minor() }
    }               


    /// The column submatrix of the mapping array indexed by matched column indices.
    /// 
    /// Returns a wrapper that implements `ViewRowAscend`.
    /// 
    /// # Design notes
    /// 
    /// It seems necessary to use the lifetime parameter 'b because (i) the matching array in the U-match decomposition
    /// is stored by value within each `Umatch` object, (ii) this is the object that contains the hashmap
    /// which we use to store the set of matched indices, however (iii) the object `OnlyKeyMinOutsideCollection`
    /// requires a *reference* to a hashmap, and (iv) said reference has to live long enough to do interesting things.  
    /// The most evident way to get that "long-enough" lifetime parameter, at present, is by passing the `Umatch`
    /// object by reference.
    pub fn mapping_matched_cols_only( &self ) -> OnlyKeyMinInsideCollection< &Mapping, &HashMap< Mapping::ColIndex, usize >, >
        // where 'a: 'b,
    {
        OnlyKeyMinInsideCollection::new( & self.mapping, self.matching.bimap_min_ref().val_to_ord_hashmap() )
    }    
    // fn mapping_matched_cols_only<'c>( &'c self ) -> OnlyKeyMinInsideCollection<'c,'c, Mapping, Mapping::ViewMajorAscend, HashMap< Mapping::ColIndex, usize >, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >
    //     where 'a: 'c, 'b: 'c,
    // {
    //     OnlyKeyMinInsideCollection::<'c, 'c,_,_,_,_,_,_>::new( self.mapping, self.matching.bimap_min().val_to_ord_hashmap() )
    // }    

    /// The column submatrix of the mapping array indexed by unmatched column indices.
    /// 
    /// Returns a wrapper that implements `ViewRowAscend`.
    /// 
    /// # Design notes
    /// 
    /// It seems necessary to use the lifetime parameter 'b because (i) the matching array in the U-match decomposition
    /// is stored by value within each `Umatch` object, (ii) this is the object that contains the hashmap
    /// which we use to store the set of matched indices, however (iii) the object `OnlyKeyMinOutsideCollection`
    /// requires a *reference* to a hashmap, and (iv) said reference has to live long enough to do interesting things.  
    /// The most evident way to get that "long-enough" lifetime parameter, at present, is by passing the `Umatch`
    /// object by reference.
    pub fn mapping_matchless_cols_only( &self ) -> OnlyKeyMinOutsideCollection< &Mapping, &HashMap< Mapping::ColIndex, usize >, >
    {
        OnlyKeyMinOutsideCollection::new( & self.mapping, self.matching.bimap_min_ref().val_to_ord_hashmap() )
    }    

    /// The row submatrix of the mapping array indexed by matched row indices.
    /// 
    /// Returns a wrapper that implements `ViewRowAscend`.
    /// 
    /// # Design notes
    /// 
    /// It seems necessary to use the lifetime parameter 'b because (i) the matching array in the U-match decomposition
    /// is stored by value within each `Umatch` object, (ii) this is the object that contains the hashmap
    /// which we use to store the set of matched indices, however (iii) the object `OnlyKeyMinOutsideCollection`
    /// requires a *reference* to a hashmap, and (iv) said reference has to live long enough to do interesting things.  
    /// The most evident way to get that "long-enough" lifetime parameter, at present, is by passing the `Umatch`
    /// object by reference.
    pub fn mapping_matched_rows_only( &self ) -> OnlyKeyMajInsideCollection< &Mapping, &HashMap< Mapping::RowIndex, usize >, >
        // where 'a: 'b,
    {
        OnlyKeyMajInsideCollection::new( & self.mapping, self.matching.bimap_maj_ref().val_to_ord_hashmap() )
    }    
    // fn mapping_matched_cols_only<'c>( &'c self ) -> OnlyKeyMinInsideCollection<'c,'c, Mapping, Mapping::ViewMajorAscend, HashMap< Mapping::ColIndex, usize >, Mapping::ColIndex, Mapping::RowIndex, Mapping::Coefficient >
    //     where 'a: 'c, 'b: 'c,
    // {
    //     OnlyKeyMinInsideCollection::<'c, 'c,_,_,_,_,_,_>::new( self.mapping, self.matching.bimap_min().val_to_ord_hashmap() )
    // }    

    /// The row submatrix of the mapping array indexed by unmatched row indices.
    /// 
    /// Returns a wrapper that implements `ViewRowAscend`.
    /// 
    /// # Design notes
    /// 
    /// It seems necessary to use the lifetime parameter 'b because (i) the matching array in the U-match decomposition
    /// is stored by value within each `Umatch` object, (ii) this is the object that contains the hashmap
    /// which we use to store the set of matched indices, however (iii) the object `OnlyKeyMinOutsideCollection`
    /// requires a *reference* to a hashmap, and (iv) said reference has to live long enough to do interesting things.  
    /// The most evident way to get that "long-enough" lifetime parameter, at present, is by passing the `Umatch`
    /// object by reference.
    pub fn mapping_matchless_rows_only( &self ) -> OnlyKeyMajOutsideCollection< &Mapping, &HashMap< Mapping::RowIndex, usize >, >
    {
        OnlyKeyMajOutsideCollection::new( & self.mapping, self.matching.bimap_maj_ref().val_to_ord_hashmap() )
    }    

    /// Returns a reference to the internally stored compressed representation of the inverse of the codomain COMB;
    /// this representation consists of a `VecOfVec` which encodes the off-diagonal entries of the square submatrix of the 
    /// inverse of the codomain COMB which is indexed by matched (i.e. pivot) indices.
    /// 
    /// Concretely, the represented matrix is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal_ref( &self ) 
        -> 
        & MatrixBimajorData<
                VecOfVec< usize, Mapping::Coefficient >,
                VecOfVec< usize, Mapping::Coefficient >,                
            >
        { & self.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal }   

        
    /// Returns a nested double reference to the internally stored compressed representation of the inverse of the codomain COMB;
    /// this representation consists of a `VecOfVec` which encodes the off-diagonal entries of the square submatrix of the 
    /// inverse of the codomain COMB which is indexed by matched (i.e. pivot) indices.
    /// 
    /// Concretely, the represented matrix is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).    
    /// 
    /// # Design note
    /// 
    /// This function exists because many parts of the `OAT` library use *references* to objects that implement
    /// matrix oracle traits.  A `VecOfVec` simple does not implement oracle traits in general, but a reference
    /// `& VecOfVec` does.  Therefore we often need to work with objects of form `&'b &'b VecOfVec`.  In
    /// practice, we find that Rust is prone to inferring the wrong lifetime if we simply write 
    /// `& self.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal_ref()` (for example, one finds errors alluding to
    /// dropped temprorary values).  This function has succeeded in sidestepping such errors in the past; please 
    /// let us know if it fails to do so successfully in future examples.
    


    /// The square, block submatrix of the inverse of the codomain COMB which is indexed by matched row indices.
    /// 
    /// Concretely, this object is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal< 'a >( &'a self ) 
        -> 
        PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend<
                &'a  MatrixBimajorData<
                            VecOfVec< usize, Mapping::Coefficient >,
                            VecOfVec< usize, Mapping::Coefficient >,                
                        >
            >
        where
            RingOperator:   Semiring< Mapping::Coefficient >,
        {   

            let prepended //: PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< &'a VecOfVec<usize, Mapping::Coefficient> >
            = PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( 
                        self.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal_ref(), 
                        RingOperator::one() 
                    );  
            
            prepended
        }  
        
    /// The square, block submatrix of the inverse of the codomain COMB which is indexed by matched row indices.
    /// 
    /// Concretely, this object is obtained from the inverse of the codomain COMB by (i) restricting to the square
    /// submatrix indexed by the row-pivot indices of the mapping array, and (ii) deleting the diagonal elements
    /// (each of which is equal to 1).
    pub fn comb_codomain_inv_matched_block_indexed_by_keymaj( &self ) 
        -> 
        ReindexSquareMatrix< 
                // PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend<
                //         &VecOfVec< usize, Mapping::Coefficient >,
                //     >,
                PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend<
                        & MatrixBimajorData<
                                VecOfVec< usize, Mapping::Coefficient >,
                                VecOfVec< usize, Mapping::Coefficient >,                
                            >
                    >,                
                &Vec< Mapping::RowIndex >,
                &HashMap< Mapping::RowIndex, usize >,
                usize,
                Mapping::RowIndex,
                Mapping::EntryMinor,
            >
        where
            RingOperator:   Semiring< Mapping::Coefficient >,
        {   

            let matrix_integer_indexed = self.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal();
            ReindexSquareMatrix::new(
                    matrix_integer_indexed,
                    self.matching.bimap_maj_ref().ord_to_val_vec(),
                    self.matching.bimap_maj_ref().val_to_ord_hashmap(),                    
                )
        }          


    /// Returns a matrix with rows indexed by `RowIndex` and columns indexed by `ColIndex`
    pub fn comb_codomain_inv_times_mapping_matched_block( &self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlock< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
        where 
            RingOperator: Clone,
            OrderOperatorRowEntries: Clone,
            OrderOperatorColEntries: Clone,           
    {
        CombCodomainInvTimesMappingMatchedBlock{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }


    pub fn comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_ordmaj( &self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
        where 
            RingOperator: Clone,
            OrderOperatorRowEntries: Clone,
            // OrderOperatorColEntries: Clone,
    {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByOrdMaj{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }    

    pub fn comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin( &self ) 
        -> 
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin< '_, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries > 
        where 
            RingOperator: Clone,
            OrderOperatorRowEntries: Clone,
            // OrderOperatorColEntries: Clone,
    {
        CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin{ umatch_ref: self } // CombCodomainInvTimesMappingMatchedBlock::new(self )
    }

    /// Solve `Rx = b`, where `R` is the codomain COMB.
    /// 
    /// Solution returns entries in strictly descending order.
    /// 
    /// `b` must iterate over entries in strictly descending order.
    pub fn solve_rx_equals_b< I >( &self, b: I ) 
        -> 
        // Vec< Mapping::EntryMinor >
        TriangularSolverMinorDescend<
                I, 
                CombCodomain< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
                Mapping::RowIndex, 
                RingOperator, 
                OrderOperatorColEntries
            >
        where 
            Mapping:               ViewColDescend,
            Mapping::Coefficient:       Clone,
            I:                          IntoIterator<Item=Mapping::EntryMinor>,
            Mapping::EntryMinor:    
                KeyValGet< Mapping::RowIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::RowIndex, Mapping::Coefficient > +
                Debug,
            Mapping::EntryMajor: 
                KeyValGet< Mapping::ColIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::ColIndex, Mapping::Coefficient > +
                Clone,
            RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            OrderOperatorRowEntries: 
                Clone + JudgePartialOrder<Mapping::EntryMajor>,
            OrderOperatorColEntries:
                Clone + JudgePartialOrder<Mapping::EntryMinor>,
            Mapping::Coefficient:       Clone + Debug,
            Mapping::ColIndex:       Clone + Debug,
            Mapping::RowIndex:       Clone + Debug,                              
    {
        TriangularSolverMinorDescend::<_,_,Mapping::RowIndex,_,_>::solve( 
                b, 
                self.comb_codomain(), 
                self.ring_operator(), 
                self.order_operator_minor()  
            )
    }




    /// Solve `xC = b`, where `C` is the domain COMB.
    /// 
    /// Solution returns entries in strictly ascending order.
    /// 
    /// `b` must iterate over entries in strictly ascending order.
    pub fn solve_xc_equals_b< I >( &self, b: I ) 
        -> 
        // Vec< Mapping::EntryMinor >
        TriangularSolverMajorAscend<
                I, 
                CombDomain< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
                Mapping::ColIndex, 
                RingOperator, 
                OrderOperatorRowEntries
            >
        where 
            Mapping:                    ViewRowAscend,
            Mapping::Coefficient:      Clone,
            I:                          IntoIterator<Item=Mapping::EntryMajor>,
            Mapping::EntryMinor:    
                KeyValGet< Mapping::RowIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::RowIndex, Mapping::Coefficient > +
                Debug,
            Mapping::EntryMajor: 
                KeyValGet< Mapping::ColIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::ColIndex, Mapping::Coefficient > +
                Clone,
            RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            OrderOperatorRowEntries: 
                Clone + JudgePartialOrder<Mapping::EntryMajor>,
            OrderOperatorColEntries:
                Clone + JudgePartialOrder<Mapping::EntryMinor>,
            Mapping::Coefficient:       Clone + Debug,
            Mapping::ColIndex:       Clone + Debug,
            Mapping::RowIndex:       Clone + Debug,                              
    {
        TriangularSolverMajorAscend::<_,_,Mapping::ColIndex,_,_>::solve( 
                b, 
                self.comb_domain(), 
                self.ring_operator(), 
                self.order_operator_major()  
            )
    }



    /// Return `x = bC` by solving `x * Inverse(C) = b` where `C` is the domain COMB.
    /// 
    /// Solution returns entries in strictly ascending order.
    /// 
    /// `b` must iterate over entries in strictly ascending order.
    /// 
    /// **Note** It is probably more efficient than computing `bC` directly, because
    /// the rows of `C^{-1}` are easier to compute than the columns of `C`.
    pub fn solve_bc_ascend< I >( &self, b: I ) 
        -> 
        // Vec< Mapping::EntryMinor >
        TriangularSolverMajorAscend<
                I, 
                CombDomainInv< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
                Mapping::ColIndex, 
                RingOperator, 
                OrderOperatorRowEntries
            >
        where 
            Mapping:                    ViewRowAscend,
            Mapping::Coefficient:      Clone,
            I:                          IntoIterator<Item=Mapping::EntryMajor>,
            Mapping::EntryMinor:    
                KeyValGet< Mapping::RowIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::RowIndex, Mapping::Coefficient > +
                Debug,
            Mapping::EntryMajor: 
                KeyValGet< Mapping::ColIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::ColIndex, Mapping::Coefficient > +
                Clone,
            RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            OrderOperatorRowEntries: 
                Clone + JudgePartialOrder<Mapping::EntryMajor>,
            OrderOperatorColEntries:
                Clone + JudgePartialOrder<Mapping::EntryMinor>,
            Mapping::Coefficient:       Clone + Debug,
            Mapping::ColIndex:       Clone + Debug,
            Mapping::RowIndex:       Clone + Debug,                              
    {
        TriangularSolverMajorAscend::<_,_,Mapping::ColIndex,_,_>::solve( 
                b, 
                self.comb_domain_inv(), 
                self.ring_operator(), 
                self.order_operator_major()  
            )
    }




    /// Solve `Dx = b`, where `D` is the mapping matrix.
    /// 
    /// Returns `None` if there is no solution.
    /// 
    /// `b` must iterate over entries in strictly descending order.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;        
    /// use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;        
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// use itertools::Itertools;
    /// 
    /// // DEFINE THE MATRIX
    /// // ===============================
    /// let matrix          =   VecOfVec::new( 
    ///                                         vec![   
    ///                                                     vec![(0,true), (1,true), (2,true)],
    ///                                                     vec![                            ], 
    ///                                                     vec![                    (2,true)], 
    ///                                         ] 
    ///                                     );
    ///                                 
    /// // COMPUTE U-MATCH
    /// // ===============================
    ///                                 
    /// let umatch
    ///     =   Umatch::factor(
    ///             & matrix,  // the matrix we wish to factor
    ///             (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
    ///             BooleanFieldOperator::new(), // the operator for the coefficient ring
    ///             OrderOperatorAuto, // order comparator for the entries in each row
    ///             OrderOperatorAuto, // order comparator for the entries in each column
    ///         );        
    /// 
    /// // SOLVE Ax = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (2,true), (0,true) ]; // note we list entries in reverse order
    /// let x   =   umatch.solve_dx_equals_b( b.clone() ).unwrap().collect_vec();
    /// let dx  =   umatch.multiply_dv(x);
    /// assert!( dx.eq( b ) );   
    /// ```
    pub fn solve_dx_equals_b< Vector >( &self, b: Vector ) 
        -> 
        Option< 
            Simplify<
                    HitMerge<
                            Scale< 
                                    CombDomainViewMinorDescend<
                                            '_, 
                                            Mapping, 
                                            RingOperator, 
                                            OrderOperatorRowEntries, 
                                            OrderOperatorColEntries,
                                        >, 
                                    Mapping::ColIndex, 
                                    RingOperator, 
                                    Mapping::Coefficient, 
                                >,
                            ReverseOrder< OrderOperatorRowEntries >,
                        >,
                    Mapping::ColIndex,
                    RingOperator,
                    Mapping::Coefficient,
                >
            >
        // TriangularSolverMinorDescend<
        //         Vector, 
        //         CombCodomain< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
        //         Mapping::RowIndex, 
        //         RingOperator, 
        //         OrderOperatorColEntries
        //     >
        where 
            Mapping:               ViewColDescend,
            Mapping::Coefficient:       Clone,
            Vector:                     IntoIterator<Item=Mapping::EntryMinor>,
            Mapping::EntryMinor:    
                KeyValGet< Mapping::RowIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::RowIndex, Mapping::Coefficient > +
                Clone + Debug,
            Mapping::EntryMajor: 
                KeyValGet< Mapping::ColIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::ColIndex, Mapping::Coefficient > +
                Clone,
            RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            OrderOperatorRowEntries: 
                Clone + JudgePartialOrder<Mapping::EntryMajor>,
            OrderOperatorColEntries:
                Clone + JudgePartialOrder<Mapping::EntryMinor>,
            Mapping::Coefficient:       Clone + Debug,
            Mapping::ColIndex:       Clone + Debug,
            Mapping::RowIndex:       Clone + Debug,                              
    {
        let matching    =   self.matching_ref();
        let comb_domain = self.comb_domain();
        let ring_operator = self.ring_operator();
        return self.solve_rx_equals_b(b)
            .map(
                |x|
                matching.keymaj_to_keymin( & x.key() ).map( 
                        |keymin| 
                        (keymin, ring_operator.invert(x.val()) ) 
                    )
                // match matching.keymaj_to_keymin( x.key() ) {
                //     Some( keymin ) => Some( ( keymin.clone(), x.val()  ) )
                // }
            )
            .collect::< Option<Vec<_>> >()
            .map(
                    |x|
                    x.multiply_matrix(
                        |i| comb_domain.view_minor_descend(i), 
                        self.ring_operator(), 
                        self.order_operator_major_reverse() 
                    )
                )
    }




    /// Solve `xD = b`, where `D` is the mapping matrix.
    /// 
    /// # Arguments
    /// 
    /// `b` must iterate over entries in strictly descending order.
    /// 
    /// # Returns
    /// 
    ///  `None` if there is no solution, otherwise `Some( iter )`, where `iter` is an iterator that runs over the entries of a solution in strictly ascendnig order
    /// 
    /// # Solution strategry
    /// 
    /// We solve `yM Cinv = b` for `y` and set `x = y Rinv`.  Then `xD = y Rinv D = y M Cinv = b`, where `RM = DC` is the U-match factorization.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// use oat_rust::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
    /// use oat_rust::utilities::order::OrderOperatorAuto;        
    /// 
    /// // define the matrix
    /// // -----------------
    /// let d = VecOfVec::new(
    ///         vec![
    ///                                 vec![  (0,true), (1,true),           ],
    ///                                 vec![            (1,true), (2,true), ],
    ///                             ]
    ///     );
    /// 
    /// // obtain a u-match factorization
    /// // ------------------------------
    /// let umatch  =   Umatch::factor( 
    ///     &d, 
    ///     0..2, 
    ///     BooleanFieldOperator::new(), 
    ///     OrderOperatorAuto, 
    ///     OrderOperatorAuto,             
    /// );
    /// 
    /// // try solving xd = b
    /// // ------------------
    /// 
    /// // Case 1: a solution exists; in this case we are gauaranteed to find one
    /// let x = umatch.solve_xd_equals_b( vec![ (0,true), (2,true), ] );        
    /// assert!( x.is_some() );
    /// assert!( x.unwrap().eq( vec![ (0,true), (1,true), ] ) );
    /// 
    /// // Case 2: no solution exists; in this case we get a certificate that no solution exists
    /// let x = umatch.solve_xd_equals_b( vec![ (0,true), (1,true), (2,true) ] );        
    /// assert!( x.is_none() );
    /// ```
    pub fn solve_xd_equals_b< Vector >( &self, b: Vector ) 
        -> 
        Option< 
            Simplify<
                    HitMerge<
                            Scale< 
                                    CombCodomainInvViewMajorAscend<
                                            Mapping, 
                                            RingOperator, 
                                        >, 
                                    Mapping::RowIndex, 
                                    RingOperator, 
                                    Mapping::Coefficient, 
                                >,
                            OrderOperatorColEntries,
                        >,
                    Mapping::RowIndex,
                    RingOperator,
                    Mapping::Coefficient,
                >
            >
        // TriangularSolverMinorDescend<
        //         Vector, 
        //         CombCodomain< Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
        //         Mapping::RowIndex, 
        //         RingOperator, 
        //         OrderOperatorColEntries
        //     >
        where 
            Mapping:                    ViewColDescend,
            Mapping::Coefficient:      Clone,
            Vector:                     IntoIterator<Item=Mapping::EntryMajor>,
            Mapping::EntryMinor:    
                KeyValGet< Mapping::RowIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::RowIndex, Mapping::Coefficient > +
                Clone + Debug,
            Mapping::EntryMajor: 
                KeyValGet< Mapping::ColIndex, Mapping::Coefficient > +
                KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                KeyValNew< Mapping::ColIndex, Mapping::Coefficient > +
                Clone,
            RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
            OrderOperatorRowEntries: 
                Clone + JudgePartialOrder<Mapping::EntryMajor>,
            OrderOperatorColEntries:
                Clone + JudgePartialOrder<Mapping::EntryMinor>,
            Mapping::Coefficient:      Clone + Debug,
            Mapping::ColIndex:          Clone + Debug,
            Mapping::RowIndex:          Clone + Debug,                              
    {

        let prod        =   self.matching_ref_packet()
                            .multiply_left( self.comb_domain_inv_packet() );
        let matching = self.matching_ref();
        let key_matching = |x| { matching.keymin_to_keymaj(&x) };
        // let key_matching    =   self.matching_ref().bijection_keymin_to_keymaj();
        EchelonSolverMajorAscendWithMajorKeys::solve(
                    b.into_iter(),
                    prod,
                    EvaluateFunctionFnMutWrapper::new( key_matching ),
                    self.ring_operator(),
                    self.order_operator_major(),
                )
                .solution()
                .map(|x| x.multiply_matrix_packet_major_ascend( self.comb_codomain_inv_packet() ) )
        // let mut tank = Vec::new();
        // while let Some(x) = y.next() { tank.push(x) }
        // match y.remainder().next().is_some() {
        //     true  => { return None }, // if there is nonzero remainder, then we have no solution
        //     false => {
        //         return Some( tank.multiply_matrix_packet_major_ascend( self.comb_codomain_inv_packet() ) )
        //     }
        // }

    }




    
    
    /// A basis for the kernel of the mapping matrix ("mapping matrix" is a term from U-match factorization)
    /// 
    /// # Inputs
    /// 
    /// `column_indices` is an iterator that runs over the columns of the mapping matrix, without repeats
    /// 
    /// # Returns
    /// 
    /// An iterator that wraps around `column_indices`.  The wrapper filters out every element
    /// of `column_indices` that indexes a nonzero column of the matching matrix; every other
    /// element is mapped to the corresponding column of the domain COMB.
    /// 
    /// # Why this works
    /// 
    /// Every vector returned by this iterator lies in the kernel of the mapping matrix, by
    /// the identity `RM = DC`.  The set of iterators is linearly independent because they form
    /// a subset of the columns of an invertible matrix.  The number of vectors equals the nullity
    /// of the matrix, again by the identity `RM = DC`.  Therefore the iterator runs over a basis for the kernel.
    pub fn kernel< ColumnIndices >( &self, column_indices: ColumnIndices ) 
        -> 
        ViewsMinorDescend<
                CombDomain
                    < Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >,
                FilterOutMembers
                    < ColumnIndices::IntoIter, & HashMap< Mapping::ColIndex, usize > >,
            >
    where   
        ColumnIndices:              IntoIterator< Item = Mapping::ColIndex >,
        Mapping:                    IndicesAndCoefficients + ViewColDescend + ViewRowAscend,
        Mapping::ColIndex:          Clone + Hash + std::cmp::Eq + Debug, 
        Mapping::RowIndex:          Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        Mapping::Coefficient:      Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        Mapping::ViewMajorAscend:   Clone, // !!! remove clone
        Mapping::EntryMinor:        Clone + Debug + 
                                    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + 
                                    KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                                    KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,
        Mapping::EntryMajor:        Clone + 
                                    KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + 
                                    KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                                    KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        OrderOperatorRowEntries:    Clone + JudgePartialOrder <  Mapping::EntryMajor >, // !!! remove clone
        OrderOperatorColEntries:    Clone + JudgePartialOrder <  Mapping::EntryMinor >,            

    {
        self.comb_domain().views_minor_descend( self.matching.filter_out_matched_minors( column_indices ) )
    }  



    /// A basis for the kernel of the mapping matrix
    /// 
    /// 
    /// # Returns
    /// 
    /// An iterator that returns a subset of the columns of the codomain COMB, which form a basis for the image.
    /// 
    /// Specifically, the iterator returns `R[:,i]` for all `i` such that `M[i,:]` is nonzero, where `R`
    /// is the codomain COMB and `M` is the generalized matching matrix.
    /// 
    /// # Why this works
    /// 
    /// Every vector returned by this iterator lies in the image of the mapping matrix, by
    /// the identity `RM = DC`.  The set of iterators in linearly independent because they form
    /// a subset of the columns of an invertible matrix.  The number of vectors equals the rank
    /// of the matrix, again by the identity `RM = DC`.  Therefore the iterator runs over a basis for the image.
    pub fn image( &self ) 
        -> 
        ViewsMinorDescend<
                CombCodomain
                    < Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >,
                Cloned< std::slice::Iter< Mapping::RowIndex > >,
            >
    where   
        Mapping:                    IndicesAndCoefficients + ViewColDescend + ViewRowAscend,
        Mapping::ColIndex:          Clone + Hash + std::cmp::Eq + Debug, 
        Mapping::RowIndex:          Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
        Mapping::Coefficient:      Clone + Debug,            // !!! remove Debug eventually       
        RingOperator:               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,
        Mapping::ViewMajorAscend:   Clone, // !!! remove clone
        Mapping::EntryMinor:        Clone + Debug + 
                                    KeyValGet< Mapping::RowIndex, Mapping::Coefficient > + 
                                    KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                                    KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,
        Mapping::EntryMajor:        Clone + 
                                    KeyValGet< Mapping::ColIndex, Mapping::Coefficient > + 
                                    KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                                    KeyValNew< Mapping::ColIndex, Mapping::Coefficient >,
        OrderOperatorRowEntries:    Clone + JudgePartialOrder <  Mapping::EntryMajor >, // !!! remove clone
        OrderOperatorColEntries:    Clone + JudgePartialOrder <  Mapping::EntryMinor >,            

    {
        self.comb_codomain().views_minor_descend( self.matching.bimap_maj_ref().ord_to_val_vec().iter().cloned() )
    }      

}





    //  MATRIX-VECTOR MULTIPLICATION
    //  ---------------------------------------------------------------------------------------------------------


impl < Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >  

    Umatch 
    < Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries, >  
    
    where   
        Mapping:                        IndicesAndCoefficients + ViewRowAscend,
        Mapping::ViewMajorAscend:       IntoIterator,     
        Mapping::EntryMinor:    
                                        KeyValGet< Mapping::RowIndex, Mapping::Coefficient > +
                                        KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + 
                                        KeyValNew< Mapping::RowIndex, Mapping::Coefficient > +
                                        Clone + Debug,
        Mapping::EntryMajor: 
                                        KeyValGet< Mapping::ColIndex, Mapping::Coefficient > +
                                        KeyValSet< Mapping::ColIndex, Mapping::Coefficient > + 
                                        KeyValNew< Mapping::ColIndex, Mapping::Coefficient > +
                                        Clone,
        Mapping::Coefficient:          Clone + Debug,
        Mapping::ColIndex:              Clone + Debug + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        Mapping::RowIndex:              Clone + Debug + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        RingOperator:                   Clone + 
                                        Semiring< Mapping::Coefficient > + 
                                        Ring< Mapping::Coefficient > + 
                                        DivisionRing< Mapping::Coefficient >,
        OrderOperatorRowEntries: 
                                        Clone + JudgePartialOrder<Mapping::EntryMajor>,
        OrderOperatorColEntries:
                                        Clone + JudgePartialOrder<Mapping::EntryMinor>,        
{

    /// Calculate the product `Dx`, where `D` is the mapping matrix.
    pub fn multiply_vd< Vector, VectorEntry >( & self, v: Vector ) 
        -> 
        Simplify<
                HitMerge<
                        Scale< 
                                Mapping::ViewMajorAscendIntoIter, 
                                Mapping::ColIndex, 
                                RingOperator, 
                                Mapping::Coefficient, 
                            >,
                        OrderOperatorRowEntries,
                    >,
                Mapping::ColIndex,
                RingOperator,
                Mapping::Coefficient,
            >
        where 
            Vector:                         IntoIterator<Item=VectorEntry>,
            VectorEntry:                    KeyValGet< Mapping::RowIndex, Mapping::Coefficient >,        
                             
    {
        let matrix = |i| self.mapping.view_major_ascend(i);
        v.multiply_matrix( matrix, self.ring_operator(), self.order_operator_major() )
    }    

    /// Calculate the product `xD`, where `D` is the mapping matrix.
    pub fn multiply_dv< Vector, VectorEntry >( &self, v: Vector ) 
        -> 
        Simplify<
                HitMerge<
                        Scale< 
                                Mapping::ViewMinorDescendIntoIter, 
                                Mapping::RowIndex, 
                                RingOperator, 
                                Mapping::Coefficient, 
                            >,
                        ReverseOrder< OrderOperatorColEntries >,
                    >,
                Mapping::RowIndex,
                RingOperator,
                Mapping::Coefficient,
            >
        where 
            Mapping:                        ViewColDescend,
            Vector:                         IntoIterator<Item=VectorEntry>,
            VectorEntry:                    KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,                             
    {
        let matrix = |i| self.mapping.view_minor_descend(i);
        v.multiply_matrix( matrix, self.ring_operator(), self.order_operator_minor_reverse() )
    }



    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    //  RETURN WHEN POSSIBLE -- THE PROBLEM IS LINE 2254
    //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    // /// Calculate the product `xC`, where `C` is the domain COMB
    // pub fn multiply_vc< Vector, VectorEntry >( & self, v: Vector ) 
    //     -> 
    //     Simplify<
    //             HitMerge<
    //                     Scale< 
    //                             CombDomainViewMajorAscend<
    //                                     Mapping, 
    //                                     RingOperator, 
    //                                     OrderOperatorRowEntries, 
    //                                 >,
    //                             Mapping::ColIndex, 
    //                             RingOperator, 
    //                             Mapping::Coefficient, 
    //                         >,
    //                     OrderOperatorRowEntries,
    //                 >,
    //             Mapping::ColIndex,
    //             RingOperator,
    //             Mapping::Coefficient,
    //         >
    //     where 
    //         Vector:                         IntoIterator<Item=VectorEntry>,
    //         VectorEntry:                    KeyValGet< Mapping::RowIndex, Mapping::Coefficient >,                            
    // {
    //     let comb_domain = self.comb_domain();        
    //     let matrix = |i| comb_domain.view_major_ascend(i);
    //     return v.multiply_matrix( matrix, self.ring_operator(), self.order_operator_major() )
    // }    



    /// Calculate the product `xD`, where `D` is the mapping matrix.
    pub fn multiply_cv< Vector, VectorEntry >( &self, v: Vector ) 
        -> 
        Simplify<
                HitMerge<
                        Scale< 
                                CombDomainViewMinorDescend<
                                        Mapping, 
                                        RingOperator, 
                                        OrderOperatorRowEntries, 
                                        OrderOperatorColEntries,
                                    >, 
                                Mapping::ColIndex, 
                                RingOperator, 
                                Mapping::Coefficient, 
                            >,
                        ReverseOrder< OrderOperatorRowEntries >,
                    >,
                Mapping::ColIndex,
                RingOperator,
                Mapping::Coefficient,
            >
        where 
            Vector:         IntoIterator<Item=VectorEntry>,
            VectorEntry:    KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
            Mapping:        ViewColDescend,                             
    {
        let comb_domain = self.comb_domain();
        let matrix = |i| comb_domain.view_minor_descend(i);
        v.multiply_matrix( matrix, self.ring_operator(), self.order_operator_major_reverse() )
    }


    /// Calculate the product `xD`, where `D` is the mapping matrix.
    pub fn multiply_vc< Vector, VectorEntry >( &self, v: Vector ) 
        -> 
        Simplify<
                HitMerge<
                        Scale< 
                                CombDomainViewMajorAscend<
                                        Mapping, 
                                        RingOperator, 
                                        OrderOperatorRowEntries, 
                                    >, 
                                Mapping::ColIndex, 
                                RingOperator, 
                                Mapping::Coefficient, 
                            >,
                        ReverseOrder< OrderOperatorRowEntries >,
                    >,
                Mapping::ColIndex,
                RingOperator,
                Mapping::Coefficient,
            >
        where 
            Vector:         IntoIterator<Item=VectorEntry>,
            VectorEntry:    KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,                          
    {
        let comb_domain = self.comb_domain();
        let matrix = |i| comb_domain.view_major_ascend(i);
        v.multiply_matrix( matrix, self.ring_operator(), self.order_operator_major_reverse() )
    }



}


//  ---------------------------------------------------------------------------------------------------------
//  U-MATCH -- IMPLEMENTATIONS FOR KEYMIN=KEMAJ
//  ---------------------------------------------------------------------------------------------------------


impl < Mapping, RingOperator, OrderOperatorRowEntries, KeyBoth, EntryBoth >  

    Umatch 
    < Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorRowEntries, >  
    
    where   
        Mapping::ColIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct
        Mapping::RowIndex:                   Clone + Hash + std::cmp::Eq, // required by the `GeneralizedMatchingArrayWithMajorOrdinals` struct 
        Mapping:                           ViewRowAscend<EntryMajor = EntryBoth> + 
                                                IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,
        Mapping::ViewMajorAscend:          IntoIterator,
        Mapping::EntryMajor:     KeyValGet< Mapping::ColIndex, Mapping::Coefficient >,
        // OrderOperatorByKey<usize, Mapping::Coefficient, (usize, Mapping::Coefficient)>: JudgePartialOrder< (usize, Mapping::Coefficient)>
{
    
     
}































//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;



    #[test]
    fn doc_test_umatchrowmajor_comprehensive_tiny() {

        // import packages
        use crate::algebra::matrices::operations::umatch::row_major::doc_test_drafts::BooleanFieldOperator;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        use crate::utilities::order::{OrderOperatorAuto};
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::operations::multiply::ProductMatrix;
        use crate::algebra::matrices::query::{ViewRowAscend};
        
        use itertools::Itertools;

        // define the coefficient ring
        let ring_operator                   =   BooleanFieldOperator::new();

        // define the matrix we wish to factor
        let num_indices_major           =   1;
        let num_indices_minor           =   1;
        let mapping_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,true),], 
                            ] );
        let mapping                               =   & mapping_data;

        // compute the U-match factorization
        let umatch
            =   Umatch::factor(
                    mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator, 
                    OrderOperatorAuto::new(), 
                    OrderOperatorAuto::new(),
                );
        
        // extract R, R^{-1}, C, C^{-1}, and M
        let matching = umatch.matching_ref();
        let comb_codomain = umatch.comb_codomain();
        let comb_codomain_inv = umatch.comb_codomain_inv();        
        let comb_domain = umatch.comb_domain();        
        let comb_domain_inv = umatch.comb_domain_inv(); 

        // get references to R, R^{-1}, C, C^{-1}, and M        
        let comb_codomain_ref         =   & comb_codomain;
        let comb_codomain_inv_ref         =   & comb_codomain_inv;
        let comb_domain_ref         =   & comb_domain;
        let comb_domain_inv_ref         =   & comb_domain_inv;   
        
        // compute some products
        let product_domain = ProductMatrix::new( comb_domain_ref, comb_domain_inv_ref, ring_operator, OrderOperatorAuto );
        let product_codomain = ProductMatrix::new( comb_codomain_ref, comb_codomain_inv_ref, ring_operator, OrderOperatorAuto );        
        let product_codomain_comb_inv_times_mapping = ProductMatrix::new( comb_codomain_inv_ref, mapping, ring_operator, OrderOperatorAuto );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrix::new( product_codomain_comb_inv_times_mapping, comb_domain_ref, ring_operator, OrderOperatorAuto );                


        // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, true) ]
            ) 
        }

        // check that the product of the codomain COMB with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, true) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }     

    }



    #[test]
    fn doc_test_umatchrowmajor_comprehensive_tiny_waist() {

        // import packages
        use crate::algebra::matrices::operations::umatch::row_major::doc_test_drafts::BooleanFieldOperator;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        use crate::utilities::order::{OrderOperatorAuto};
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::operations::multiply::ProductMatrix;
        use crate::algebra::matrices::query::{ViewRowAscend};
        
        use itertools::Itertools;

        // define the coefficient ring
        let ring_operator                   =   BooleanFieldOperator::new();

        // define the matrix we wish to factor
        let num_indices_major           =   2;
        let num_indices_minor           =   1;
        let mapping_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,true),], 
                                vec![(0,true),],                                 
                            ] );
        let mapping                               =   & mapping_data;

        // compute the U-match factorization
        let umatch
            =   Umatch::factor(
                    mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator, 
                    OrderOperatorAuto::new(), 
                    OrderOperatorAuto::new(),
                );
        
        // extract R, R^{-1}, C, C^{-1}, and M
        let matching = umatch.matching_ref();
        let comb_codomain = umatch.comb_codomain();
        let comb_codomain_inv = umatch.comb_codomain_inv();        
        let comb_domain = umatch.comb_domain();        
        let comb_domain_inv = umatch.comb_domain_inv(); 

        // get references to R, R^{-1}, C, C^{-1}, and M        
        let comb_codomain_ref         =   & comb_codomain;
        let comb_codomain_inv_ref         =   & comb_codomain_inv;
        let comb_domain_ref         =   & comb_domain;
        let comb_domain_inv_ref         =   & comb_domain_inv;   
        
        // compute some products
        let product_domain = ProductMatrix::new( comb_domain_ref, comb_domain_inv_ref, ring_operator, OrderOperatorAuto );
        let product_codomain = ProductMatrix::new( comb_codomain_ref, comb_codomain_inv_ref, ring_operator, OrderOperatorAuto );        
        let product_codomain_comb_inv_times_mapping = ProductMatrix::new( comb_codomain_inv_ref, mapping, ring_operator, OrderOperatorAuto );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrix::new( product_codomain_comb_inv_times_mapping, comb_domain_ref, ring_operator, OrderOperatorAuto );                


        // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, true) ]
            ) 
        }

        // check that the product of the codomain COMB with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, true) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }     

    }    





    #[test]
    fn doc_test_umatchrowmajor_comprehensive_tiny_height() {

        // import packages
        use crate::algebra::matrices::operations::umatch::row_major::doc_test_drafts::BooleanFieldOperator;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        use crate::utilities::order::{OrderOperatorAuto};
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::operations::multiply::ProductMatrix;
        use crate::algebra::matrices::query::{ViewRowAscend};
        
        use itertools::Itertools;

        // define the coefficient ring
        let ring_operator                   =   BooleanFieldOperator::new();

        // define the matrix we wish to factor
        let num_indices_major           =   1;
        let num_indices_minor           =   2;
        let mapping_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,true), (1,true)], 
                            ] );
        let mapping                               =   & mapping_data;

        // compute the U-match factorization
        let umatch
            =   Umatch::factor(
                    mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator, 
                    OrderOperatorAuto::new(), 
                    OrderOperatorAuto::new(),
                );
        
        // extract R, R^{-1}, C, C^{-1}, and M
        let matching = umatch.matching_ref();
        let comb_codomain = umatch.comb_codomain();
        let comb_codomain_inv = umatch.comb_codomain_inv();        
        let comb_domain = umatch.comb_domain();        
        let comb_domain_inv = umatch.comb_domain_inv(); 

        // get references to R, R^{-1}, C, C^{-1}, and M        
        let comb_codomain_ref         =   & comb_codomain;
        let comb_codomain_inv_ref         =   & comb_codomain_inv;
        let comb_domain_ref         =   & comb_domain;
        let comb_domain_inv_ref         =   & comb_domain_inv;   
        
        // compute some products
        let product_domain = ProductMatrix::new( comb_domain_ref, comb_domain_inv_ref, ring_operator, OrderOperatorAuto );
        let product_codomain = ProductMatrix::new( comb_codomain_ref, comb_codomain_inv_ref, ring_operator, OrderOperatorAuto );        
        let product_codomain_comb_inv_times_mapping = ProductMatrix::new( comb_codomain_inv_ref, mapping, ring_operator, OrderOperatorAuto );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrix::new( product_codomain_comb_inv_times_mapping, comb_domain_ref, ring_operator, OrderOperatorAuto );                


        // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, true) ]
            ) 
        }

        // check that the product of the codomain COMB with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, true) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }     

    }





    #[test]
    fn doc_test_umatchrowmajor_comprehensive_small() {

        // import packages
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        use crate::utilities::order::{OrderOperatorAuto};
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::operations::multiply::ProductMatrix;
        use crate::algebra::matrices::query::{ViewRowAscend};
        
        use itertools::Itertools;

        // define the coefficient ring
        let modulus                     =   5;
        let ring_operator                   =   PrimeOrderFieldOperator::new( modulus );        

        // define the matrix we wish to factor
        let num_indices_major           =   2;
        let num_indices_minor           =   3;
        let mapping_data              =   VecOfVec::new( 
            vec![   
                                vec![(0,1), (1,2), (2,3)], 
                                vec![              (2,1)]  
                            ] );
        let mapping                               =   & mapping_data;

        // compute the U-match factorization
        let umatch
            =   Umatch::factor(
                    mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator, 
                    OrderOperatorAuto::new(), 
                    OrderOperatorAuto::new(),
                );
        
        // extract R, R^{-1}, C, C^{-1}, and M
        let matching = umatch.matching_ref();
        let comb_codomain = umatch.comb_codomain();
        let comb_codomain_inv = umatch.comb_codomain_inv();        
        let comb_domain = umatch.comb_domain();        
        let comb_domain_inv = umatch.comb_domain_inv(); 

        // get references to R, R^{-1}, C, C^{-1}, and M        
        let comb_codomain_ref         =   & comb_codomain;
        let comb_codomain_inv_ref         =   & comb_codomain_inv;
        let comb_domain_ref         =   & comb_domain;
        let comb_domain_inv_ref         =   & comb_domain_inv;   
        
        // compute some products
        let product_domain = ProductMatrix::new( comb_domain_ref, comb_domain_inv_ref, ring_operator, OrderOperatorAuto );
        let product_codomain = ProductMatrix::new( comb_codomain_ref, comb_codomain_inv_ref, ring_operator, OrderOperatorAuto );        
        let product_codomain_comb_inv_times_mapping = ProductMatrix::new( comb_codomain_inv_ref, mapping, ring_operator, OrderOperatorAuto );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrix::new( product_codomain_comb_inv_times_mapping, comb_domain_ref, ring_operator, OrderOperatorAuto );                


        // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, 1) ]
            ) 
        }

        // check that the product of the codomain COMB with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, 1) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }     

    }


}    


//  ---------------------------------------------------------------------
//  Unit tests
//  ---------------------------------------------------------------------

#[cfg(test)]
mod unit_tests {
    use std::{iter::Cloned};

    use itertools::{assert_equal, Itertools};

    use crate::{algebra::{matrices::{display::print_indexed_major_views, operations::{multiply::vector_matrix_multiply_major_ascend_simplified}, }}, utilities::order::{is_sorted_strictly}};
    use crate::algebra::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
    
    use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
     
    use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;



    //  ===================================================================================
    //  CONSTRUCTION OF PIVOT BLOCK OF INVERSE OF THE CODOMAIN COMB
    //  ===================================================================================    

    

    /// This test targets the initial computation of the pivot block of the inverse of the codomain COMB.
    #[cfg(test)]
    fn test_initial_decomposition() {
        use crate::algebra::matrices::operations::umatch::row_major::get_codomain_comb_inv_off_diag_pivot_block;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::order::OrderOperatorAuto;

        let mapping       =   VecOfVec::new(
                                            vec![
                                                vec![ (0, 1), (1, 1) ],
                                                vec![ (0, 1), (1, 1) ],
                                                vec![ (0, 1), (1, 1) ],                                                
                                            ]
                                        );
        let iter_keymaj         = (0..2).rev();
        let ring_operator       =   PrimeOrderFieldOperator::new(5);
        
        let _codomain_comb_inv_pivot_block = get_codomain_comb_inv_off_diag_pivot_block(
                                            & (& mapping),
                                            iter_keymaj,
                                            ring_operator,
                                            OrderOperatorAuto,
                                            // OrderOperatorAuto,
                                        );       
                                        
        // println!("{:?}", & codomain_comb_inv_pivot_block.0);
        // println!("{:?}", & codomain_comb_inv_pivot_block.1);        
    }


    /// This test targets the initial computation of the pivot block of the inverse of the codomain COMB.    
    #[test]
    fn test_initial_decomposition_another_example() {
        use itertools::Itertools;

        use crate::algebra::matrices::operations::umatch::row_major::get_codomain_comb_inv_off_diag_pivot_block;
                
        use crate::algebra::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
        use crate::algebra::matrices::operations::transform_vector_wise::VecWiseTransformed;
        
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewRowAscend;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::order::{OrderOperatorByKey, OrderOperatorAuto};
        

        use std::slice::Iter;

        let matrix_size = 5;
        let modulus = 7;
        let mapping = VecOfVec::new(
                                    vec![
                                        vec![(0, 1), (1, 2),                 (4, 0)],
                                        vec![        (1, 1), (2, 0),         (4, 1)],
                                        vec![                (2, 1), (3, 0), (4, 0)],
                                        vec![                        (3, 1), (4, 0)],
                                        vec![                                (4, 1)],
                                    ]
                                );
        let array_mapping_ref = & mapping;
        let iter_keymaj         = (0 .. matrix_size).rev();
        let ring_operator       =   PrimeOrderFieldOperator::new( modulus );

        // compute the inverse
        let inverse =   InverseOfTriangularArrayLazyAscending::new( 
            & array_mapping_ref,
            ring_operator,
            OrderOperatorByKey::new(),
        );    
        
        // let `mapping_transformed` be the matrix obtained by reversing the order of columns of `mapping`
        let vector_transformer = |x: Cloned< Iter< '_, (usize,usize) > >| 
                                                                    x
                                                                        .rev()
                                                                        .map(   |(a,b)| 
                                                                                ( matrix_size - 1 - a, b) 
                                                                            )
                                                                        .collect_vec() ;

        let mapping_transformed =   VecWiseTransformed::new( 
                                            array_mapping_ref,
                                            vector_transformer,
                                        );

        // compute the codomain COMB of `mapping_transformed`
        //
        // NOTE: the codomain COMB of `mapping_transformed` is the inverse of `mapping`, where -- by contrast -- 
        //       the codomain COMB of `mapping` is the identity matrix (which is less interesting)
        let (array_comb, matching) = get_codomain_comb_inv_off_diag_pivot_block(
            & mapping_transformed, // the matrix mapping_transformed ought to implement Copy
            iter_keymaj,
            ring_operator,
            OrderOperatorAuto,
            // OrderOperatorAuto,
        );          


        for keymaj in 0 .. matrix_size {

            // obtain a row of the inverse
            let inv_row     =    inverse.view_major_ascend(keymaj).collect_vec();

            // obtain a row of the codomain COMB
            let ordmaj    =   matching.keymaj_to_ord( &keymaj ).unwrap();  

            let comb_off_diag_view =    (& array_comb)
                                        .view_major_ascend( ordmaj )
                                        // .iter()
                                        .map(   | (x,y)| 
                                                ( matching.ord_to_keymaj( x ), y ) // reindex the row from ordinals to major keys
                                            );
            let mut comb_view   =   comb_off_diag_view.collect_vec();
            comb_view.push( (keymaj, 1) );
            comb_view.sort();       
            
            // compute products
            let product_inv=   vector_matrix_multiply_major_ascend_simplified(
                                                                        inv_row.clone(),
                                                                        & mapping_transformed,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );
            let product_umatch=   vector_matrix_multiply_major_ascend_simplified(
                                                                        comb_view.clone(),
                                                                        & mapping_transformed,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );                                                                    
                                                                

            // compare
            // println!("MATRIX:");
            for k in 0 .. matrix_size { println!("{:?}", array_mapping_ref.view_major_ascend(k).collect_vec()  ) }
            // println!("MATRIX TRANSFORMED:");
            for k in 0 .. matrix_size { println!("{:?}", mapping_transformed.view_major_ascend(k)  ) }            
            // println!("product_inv: {:?}", product_inv.clone().collect_vec() );
            // println!("product_umatch: {:?}", product_umatch.clone().collect_vec());            
            // println!("COMB row len: {:?}, inv row len{:?}", comb_view.len(), inv_row.len());
            assert_eq!(comb_view, inv_row);
                            
        } 
        
        
    }    


    /// This test targets the initial computation of the pivot block of the inverse of the codomain COMB.
    ///
    /// The key idea of this test is the following fact: let M be a square upper-unitriangular matrix,
    /// and let N be the matrix obtained by reversing the order of columns of M.  Then the standard
    /// "cohomology algorithm," applied to N, produces a codomain COMB equal to M^{-1}.
    /// 
    /// This test applies the standard cohomology algorithm to compute a codomain COMB of N.  We 
    /// check to ensure that this codomain COMB equals M^{-1}.
    #[test]
    fn test_initial_decomposition_larger() {
        use itertools::Itertools;

        use crate::algebra::matrices::operations::umatch::row_major::get_codomain_comb_inv_off_diag_pivot_block;
                
        use crate::algebra::matrices::operations::invert::InverseOfTriangularArrayLazyAscending;
        use crate::algebra::matrices::operations::transform_vector_wise::VecWiseTransformed;
        use crate::algebra::matrices::query::ViewRowAscend;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::order::{OrderOperatorByKey, OrderOperatorAuto};
        

        use std::slice::Iter;

        let matrix_size     =   10;
        let modulus             =   7;
        
        let mapping = VecOfVec::random_mod_p_upper_unitriangular( matrix_size, modulus );
        let array_mapping_ref = & mapping;
        let iter_keymaj         = (0 .. matrix_size).rev();
        let ring_operator       =   PrimeOrderFieldOperator::new( modulus );

        // compute the inverse
        let inverse =   InverseOfTriangularArrayLazyAscending::new( 
            & array_mapping_ref,
            ring_operator,
            OrderOperatorByKey::new(),
        );    
        
        // let `mapping_transformed` be the matrix obtained by reversing the order of columns of `mapping`        
        let vector_transformer = |x: Cloned< Iter< '_, (usize,usize) > >| 
                                                                    x
                                                                        .rev()
                                                                        .map(   |(a,b)| 
                                                                                ( matrix_size - 1 - a, b) 
                                                                            )
                                                                        .collect_vec() ;

        let mapping_transformed  =   VecWiseTransformed::new( 
                                            array_mapping_ref,
                                            vector_transformer,
                                        );

        // compute the codomain COMB of `mapping_transformed`
        //
        // NOTE: the codomain COMB of `mapping_transformed` is the inverse of `mapping`, where -- by contrast -- 
        //       the codomain COMB of `mapping` is the identity matrix (which is less interesting)        
        let (array_comb, matching) = get_codomain_comb_inv_off_diag_pivot_block(
            & mapping_transformed,
            iter_keymaj,
            ring_operator,
            OrderOperatorAuto,
            // OrderOperatorAuto,
        );          


        for keymaj in 0 .. matrix_size {

            // obtain a row of the inverse
            let inv_row     =    inverse.view_major_ascend(keymaj).collect_vec();

            // obtain a row of the codomain COMB
            let ordmaj    =   matching.keymaj_to_ord( &keymaj ).unwrap();

            let comb_off_diag_view =    (& array_comb)
                                        .view_major_ascend( ordmaj )
                                        // .iter()
                                        .map(   | (x,y)| 
                                                ( matching.ord_to_keymaj( x ), y ) // reindex the row from ordinals to major keys
                                            );
            let mut comb_view   =   comb_off_diag_view.collect_vec();
            comb_view.push( (keymaj, 1) );
            comb_view.sort();       
            
            // compute products
            let product_inv=   vector_matrix_multiply_major_ascend_simplified(
                                                                        inv_row.clone(),
                                                                        & mapping_transformed,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );
            let product_umatch=   vector_matrix_multiply_major_ascend_simplified(
                                                                        comb_view.clone(),
                                                                        & mapping_transformed,
                                                                        ring_operator,
                                                                        OrderOperatorByKey::new()
                                                                    );                                                                    
                                                                

            // compare
            // println!("MATRIX:");
            for k in 0 .. matrix_size { println!("{:?}", array_mapping_ref.view_major_ascend(k).collect_vec()  ) }
            // println!("MATRIX TRANSFORMED:");
            for k in 0 .. matrix_size { println!("{:?}", mapping_transformed.view_major_ascend(k)  ) }            
            // println!("product_inv: {:?}", product_inv.clone().collect_vec() );
            // println!("product_umatch: {:?}", product_umatch.clone().collect_vec());            
            // println!("COMB row len: {:?}, inv row len{:?}", comb_view.len(), inv_row.len());
            assert_eq!(comb_view, inv_row);
                            
        } 
        
        
    }    



    //  ===================================================================================
    //  RECOVERY OF MAJOR VIEWS OF COMBS -- COMPREHENSIVE
    //  ===================================================================================    


    /// Checks that Umatch decomposition is correct (using a small example matrix, D) in the following sense:
    /// R^{-1} * R = I
    /// C^{-1} * C = I
    /// R^{-1} * D * C = M
    /// And the rows of R, R^{-1}, C, and C^{-1} appear in strictly ascending order
    #[test]
    fn test_umatchrowmajor_comprehensive_small() {

        
        use crate::utilities::order::{OrderOperatorAuto, OrderOperatorByKeyCutsom, };
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::operations::multiply::ProductMatrix;
        use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend};
        use crate::utilities::iterators::is_sorted::IsSortedBy;

        let num_indices_major           =   2;
        let num_indices_minor           =   3;
        // let approximate_density         =   0.3;
        let modulus                     =   3;
        // let allow_nonstructural_zero         =   true;

        let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
        let mapping_data       =   VecOfVec::new( 
            vec![   
                                vec![(0,1), (1,2), (2,0)], 
                                vec![              (2,0)]  
                            ] );
        let mapping                               =   & mapping_data;

        // !!! CONSIDER CHANGING OrderOperatorAuto TO INCLUDE PHANTOM DATA WHICH FIXES THE TYPE OF THE COMPARED OBJECTS

        
        let umatch: Umatch<
                        &VecOfVec<usize, usize>, 
                        PrimeOrderFieldOperator, 
                        OrderOperatorByKeyCutsom<usize, usize, (usize, usize), OrderOperatorAuto>, 
                        OrderOperatorByKeyCutsom<usize, usize, (usize, usize), OrderOperatorAuto>,
                    > 
            =   Umatch::factor(
                    mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator, 
                    OrderOperatorAuto, 
                    OrderOperatorAuto,
                );


        


        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
        let matching = umatch.matching_ref();

        let comb_codomain = umatch.comb_codomain();
        let comb_codomain_inv = umatch.comb_codomain_inv();        
        let comb_domain = umatch.comb_domain();        
        let comb_domain_inv = umatch.comb_domain_inv(); 

        let comb_codomain_ref         =   & comb_codomain;
        let comb_codomain_inv_ref         =   & comb_codomain_inv;
        let comb_domain_ref         =   & comb_domain;
        let comb_domain_inv_ref         =   & comb_domain_inv;                        
        
        
        let product_domain = ProductMatrix::new( comb_domain_ref, comb_domain_inv_ref, ring_operator, OrderOperatorAuto );
        let product_codomain = ProductMatrix::new( comb_codomain_ref, comb_codomain_inv_ref, ring_operator, OrderOperatorAuto );        
        let product_codomain_comb_inv_times_mapping = ProductMatrix::new( comb_codomain_inv_ref, mapping, ring_operator, OrderOperatorAuto );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrix::new( product_codomain_comb_inv_times_mapping, comb_domain_ref, ring_operator, OrderOperatorAuto );                


        // println!("mapping:");
        // print_indexed_major_views( & mapping, 0 .. num_indices_major );
        // println!("matching:");
        // print_indexed_major_views( & matching, 0 .. num_indices_major );        
        // println!("comb_domain:");
        // print_indexed_major_views( & comb_domain, 0 .. num_indices_minor );        
        // println!("comb_domain_inv:");
        // print_indexed_major_views( & comb_domain_inv, 0 .. num_indices_minor );     
        // println!("comb_codomain:");
        // print_indexed_major_views( & comb_codomain, 0 .. num_indices_major );        
        // println!("comb_codomain_inv:");
        // print_indexed_major_views( & comb_codomain_inv, 0 .. num_indices_major );    
        // println!("comb_codomain_inv * mapping * comb_domain:");
        // print_indexed_major_views( & product_codomain_comb_inv_times_mapping_times_domain_comb, 0 .. num_indices_major );                                
        for column_index in 0 .. num_indices_minor {
            // println!("{:?}", product_domain.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( product_domain.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }


        // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, 1) ]
            ) 
        }

        // check that the product of the codomain COMB with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, 1) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }     
        
        // check that rows are sorted in strictly ascending order
        for keymaj in 0 .. num_indices_major { 
            assert!(    mapping.view_major_ascend( keymaj             ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_codomain.view_major_ascend( keymaj       ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_codomain_inv.view_major_ascend( keymaj   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_domain.view_major_ascend( keymaj         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_domain_inv.view_major_ascend( keymaj     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }        
        
        // check that the major and minor views of the inverse of the codomain COMB agree
        let comb_codomain_vec_of_vec_simple     
            =   VecOfVec::from_iterable( (0..num_indices_major).map( |k| comb_codomain_inv.view_major_ascend(k) ) );
        for row_index in 0..num_indices_major {
            itertools::assert_equal( 
                    comb_codomain.view_minor_descend( row_index ),
                    (& comb_codomain_vec_of_vec_simple).view_minor_descend( row_index )
                )
        }

        // check that the major and minor views of `CombCodomainInv` agree
        let comb_codomain_inv_vec_of_vec_simple     
            =   VecOfVec::from_iterable( (0..num_indices_major).map( |k| comb_codomain_inv.view_major_ascend(k) ) );
        for row_index in 0..num_indices_major {
            // println!("PRINTING HERE: see below");
            // println!("ROW INDEX = {:?}", row_index );
            // println!("{:?}", comb_codomain_inv.view_minor_descend( row_index ).collect_vec());
            // println!("{:?}", (& comb_codomain_inv_vec_of_vec_simple).view_minor_descend( row_index ).collect_vec());            
            assert_equal( 
                    comb_codomain_inv.view_minor_descend( row_index ).collect_vec(),
                    (& comb_codomain_inv_vec_of_vec_simple).view_minor_descend( row_index ).collect_vec()
                )
        }

    }



    /// Checks that Umatch decomposition is correct (using a random example matrix, D) in the following sense:
    /// R^{-1} * R = I
    /// C^{-1} * C = I
    /// R^{-1} * D * C = M   
    /// And the rows of R, R^{-1}, C, and C^{-1} appear in strictly ascending order 
    #[test]
    fn test_umatchrowmajor_comprehensive() {

        use crate::utilities::order::{OrderOperatorAuto, };
        
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::operations::multiply::ProductMatrix;
        use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend};
        use crate::utilities::iterators::is_sorted::IsSortedBy;        

        let num_indices_major           =   10;
        let num_indices_minor           =   20;
        let approximate_density           =   0.2;
        let modulus                     =   17;
        let allow_nonstructural_zero     =   true;

        let ring_operator           =   PrimeOrderFieldOperator::new( modulus );
        let mapping_data       =   VecOfVec::random_mod_p_with_density( num_indices_major, num_indices_minor, approximate_density, modulus, allow_nonstructural_zero );
        let mapping                               =   & mapping_data;

        let umatch 
            =   Umatch::factor( 
                    mapping, 
                    (0..num_indices_major).rev(), 
                    ring_operator, 
                    OrderOperatorAuto, 
                    OrderOperatorAuto,
                );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch );
        let matching = umatch.matching_ref();

        let comb_codomain = umatch.comb_codomain();
        let comb_codomain_inv = umatch.comb_codomain_inv();        
        let comb_domain = umatch.comb_domain();        
        let comb_domain_inv = umatch.comb_domain_inv();  
        let comb_codomain_inv_times_mapping_matched_block = umatch.comb_codomain_inv_times_mapping_matched_block();  
        let comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin = umatch.comb_codomain_inv_times_mapping_matched_block_with_rows_indexed_by_matched_keymin();

        let comb_codomain_ref         =   & comb_codomain;
        let comb_codomain_inv_ref         =   & comb_codomain_inv;
        let comb_domain_ref         =   & comb_domain;
        let comb_domain_inv_ref         =   & comb_domain_inv;            
        let _comb_codomain_inv_times_mapping_matched_block_ref     =   & comb_codomain_inv_times_mapping_matched_block;
        let comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref = & comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin;
        

        let product_domain = ProductMatrix::new( comb_domain_ref, comb_domain_inv_ref, ring_operator, OrderOperatorAuto );
        let product_codomain = ProductMatrix::new( comb_codomain_ref, comb_codomain_inv_ref, ring_operator, OrderOperatorAuto );        
        let product_codomain_comb_inv_times_mapping = ProductMatrix::new( comb_codomain_inv_ref, mapping, ring_operator, OrderOperatorAuto );      
        let product_codomain_comb_inv_times_mapping_times_domain_comb = ProductMatrix::new( product_codomain_comb_inv_times_mapping, comb_domain_ref, ring_operator, OrderOperatorAuto );                        


        // println!("mapping:");
        // print_indexed_major_views( & mapping, 0 .. num_indices_major );
        // println!("matching:");
        // print_indexed_major_views( & matching, 0 .. num_indices_major );        
        // println!("comb_domain:");
        // print_indexed_major_views( & comb_domain, 0 .. num_indices_minor );        
        // println!("comb_domain_inv:");
        // print_indexed_major_views( & comb_domain_inv, 0 .. num_indices_minor );     
        // println!("comb_codomain:");
        // print_indexed_major_views( & comb_codomain, 0 .. num_indices_major );        
        // println!("comb_codomain_inv:");
        // print_indexed_major_views( & comb_codomain_inv, 0 .. num_indices_major );                        
        // println!("comb_codomain_inv * mapping * comb_domain:");
        // print_indexed_major_views( & product_codomain_comb_inv_times_mapping_times_domain_comb, 0 .. num_indices_major );                                        
        for column_index in 0 .. num_indices_minor {
            // println!("{:?}", product_domain.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( product_domain.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }



        // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
        for keymin in 0 .. num_indices_minor { 
            assert_eq!(
                product_domain.view_major_ascend( keymin ).collect_vec(),
                vec![ (keymin, 1) ]
            ) 
        }

        // check that the product of the codomain COMB with its inverse is identity R * R^{-1} = I
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain.view_major_ascend( keymaj ).collect_vec(),
                vec![ (keymaj, 1) ]
            ) 
        }    
        
        // check the factorization R^{-1} * D * C = M
        for keymaj in 0 .. num_indices_major { 
            assert_eq!(
                product_codomain_comb_inv_times_mapping_times_domain_comb.view_major_ascend( keymaj ).collect_vec(),
                matching.view_major_ascend( keymaj ).collect_vec()
            ) 
        }    

        // check that the major and minor views of `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` agree (i.e. that, taken all together, they run over the same entries)     
        verify_viewmajorascend_compatible_with_viewminordescend(
                comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref,
                umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned(),
                umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned(),                
            );
        
        // check that the major and minor views of `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` are sorted
        for keymin in umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned() {
            assert!(    comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_minor_descend( keymin ).is_sorted_by( |x, y| x.0 > y.0 )     );
        }
        for keymin in umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned() {
            assert!(    comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_major_ascend( keymin ).is_sorted_by( |x, y| x.0 < y.0 )     );
        }     



        // ----------------------------------------------------------------------------------------------------------------
        // check that the major and minor views of `CombCodomain` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                comb_codomain_ref,
                0..num_indices_major,
                0..num_indices_major,
            );

        // check that the minor views of `CombCodomain` are sorted
        for keymin in 0..num_indices_major {
            let view_minor_descend = comb_codomain_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend.collect_vec(), 
                    & umatch.order_operator_minor_reverse() 
                ) );
        }

        // check that the major views of `CombCodomain` are sorted
        for keymin in 0..num_indices_major {
            let view_major_ascend = comb_codomain_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_operator_major() 
                ) );
        }        


        // ----------------------------------------------------------------------------------------------------------------
        // check that the major and minor views of `CombCodomainInv` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                comb_codomain_inv_ref,
                0..num_indices_major,
                0..num_indices_major,
            );

        // check that the minor views of `CombCodomainInv` are sorted
        for keymin in 0..num_indices_major {
            let view_minor_descend = comb_codomain_inv_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend.collect_vec(), 
                    & umatch.order_operator_minor_reverse() 
                ) );
        }

        // check that the major views of `CombCodomainInv` are sorted
        for keymin in 0..num_indices_major {
            let view_major_ascend = comb_codomain_inv_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_operator_major() 
                ) );
        }        

        // ----------------------------------------------------------------------------------------------------------------        
        // check that the major and minor views of `CombDomain` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                comb_domain_ref,
                0..num_indices_minor,
                0..num_indices_minor,
            );

        // check that the minor views of `CombDomain` are sorted
        for keymin in 0..num_indices_minor {
            let view_minor_descend = comb_domain_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend.collect_vec(), 
                    & umatch.order_operator_minor_reverse() 
                ) );
        }  

        // check that the major views of `CombDomain` are sorted
        for keymin in 0..num_indices_minor {
            let view_major_ascend = comb_domain_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_operator_major() 
                ) );
        }          
        
        // ----------------------------------------------------------------------------------------------------------------        
        // check that the major and minor views of `CombDomainInv` agree
        verify_viewmajorascend_compatible_with_viewminordescend(
                comb_domain_inv_ref,
                0..num_indices_minor,
                0..num_indices_minor,
            );

        // check that the minor views of `CombDomainInv` are sorted
        for keymin in 0..num_indices_minor {
            let view_minor_descend = comb_domain_inv_ref.view_minor_descend( keymin );
            assert!( is_sorted_strictly(  
                    & view_minor_descend, 
                    & umatch.order_operator_minor_reverse() 
                ) );
        }   
        
        // check that the major views of `CombDomainInv` are sorted
        for keymin in 0..num_indices_minor {
            let view_major_ascend = comb_domain_inv_ref.view_major_ascend( keymin );
            assert!( is_sorted_strictly(  
                    & view_major_ascend.collect_vec(), 
                    & umatch.order_operator_major() 
                ) );
        }           






        // check that `CombCodomainInvTimesMappingMatchedBlockRowsIndexedByKeyMin` is upper triangular
        for keymin in umatch.matching_ref().bimap_min_ref().ord_to_val_vec().iter().cloned() {
            assert!(    comb_codomain_inv_times_mapping_matched_block_rows_indexed_by_keymin_ref.view_minor_descend( keymin ).next().unwrap().0 == keymin     );
        }        


// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:BEGIN        
        // check that columns are sorted in strictly descending order
        //  NOTE: THIS IS UNNECESSARY FOR THE COMBS, SINCE WE TEST THAT THEIR MINOR VIEWS EQUAL THOSE OF VecOfVec objects, WHOSE MINOR DESCENDING VIEWS ARE *ALWAYS* STRICTLY DECREASING IN INDEX
        for keymaj in 0 .. num_indices_minor { 
            assert!(    mapping.view_minor_descend( keymaj             ).is_sorted_by( |x, y| x.0 > y.0 )     );
            assert!(    comb_codomain.view_minor_descend( keymaj       ).is_sorted_by( |x, y| x.0 > y.0 )     );
            // assert!(    comb_codomain_inv.view_minor_descend( keymaj.clone()   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            // assert!(    comb_domain.view_minor_descend( keymaj.clone()         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            // assert!(    comb_domain_inv.view_minor_descend( keymaj.clone()     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }          
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END  

        
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG: BEGIN   
        // check that the major and minor views of the inverse of the codomain COMB agree
        let comb_codomain_vec_of_vec_simple     
            =   VecOfVec::from_iterable( (0..num_indices_major).map( |k| comb_codomain.view_major_ascend(k) ) );
        for keymaj in 0..num_indices_major {
            // println!("VIEW MAJOR DESCEND IS STARTING FOR THIS ROUND: keymaj = {:?}", keymaj);
            // println!("VIEW MAJOR DESCEND LAZY CONSTRUCTION: {:?}", comb_codomain.view_minor_descend( keymaj ).collect_vec());
            // println!("VIEW MAJOR DESCEND FROM MAJOR VIEW: {:?}", (& comb_codomain_vec_of_vec_simple).view_minor_descend( keymaj ).collect_vec());            
            // println!("VIEW MAJOR DESCEND IS FINISHED FOR THIS ROUND");
            itertools::assert_equal( 
                    comb_codomain.view_minor_descend( keymaj ),
                    (& comb_codomain_vec_of_vec_simple).view_minor_descend( keymaj )
                )
        }      
// ---------------- UNCOMMENT THIS WHEN READY TO RESUME DEBUG:END      
        
        // check that rows are sorted in strictly ascending order
        for keymaj in 0 .. num_indices_major { 
            assert!(    mapping.view_major_ascend( keymaj             ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_codomain.view_major_ascend( keymaj       ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_codomain_inv.view_major_ascend( keymaj   ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_domain.view_major_ascend( keymaj         ).is_sorted_by( |x, y| x.0 < y.0 )     );
            assert!(    comb_domain_inv.view_major_ascend( keymaj     ).is_sorted_by( |x, y| x.0 < y.0 )     );                                
        }  
              

    }



    //  ===================================================================================
    //  RECOVERY OF COMBS -- TARGETTED AT SPECIFIC POINTS IN THE PROCESS
    //  ===================================================================================

    
    //  COMB DOMAIN INV (SMALL + LARGE)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #CombDomainInvImplementViewRowAscend    

    #[test]
    fn test_retreival() {
        use itertools::Itertools;

        use crate::algebra::matrices::operations::umatch::row_major::{Umatch, CombCodomainInvTimesMappingMatchedBlock};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewRowAscend;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        
        use crate::utilities::order::{OrderOperatorAuto};

        let mapping: VecOfVec< usize, usize >   =   
        VecOfVec::new(  
                vec![
                    vec![ (0, 1), (1, 1), (2, 2) ],
                    vec![ (0, 1),         (2, 1) ],
                    vec![         (1, 1), (2, 1) ],
                ]
            );
        let _mapping_ref   =   & mapping; // the oracle trait is only implemented for *references* to a `VecOfVec` object
        
        let order_operator_major                =   OrderOperatorAuto;     
        let order_operator_minor                =   OrderOperatorAuto;                
        let ring_operator       =   PrimeOrderFieldOperator::new( 13 );

        let mapping_ref   =   & mapping;
        let umatch  
            =   Umatch::factor( 
                        mapping_ref, 
                        (0..3).rev(), 
                        ring_operator,
                        order_operator_major,
                        order_operator_minor,                        
                    );
        let _umatch_ref = & umatch;
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( umatch_ref);
        // let umatch_with_refs_ref = & umatch_with_refs;
        
        //  the "seed" matrix A (equal to the pivot block of the inverse of the codomain COMB times the pivot block of the matching array)
        let A   =   CombCodomainInvTimesMappingMatchedBlock::new( & umatch );
        let A_ref = &A;
        // println!("{:?}", umatch_with_refs.matching_ref());

        // for keymaj in 1..3 { println!( "{:?}", A_ref.view_major_ascend(keymaj).collect_vec() ) }

        //  the domain COMB
        let comb_domain_inv = umatch.comb_domain_inv();
        let comb_domain_inv_ground_truth
                =   VecOfVec::new(
                            vec![
                                vec![ (0, 1),         (2, 1) ],
                                vec![         (1, 1), (2, 1) ],
                                vec![                 (2, 1 )                        ],
                            ]
                        );
        let comb_domain_inv_ground_truth_ref = & comb_domain_inv_ground_truth;
        for keymaj in 0 .. 3 {
            // println!("GROUND TRUTH  : {:?}", comb_domain_inv_ground_truth_ref.view_major_ascend( keymaj ).collect_vec() );
            // println!("UNPACKED      : {:?}", comb_domain_inv.view_major_ascend( keymaj ).into_iter().collect_vec() );   
            // println!("SCALE FACTORS : {:?}", umatch_with_refs.matching.vec_snzval_ref() );    
            // println!("keymaj        : {:?}", keymaj );                                
            itertools::assert_equal(
                    comb_domain_inv_ground_truth_ref.view_major_ascend( keymaj ),
                    comb_domain_inv.view_major_ascend( keymaj ),
                )
        }
    }





    //  COMB DOMAIN (SMALL)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #CombDomainImplementViewRowAscend

    #[test]
    fn test_umatchrowmajor_comb_domain_small_example() {

        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewRowAscend;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        
        use crate::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};   
        use crate::algebra::matrices::operations::multiply::ProductMatrix;     

        // let num_rows = 2; let num_cols = 2; let modulus = 7;
        // let mapping = VecOfVec::new( vec![ vec![ (0usize,5usize), (1,5)], vec![ (1,6)]] );
        let num_rows = 1; let num_cols = 4; let modulus = 7;
        let mapping = VecOfVec::new( vec![ vec![ (2usize, 6usize), (3,1)], ]  );        
        //  NOTE: mapping can be regarded as     [  0  0  6  1  ]
        let mapping_ref = & mapping;
        let ring_operator = PrimeOrderFieldOperator::new( modulus );
        let order_operator_major = OrderOperatorAuto;
        let order_operator_minor = OrderOperatorAuto;        

        let umatch_root = 
                Umatch::factor( 
                        mapping_ref, 
                        (0 .. num_rows).rev(), 
                        ring_operator,
                        order_operator_major,
                        order_operator_minor,
                    );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch_root );
        // let umatch_ref = & umatch_with_refs;

        let comb_domain = umatch_root.comb_domain();
        let comb_domain_inv = umatch_root.comb_domain_inv();     
        let mapping_matched_cols_only = umatch_root.mapping_matched_cols_only();          
                
        // check that C * C^{-1} = identity
        let c_times_c_inv = 
            ProductMatrix::new( 
                    &comb_domain, 
                    &comb_domain_inv, 
                    ring_operator, 
                    OrderOperatorByKey::new() 
                );

        // println!("mapping:");
        // print_indexed_major_views( & mapping_ref, 0 .. num_rows );
        // println!("mapping_matched_cols_only:");
        // print_indexed_major_views( & mapping_matched_cols_only, 0 .. num_rows );        
        // println!("matching:");
        // print_indexed_major_views( & umatch_root.matching_ref(), 0 .. num_rows );    
        // println!("comb_codomain_inv_times_mapping_matched_block (recall that num_rows = {:?}):", num_rows);        
        // print_indexed_major_views( && umatch_root.comb_codomain_inv_times_mapping_matched_block(), 0 .. num_rows );
        // println!("comb_domain (recall that num_cols = {:?}) (THIS FUNCTION CALL SEEMS TO BREAK DOWN INTERNALLY WHERE comb_codomain_inv_times_mapping_matched_block IS CALLED):", num_cols);
        // print_indexed_major_views( & comb_domain, 0 .. num_cols );        
        // println!("comb_domain_inv:");
        // print_indexed_major_views( & comb_domain_inv, 0 .. num_cols );                
        for column_index in 0 .. num_cols {
            // println!("{:?}", c_times_c_inv.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( c_times_c_inv.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }

        // check D * C is right-reduced

    }


    //  COMB DOMAIN (LARGER + RANDOM)
    //  ------------------------------------------------------------------------------------
    //  TAG USED TO IDENTIFY THESE TESTS: #CombDomainImplementViewRowAscend)    


    #[test]
    fn test_umatchrowmajor_comb_domain() {

        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};
        
        use crate::algebra::matrices::query::ViewRowAscend;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        
        use crate::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};   
        use crate::algebra::matrices::operations::multiply::ProductMatrix;     

        let num_rows = 10; let num_cols = 10; let modulus = 7;
        let mapping = VecOfVec::random_mod_p(num_rows, num_cols, modulus);
        let mapping_ref = & mapping;
        let ring_operator = PrimeOrderFieldOperator::new( modulus );
        let order_operator_major = OrderOperatorAuto;
        let order_operator_minor = OrderOperatorAuto;        

        let umatch_root = 
                Umatch::factor( 
                        mapping_ref, 
                        (0 .. num_rows).rev(), 
                        ring_operator,
                        order_operator_major,
                        order_operator_minor,                        
                    );
        // let umatch_with_refs = UmatchRowMajorWithRefs::new( &umatch_root );
        // let umatch_ref = & umatch_with_refs;

        let comb_domain = umatch_root.comb_domain();
        let comb_domain_inv = umatch_root.comb_domain_inv();     
                
        // check that C * C^{-1} = identity
        let c_times_c_inv = 
            ProductMatrix::new( 
                    & comb_domain, 
                    & comb_domain_inv, 
                    ring_operator, 
                    OrderOperatorByKey::new() 
                );

        // println!("mapping:");
        // print_indexed_major_views( & mapping_ref, 0 .. num_rows );
        // println!("matching:");
        // print_indexed_major_views( & umatch_root.matching_ref(), 0 .. num_rows );        
        // println!("comb_domain:");
        // print_indexed_major_views( & comb_domain, 0 .. num_cols );        
        // println!("comb_domain_inv:");
        // print_indexed_major_views( & comb_domain_inv, 0 .. num_cols );                
        for column_index in 0 .. num_cols {
            // println!("{:?}", c_times_c_inv.view_major_ascend( column_index ).collect_vec() );
            itertools::assert_equal( c_times_c_inv.view_major_ascend( column_index ), std::iter::once( (column_index, 1) )   );
        }

        // check D * C is right-reduced

    }


    #[test]
    fn doc_test() {
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec};
        use crate::algebra::matrices::operations::{umatch::row_major::Umatch, multiply::ProductMatrix};
        use crate::algebra::matrices::debug::verify_that_product_is_identity;
        use crate::algebra::matrices::query::{ViewRowAscend};
        use crate::algebra::matrices::display::print_indexed_minor_views;
        
        use crate::utilities::order::{OrderOperatorAuto};
        use itertools::Itertools;

        // DEFINE INPUTS
        // ===============================

        // define the ring operator and order operator
        let modulus               =   5;
        let ring_operator         =   PrimeOrderFieldOperator::new( modulus );        
        let order_operator              =   OrderOperatorAuto;

        // define the matrix we wish to factor
        let mapping_data          =   VecOfVec::new( 
                                                vec![   
                                                            vec![(0,1), (1,1), (2,1)],
                                                            vec![                   ], 
                                                            vec![              (2,1)], 
                                                ] 
                                            );
        let mapping               =   & mapping_data;
                                        
        // COMPUTE U-MATCH
        // ===============================
                                        
        let umatch
            =   Umatch::factor(
                    mapping,  // the matrix we wish to factor
                    (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
                    ring_operator, // the operator for the coefficient ring
                    order_operator.clone(), // order comparator for the entries in each row
                    order_operator, // order comparator for the entries in each column
                );

        println!("matrix_major");
        umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal.matrix_major_data.print_dense(0);
        println!("{:?}", umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal.matrix_major_data.vec_of_vec() );
        
        println!("matrix_minor");
        umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal.matrix_minor_data.print_dense(0);        
        println!("{:?}", umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal.matrix_minor_data.vec_of_vec() );        
        
        println!("number of pairs: {:?}", umatch.matching.num_pairs() );        
        
        println!(
            "result of transpose: {:?}", 
            umatch.comb_codomain_inv_matched_block_indexed_by_matched_keymaj_ordinal_off_diagonal.matrix_major_data.transpose_deep( umatch.matching_ref().num_pairs() ).unwrap().vec_of_vec()
        );
            
            
        // INSPECT FACTORIZATION
        // ===============================
            
        // extract R, R^{-1}, C, C^{-1}, and M
        let r           =   umatch.comb_codomain();        // the codomain COMB
        let rinv        =   umatch.comb_codomain_inv();    // inverse of the the codomain COMB
        let c           =   umatch.comb_domain();          // the domain COMB
        let cinv        =   umatch.comb_domain_inv();      // inverse of the domain COMB
        let m           =   umatch.matching_ref();         // the generalized matching matrix
            
            
        println!("\nMinor views of the codomain COMB");   print_indexed_minor_views( &r, 0..3 ); 
        println!("\nMinor views of the   domain COMB");   print_indexed_minor_views( &c, 0..3 ); 
        println!("\nMinor views of the matching matrix"); print_indexed_minor_views( &m, 0..3 ); 
            
        // this will print the following:
        //
        // Minor views of the codomain COMB
        // minor_view 0: [(0, 1)]
        // minor_view 1: [(1, 1)]
        // 
        // Minor views of the   domain COMB
        // minor_view 0: [(0, 1)]
        // minor_view 1: [(1, 1), (0, 3)]
        // minor_view 2: [(2, 1), (0, 2)]
        // 
        // Minor views of the matching matrix
        // minor_view 0: [(0, 1)]
        // minor_view 1: []
        // minor_view 2: [(1, 1)]

        // SOLVE Ax = b FOR x
        // ===============================
        
        let b   =   [ (2,1), (0,1) ]; // note we list entries in reverse order
        let x   =   umatch.solve_dx_equals_b( b ).unwrap().collect_vec();
        let dx  =   umatch.multiply_dv(x);
        assert!( dx.eq( b ) );
            
            
        // VERIFY THE CALCULATION
        // ===============================
            
        // check that the product of the domain COMB with its inverse is identity: C * C^{-1} = I
        verify_that_product_is_identity( &c, &cinv, 0..3, ring_operator, OrderOperatorAuto );
            
        // check that the product of the codomain COMB with its inverse is identity: R * R^{-1} = I
        verify_that_product_is_identity( &r, &rinv, 0..3, ring_operator, OrderOperatorAuto );
            
        // check the factorization: R^{-1} * D * C = M
        let rinv_d   = ProductMatrix::new( &rinv,   &mapping, ring_operator, OrderOperatorAuto );      
        let rinv_d_c = ProductMatrix::new( &rinv_d, &c,       ring_operator, OrderOperatorAuto );                
        for keymaj in 0 .. 3 { 
            assert_eq!(
                rinv_d_c.view_major_ascend( keymaj ).collect_vec(),
                m.view_major_ascend( keymaj ).collect_vec()
            ) 
        }            
    }


    
 
}    


#[cfg(test)]
mod doc_test_solvers {
    

    

    


    #[test]
    fn doc_test_solve_dx_equals_b() {
        
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;        
        use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;        
        use crate::utilities::order::OrderOperatorAuto;
        use itertools::Itertools;

        // DEFINE THE MATRIX
        // ===============================
        let matrix          =   VecOfVec::new( 
                                                vec![   
                                                            vec![(0,true), (1,true), (2,true)],
                                                            vec![                            ], 
                                                            vec![                    (2,true)], 
                                                ] 
                                            );
                                        
        // COMPUTE U-MATCH
        // ===============================
                                        
        let umatch
            =   Umatch::factor(
                    & matrix,  // the matrix we wish to factor
                    (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
                    BooleanFieldOperator::new(), // the operator for the coefficient ring
                    OrderOperatorAuto, // order comparator for the entries in each row
                    OrderOperatorAuto, // order comparator for the entries in each column
                );        

        // SOLVE Ax = b FOR x
        // ===============================
        
        let b   =   [ (2,true), (0,true) ]; // note we list entries in reverse order
        let x   =   umatch.solve_dx_equals_b( b ).unwrap().collect_vec();
        let dx  =   umatch.multiply_dv(x);
        assert!( dx.eq( b ) );        
    }






    #[test]
    fn doc_test_solve_xd_equals_b() {

        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::rings::operator_structs::field_prime_order::BooleanFieldOperator;
        use crate::utilities::order::OrderOperatorAuto;        
        
        // define the matrix
        // -----------------
        let d = VecOfVec::new(
                vec![
                                        vec![  (0,true), (1,true),           ],
                                        vec![            (1,true), (2,true), ],
                                    ]
            );

        // obtain a u-match factorization
        // ------------------------------
        let umatch  =   Umatch::factor( 
            &d, 
            0..2, 
            BooleanFieldOperator::new(), 
            OrderOperatorAuto, 
            OrderOperatorAuto,             
        );
        
        // try solving xd = b
        // ------------------
        
        // Case 1: a solution exists; in this case we are gauaranteed to find one
        let x = umatch.solve_xd_equals_b( vec![ (0,true), (2,true), ] );        
        assert!( x.is_some() );
        assert!( x.unwrap().eq( vec![ (0,true), (1,true), ] ) );

        // Case 2: no solution exists; in this case we get a certificate that no solution exists
        let x = umatch.solve_xd_equals_b( vec![ (0,true), (1,true), (2,true) ] );        
        assert!( x.is_none() );
    }


    #[test]
    fn doc_test_solve_xd_equals_b__withfloatcoefficients() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::rings::operator_structs::ring_native::FieldFloat64;
        use crate::utilities::order::OrderOperatorAuto;            
        
        let d = VecOfVec::new(
                vec![
                                        vec![  (0,3.), (1,3.),         ],
                                        vec![          (1,3.), (2,3.), ],
                                    ]
            );
        let umatch  =   Umatch::factor( 
            &d, 
            0..2, 
            FieldFloat64::new(), 
            OrderOperatorAuto, 
            OrderOperatorAuto,             
        );
        
        // instance where there exists a solution
        let x = umatch.solve_xd_equals_b( vec![ (0,3.), (1,6.), (2,3.), ] );        
        assert!( x.is_some() );
        assert!( x.unwrap().eq( vec![ (0,1.), (1,1.) ] ) );

        // instance where there does not exist a solution
        let x = umatch.solve_xd_equals_b( vec![ (0,1.), (1,-1.) ] );        
        assert!( x.is_none() );
    }    


}    














