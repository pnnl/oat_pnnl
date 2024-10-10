//! Matrix multiplication, inversion, factorization, etc.
//! 
//!
//! - [accessing entries in rows and columns](crate::algebra::matrices::query)
//! - [matrix and vector multiplication](crate::algebra::matrices::operations::multiply)
//! - [solving systems of equations](crate::algebra::matrices::operations::solve)
//! - [matrix factorization](crate::algebra::matrices::operations::umatch)
//! - [matrix inversion](crate::algebra::matrices::operations::invert)
//! - [kernels](crate::algebra::matrices::operations::umatch::row_major::Umatch::kernel)
//! - [images](crate::algebra::matrices::operations::umatch::row_major::Umatch::image)
//! - [transforming and combining vectors (scale, add, drop zeros, reindex, etc.)](crate::algebra::vectors::operations)

pub mod multiply; 
pub mod invert;
pub mod vec_of_vec_reduction;
pub mod solve;
pub mod transform_entry_wise;
pub mod transform_vector_wise;
pub mod umatch;
// pub mod display;
pub mod transpose;


// =========================================================================


use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};
use crate::algebra::vectors::entries::{KeyValSet, KeyValGet};

use crate::utilities::order::{OrderOperatorByKeyCutsom, JudgePartialOrder};

use self::{multiply::ProductMatrix, umatch::row_major::Umatch};

use super::query::{ViewRowAscend, IndicesAndCoefficients};
use super::types::packet::MatrixAlgebraPacket;
use super::types::transpose::{Transpose, AntiTranspose};
use super::types::reverse::Reverse;

use std::hash::Hash;
use std::fmt::Debug;


/// Convenient methods for matrix operations
pub trait MatrixOperations {



    // fn map_left< Vector, RingOperator, OrderOperator, >(
    //     self,
    //     vector: Vector,
    //     ring_operator: RingOperator, 
    //     order_operator: OrderOperator 
    // )
    // ->  LinearCombinationSimplified< MatrixView::IntoIter, MatrixIndex, RingElement, RingOperator >
    
    // Simplify<
    //             HitMerge<
    //                     Scale< MatrixView::IntoIter, MatrixIndex, RingOperator, RingElement >,
    //                     OrderOperator,
    //                 >,
    //             MatrixIndex,
    //             RingOperator,
    //             RingElement,
    //         >
    //     where 
    //         Self:                                   Sized,         
    //         Matrix:                                 FnMut( Index ) -> MatrixView,
    //         MatrixView:                             IntoIterator,
    //         < MatrixView as IntoIterator >::Item:   KeyValGet< MatrixIndex, RingElement > +  KeyValSet< MatrixIndex, RingElement >,
    //         OrderOperator:                          JudgePartialOrder< < MatrixView as IntoIterator >::Item >,
    //         MatrixIndex:                            PartialEq,
    //         RingElement:                            Clone,
    //         RingOperator:                           Clone + Semiring< RingElement >,
    // {
    //     return self.matrix_multiply_unsimplified(matrix, ring_operator.clone(), order_operator).simplify( ring_operator )
    // }      



    /// Lefthand matrix multiplication
    /// 
    /// Returns `self * other`
    fn multiply_left< Othr, RingOperator, OrderOperator >( 
                self, 
                other: Othr, 
                ring_operator: RingOperator, 
                order_operator: OrderOperator 
            )
        ->
        ProductMatrix< Self, Othr, RingOperator, OrderOperator > 

        where 
            Self:                            ViewRowAscend + IndicesAndCoefficients + Sized,
            Othr:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Self::Coefficient, RowIndex = Self::ColIndex >, 
            Self::ViewMajorAscend:           IntoIterator,
            Othr::ViewMajorAscend:           IntoIterator,
            Self::EntryMajor:                KeyValGet < Self::ColIndex, Self::Coefficient >,
            Othr::EntryMajor:                KeyValGet < Othr::ColIndex, Othr::Coefficient >,  
            Othr::ColIndex:                  Clone,
            Othr::Coefficient:              Clone,                          
            RingOperator:                    Clone + Semiring< Self::Coefficient >,
            OrderOperator:                   Clone + JudgePartialOrder<  Othr::EntryMajor >,        
    {
        ProductMatrix::new( self, other, ring_operator, order_operator)
    }

    /// Righthand matrix multiplication
    /// 
    /// Returns `other * self`
    fn multiply_right< Othr, RingOperator, OrderOperator >( 
                self, 
                other: Othr, 
                ring_operator: RingOperator, 
                order_operator: OrderOperator 
            )
        ->
        ProductMatrix< Othr, Self, RingOperator, OrderOperator > 
        
        where 
            Self:                            ViewRowAscend + IndicesAndCoefficients + Sized,
            Othr:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Self::Coefficient, ColIndex = Self::RowIndex >, 
            Self::ViewMajorAscend:           IntoIterator,
            Othr::ViewMajorAscend:           IntoIterator,
            Self::EntryMajor:                KeyValGet < Self::ColIndex, Self::Coefficient >,
            Othr::EntryMajor:                KeyValGet < Othr::ColIndex, Othr::Coefficient >,  
            Self::ColIndex:                  Clone,
            Self::Coefficient:              Clone,                          
            RingOperator:                    Clone + Semiring< Self::Coefficient >,
            OrderOperator:                   Clone + JudgePartialOrder< Self::EntryMajor >, 
    {
        ProductMatrix::new( other, self, ring_operator, order_operator)
    }  


    /// Lefthand matrix multiplication
    /// 
    /// Returns `self * other`, where `other` is a matrix "packet," containing a matrix, a ring operator, and order operators for row and column entries.
    /// 
    /// **Note** There is no corresponding function `multiply_right_packet`, because we need an order operator from the packet to combine entries
    fn multiply_left_packet< Othr, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >( 
                self, 
                other: MatrixAlgebraPacket< Othr, RingOperator, OrderOperatorRowEntries, OrderOperatorColEntries >, 
            )
        ->
        ProductMatrix< Self, Othr, RingOperator, OrderOperatorRowEntries > 

        where 
            Self:                            ViewRowAscend + IndicesAndCoefficients + Sized,
            Othr:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Self::Coefficient, RowIndex = Self::ColIndex >, 
            Self::ViewMajorAscend:           IntoIterator,
            Othr::ViewMajorAscend:           IntoIterator,
            Self::EntryMajor:                KeyValGet < Self::ColIndex, Self::Coefficient >,
            Othr::EntryMajor:                KeyValGet < Othr::ColIndex, Othr::Coefficient >,  
            Othr::ColIndex:                  Clone,
            Othr::Coefficient:              Clone,                          
            RingOperator:                    Clone + Semiring< Self::Coefficient >,
            OrderOperatorRowEntries:         Clone + JudgePartialOrder<  Othr::EntryMajor >,        
    {
        ProductMatrix::new( self, other.matrix, other.ring, other.row_entry_order )
    }


    /// Returns a U-match factorization
    /// 
    /// The argument `iter_keymaj` must run over row indices in *strictly descending order*, as determined by `order_operator_RowIndex`.
    fn umatch< Other, RingOperator, OrderOperatorRowIndex, OrderOperatorColIndex, IterRowIndex >(
                self, 
                iter_keymaj:                IterRowIndex, 
                ring_operator:              RingOperator, 
                order_operator_RowIndex:    OrderOperatorRowIndex,                                
                order_operator_ColIndex:    OrderOperatorColIndex,
            )
        ->
        Umatch 
            < 
                Self, 
                RingOperator, 
                OrderOperatorByKeyCutsom< Self::ColIndex, Self::Coefficient, Self::EntryMajor, OrderOperatorColIndex, >,  // recall that entries in the major views have minor keys
                OrderOperatorByKeyCutsom< Self::RowIndex, Self::Coefficient, Self::EntryMinor, OrderOperatorRowIndex, >, // recall that entries in the minor views have major keys
            >  
        where   Self:                    Sized + ViewRowAscend + IndicesAndCoefficients, 
                IterRowIndex:                 Iterator < Item = Self::RowIndex >,
                Self::ColIndex:          Clone + Hash + std::cmp::Eq + Debug, 
                Self::RowIndex:          Clone + Hash + std::cmp::Eq + Debug,     // !!! remove Debug eventually       
                Self::Coefficient:      Clone + Debug,            // !!! remove Debug eventually       
                RingOperator:               Clone + Semiring< Self::Coefficient > + Ring< Self::Coefficient > + DivisionRing< Self::Coefficient >,
                Self::ViewMajorAscend:   Clone, // !!! remove clone
                Self::ViewMajorAscendIntoIter: Clone,                
                Self::EntryMajor:        KeyValSet< Self::ColIndex, Self::Coefficient > + Debug + Clone,    // !!!!!!!!!!!!!! REMOVE THE DEBUG REQUIREMENT AFTER DEBUGGING!!!!!            
                OrderOperatorColIndex:        Clone + JudgePartialOrder <  Self::ColIndex >, // !!! remove clone
                OrderOperatorRowIndex:        Clone + JudgePartialOrder <  Self::RowIndex >,    

        {
            Umatch::factor( self, iter_keymaj, ring_operator, order_operator_RowIndex, order_operator_ColIndex )
        }   

    /// Wraps `self` in a [crate::algebra::matrices::types::transpose::Transpose] struct.
    fn transpose( self ) -> Transpose<Self> 
        where
            Self:   Sized,
    { Transpose::new( self )  } 

    /// Wraps `& self` in a [crate::algebra::matrices::types::transpose::Transpose] struct.
    fn transpose_ref( &self ) -> Transpose<&Self> 
        where
    { Transpose::new( & self )  }         

    /// Wraps `self` in a [crate::algebra::matrices::types::transpose::Transpose] struct.
    /// 
    /// # Caution
    /// 
    /// There are three important differences between [AntiTranspose] and the matrix returned by [antitranspose_deep](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep):
    /// 
    /// - [AntiTranspose] is a lazy object that does not generate new data
    /// - The set of (key,val) pairs that appear in a major (respectively, minor) view of `AntiTranspose::new(matrix)`
    ///   are the *same* as the entries in a minor (respectively, major) view of `matrix`; only the sequence in which those entries appear is different.
    ///   By contrast, the keys in the (key,val) pairs of [matrix.antitranspose_deep()](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) are different;
    ///   they are obtained by subtracting the original keys from (# rows in the antitransposed matrix - 1).
    /// - For this reason, [matrix.antitranspose_deep()](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) is only available for
    ///   very specific types of matrices; [AntiTranspose] is available for a much broader class.
    fn antitranspose( self ) -> AntiTranspose<Self> 
        where
            Self:   Sized,
    { AntiTranspose::new( self )  } 

    /// Wraps `& self` in a [crate::algebra::matrices::types::transpose::AntiTranspose] struct.
    fn antitranspose_ref( &self ) -> AntiTranspose<&Self> 
        where
            Self:   Sized,    
    { AntiTranspose::new( &self ) }    

    /// Wraps `self` in a [crate::algebra::matrices::types::transpose::Transpose] struct.
    fn reverse( self ) -> Reverse<Self> 
        where
            Self:   Sized,
    { Reverse::new( self )  } 

    /// Wraps `& self` in a [crate::algebra::matrices::types::transpose::Transpose] struct.
    fn reverse_ref( &self ) -> Reverse<&Self> 
        where
    { Reverse::new( & self )  } 
}



//  Auto-implement this trait on all types
//  --------------------------------------------------------------
impl < Matrix > 

    MatrixOperations for 
    
    Matrix

{}    