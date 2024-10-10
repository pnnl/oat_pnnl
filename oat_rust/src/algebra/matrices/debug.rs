//! Analyze and debug matrices.


use crate::{utilities::order::{JudgePartialOrder, is_sorted_strictly}, algebra::rings::operator_traits::Semiring};
use crate::algebra::vectors::entries::{KeyValGet, KeyValNew, KeyValSet};
use std::fmt::Debug;
use itertools::Itertools;
use crate::algebra::matrices::operations::multiply::ProductMatrix;

use super::{query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients, ViewRowDescend, ViewColAscend}, types::{transpose::AntiTranspose}};


//  VERIFY THAT VIEWS ARE SORTED
//  -----------------------------------------------------------------------------------------------

/// Check that entries of major ascending views appearin in strictly ascending order; only checks views indexed by `iter_keymaj`
pub fn verify_view_major_ascend_is_sorted_strictly< Matrix, IterKeyMaj, OrderOperatorKeyMin, PrintFunction >(
                matrix:                     Matrix,
                iter_keymaj:                IterKeyMaj,   
                order_operator_keymin:    OrderOperatorKeyMin,
                mut print_function:         PrintFunction,
            )
    where
        Matrix:                             ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::EntryMajor:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::ColIndex, Matrix::Coefficient > + KeyValNew< Matrix::ColIndex, Matrix::Coefficient >,
        Matrix::ColIndex:                     Clone + Debug,
        Matrix::RowIndex:                     Clone + Debug,
        Matrix::Coefficient:                     Debug,
        IterKeyMaj:                         Clone + Iterator< Item = Matrix::RowIndex >,
        OrderOperatorKeyMin:              JudgePartialOrder< Matrix::ColIndex  >,
        // IterKeyMin:                         Clone + Iterator< Item = Matrix::ColIndex >,  
        PrintFunction:                      FnMut( Vec< Matrix::ColIndex > ),      
{
    for keymaj in iter_keymaj {
        let viewmajascend = matrix.view_major_ascend( keymaj.clone() )
            .into_iter()
            .map(|x| x.key() )
            .collect_vec();
        
        if ! is_sorted_strictly( &viewmajascend, &order_operator_keymin ) { 
            print_function( viewmajascend );
            panic!("Major views are not strictly sorted") 
        }
        
    }  
}

/// Check that entries of minor descending views appearin in strictly ascending order; only checks views indexed by `iter_keymin`
pub fn verify_view_minor_descend_is_sorted_strictly< Matrix, IterKeyMin, OrderOperatorKeyMaj, PrintFunction >(
    matrix:                     Matrix,
    iter_keymin:                IterKeyMin,   
    order_operator_keymaj:    OrderOperatorKeyMaj,
    mut print_function:         PrintFunction,
)
    where
        Matrix:                             ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Matrix::ViewMinorDescend:           IntoIterator,
        Matrix::EntryMinor:      Debug + std::cmp::PartialEq + KeyValGet< Matrix::RowIndex, Matrix::Coefficient > + KeyValNew< Matrix::RowIndex, Matrix::Coefficient >,
        Matrix::ColIndex:                     Clone + Debug,
        Matrix::RowIndex:                     Clone + Debug,
        Matrix::Coefficient:                     Debug,
        IterKeyMin:                         Clone + Iterator< Item = Matrix::ColIndex >,
        OrderOperatorKeyMaj:              JudgePartialOrder< Matrix::RowIndex  >,
        // IterKeyMin:                         Clone + Iterator< Item = Matrix::ColIndex >,  
        PrintFunction:                      FnMut( Vec< Matrix::RowIndex > ),                    
{
    for keymin in iter_keymin {
        let viewmindescend = matrix.view_minor_descend( keymin.clone() )
            .into_iter()
            .map(|x| x.key() )
            .collect_vec();

        if ! is_sorted_strictly( &viewmindescend, &order_operator_keymaj ) { 
            print_function( viewmindescend );
            panic!("Minor views are not strictly sorted") 
        }

    }  
}



//  VERIFY THAT MAJOR VIEWS ARE COMPATIBLE WITH MINOR VIEWS
//  -----------------------------------------------------------------------------------------------



/// Ensures that each entry in each major ascending view also appears in the corresponding minor descending view; panics otherwise.
/// 
/// Only checks the major ascending views returned by `iter_keymaj`.
/// 
/// See comments in source code for details.
pub fn verify_viewmajorascend_compatible_with_viewminordescend_helper_function< Matrix, IterKeyMaj >(
                matrix:             Matrix,
                iter_keymaj:        IterKeyMaj,
                interp_of_viewmaj:  &str,
                interp_of_viewmin:  &str,                
            )
    where
        Matrix:                             ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::EntryMajor:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::ColIndex, Matrix::Coefficient > + KeyValNew< Matrix::ColIndex, Matrix::Coefficient >,
        Matrix::ViewMinorDescend:           IntoIterator,
        Matrix::EntryMinor:      Debug + std::cmp::PartialEq + KeyValGet< Matrix::RowIndex, Matrix::Coefficient > + KeyValNew< Matrix::RowIndex, Matrix::Coefficient >,
        Matrix::ColIndex:                     Clone + Debug,
        Matrix::RowIndex:                     Clone + Debug,
        Matrix::Coefficient:                     Debug,
        IterKeyMaj:                         Clone + Iterator< Item = Matrix::RowIndex >,
        // IterKeyMin:                         Clone + Iterator< Item = Matrix::ColIndex >,        
{
    // Make sure all the entries in the major ascending view also appear in the minor descending view
    for keymaj in iter_keymaj {
        let viewmajascend = matrix.view_major_ascend( keymaj.clone() );
        // for each entry in the major ascending view ..
        for viewmajascend_entry in viewmajascend {
            // construct the corresponding entry of the corresponding minor descending view
            let viewmindescend_entry 
                    = Matrix::EntryMinor::new( keymaj.clone(), viewmajascend_entry.val() );
            // check that this entry does in fact lie in the minor descending view
            let viewminordescend = matrix.view_minor_descend( viewmajascend_entry.key() );
            let exists_in_minor_view = viewminordescend.into_iter().any( |x| x == viewmindescend_entry );
            if ! exists_in_minor_view {
                println!();
                println!("error: entry {:?} appears in the {:?} {:?} view of the matrix, but entry {:?} does not appear in the {:?} {:?} view", 
                    (viewmajascend_entry.key(), viewmajascend_entry.val()), 
                    keymaj.clone(),
                    interp_of_viewmaj,        
                    (viewmindescend_entry.key(), viewmindescend_entry.val()),            
                    viewmajascend_entry.key(),
                    interp_of_viewmin,
                );      
                let first_10_of_viewmaj = matrix.view_major_ascend( keymaj.clone() ).into_iter().take(10).collect_vec();
                let first_10_of_viewmin = matrix.view_minor_descend( viewmajascend_entry.key() ).into_iter().take(10).collect_vec();                
                println!();
                println!("first 10 entries of the {:?} {:?} view:", keymaj.clone(), interp_of_viewmaj, );
                println!("{:?}", first_10_of_viewmaj);
                println!();
                println!("first 10 entries of the {:?} {:?} view:", viewmajascend_entry.key(), interp_of_viewmin, );                
                println!("{:?}", first_10_of_viewmin);                
            }
            assert!( exists_in_minor_view );
        }
    }  
}


/// Ensures that each entry in each major ascending view also appears in the corresponding minor descending view, and vice versa; panics otherwise.
pub fn verify_viewmajorascend_compatible_with_viewminordescend< Matrix, IterKeyMaj, IterKeyMin >(
                matrix:         Matrix,
                iter_keymin:    IterKeyMin,                
                iter_keymaj:    IterKeyMaj,
            )
    where
        Matrix:                             ViewRowAscend + ViewColDescend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::EntryMajor:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::ColIndex, Matrix::Coefficient > + KeyValNew< Matrix::ColIndex, Matrix::Coefficient >,
        Matrix::ViewMinorDescend:           IntoIterator,
        Matrix::EntryMinor:      Debug + std::cmp::PartialEq + KeyValGet< Matrix::RowIndex, Matrix::Coefficient > + KeyValNew< Matrix::RowIndex, Matrix::Coefficient >,
        Matrix::ColIndex:                     Clone + Debug,
        Matrix::RowIndex:                     Clone + Debug,
        Matrix::Coefficient:                     Debug,        
        IterKeyMaj:                         Clone + Iterator< Item = Matrix::RowIndex >,
        IterKeyMin:                         Clone + Iterator< Item = Matrix::ColIndex >,        
{
    verify_viewmajorascend_compatible_with_viewminordescend_helper_function( 
            & matrix, 
            // iter_keymin.clone(), 
            iter_keymaj,      
            "major",
            "minor",            
        );

    verify_viewmajorascend_compatible_with_viewminordescend_helper_function( 
            AntiTranspose::new( matrix ), 
            // iter_keymaj.clone(),
            iter_keymin,              
            "minor",            
            "major",                        
        );        
}

/// Panics if a major ascending view has entries different from the major descending view.
pub fn verify_viewmajorascend_compatible_with_viewmajordescend< Matrix, IterKeyMaj >(
    matrix:         Matrix,
    iter_keymaj:    IterKeyMaj,
)
    where
        Matrix:                             ViewRowAscend + ViewRowDescend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::EntryMajor:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::ColIndex, Matrix::Coefficient > + KeyValNew< Matrix::ColIndex, Matrix::Coefficient >,
        Matrix::RowIndex:                     Clone + Debug,
        IterKeyMaj:                         Clone + Iterator< Item = Matrix::RowIndex >,
{
    for keymaj in iter_keymaj {
        let ascend  =   matrix.view_major_ascend( keymaj.clone() );
        let descend  =   matrix.view_major_descend( keymaj.clone() );        
        let mut descend_vec: Vec<_> = descend.into_iter().collect();
        descend_vec.reverse();
        itertools::assert_equal(ascend, descend_vec);
    }
}  


/// Panics if a minor ascending view has entries different from the minor descending view.
pub fn verify_viewminorascend_compatible_with_viewminordescend< Matrix, IterKeyMin >(
    matrix:         Matrix,
    iter_keymin:    IterKeyMin,
)
    where
        Matrix:                             ViewColAscend + ViewColDescend + IndicesAndCoefficients,
        Matrix::ViewMinorAscend:            IntoIterator,
        Matrix::EntryMinor:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::ColIndex, Matrix::Coefficient > + KeyValNew< Matrix::ColIndex, Matrix::Coefficient >,
        Matrix::ColIndex:                     Clone + Debug,
        IterKeyMin:                         Clone + Iterator< Item = Matrix::ColIndex >,
{
    for keymin in iter_keymin {
        let ascend   =   matrix.view_minor_ascend( keymin.clone() );
        let descend  =   matrix.view_minor_descend( keymin.clone() );        
        let mut descend_vec: Vec<_> = descend.into_iter().collect();
        descend_vec.reverse();
        itertools::assert_equal(ascend, descend_vec);
    }
}  





//  UTILITY FUNCTION FOR TESTING: CHECK THAT A PRODUCT OF TWO MATRICES IS IDENTITY
//  ==============================================================================

/// Panics if view `k` of the product of two matrices does not equal the `k`th standard unit vector, for any `k` in `iter_keymaj`.
pub fn verify_that_product_is_identity<
            Matrix1, 
            Matrix2,               
            RingOperator,
            OrderOperator,
            IterKeyMaj,
        > 
        (
            matrix_1:           Matrix1,
            matrix_2:           Matrix2,
            iter_keymaj:        IterKeyMaj,
            ring_operator:      RingOperator,
            order_operator:   OrderOperator
        )
        ->
        bool
    where
        Matrix1:                            ViewRowAscend + IndicesAndCoefficients,
        Matrix2:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = Matrix1::Coefficient, RowIndex = Matrix1::ColIndex, ColIndex = Matrix1::RowIndex >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::EntryMajor:      KeyValGet < Matrix1::ColIndex, Matrix1::Coefficient >,
        Matrix2::EntryMajor:      KeyValGet < Matrix2::ColIndex, Matrix2::Coefficient > + KeyValSet < Matrix2::ColIndex, Matrix2::Coefficient >,  
        Matrix1::RowIndex:                    Clone + std::fmt::Debug, // Debug is required by  itertools::assert_equal
        Matrix2::ColIndex:                    Clone + PartialEq, // PartialEq is required by the struct that simplifies linear combinations of iterators
        Matrix2::Coefficient:                    Clone + PartialEq + std::fmt::Debug, // PartialEq and Debug are required by  itertools::assert_equal
        IterKeyMaj:                         Iterator< Item = Matrix1::RowIndex >,                    
        RingOperator:                       Clone + Semiring< Matrix1::Coefficient >,
        OrderOperator:                    Clone + JudgePartialOrder<  Matrix2::EntryMajor >,     
{

    let product = ProductMatrix::new( matrix_1, matrix_2, ring_operator, order_operator );
    let one = RingOperator::one();

    for keymaj in iter_keymaj {
        let view = product.view_major_ascend( keymaj.clone() );
        itertools::assert_equal(
            view.map( |x| (x.key(), x.val()) ),
            std::iter::once( ( keymaj, one.clone() ) )
        );
        // match equals_standard_unit_vector { true => {continue}, false=> {return false } }
    }
    true
}