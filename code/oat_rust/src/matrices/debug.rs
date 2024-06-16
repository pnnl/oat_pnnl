//! Tools to analyze / debug matrices.

//  VERIFY THAT MAJOR VIEWS ARE COMPATIBLE WITH MINOR VIEWS
//  -----------------------------------------------------------------------------------------------

use crate::{entries::{KeyValGet, KeyValNew, KeyValSet}, utilities::partial_order::StrictlyLess, rings::operator_traits::Semiring};
use std::fmt::Debug;
use itertools::Itertools;
use crate::matrices::operations::multiply::ProductMatrixLazyMajorAscendSimplified;

use super::{matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend, IndicesAndCoefficients, OracleMajorDescend, OracleMinorAscend}, matrix_types::{oracle_ref::OracleRef, transpose::AntitransposeLazy}};

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
        Matrix:                             OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::ViewMajorAscendEntry:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::KeyMin, Matrix::SnzVal > + KeyValNew< Matrix::KeyMin, Matrix::SnzVal >,
        Matrix::ViewMinorDescend:           IntoIterator,
        Matrix::ViewMinorDescendEntry:      Debug + std::cmp::PartialEq + KeyValGet< Matrix::KeyMaj, Matrix::SnzVal > + KeyValNew< Matrix::KeyMaj, Matrix::SnzVal >,
        Matrix::KeyMin:                     Clone + Debug,
        Matrix::KeyMaj:                     Clone + Debug,
        Matrix::SnzVal:                     Debug,
        IterKeyMaj:                         Clone + Iterator< Item = Matrix::KeyMaj >,
        // IterKeyMin:                         Clone + Iterator< Item = Matrix::KeyMin >,        
{
    // Make sure all the entries in the major ascending view also appear in the minor descending view
    for keymaj in iter_keymaj {
        let viewmajascend = matrix.view_major_ascend( keymaj.clone() );
        // for each entry in the major ascending view ..
        for viewmajascend_entry in viewmajascend {
            // construct the corresponding entry of the corresponding minor descending view
            let viewmindescend_entry 
                    = Matrix::ViewMinorDescendEntry::new( keymaj.clone(), viewmajascend_entry.val() );
            // check that this entry does in fact lie in the minor descending view
            let viewminordescend = matrix.view_minor_descend( viewmajascend_entry.key() );
            let exists_in_minor_view = viewminordescend.into_iter().any( |x| x == viewmindescend_entry );
            if ! exists_in_minor_view {
                println!("");
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
                println!("");
                println!("first 10 entries of the {:?} {:?} view:", keymaj.clone(), interp_of_viewmaj, );
                println!("{:?}", first_10_of_viewmaj);
                println!("");
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
        Matrix:                             OracleMajorAscend + OracleMinorDescend + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::ViewMajorAscendEntry:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::KeyMin, Matrix::SnzVal > + KeyValNew< Matrix::KeyMin, Matrix::SnzVal >,
        Matrix::ViewMinorDescend:           IntoIterator,
        Matrix::ViewMinorDescendEntry:      Debug + std::cmp::PartialEq + KeyValGet< Matrix::KeyMaj, Matrix::SnzVal > + KeyValNew< Matrix::KeyMaj, Matrix::SnzVal >,
        Matrix::KeyMin:                     Clone + Debug,
        Matrix::KeyMaj:                     Clone + Debug,
        Matrix::SnzVal:                     Debug,        
        IterKeyMaj:                         Clone + Iterator< Item = Matrix::KeyMaj >,
        IterKeyMin:                         Clone + Iterator< Item = Matrix::KeyMin >,        
{
    verify_viewmajorascend_compatible_with_viewminordescend_helper_function( 
            OracleRef::new( & matrix), 
            // iter_keymin.clone(), 
            iter_keymaj.clone(),      
            "major",
            "minor",            
        );

    verify_viewmajorascend_compatible_with_viewminordescend_helper_function( 
            AntitransposeLazy::new( matrix ), 
            // iter_keymaj.clone(),
            iter_keymin.clone(),              
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
        Matrix:                             OracleMajorAscend + OracleMajorDescend< ViewMajorDescendEntry = Matrix::ViewMajorAscendEntry > + IndicesAndCoefficients,
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::ViewMajorAscendEntry:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::KeyMin, Matrix::SnzVal > + KeyValNew< Matrix::KeyMin, Matrix::SnzVal >,
        Matrix::KeyMaj:                     Clone + Debug,
        IterKeyMaj:                         Clone + Iterator< Item = Matrix::KeyMaj >,
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
        Matrix:                             OracleMinorAscend + OracleMinorDescend< ViewMinorDescendEntry = Matrix::ViewMinorAscendEntry > + IndicesAndCoefficients,
        Matrix::ViewMinorAscend:            IntoIterator,
        Matrix::ViewMinorAscendEntry:       Debug + std::cmp::PartialEq + KeyValGet< Matrix::KeyMin, Matrix::SnzVal > + KeyValNew< Matrix::KeyMin, Matrix::SnzVal >,
        Matrix::KeyMin:                     Clone + Debug,
        IterKeyMin:                         Clone + Iterator< Item = Matrix::KeyMin >,
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
            OrderComparator,
            IterKeyMaj,
        > 
        (
            matrix_1:           Matrix1,
            matrix_2:           Matrix2,
            iter_keymaj:        IterKeyMaj,
            ring_operator:      RingOperator,
            order_comparator:   OrderComparator
        )
        ->
        bool
    where
        Matrix1:                            OracleMajorAscend + IndicesAndCoefficients,
        Matrix2:                            OracleMajorAscend + IndicesAndCoefficients< SnzVal = Matrix1::SnzVal, KeyMaj = Matrix1::KeyMin, KeyMin = Matrix1::KeyMaj >, 
        Matrix1::ViewMajorAscend:           IntoIterator,
        Matrix2::ViewMajorAscend:           IntoIterator,
        Matrix1::ViewMajorAscendEntry:      KeyValGet < Matrix1::KeyMin, Matrix1::SnzVal >,
        Matrix2::ViewMajorAscendEntry:      KeyValGet < Matrix2::KeyMin, Matrix2::SnzVal > + KeyValSet < Matrix2::KeyMin, Matrix2::SnzVal >,  
        Matrix1::KeyMaj:                    Clone + std::fmt::Debug, // Debug is required by  itertools::assert_equal
        Matrix2::KeyMin:                    Clone + PartialEq, // PartialEq is required by the struct that simplifies linear combinations of iterators
        Matrix2::SnzVal:                    Clone + PartialEq + std::fmt::Debug, // PartialEq and Debug are required by  itertools::assert_equal
        IterKeyMaj:                         Iterator< Item = Matrix1::KeyMaj >,                    
        RingOperator:                       Clone + Semiring< Matrix1::SnzVal >,
        OrderComparator:                    Clone + StrictlyLess<  Matrix2::ViewMajorAscendEntry >,     
{

    let product = ProductMatrixLazyMajorAscendSimplified::new( matrix_1, matrix_2, ring_operator.clone(), order_comparator.clone() );
    let one = RingOperator::one();

    for keymaj in iter_keymaj {
        let view = product.view_major_ascend( keymaj.clone() );
        itertools::assert_equal(
            view.map( |x| (x.key(), x.val()) ),
            std::iter::once( ( keymaj, one.clone() ) )
        );
        // match equals_standard_unit_vector { true => {continue}, false=> {return false } }
    }
    return true
}