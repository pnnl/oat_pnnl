//! Matrix wrapper that prepends an entry to each ascending major view and descending minor view of the wrapped matrix.


use std::{iter::{Once, once}};

use crate::{algebra::matrices::query::{ViewRowAscend, IndicesAndCoefficients, ViewColDescend}, };
use crate::algebra::vectors::entries::KeyValNew;



/// Prepends a entry equivalent to `( keymaj, diagonal_scalar )` to each ascending major
/// view of form `self.view_major_ascend( keymaj )`. This operation is unsafe, in the 
/// sense that it does not check that the resulting matrix is actually upper triangular.
/// 
/// # Design notes
/// 
/// The developers originally attempted to write this struct such that the `matrix_unprepended`
/// field *owned* `MatrixUnprepended`, not just a *reference* to `MatrixUnprepended`.  However,
/// this lead to errors which we didn't fully understand.  We hypothesize that this may have to
/// do with lifetime parameters.  In particular, in the original version where `matrix_unprepended`
/// owned its matrix and not a reference, the compiler asked us to give an explicit type
/// annotation for the `ViewMajorAscend` type of `MatrixUnprepended`; that ran us into difficulty,
/// becasue the `ViewMajorAscend` type involved a lifetime parameter.
pub struct PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend < MatrixUnprepended: IndicesAndCoefficients >
{
    matrix_unprepended:         MatrixUnprepended,
    diagonal_scalar:            MatrixUnprepended::Coefficient,
}

// Implement the struct
impl < MatrixUnprepended: IndicesAndCoefficients >

    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
    < MatrixUnprepended >
{
    pub fn new( matrix_unprepended: MatrixUnprepended, diagonal_scalar: MatrixUnprepended::Coefficient ) -> Self {  
        PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend{
                matrix_unprepended,
                diagonal_scalar, 
            }
    }

}

//  IndicesAndCoefficients

impl < MatrixUnprepended: IndicesAndCoefficients > 

    IndicesAndCoefficients for 
    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend < MatrixUnprepended >

{   type ColIndex = MatrixUnprepended::ColIndex; 
    type RowIndex = MatrixUnprepended::RowIndex; 
    type Coefficient = MatrixUnprepended::Coefficient;
    type EntryMajor = MatrixUnprepended::EntryMajor;
    type EntryMinor = MatrixUnprepended::EntryMinor;    
 }    


//  ViewRowAscend
impl < 'a, MatrixUnprepended > 

    ViewRowAscend  for 

    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
        < MatrixUnprepended > 
    
    where
        MatrixUnprepended:                          ViewRowAscend + IndicesAndCoefficients< ColIndex = Self::RowIndex >, // this construction only works when minor keys have the same type as major keys
        MatrixUnprepended::ViewMajorAscend:         IntoIterator,
        MatrixUnprepended::EntryMajor:    KeyValNew< MatrixUnprepended::ColIndex, MatrixUnprepended::Coefficient >,
        Self::RowIndex:                               Clone,
        Self::Coefficient:                               Clone,
{
    type ViewMajorAscend          =   core::iter::Chain < Once < MatrixUnprepended::EntryMajor >, MatrixUnprepended::ViewMajorAscendIntoIter >;
    type ViewMajorAscendIntoIter  =   Self::ViewMajorAscend;
    
     fn view_major_ascend( &self, keymaj: Self::RowIndex ) -> Self::ViewMajorAscend {
        std::iter::Iterator::chain(
                once( Self::EntryMajor::new( keymaj.clone(), self.diagonal_scalar.clone() )  ),
                self.matrix_unprepended.view_major_ascend( keymaj ),
            )
     }
}


//  ViewRowAscend
impl < 'a, MatrixUnprepended > 

    ViewColDescend  for 

    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
        < MatrixUnprepended > 
    
    where
        MatrixUnprepended:                          ViewColDescend + IndicesAndCoefficients< ColIndex = Self::RowIndex >, // this construction only works when minor keys have the same type as major keys
        MatrixUnprepended::ViewMinorDescend:        IntoIterator,
        MatrixUnprepended::EntryMinor:   KeyValNew< MatrixUnprepended::ColIndex, MatrixUnprepended::Coefficient >,
        Self::RowIndex:                               Clone,
        Self::Coefficient:                               Clone,
{
    type ViewMinorDescend          =   core::iter::Chain < Once < MatrixUnprepended::EntryMinor >, MatrixUnprepended::ViewMinorDescendIntoIter >;
    type ViewMinorDescendIntoIter  =   Self::ViewMinorDescend;
    
     fn view_minor_descend( &self, keymaj: Self::RowIndex ) -> Self::ViewMinorDescend {
        std::iter::Iterator::chain(
                once( Self::EntryMinor::new( keymaj.clone(), self.diagonal_scalar.clone() )  ),
                self.matrix_unprepended.view_minor_descend( keymaj ),
            )
     }
}








// struct PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend < MatrixUnprepended, Scalar >
// {
//     matrix_unprepended:     MatrixUnprepended,
//     diagonal_scalar:        Scalar
// }

// impl < MatrixUnprepended, ViewMajorAscend, RowIndex, Coefficient > 

//     ViewRowAscend
//     < RowIndex, core::iter::Chain < Once < ViewMajorAscend::Item >, ViewMajorAscend::IntoIter > > for 

//     PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
//     < MatrixUnprepended, Coefficient > where

//     MatrixUnprepended:      ViewRowAscend < RowIndex, ViewMajorAscend >,
//     ViewMajorAscend:        IntoIterator,
//     ViewMajorAscend::Item:  Clone + KeyValNew < RowIndex, Coefficient >,
// {
//      fn view_major_ascend( &self, keymaj: RowIndex ) -> Chain < Once < ViewMajorAscend::Item >, ViewMajorAscend > {
//         std::iter::Iterator::chain(
//                 once( ViewMajorAscend::Item::KeyValMake( keymaj.clone(), self.diagonal_scalar.clone() )  ),
//                 self.matrix_unprepended.view_major_ascend( keymaj ).into_iter(),
//             )
//      }
// }


//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    use std::iter::Cloned;

    

    


    
    #[test] 
    fn test_PrependEntryToViewMajorAscend() {

        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::types::prepend_viewmaj::PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend;
        use crate::algebra::matrices::query::ViewRowAscend;
        use std::slice::Iter;

        type ViewMajorAscend<'a> = Cloned< Iter< 'a, (i32, i32) > >;
        let matrix_unprepended  =   VecOfVec::new(
                                            vec![
                                                vec![ (1, 1), (2, 1) ],
                                                vec![         (2, 1) ],
                                                vec![                ],
                                            ]
                                        );

        let matrix_prepended
        // :   PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< 
        //                                 & VecOfVec<i32, i32>, 
        //                                 Cloned< Iter< (i32, i32) > >, 
        //                                 i32, 
        //                                 i32
        //                             >
            =   PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( & matrix_unprepended, 1 );

        let intended_result  =   VecOfVec::new(
                                            vec![
                                                vec![ (0, 1), (1, 1), (2, 1) ],
                                                vec![         (1, 1), (2, 1) ],
                                                vec![                 (2, 1) ],
                                            ]
                                        );     
        // let intended_result_ref: &'a VecOfVec< i32, i32 > = & &intended_result;
                                        
        for keymaj in 0 .. 3 {
            itertools::assert_equal( 
                matrix_prepended.view_major_ascend( keymaj ),
                (& intended_result).view_major_ascend( keymaj ),
            )
        }                                              

    }

}