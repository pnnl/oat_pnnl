//! Matrix wrapper that prepends an entry to each ascending major view and descending minor view of the wrapped matrix.


use std::{iter::{Chain, Once, once}, marker::PhantomData};

use crate::{matrices::matrix_oracle_traits::{OracleMajorAscend, IndicesAndCoefficients, OracleMinorDescend}, entries::KeyValNew};



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
    diagonal_scalar:            MatrixUnprepended::SnzVal,
}

// Implement the struct
impl < MatrixUnprepended: IndicesAndCoefficients >

    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
    < MatrixUnprepended >
{
    pub fn new( matrix_unprepended: MatrixUnprepended, diagonal_scalar: MatrixUnprepended::SnzVal ) -> Self {  
        PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend{
                matrix_unprepended:         matrix_unprepended,
                diagonal_scalar:            diagonal_scalar, 
            }
    }

}

//  IndicesAndCoefficients

impl < MatrixUnprepended: IndicesAndCoefficients > 

    IndicesAndCoefficients for 
    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend < MatrixUnprepended >

{ type KeyMin = MatrixUnprepended::KeyMin; type KeyMaj = MatrixUnprepended::KeyMaj; type SnzVal = MatrixUnprepended::SnzVal; }    


//  OracleMajorAscend
impl < 'a, MatrixUnprepended > 

    OracleMajorAscend  for 

    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
        < MatrixUnprepended > 
    
    where
        MatrixUnprepended:                          OracleMajorAscend + IndicesAndCoefficients< KeyMin = Self::KeyMaj >, // this construction only works when minor keys have the same type as major keys
        MatrixUnprepended::ViewMajorAscend:         IntoIterator,
        MatrixUnprepended::ViewMajorAscendEntry:    KeyValNew< MatrixUnprepended::KeyMin, MatrixUnprepended::SnzVal >,
        Self::KeyMaj:                               Clone,
        Self::SnzVal:                               Clone,
{
    type ViewMajorAscend          =   core::iter::Chain < Once < MatrixUnprepended::ViewMajorAscendEntry >, MatrixUnprepended::ViewMajorAscendIntoIter >;
    type ViewMajorAscendEntry     =   < MatrixUnprepended::ViewMajorAscend as IntoIterator >::Item;
    type ViewMajorAscendIntoIter  =   Self::ViewMajorAscend;
    
     fn view_major_ascend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMajorAscend {
        std::iter::Iterator::chain(
                once( Self::ViewMajorAscendEntry::new( keymaj.clone(), self.diagonal_scalar.clone() )  ),
                self.matrix_unprepended.view_major_ascend( keymaj ),
            )
     }
}


//  OracleMajorAscend
impl < 'a, MatrixUnprepended > 

    OracleMinorDescend  for 

    PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
        < MatrixUnprepended > 
    
    where
        MatrixUnprepended:                          OracleMinorDescend + IndicesAndCoefficients< KeyMin = Self::KeyMaj >, // this construction only works when minor keys have the same type as major keys
        MatrixUnprepended::ViewMinorDescend:        IntoIterator,
        MatrixUnprepended::ViewMinorDescendEntry:   KeyValNew< MatrixUnprepended::KeyMin, MatrixUnprepended::SnzVal >,
        Self::KeyMaj:                               Clone,
        Self::SnzVal:                               Clone,
{
    type ViewMinorDescend          =   core::iter::Chain < Once < MatrixUnprepended::ViewMinorDescendEntry >, MatrixUnprepended::ViewMinorDescendIntoIter >;
    type ViewMinorDescendEntry     =   < MatrixUnprepended::ViewMinorDescend as IntoIterator >::Item;
    type ViewMinorDescendIntoIter  =   Self::ViewMinorDescend;
    
     fn view_minor_descend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMinorDescend {
        std::iter::Iterator::chain(
                once( Self::ViewMinorDescendEntry::new( keymaj.clone(), self.diagonal_scalar.clone() )  ),
                self.matrix_unprepended.view_minor_descend( keymaj ),
            )
     }
}








// struct PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend < MatrixUnprepended, Scalar >
// {
//     matrix_unprepended:     MatrixUnprepended,
//     diagonal_scalar:        Scalar
// }

// impl < MatrixUnprepended, ViewMajorAscend, KeyMaj, SnzVal > 

//     OracleMajorAscend
//     < KeyMaj, core::iter::Chain < Once < ViewMajorAscend::Item >, ViewMajorAscend::IntoIter > > for 

//     PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend 
//     < MatrixUnprepended, SnzVal > where

//     MatrixUnprepended:      OracleMajorAscend < KeyMaj, ViewMajorAscend >,
//     ViewMajorAscend:        IntoIterator,
//     ViewMajorAscend::Item:  Clone + KeyValNew < KeyMaj, SnzVal >,
// {
//      fn view_major_ascend( &self, keymaj: KeyMaj ) -> Chain < Once < ViewMajorAscend::Item >, ViewMajorAscend > {
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

        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_types::prepend_viewmaj::PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use std::slice::Iter;

        type ViewMajorAscend<'a> = Cloned< Iter< 'a, (i32, i32) > >;
        let matrix_unprepended  =   VecOfVecSimple::new(
                                            vec![
                                                vec![ (1, 1), (2, 1) ],
                                                vec![         (2, 1) ],
                                                vec![                ],
                                            ]
                                        );

        let matrix_prepended
        // :   PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend< 
        //                                 & VecOfVecSimple<i32, i32>, 
        //                                 Cloned< Iter< (i32, i32) > >, 
        //                                 i32, 
        //                                 i32
        //                             >
            =   PrependDiagonalEntryToViewMajorAscendAndViewMinorDescend::new( & matrix_unprepended, 1 );

        let intended_result  =   VecOfVecSimple::new(
                                            vec![
                                                vec![ (0, 1), (1, 1), (2, 1) ],
                                                vec![         (1, 1), (2, 1) ],
                                                vec![                 (2, 1) ],
                                            ]
                                        );     
        // let intended_result_ref: &'a VecOfVecSimple< i32, i32 > = & &intended_result;
                                        
        for keymaj in 0 .. 3 {
            itertools::assert_equal( 
                matrix_prepended.view_major_ascend( keymaj ),
                (& intended_result).view_major_ascend( keymaj ),
            )
        }                                              

    }

}