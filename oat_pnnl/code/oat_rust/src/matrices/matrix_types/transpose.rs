//! Lazy transpose and anti-transpose; wraps around another matrix, and swaps order and/or major vs. minor views.

use crate::{matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMinorDescend, OracleMajorDescend, OracleMinorAscend, IndicesAndCoefficients}, entries::{KeyValGet, KeyValNew}};

use super::oracle_ref::OracleRef;



//  ANTITRANSPOSE
//  -----------------------------------------------------------------------------------------------

/// Antitranpose of a matrix, evaluated in a lazy fashion. (warning: doesn't work quite the same as for ordinary dense matrices)
/// 
/// Concretely, the antitranpose is obtained by (i) transposing, and then (ii) reversing the order of rows and reversing the order of columns.
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  For example, to execute the method `view_major_ascend`, an `AntitransposeLazy` struct simply calls `view_minor_descend` on
/// the underlying matrix.   **Note that doesn't work quite the same as for ordinary dense matrices when indices happen to be integers.**
pub struct AntitransposeLazy< Matrix > { unantitransposed: Matrix }

impl < Matrix >

    AntitransposeLazy
        < Matrix >
{
    pub fn new( unantitransposed: Matrix ) -> AntitransposeLazy< Matrix > { AntitransposeLazy { unantitransposed } }
}


    //  CLONE
impl < Matrix: Clone > 

    Clone for

    AntitransposeLazy< Matrix >

{ fn clone(&self) -> Self { AntitransposeLazy { unantitransposed: self.unantitransposed.clone() } } }    


    //  INDICES AND COEFFICIENTS
impl < Matrix: IndicesAndCoefficients > 

    IndicesAndCoefficients for

    AntitransposeLazy< Matrix >

{ type KeyMin = Matrix::KeyMaj; type KeyMaj = Matrix::KeyMin; type SnzVal = Matrix::SnzVal; }    


    //  MAJOR ASCEND
impl < Matrix > 

    OracleMajorAscend for 

    AntitransposeLazy< Matrix >

    where
        Matrix:                     OracleMinorDescend + IndicesAndCoefficients,
{
    type ViewMajorAscend            = Matrix::ViewMinorDescend;
    type ViewMajorAscendEntry       = Matrix::ViewMinorDescendEntry;
    type ViewMajorAscendIntoIter    = Matrix::ViewMinorDescendIntoIter;

    fn   view_major_ascend( &self, index: Self::KeyMaj ) -> Self::ViewMajorAscend {  self.unantitransposed.view_minor_descend( index ) } 
}


//  MAJOR DESCEND
impl < Matrix > 

    OracleMajorDescend for 

    AntitransposeLazy< Matrix >

    where
        Matrix:                     OracleMinorAscend + IndicesAndCoefficients,
{
    type ViewMajorDescend           =   Matrix::ViewMinorAscend;
    type ViewMajorDescendEntry      =   Matrix::ViewMinorAscendEntry;
    type ViewMajorDescendIntoIter   =   Matrix::ViewMinorAscendIntoIter;

    fn   view_major_descend( &self, index: Self::KeyMaj ) -> Self::ViewMajorDescend {  self.unantitransposed.view_minor_ascend( index ) } 
}


//  MINOR ASCEND
impl < Matrix > 

    OracleMinorAscend for 

    AntitransposeLazy< Matrix >

    where
        Matrix:                     OracleMajorDescend + IndicesAndCoefficients,
{
    type ViewMinorAscend            =   Matrix::ViewMajorDescend;
    type ViewMinorAscendEntry       =   Matrix::ViewMajorDescendEntry;
    type ViewMinorAscendIntoIter    =   Matrix::ViewMajorDescendIntoIter;

    fn   view_minor_ascend( &self, index: Self::KeyMin ) -> Self::ViewMinorAscend {  self.unantitransposed.view_major_descend( index ) } 
}


//  MINOR DESCEND
impl < Matrix > 

    OracleMinorDescend for 

    AntitransposeLazy< Matrix >

    where
        Matrix:                     OracleMajorAscend + IndicesAndCoefficients,
{
    type ViewMinorDescend           =   Matrix::ViewMajorAscend;
    type ViewMinorDescendIntoIter   =   Matrix::ViewMajorAscendIntoIter;
    type ViewMinorDescendEntry      =   Matrix::ViewMajorAscendEntry;
    
    fn   view_minor_descend( &self, index: Self::KeyMin ) -> Self::ViewMinorDescend {  self.unantitransposed.view_major_ascend( index ) } 
}




// //  ANTITRANSPOSE
// //  -----------------------------------------------------------------------------------------------

// /// Antitranpose of a matrix.
// /// 
// /// Concretely, the antitranpose is obtained by (i) transposing, and then (ii) reversing the order of rows and reversing the order of columns.
// pub struct AntitransposeLazy< 'antitranspose, Matrix > { unantitransposed: &'antitranspose Matrix }

// impl < 'antitranspose, Matrix >

//     AntitransposeLazy
//         < 'antitranspose, Matrix >
//     {
//         pub fn new( unantitransposed: &'antitranspose Matrix ) -> AntitransposeLazy< 'antitranspose, Matrix > { AntitransposeLazy { unantitransposed } }
//     }

// //  MAJOR ASCEND
// impl < 'antitranspose, Matrix, KeyMin, ViewMinorDescend > 

//     OracleMajorAscend
//         < KeyMin, ViewMinorDescend > for 

//     AntitransposeLazy< 'antitranspose, Matrix >

//     where
//         Matrix:                     OracleMinorDescend< KeyMin, ViewMinorDescend >,
//         ViewMinorDescend:           IntoIterator,
// {
//     fn   view_major_ascend( &self, index: KeyMin ) -> ViewMinorDescend {  self.unantitransposed.view_minor_descend( index ) } 
// }

// //  MAJOR DESCEND
// impl < 'antitranspose, Matrix, KeyMin, ViewMinorAscend > 

//     OracleMajorDescend
//         < KeyMin, ViewMinorAscend > for 

//     AntitransposeLazy< 'antitranspose, Matrix >

//     where
//         Matrix:                     OracleMinorAscend< KeyMin, ViewMinorAscend >,
//         ViewMinorAscend:            IntoIterator,
// {
//     fn   view_major_descend( &self, index: KeyMin ) -> ViewMinorAscend {  self.unantitransposed.view_minor_ascend( index ) } 
// }

// //  MINOR ASCEND
// impl < 'antitranspose, Matrix, KeyMaj, ViewMajorDescend > 

//     OracleMinorAscend
//         < KeyMaj, ViewMajorDescend > for 

//     AntitransposeLazy< 'antitranspose, Matrix >

//     where
//         Matrix:                     OracleMajorAscend< KeyMaj, ViewMajorDescend >,
//         ViewMajorDescend:           IntoIterator,
// {
//     fn   view_minor_ascend( &self, index: KeyMaj ) -> ViewMajorDescend {  self.unantitransposed.view_major_ascend( index ) } 
// }

// //  MINOR DESCEND
// impl < 'antitranspose, Matrix, KeyMaj, ViewMajorAscend > 

//     OracleMinorDescend
//         < KeyMaj, ViewMajorAscend > for 

//     AntitransposeLazy< 'antitranspose, Matrix >

//     where
//         Matrix:                     OracleMajorAscend< KeyMaj, ViewMajorAscend >,
//         ViewMajorAscend:            IntoIterator,
// {
//     fn   view_minor_descend( &self, index: KeyMaj ) -> ViewMajorAscend {  self.unantitransposed.view_major_ascend( index ) } 
// }














//  TRANSPOSE
//  -----------------------------------------------------------------------------------------------

/// Transpose of a matrix, evaluated in a lazy fashion.
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  For example, to execute the method `view_major_ascend`, a `TransposeLazy` struct simply calls `view_minor_ascend` on
/// the underlying matrix.
pub struct TransposeLazy< Matrix > { untransposed: Matrix }

impl < Matrix >

    TransposeLazy
        < Matrix >
{
    pub fn new( untransposed: Matrix ) -> TransposeLazy< Matrix > { TransposeLazy { untransposed } }
}

    //  CLONE
impl < Matrix: Clone > 

    Clone for

    TransposeLazy< Matrix >

{ fn clone(&self) -> Self { TransposeLazy { untransposed: self.untransposed.clone() } } }    


    //  INDICES AND COEFFICIENTS
impl < Matrix: IndicesAndCoefficients > 

    IndicesAndCoefficients for

    TransposeLazy< Matrix >

{ type KeyMin = Matrix::KeyMaj; type KeyMaj = Matrix::KeyMin; type SnzVal = Matrix::SnzVal; }    


    //  MAJOR ASCEND
impl < Matrix > 

    OracleMajorAscend for 

    TransposeLazy< Matrix >

    where
        Matrix:                     OracleMinorAscend + IndicesAndCoefficients,
{
    type ViewMajorAscend            = Matrix::ViewMinorAscend;
    type ViewMajorAscendEntry       = Matrix::ViewMinorAscendEntry;
    type ViewMajorAscendIntoIter    = Matrix::ViewMinorAscendIntoIter;

    fn   view_major_ascend( &self, index: Self::KeyMaj ) -> Self::ViewMajorAscend {  self.untransposed.view_minor_ascend( index ) } 
}


//  MAJOR DESCEND
impl < Matrix > 

    OracleMajorDescend for 

    TransposeLazy< Matrix >

    where
        Matrix:                     OracleMinorDescend + IndicesAndCoefficients,
{
    type ViewMajorDescend           =   Matrix::ViewMinorDescend;
    type ViewMajorDescendEntry      =   Matrix::ViewMinorDescendEntry;
    type ViewMajorDescendIntoIter   =   Matrix::ViewMinorDescendIntoIter;

    fn   view_major_descend( &self, index: Self::KeyMaj ) -> Self::ViewMajorDescend {  self.untransposed.view_minor_descend( index ) } 
}


//  MINOR ASCEND
impl < Matrix > 

    OracleMinorAscend for 

    TransposeLazy< Matrix >

    where
        Matrix:                     OracleMajorAscend + IndicesAndCoefficients,
{
    type ViewMinorAscend            =   Matrix::ViewMajorAscend;
    type ViewMinorAscendEntry       =   Matrix::ViewMajorAscendEntry;
    type ViewMinorAscendIntoIter    =   Matrix::ViewMajorAscendIntoIter;

    fn   view_minor_ascend( &self, index: Self::KeyMin ) -> Self::ViewMinorAscend {  self.untransposed.view_major_ascend( index ) } 
}


//  MINOR DESCEND
impl < Matrix > 

    OracleMinorDescend for 

    TransposeLazy< Matrix >

    where
        Matrix:                     OracleMajorDescend + IndicesAndCoefficients,
{
    type ViewMinorDescend           =   Matrix::ViewMajorDescend;
    type ViewMinorDescendIntoIter   =   Matrix::ViewMajorDescendIntoIter;
    type ViewMinorDescendEntry      =   Matrix::ViewMajorDescendEntry;
    
    fn   view_minor_descend( &self, index: Self::KeyMin ) -> Self::ViewMinorDescend {  self.untransposed.view_major_descend( index ) } 
}








