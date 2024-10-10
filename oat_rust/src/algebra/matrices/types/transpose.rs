//! Lazy transpose and anti-transpose; wraps around another matrix, and swaps order and/or major vs. minor views.

use crate::{algebra::matrices::query::{ViewRowAscend, ViewColDescend, ViewRowDescend, ViewColAscend, IndicesAndCoefficients, ViewRow, ViewCol, MatrixOracle}, };






//  ANTITRANSPOSE
//  -----------------------------------------------------------------------------------------------

/// Antitranpose of a matrix, evaluated in a lazy fashion. (warning: doesn't work quite the same as for ordinary dense matrices)
/// 
/// Concretely, the antitranpose is obtained by (i) transposing, and then (ii) reversing the order of rows and reversing the order of columns.
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  For example, to execute the method `view_major_ascend`, an `AntiTranspose` struct simply calls `view_minor_descend` on
/// the underlying matrix.   **Note that doesn't work quite the same as for ordinary dense matrices when indices happen to be integers.**
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
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::transpose::{Transpose, AntiTranspose};
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use itertools;
/// 
/// // matrix
/// let a     =   VecOfVec::new( vec![
///                                 vec![ (0,0), (1,1), (2,2) ],
///                                 vec![ (0,3), (1,4), (2,5) ],
///                              ] );
/// 
/// // its transpose
/// let tran    =   Transpose::new( &a );
/// 
/// // its antitranspose
/// let anti    =   AntiTranspose::new( &a );                                                                                                                                 
/// 
/// 
/// for row in 0 .. 2 {
///     
///     assert!( itertools::equal( (&a).row(row),                           tran.column( row ) )                        );
///     assert!( itertools::equal( (&a).row_opt(row).unwrap(),              tran.column_opt( row ).unwrap() )           );
///     assert!( itertools::equal( (&a).row_reverse(row),                   tran.column_reverse( row))                  );
///     assert!( itertools::equal( (&a).row_reverse_opt(row).unwrap(),      tran.column_reverse_opt( row).unwrap() )    );
/// 
///     assert!( itertools::equal( (&a).row_reverse(row),                   anti.column( row) )                         );
///     assert!( itertools::equal( (&a).row_reverse_opt(row).unwrap(),      anti.column_opt( row).unwrap() )            ); 
///     assert!( itertools::equal( (&a).row(row),                           anti.column_reverse( row) )                 );
///     assert!( itertools::equal( (&a).row_opt(row).unwrap(),              anti.column_reverse_opt( row).unwrap() )    );
/// }
/// 
/// for col in 0 .. 3 {
/// 
///     assert!( itertools::equal( (&a).column(col),                        tran.row(col) )                             );
///     assert!( itertools::equal( (&a).column_opt(col).unwrap(),           tran.row_opt(col).unwrap() )                );
///     assert!( itertools::equal( (&a).column_reverse(col),                tran.row_reverse(col) )                     );
///     assert!( itertools::equal( (&a).column_reverse_opt(col).unwrap(),   tran.row_reverse_opt(col).unwrap() )        );
///     
///     assert!( itertools::equal( (&a).column(col),                        anti.row_reverse(col) )                     );
///     assert!( itertools::equal( (&a).column_opt(col).unwrap(),           anti.row_reverse_opt(col).unwrap() )        );
///     assert!( itertools::equal( (&a).column_reverse(col),                anti.row(col) )                             );
///     assert!( itertools::equal( (&a).column_reverse_opt(col).unwrap(),   anti.row_opt(col).unwrap() )                );
/// }  
/// ```
pub struct AntiTranspose< Matrix > { unantitransposed: Matrix }

impl < Matrix >

    AntiTranspose
        < Matrix >
{
    pub fn new( unantitransposed: Matrix ) -> AntiTranspose< Matrix > { AntiTranspose { unantitransposed } }
}


    //  CLONE
impl < Matrix: Clone > 

    Clone for

    AntiTranspose< Matrix >

{ fn clone(&self) -> Self { AntiTranspose { unantitransposed: self.unantitransposed.clone() } } }    

    //  MATRIX ORACLE
impl< Matrix > 

    MatrixOracle for 
    
    AntiTranspose< Matrix >

    where
        Matrix:     MatrixOracle
{
    type Coefficient        =   Matrix::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   Matrix::ColumnIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   Matrix::RowIndex;       // The type of column indices
    
    type RowEntry           =   Matrix::ColumnEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   Matrix::RowEntry;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Matrix::ColumnReverse;               // What you get when you ask for a row.
    type RowIter            =   Matrix::ColumnReverseIter;           // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse         =   Matrix::Column;        // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter     =   Matrix::ColumnIter;    // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    
    type Column             =   Matrix::RowReverse;            // What you get when you ask for a column
    type ColumnIter         =   Matrix::RowReverseIter;        // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse      =   Matrix::Row;     // What you get when you ask for a column with the order of entries reversed 
    type ColumnReverseIter  =   Matrix::RowIter; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    fn entry(                   & self,  row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient > {
        self.unantitransposed.entry( column, row )
    }

    fn row(                     & self,  index: Self::RowIndex    )       -> Self::Row 
        { self.unantitransposed.column_reverse(index) }
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>
        { self.unantitransposed.column_reverse_opt(index) }     
    fn row_reverse(             & self,  index: Self::RowIndex    )       -> Self::RowReverse
        { self.unantitransposed.column(index) }    
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>
        { self.unantitransposed.column_opt(index) }
    
    fn column(                  & self,  index: Self::ColumnIndex )       -> Self::Column
        { self.unantitransposed.row_reverse(index) }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column>
        { self.unantitransposed.row_reverse_opt(index) }    
    fn column_reverse(          & self,  index: Self::ColumnIndex )       -> Self::ColumnReverse
        { self.unantitransposed.row(index) }
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse>
        { self.unantitransposed.row_opt(index) }
} 


    //  INDICES AND COEFFICIENTS
impl < Matrix: IndicesAndCoefficients > 

    IndicesAndCoefficients for

    AntiTranspose< Matrix >

{ 
    type ColIndex = Matrix::RowIndex; 
    type RowIndex = Matrix::ColIndex; 
    type Coefficient = Matrix::Coefficient;
    type EntryMajor = Matrix::EntryMinor;
    type EntryMinor = Matrix::EntryMajor;
}    


    //  MAJOR ASCEND
impl < Matrix > 

    ViewRowAscend for 

    AntiTranspose< Matrix >

    where
        Matrix:                     ViewColDescend + IndicesAndCoefficients,
{
    type ViewMajorAscend            = Matrix::ViewMinorDescend;
    type ViewMajorAscendIntoIter    = Matrix::ViewMinorDescendIntoIter;

    fn   view_major_ascend( &self, index: Self::RowIndex ) -> Self::ViewMajorAscend {  self.unantitransposed.view_minor_descend( index ) } 
}


//  MAJOR DESCEND
impl < Matrix > 

    ViewRowDescend for 

    AntiTranspose< Matrix >

    where
        Matrix:                     ViewColAscend + IndicesAndCoefficients,
{
    type ViewMajorDescend           =   Matrix::ViewMinorAscend;
    type ViewMajorDescendIntoIter   =   Matrix::ViewMinorAscendIntoIter;

    fn   view_major_descend( &self, index: Self::RowIndex ) -> Self::ViewMajorDescend {  self.unantitransposed.view_minor_ascend( index ) } 
}


//  MINOR ASCEND
impl < Matrix > 

    ViewColAscend for 

    AntiTranspose< Matrix >

    where
        Matrix:                     ViewRowDescend + IndicesAndCoefficients,
{
    type ViewMinorAscend            =   Matrix::ViewMajorDescend;
    type ViewMinorAscendIntoIter    =   Matrix::ViewMajorDescendIntoIter;

    fn   view_minor_ascend( &self, index: Self::ColIndex ) -> Self::ViewMinorAscend {  self.unantitransposed.view_major_descend( index ) } 
}


//  MINOR DESCEND
impl < Matrix > 

    ViewColDescend for 

    AntiTranspose< Matrix >

    where
        Matrix:                     ViewRowAscend + IndicesAndCoefficients,
{
    type ViewMinorDescend           =   Matrix::ViewMajorAscend;
    type ViewMinorDescendIntoIter   =   Matrix::ViewMajorAscendIntoIter;
    
    fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend {  self.unantitransposed.view_major_ascend( index ) } 
}




// //  ANTITRANSPOSE
// //  -----------------------------------------------------------------------------------------------

// /// Antitranpose of a matrix.
// /// 
// /// Concretely, the antitranpose is obtained by (i) transposing, and then (ii) reversing the order of rows and reversing the order of columns.
// pub struct AntiTranspose< 'antitranspose, Matrix > { unantitransposed: &'antitranspose Matrix }

// impl < 'antitranspose, Matrix >

//     AntiTranspose
//         < 'antitranspose, Matrix >
//     {
//         pub fn new( unantitransposed: &'antitranspose Matrix ) -> AntiTranspose< 'antitranspose, Matrix > { AntiTranspose { unantitransposed } }
//     }

// //  MAJOR ASCEND
// impl < 'antitranspose, Matrix, ColIndex, ViewMinorDescend > 

//     ViewRowAscend
//         < ColIndex, ViewMinorDescend > for 

//     AntiTranspose< 'antitranspose, Matrix >

//     where
//         Matrix:                     ViewColDescend< ColIndex, ViewMinorDescend >,
//         ViewMinorDescend:           IntoIterator,
// {
//     fn   view_major_ascend( &self, index: ColIndex ) -> ViewMinorDescend {  self.unantitransposed.view_minor_descend( index ) } 
// }

// //  MAJOR DESCEND
// impl < 'antitranspose, Matrix, ColIndex, ViewMinorAscend > 

//     ViewRowDescend
//         < ColIndex, ViewMinorAscend > for 

//     AntiTranspose< 'antitranspose, Matrix >

//     where
//         Matrix:                     ViewColAscend< ColIndex, ViewMinorAscend >,
//         ViewMinorAscend:            IntoIterator,
// {
//     fn   view_major_descend( &self, index: ColIndex ) -> ViewMinorAscend {  self.unantitransposed.view_minor_ascend( index ) } 
// }

// //  MINOR ASCEND
// impl < 'antitranspose, Matrix, RowIndex, ViewMajorDescend > 

//     ViewColAscend
//         < RowIndex, ViewMajorDescend > for 

//     AntiTranspose< 'antitranspose, Matrix >

//     where
//         Matrix:                     ViewRowAscend< RowIndex, ViewMajorDescend >,
//         ViewMajorDescend:           IntoIterator,
// {
//     fn   view_minor_ascend( &self, index: RowIndex ) -> ViewMajorDescend {  self.unantitransposed.view_major_ascend( index ) } 
// }

// //  MINOR DESCEND
// impl < 'antitranspose, Matrix, RowIndex, ViewMajorAscend > 

//     ViewColDescend
//         < RowIndex, ViewMajorAscend > for 

//     AntiTranspose< 'antitranspose, Matrix >

//     where
//         Matrix:                     ViewRowAscend< RowIndex, ViewMajorAscend >,
//         ViewMajorAscend:            IntoIterator,
// {
//     fn   view_minor_descend( &self, index: RowIndex ) -> ViewMajorAscend {  self.unantitransposed.view_major_ascend( index ) } 
// }














//  TRANSPOSE
//  -----------------------------------------------------------------------------------------------

/// Transpose of a matrix, evaluated in a lazy fashion.
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  For example, to execute the method `view_major_ascend`, a `Transpose` struct simply calls `view_minor_ascend` on
/// the underlying matrix.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::transpose::{Transpose, AntiTranspose};
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use itertools;
/// 
/// // matrix
/// let a     =   VecOfVec::new( vec![
///                                 vec![ (0,0), (1,1), (2,2) ],
///                                 vec![ (0,3), (1,4), (2,5) ],
///                              ] );
/// 
/// // its transpose
/// let tran    =   Transpose::new( &a );
/// 
/// // its antitranspose
/// let anti    =   AntiTranspose::new( &a );                                                                                                                                 
/// 
/// 
/// for row in 0 .. 2 {
///     
///     assert!( itertools::equal( (&a).row(row),                           tran.column( row ) )                        );
///     assert!( itertools::equal( (&a).row_opt(row).unwrap(),              tran.column_opt( row ).unwrap() )           );
///     assert!( itertools::equal( (&a).row_reverse(row),                   tran.column_reverse( row))                  );
///     assert!( itertools::equal( (&a).row_reverse_opt(row).unwrap(),      tran.column_reverse_opt( row).unwrap() )    );
/// 
///     assert!( itertools::equal( (&a).row_reverse(row),                   anti.column( row) )                         );
///     assert!( itertools::equal( (&a).row_reverse_opt(row).unwrap(),      anti.column_opt( row).unwrap() )            ); 
///     assert!( itertools::equal( (&a).row(row),                           anti.column_reverse( row) )                 );
///     assert!( itertools::equal( (&a).row_opt(row).unwrap(),              anti.column_reverse_opt( row).unwrap() )    );
/// }
/// 
/// for col in 0 .. 3 {
/// 
///     assert!( itertools::equal( (&a).column(col),                        tran.row(col) )                             );
///     assert!( itertools::equal( (&a).column_opt(col).unwrap(),           tran.row_opt(col).unwrap() )                );
///     assert!( itertools::equal( (&a).column_reverse(col),                tran.row_reverse(col) )                     );
///     assert!( itertools::equal( (&a).column_reverse_opt(col).unwrap(),   tran.row_reverse_opt(col).unwrap() )        );
///     
///     assert!( itertools::equal( (&a).column(col),                        anti.row_reverse(col) )                     );
///     assert!( itertools::equal( (&a).column_opt(col).unwrap(),           anti.row_reverse_opt(col).unwrap() )        );
///     assert!( itertools::equal( (&a).column_reverse(col),                anti.row(col) )                             );
///     assert!( itertools::equal( (&a).column_reverse_opt(col).unwrap(),   anti.row_opt(col).unwrap() )                );
/// }  
/// ```
/// 
pub struct Transpose< Matrix > { untransposed: Matrix }

impl < Matrix >

    Transpose
        < Matrix >
{
    pub fn new( untransposed: Matrix ) -> Transpose< Matrix > { Transpose { untransposed } }
}

    //  CLONE
impl < Matrix: Clone > 

    Clone for

    Transpose< Matrix >

{ fn clone(&self) -> Self { Transpose { untransposed: self.untransposed.clone() } } }    


    //  MATRIX ORACLE
    impl< Matrix > 

    MatrixOracle for 
    
    Transpose< Matrix >

    where
        Matrix:     MatrixOracle

{
    type Coefficient        =   Matrix::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   Matrix::ColumnIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   Matrix::RowIndex;       // The type of column indices
    
    type RowEntry           =   Matrix::ColumnEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   Matrix::RowEntry;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Matrix::Column;               // What you get when you ask for a row.
    type RowIter            =   Matrix::ColumnIter;           // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse         =   Matrix::ColumnReverse;        // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter     =   Matrix::ColumnReverseIter;    // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    
    type Column             =   Matrix::Row;            // What you get when you ask for a column
    type ColumnIter         =   Matrix::RowIter;        // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse      =   Matrix::RowReverse;     // What you get when you ask for a column with the order of entries reversed 
    type ColumnReverseIter  =   Matrix::RowReverseIter; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    fn entry(                   & self,  row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient > {
        self.untransposed.entry( column, row )
    }

    fn row(                     & self,  index: Self::RowIndex    )       -> Self::Row 
        { self.untransposed.column(index) }
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>
        { self.untransposed.column_opt(index) }     
    fn row_reverse(             & self,  index: Self::RowIndex    )       -> Self::RowReverse
        { self.untransposed.column_reverse(index) }    
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>
        { self.untransposed.column_reverse_opt(index) }
    
    fn column(                  & self,  index: Self::ColumnIndex )       -> Self::Column
        { self.untransposed.row(index) }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column>
        { self.untransposed.row_opt(index) }    
    fn column_reverse(          & self,  index: Self::ColumnIndex )       -> Self::ColumnReverse
        { self.untransposed.row_reverse(index) }
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse>
        { self.untransposed.row_reverse_opt(index) }
} 



    //  INDICES AND COEFFICIENTS
impl < Matrix: IndicesAndCoefficients > 

    IndicesAndCoefficients for

    Transpose< Matrix >

{ 
    type EntryMajor = Matrix::EntryMinor;
    type EntryMinor = Matrix::EntryMajor;
    type ColIndex = Matrix::RowIndex; 
    type RowIndex = Matrix::ColIndex; 
    type Coefficient = Matrix::Coefficient; 
}    


//  MAJOR ASCEND
impl < Matrix > 

    ViewRow for 

    Transpose< Matrix >

    where
        Matrix:                     ViewCol + IndicesAndCoefficients,
{
    type ViewMajor            = Matrix::ViewMinor;
    type ViewMajorIntoIter    = Matrix::ViewMinorIntoIter;

    fn   view_major( &self, index: Self::RowIndex ) -> Self::ViewMajor {  self.untransposed.view_minor( index ) } 
}


//  MAJOR ASCEND
impl < Matrix > 

    ViewRowAscend for 

    Transpose< Matrix >

    where
        Matrix:                     ViewColAscend + IndicesAndCoefficients,
{
    type ViewMajorAscend            = Matrix::ViewMinorAscend;
    type ViewMajorAscendIntoIter    = Matrix::ViewMinorAscendIntoIter;

    fn   view_major_ascend( &self, index: Self::RowIndex ) -> Self::ViewMajorAscend {  self.untransposed.view_minor_ascend( index ) } 
}


//  MAJOR DESCEND
impl < Matrix > 

    ViewRowDescend for 

    Transpose< Matrix >

    where
        Matrix:                     ViewColDescend + IndicesAndCoefficients,
{
    type ViewMajorDescend           =   Matrix::ViewMinorDescend;
    type ViewMajorDescendIntoIter   =   Matrix::ViewMinorDescendIntoIter;

    fn   view_major_descend( &self, index: Self::RowIndex ) -> Self::ViewMajorDescend {  self.untransposed.view_minor_descend( index ) } 
}


//  MAJOR ASCEND
impl < Matrix > 

    ViewCol for 

    Transpose< Matrix >

    where
        Matrix:                     ViewRow + IndicesAndCoefficients,
{
    type ViewMinor            = Matrix::ViewMajor;
    type ViewMinorIntoIter    = Matrix::ViewMajorIntoIter;

    fn   view_minor( &self, index: Self::ColIndex ) -> Self::ViewMinor {  self.untransposed.view_major( index ) } 
}


//  MINOR ASCEND
impl < Matrix > 

    ViewColAscend for 

    Transpose< Matrix >

    where
        Matrix:                     ViewRowAscend + IndicesAndCoefficients,
{
    type ViewMinorAscend            =   Matrix::ViewMajorAscend;
    type ViewMinorAscendIntoIter    =   Matrix::ViewMajorAscendIntoIter;

    fn   view_minor_ascend( &self, index: Self::ColIndex ) -> Self::ViewMinorAscend {  self.untransposed.view_major_ascend( index ) } 
}


//  MINOR DESCEND
impl < Matrix > 

    ViewColDescend for 

    Transpose< Matrix >

    where
        Matrix:                     ViewRowDescend + IndicesAndCoefficients,
{
    type ViewMinorDescend           =   Matrix::ViewMajorDescend;
    type ViewMinorDescendIntoIter   =   Matrix::ViewMajorDescendIntoIter;
    
    fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend {  self.untransposed.view_major_descend( index ) } 
}







//  ---------------------------------------------------------------------------
//  TEST
//  ---------------------------------------------------------------------------


#[cfg(test)]
mod tests {
    use itertools::Itertools;

    


    #[test]
    fn test_transpose_and_antitranspose() {
        
        use crate::algebra::matrices::types::transpose::{Transpose, AntiTranspose};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        use itertools;
        
        // matrix
        let a     =   VecOfVec::new( vec![
                                                vec![ (0,0), (1,1), (2,2) ],
                                                vec![ (0,3), (1,4), (2,5) ],
                                            ] );
        
        // its transpose
        let tran    =   Transpose::new( &a );

        // its antitranspose
        let anti    =   AntiTranspose::new( &a );                                                                                                                                 
        

        for row in 0 .. 2 {
            
            assert!( itertools::equal( (&a).row(row),                           tran.column( row ) )                        );
            assert!( itertools::equal( (&a).row_opt(row).unwrap(),              tran.column_opt( row ).unwrap() )           );
            assert!( itertools::equal( (&a).row_reverse(row),                   tran.column_reverse( row))                  );
            assert!( itertools::equal( (&a).row_reverse_opt(row).unwrap(),      tran.column_reverse_opt( row).unwrap() )    );

            assert!( itertools::equal( (&a).row_reverse(row),                   anti.column( row) )                         );
            assert!( itertools::equal( (&a).row_reverse_opt(row).unwrap(),      anti.column_opt( row).unwrap() )            ); 
            assert!( itertools::equal( (&a).row(row),                           anti.column_reverse( row) )                 );
            assert!( itertools::equal( (&a).row_opt(row).unwrap(),              anti.column_reverse_opt( row).unwrap() )    );
        }

        for col in 0 .. 3 {

            assert!( itertools::equal( (&a).column(col),                        tran.row(col) )                             );
            assert!( itertools::equal( (&a).column_opt(col).unwrap(),           tran.row_opt(col).unwrap() )                );
            assert!( itertools::equal( (&a).column_reverse(col),                tran.row_reverse(col) )                     );
            assert!( itertools::equal( (&a).column_reverse_opt(col).unwrap(),   tran.row_reverse_opt(col).unwrap() )        );
            
            assert!( itertools::equal( (&a).column(col),                        anti.row_reverse(col) )                     );
            assert!( itertools::equal( (&a).column_opt(col).unwrap(),           anti.row_reverse_opt(col).unwrap() )        );
            assert!( itertools::equal( (&a).column_reverse(col),                anti.row(col) )                             );
            assert!( itertools::equal( (&a).column_reverse_opt(col).unwrap(),   anti.row_opt(col).unwrap() )                );
            
        }        


    }



    // Simplify entry-iterators
    // =====================================================================================================================

    // * SEE DOC TESTS    

}




