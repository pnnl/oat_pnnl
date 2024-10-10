//! Lazy transpose and anti-transpose; wraps around another matrix, and swaps order and/or major vs. minor views.

use crate::{algebra::matrices::query::{ViewRowAscend, ViewColDescend, ViewRowDescend, ViewColAscend, IndicesAndCoefficients, ViewRow, ViewCol, MatrixOracle}, };






//  REVERSE
//  -----------------------------------------------------------------------------------------------

/// Wraps around a matrix; returns entries of each row (respectively, column) in reverse order (warning: doesn't work quite the same as for ordinary dense matrices)
/// 
/// Concretely, the antitranpose is obtained by (i) transposing, and then (ii) reversing the order of rows and reversing the order of columns.
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  For example, to execute the method `row`, an `Reverse` struct simply calls `row_reverse` on
/// the underlying matrix oracle.   **Note that this doesn't work quite the same as for ordinary dense matrices when indices happen to be integers.**
/// 
/// # Caution
/// 
/// There are three important differences between [Reverse] and the matrix returned by [antitranspose_deep](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep):
/// 
/// - [Reverse] is a lazy object that does not generate new data
/// - The set of (key,val) pairs that appear in a major (respectively, minor) view of `Reverse::new(matrix)`
///   are the *same* as the entries in a minor (respectively, major) view of `matrix`; only the sequence in which those entries appear is different.
///   By contrast, the keys in the (key,val) pairs of [matrix.antitranspose_deep()](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) are different;
///   they are obtained by subtracting the original keys from (# rows in the antitransposed matrix - 1).
/// - For this reason, [matrix.antitranspose_deep()](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) is only available for
///   very specific types of matrices; [Reverse] is available for a much broader class.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{reverse::Reverse, transpose::Transpose};
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
/// // its reverse
/// let rev    =   Reverse::new( &a );
/// 
/// // check that rows are reversed correctly
/// for row in 0 .. 2 {
///     assert!( itertools::equal( (&a).row_reverse(row),                   rev.row( row) )                         );
///     assert!( itertools::equal( (&a).row_reverse_opt(row).unwrap(),      rev.row_opt( row).unwrap() )            ); 
///     assert!( itertools::equal( (&a).row(row),                           rev.row_reverse( row) )                 );
///     assert!( itertools::equal( (&a).row_opt(row).unwrap(),              rev.row_reverse_opt( row).unwrap() )    );
/// }
/// 
/// // check that columns are reversed correctly
/// for col in 0 .. 3 {
///     assert!( itertools::equal( (&a).column(col),                        rev.column_reverse(col) )                     );
///     assert!( itertools::equal( (&a).column_opt(col).unwrap(),           rev.column_reverse_opt(col).unwrap() )        );
///     assert!( itertools::equal( (&a).column_reverse(col),                rev.column(col) )                             );
///     assert!( itertools::equal( (&a).column_reverse_opt(col).unwrap(),   rev.column_opt(col).unwrap() )                );
/// }  
/// ```
pub struct Reverse< Matrix > { unreversed: Matrix }

impl < Matrix >

    Reverse
        < Matrix >
{
    pub fn new( unreversed: Matrix ) -> Reverse< Matrix > { Reverse { unreversed } }
}


    //  CLONE
impl < Matrix: Clone > 

    Clone for

    Reverse< Matrix >

{ fn clone(&self) -> Self { Reverse { unreversed: self.unreversed.clone() } } }    

    //  MATRIX ORACLE
impl< Matrix > 

    MatrixOracle for 
    
    Reverse< Matrix >

    where
        Matrix:     MatrixOracle
{
    type Coefficient        =   Matrix::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   Matrix::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   Matrix::ColumnIndex;       // The type of column indices
    
    type RowEntry           =   Matrix::RowEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   Matrix::ColumnEntry;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Matrix::RowReverse;               // What you get when you ask for a row.
    type RowIter            =   Matrix::RowReverseIter;           // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse         =   Matrix::Row;        // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter     =   Matrix::RowIter;    // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    
    type Column             =   Matrix::ColumnReverse;            // What you get when you ask for a column
    type ColumnIter         =   Matrix::ColumnReverseIter;        // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse      =   Matrix::Column;     // What you get when you ask for a column with the order of entries reversed 
    type ColumnReverseIter  =   Matrix::ColumnIter; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    fn entry(                   & self,  row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient > {
        self.unreversed.entry( row, column )
    }

    fn row(                     & self,  index: Self::RowIndex    )       -> Self::Row 
        { self.unreversed.row_reverse(index) }
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>
        { self.unreversed.row_reverse_opt(index) }     
    fn row_reverse(             & self,  index: Self::RowIndex    )       -> Self::RowReverse
        { self.unreversed.row(index) }    
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>
        { self.unreversed.row_opt(index) }
    
    fn column(                  & self,  index: Self::ColumnIndex )       -> Self::Column
        { self.unreversed.column_reverse(index) }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column>
        { self.unreversed.column_reverse_opt(index) }    
    fn column_reverse(          & self,  index: Self::ColumnIndex )       -> Self::ColumnReverse
        { self.unreversed.column(index) }
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse>
        { self.unreversed.column_opt(index) }
} 