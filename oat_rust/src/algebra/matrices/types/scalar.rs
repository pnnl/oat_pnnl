//! Memory-efficient scalar matrices (ie diagonal matrices where all entries on the diagonal are equal)

use crate::algebra::matrices::query::{   ViewRow,
                                        ViewRowAscend,
                                        ViewRowDescend,
                                        ViewCol, 
                                        ViewColAscend,
                                        ViewColDescend, IndicesAndCoefficients, MatrixOracle,
                                    };
use std::{iter::{once, Once}, marker::PhantomData};


//  ---------------------------------------------------------------------------
//  SCALAR MATRICES (INDICES CAN BE OF ANY TYPE)
//  ---------------------------------------------------------------------------


//  STRUCT
//  ------

/// Represents a scalar matrix.
///
/// Concretely, for any `index` (whether major or minor), each major/minor view
/// of the matrix returns an interator of form `Once< (index, scalar) >`, 
/// where `scalar` is the given scalar. 
/// 
/// Memory efficient, in the sense that this struct stores only one scalar value (rather than one scalar value for each row of the matrix).
///
/// # Examples
///
/// ```
/// use oat_rust::algebra::matrices::types::scalar::ScalarMatrix;
/// use oat_rust::algebra::matrices::query::{ViewRow};   
/// 
/// // create a scalar matrix indexed by `isize` integers
/// let two = ScalarMatrix::new( 2. );
/// assert!(itertools::equal( 
///         two.view_major( 0 ),   
///         std::iter::once( (0,2.) ),
///     ));  
/// ```
pub struct ScalarMatrix < Key, Val >
{
    scalar:         Val,
    phantom_key:    PhantomData< Key >
}

impl    < Key, Val >
        
        ScalarMatrix 
            < Key, Val > 
{
    /// Create new scalar matrix.
    pub fn new( scalar: Val ) 
            -> 
            Self  
    {
        ScalarMatrix { 
                scalar,
                phantom_key:        PhantomData,
            }
    }
}


//  ---------------------
//  TRAIT IMPLEMENTATIONS
//  ---------------------


//  MATRIX ORACLE
//  ---------------------------------------------------------------------------

impl    < Index, Scalar >

    MatrixOracle      for 

    ScalarMatrix < Index, Scalar >

    where
        Scalar:     Clone,
        Index:      Clone + Eq,     // Clone is needed for (Index,Scalar) to implement KeyValGet
                                    // Eq is needed for the `entry(row,column)` method, to determine if the row and column indices are equal

{
    type Coefficient            =      Scalar; // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =      Index; // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =      Index; // The type of column indices    
    type RowEntry               =      (Index,Scalar);  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =      (Index,Scalar);  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =      Once< (Index, Scalar) >;  // What you get when you ask for a row.
    type RowIter                =      Once< (Index, Scalar) >;  // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse             =      Once< (Index, Scalar) >;  // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter         =      Once< (Index, Scalar) >;  // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    type Column                 =      Once< (Index, Scalar) >;  // What you get when you ask for a column
    type ColumnIter             =      Once< (Index, Scalar) >;  // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse          =      Once< (Index, Scalar) >;  // What you get when you ask for a column with the order of entries reversed                             
    type ColumnReverseIter      =      Once< (Index, Scalar) >;  // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    fn entry(                   &   self, row: Self::RowIndex, column: Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        if row == column { 
            Some( self.scalar.clone() ) 
        } else { 
            None 
        }
    }
    fn row(                     &   self, index: Self::RowIndex   )   -> Self::Row                  
        { once( ( index, self.scalar.clone() ) ) } 
    fn row_opt(                 &   self, index: Self::RowIndex   )   -> Option<Self::Row>
        { Some( once( ( index, self.scalar.clone() ) ) ) }  
    fn row_reverse(             &   self, index: Self::RowIndex   )   -> Self::RowReverse
        { once( ( index, self.scalar.clone() ) ) } 
    fn row_reverse_opt(         &   self, index: Self::RowIndex   )   -> Option<Self::RowReverse>
        { Some( once( ( index, self.scalar.clone() ) ) ) } 
    fn column(                  &   self, index: Self::ColumnIndex)   -> Self::Column
        { once( ( index, self.scalar.clone() ) ) }     
    fn column_opt(              &   self, index: Self::ColumnIndex)   -> Option<Self::Column>
        { Some( once( ( index, self.scalar.clone() ) ) ) } 
    fn column_reverse(          &   self, index: Self::ColumnIndex)   -> Self::ColumnReverse
        { once( ( index, self.scalar.clone() ) ) } 
    fn column_reverse_opt(      &   self, index: Self::ColumnIndex)   -> Option<Self::ColumnReverse>
        { Some( once( ( index, self.scalar.clone() ) ) ) } 

} 




//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------

impl    < Key, Val >

    IndicesAndCoefficients      for 
    ScalarMatrix                < Key, Val >
{   type EntryMajor = (Key, Val);
    type EntryMinor = (Key, Val);
    type ColIndex = Key; 
    type RowIndex = Key; 
    type Coefficient = Val; }


//  MAJOR VIEWS
//  ---------------------------------------------------------------------------


//  ViewRow
//  
impl     < Key, Val >

        ViewRow       for        
        ScalarMatrix      < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMajor = Once < (Key, Val) >;
    type ViewMajorIntoIter = < Self::ViewMajor as IntoIterator >::IntoIter;

    fn view_major( & self, index: Key ) -> Once < (Key, Val) > 
    { 
        once( ( index, self.scalar.clone() ) )
    }
}

//  ViewRowAscend
//  
impl     < Key, Val >

        ViewRowAscend   for
        ScalarMatrix        < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMajorAscend = Once < (Key, Val) >;
    type ViewMajorAscendIntoIter = < Self::ViewMajorAscend as IntoIterator >::IntoIter;

    fn view_major_ascend( & self, index: Key ) -> Once < (Key, Val) >
    { 
        once( ( index, self.scalar.clone() ) )
    }
}


//  ViewRowDescend
//  
impl     < Key, Val >

        ViewRowDescend  for
        ScalarMatrix  < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMajorDescend           =   Once < (Key, Val) >;
    type ViewMajorDescendIntoIter   =   < Self::ViewMajorDescend as IntoIterator >::IntoIter; 

    fn view_major_descend( & self, index: Key ) -> Once < (Key, Val) >
    { 
        once( ( index, self.scalar.clone() ) )
    }
}


//  MINOR VIEWS
//  ---------------------------------------------------------------------------


//  ViewCol
//  
impl     < 'a, Key, Val >

        ViewCol       for
        ScalarMatrix      < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMinor          =   Once < (Key, Val) >;
    type ViewMinorIntoIter  =   Once < (Key, Val) >;  

    fn view_minor( & self, index: Key ) -> Once < (Key, Val) > 
    { 
        once( ( index, self.scalar.clone() ) )
    }
}

//  ViewColAscend
//  
impl     < 'a, Key, Val >

        ViewColAscend   for
        ScalarMatrix        < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMinorAscend           =   Once < (Key, Val) >;
    type ViewMinorAscendIntoIter   =   < Self::ViewMinorAscend as IntoIterator >::IntoIter;

    fn view_minor_ascend( & self, index: Key ) -> Once < (Key, Val) >
    { 
        once( ( index, self.scalar.clone() ) )
    }
}


//  ViewColDescend
//  
impl     < 'a, Key, Val >

        ViewColDescend      for
        ScalarMatrix            < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMinorDescend           =   Once < (Key, Val) >;
    type ViewMinorDescendIntoIter   =   Once < (Key, Val) >;

    fn view_minor_descend( & self, index: Key ) -> Once < (Key, Val) >
    { 
        once( ( index, self.scalar.clone() ) )
    }
}


//  ORACLE INHERITANCE (POSSIBLY DEPRECATED?)
//  ---------------------------------------------------------------------------

// impl    < 'a, Val: Clone >

//         OracleRefInherit    for
//         ScalarMatrix        < Key, Val >
// {}        


//  MAJOR DIMENSION
//  ---------------------------------------------------------------------------


//  WHICH MAJOR 
//  

// impl    < Key, Val >
        
//         WhichMajor  for 
//         ScalarMatrix      < Key, Val > 

// { fn major_dimension( &self ) -> MajorDimension { self.major_dimension.clone() } }









//  TESTS
//  =========================================================================================================


//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    

    #[test] 
    fn test_scalar_array() {
        use crate::algebra::matrices::types::scalar::ScalarMatrix;
        use crate::algebra::matrices::query::{ViewRow};   

        // create a scalar matrix indexed by `isize` integers
        let two = ScalarMatrix::new( 2. );
        assert!(itertools::equal( 
                two.view_major( 0 ),   
                std::iter::once( (0,2.) ),
            ));         
    }    
}


// DEPRECATED -- OK TO DELETE
// ================================================================================================

//  THE FOLLOWING IS DEPRECATED AND OK TO DELETE


// //  ---------------------------------------------------------------------------
// //  SCALAR MATRICES (INDEXED ONLY BY INTEGERS; SEE BELOW FOR A GENERALIZATION)
// //  ---------------------------------------------------------------------------
// 
// 
// //  STRUCT
// //  ------
// 
// /// Represents a scalar matrix indexed by integers (see below for a more general struct).
// ///
// /// Concretely, for any `index` (whether major or minor), each major/minor view
// /// of the matrix returns an interator of form `Once< (index, scalar) >` out, 
// /// where `scalar` is the given scalar. 
// ///
// /// # Examples
// ///
// /// ```
// /// use oat_rust::algebra::matrices::types::scalar_algebra::matrices::ScalarMatrixOracleUsize;
// /// use oat_rust::algebra::matrices::query::{ViewRow};   
// /// 
// /// // create a scalar matrix indexed by `usize` integers
// /// let two = ScalarMatrixOracleUsize::new( 2. );
// /// assert!(itertools::equal( 
// ///         two.view_major(0),   
// ///         std::iter::once( (0,2.) ),
// ///     ));
// /// ```
// pub struct ScalarMatrixOracleUsize < Key, Val >
// {
//     scalar: Val,
//     // major_dimension: MajorDimension,
// }
// 
// impl    < Key, Val >
//         ScalarMatrixOracleUsize
//         < Key, Val > 
// {
//     /// Create new scalar matrix.
//     pub fn new( 
//                     scalar: Val, 
//                     // major_dimension: MajorDimension 
//                 ) 
//             -> 
//             Self  
//     {
//         ScalarMatrixOracleUsize { 
//                 scalar:             scalar,
//                 // major_dimension:    major_dimension,
//             }
//     }
// }
// 
// 
// //  ---------------------
// //  TRAIT IMPLEMENTATIONS
// //  ---------------------
// 
// 
// //  MAJOR DIMENSION
// //  ---------------------------------------------------------------------------
// 
// 
// //  WHICH MAJOR 
// //  
// 
// // impl     < Key, Val >
//         
// //         WhichMajor                  for 
// //         ScalarMatrixOracleUsize     < Key, Val > 
// 
// // { fn major_dimension( &self ) -> MajorDimension { self.major_dimension.clone() } }
// 
// 
// //  MAJORS
// //  ---------------------------------------------------------------------------
// 
// 
// //  ViewRow
// //  
// impl     < 'a, Val >
// 
//         ViewRow < 'a, usize, Once< (usize, Val) > > for 
//         
//         ScalarMatrixOracleUsize < Key, Val > 
//         
//         where   Val: 'a + Clone, // hard to drop this requirement (tuples give move errors if no clone) 
// {
// 
//     fn view_major<'b: 'a>( &'b self, index: usize ) -> Once< (usize, Val) > 
//     { 
//         once( ( index, self.scalar.clone() ) )
//     }
// }
// 
// //  ViewRowAscend
// //  
// impl    < 'a, Val >
//         
//         ViewRowAscend < 'a, usize, Once< (usize, Val) > > for 
//         ScalarMatrixOracleUsize < Key, Val > 
//         
//         where   Val:    Clone, // hard to drop this requirement (tuples give move errors if no clone) 
// {
// 
// 
//     fn view_major_ascend<'b: 'a>( &'b self, index: usize ) -> Once< (usize, Val) >
//     { 
//         once( ( index, self.scalar.clone() ) )
//     }
// }
// 
// 
// //  ViewRowDescend
// //  
// impl     < 'a, Val >
// 
//         ViewRowDescend          < 'a, usize, Once < (usize, Val) > > for         
//         ScalarMatrixOracleUsize     < Key, Val > 
//         
//         where   Val:    Clone, // hard to drop this requirement (tuples give move errors if no clone) 
// {
// 
//     fn view_major_descend<'b: 'a>( &'b self, index: usize ) -> Once < (usize, Val) >
//     { 
//         once( ( index, self.scalar.clone() ) )
//     }
// }
// 
// 
// //  MINORS
// //  ---------------------------------------------------------------------------
// 
// 
// //  ViewCol
// //  
// impl     < 'a, Val >
//         
//         ViewCol                 < 'a, usize, Once < (usize, Val) > > for
//         ScalarMatrixOracleUsize     < Key, Val > 
//         
//         where   Val:    Clone, // hard to drop this requirement (tuples give move errors if no clone) 
// {
//     fn view_minor<'b: 'a>( &'b self, index: usize ) -> Once < (usize, Val) > 
//     { 
//         once( ( index, self.scalar.clone() ) )
//     }
// }
// 
// //  ViewColAscend
// //  
// impl     < 'a, Val >
//         
//         ViewColAscend           < 'a, usize, Once < (usize, Val) > > for
//         ScalarMatrixOracleUsize     < Key, Val > 
//         
//         where   Val:    Clone, // hard to drop this requirement (tuples give move errors if no clone) 
// {
//     fn view_minor_ascend<'b: 'a>( &'b self, index: usize ) -> Once < (usize, Val) >
//     { 
//         once( ( index, self.scalar.clone() ) )
//     }
// }
// 
// 
// //  ViewColDescend
// //  
// impl     < 'a, Val >
//         
//         ViewColDescend          < 'a, usize, Once < (usize, Val) > > for
//         ScalarMatrixOracleUsize     < Key, Val > 
//         
//         where   Val:    Clone, // hard to drop this requirement (tuples give move errors if no clone) 
// {
//     fn view_minor_descend<'b: 'a>( &'b self, index: usize ) -> Once < (usize, Val) >
//     { 
//         once( ( index, self.scalar.clone() ) )
//     }
// }











