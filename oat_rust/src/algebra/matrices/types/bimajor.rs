//! Stores two matrices, A and B; returns the major views of B as minor views.
//! 
//! For example, calling `view_minor_ascend` on this struct will return `B.view_major_ascend`.
//! 
//! This module offers two data structures:
//! 
//! - `MatrixBimajor` holds a pair of oracles
//! - `MatrixBimajorData` holds a pair of objects that become oracles when references (such as [VecOfVec](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec).)

use derive_getters::{Getters, Dissolve};
use derive_new::new;

use crate::algebra::matrices::query::{IndicesAndCoefficients, ViewRowAscend, ViewRowDescend, ViewColAscend, ViewColDescend, ViewRow, ViewCol, MatrixOracle};



//  ==========================================================
//  BIMAJOR (REQUIRES TWO ORACLES)
//  ==========================================================

/// Stores two matrices, A and B; returns the major views of B as minor views.
/// 
/// For example, calling `view_minor_ascend` on this struct will return `B.view_major_ascend`.
/// 
/// # See also
/// 
/// Also check out [MatrixBimajorData], which serves a similar function.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::MatrixOperations;   
/// use oat_rust::algebra::matrices::query::{ViewRowAscend, ViewRowDescend, ViewColDescend};
/// use oat_rust::algebra::matrices::query::{ViewRow, ViewColAscend, ViewCol};        
/// use oat_rust::algebra::matrices::types::bimajor::MatrixBimajor;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// 
/// use itertools::Itertools;        
///   
/// let matrix_data = vec![ 
///         vec![ (0, 1.),                          ], 
///         vec![ (0,-1.),  (1 ,1.),                ], 
///         vec![           (1,-1.),    (2, 1.),    ], 
///         vec![                       (2,-1.),    ],                 
///     ];
///      
/// // define a matrix of form
/// // |  1   0   0 |
/// // | -1   1   0 |
/// // |  0  -1   1 |
/// // |  0   0  -1 |        
/// let matrix_data = VecOfVec::new(matrix_data);
/// let matrix = & matrix_data;
/// let transpose   = matrix.transpose();
/// let bimajor     =   MatrixBimajor{ matrix_major: &matrix, matrix_minor: & transpose };
///    
/// // check that major views agree
/// for index in 0..4 {
///     assert_eq!(     matrix.view_major_ascend(index).collect_vec(), 
///                     bimajor.view_major_ascend(index).collect_vec(),     );
///     assert_eq!(     matrix.view_major_descend(index).collect_vec(), 
///                     bimajor.view_major_descend(index).collect_vec(),    );                            
///     assert_eq!(     matrix.view_major(index).collect_vec(), 
///                     bimajor.view_major(index).collect_vec(),            ); 
/// }
/// 
/// // check that minor  views agree
/// for index in 0..3 {
///     assert_eq!(     matrix.view_minor_ascend(index).collect_vec(), 
///                     bimajor.view_minor_ascend(index).collect_vec(),     );
///     assert_eq!(     matrix.view_minor_descend(index).collect_vec(), 
///                     bimajor.view_minor_descend(index).collect_vec(),    );                            
///     assert_eq!(     matrix.view_minor(index).collect_vec(), 
///                     bimajor.view_minor(index).collect_vec(),            );
/// }   
/// ```        
#[derive(new,Clone,Copy,Debug,Getters,Dissolve,Eq,PartialEq)]
pub struct MatrixBimajor< MatrixMajor, MatrixMinor >
{
    pub matrix_major:  MatrixMajor,
    pub matrix_minor:  MatrixMinor,
}


//  MATRIX ORACLE
//  ---------------------------------------------------------------------

impl < MatrixMajor, MatrixMinor >

    MatrixOracle for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    MatrixOracle,
        MatrixMinor:    MatrixOracle< 
                                RowIndex        =   MatrixMajor::ColumnIndex, 
                                ColumnIndex     =   MatrixMajor::RowIndex,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >,    
{
    type Coefficient    =   MatrixMajor::Coefficient;
    
    type RowIndex       =   MatrixMajor::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex    =   MatrixMajor::ColumnIndex;       // The type of column indices
    
    type RowEntry       =   MatrixMajor::RowEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry    =   MatrixMinor::RowEntry;          // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    

    
    type Row =              // What you get when you ask for a row.
                            MatrixMajor::Row;
    type RowIter =          // What you get when you call `row.into_iter()`, where `row` is a row
                            MatrixMajor::RowIter;
    type RowReverse =       // What you get when you ask for a row with the order of entries reversed
                            MatrixMajor::RowReverse;
    type RowReverseIter =   // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
                            MatrixMajor::RowReverseIter;
    type Column =           // What you get when you ask for a column
                            MatrixMinor::Row;
    type ColumnIter =       // What you get when you call `column.into_iter()`, where `column` is a column
                            MatrixMinor::RowIter;
    type ColumnReverse =    // What you get when you ask for a column with the order of entries reversed                             
                            MatrixMinor::RowReverse;
    type ColumnReverseIter =// What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)
                            MatrixMinor::RowReverseIter;

    fn entry(                   & self,  row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient >
        { self.matrix_major.entry( row, column ) }

    fn row(                     & self,  index: Self::RowIndex    )   -> Self::Row
        { self.matrix_major.row( index ) }
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>
        { self.matrix_major.row_opt( index ) }
    fn row_reverse(             & self,  index: Self::RowIndex    )   -> Self::RowReverse
        { self.matrix_major.row_reverse( index ) }
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>
        { self.matrix_major.row_reverse_opt( index ) }    
    
    fn column(                  & self,  index: Self::ColumnIndex )   -> Self::Column
        { self.matrix_minor.row(index) }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column>
        { self.matrix_minor.row_opt(index) }    
    fn column_reverse(          & self,  index: Self::ColumnIndex )   -> Self::ColumnReverse
        { self.matrix_minor.row_reverse(index) }
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse>
        { self.matrix_minor.row_reverse_opt(index) }    
        
}



//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------

impl < MatrixMajor, MatrixMinor >

    IndicesAndCoefficients for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    IndicesAndCoefficients,
        MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   MatrixMajor::ColIndex, 
                                ColIndex        =   MatrixMajor::RowIndex,
                                EntryMajor      =   MatrixMajor::EntryMinor,
                                EntryMinor      =   MatrixMajor::EntryMajor,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >,    
{
    type Coefficient    =   MatrixMajor::Coefficient;
    type EntryMajor     =   MatrixMajor::EntryMajor;
    type EntryMinor     =   MatrixMajor::EntryMinor;
    type RowIndex       =   MatrixMajor::RowIndex;
    type ColIndex       =   MatrixMajor::ColIndex;
}


//  ASCENDING MAJOR VIEWS
//  ---------------------------------------------------------------------


impl < MatrixMajor, MatrixMinor >

    ViewRow for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    IndicesAndCoefficients + ViewRow,
        MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   MatrixMajor::ColIndex, 
                                ColIndex        =   MatrixMajor::RowIndex,
                                EntryMajor      =   MatrixMajor::EntryMinor,
                                EntryMinor      =   MatrixMajor::EntryMajor,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >,    
{
    type ViewMajor            =   MatrixMajor::ViewMajor;
    type ViewMajorIntoIter    =   MatrixMajor::ViewMajorIntoIter;    

    fn   view_major( &self, index: Self::RowIndex ) -> Self::ViewMajor {
        self.matrix_major.view_major(index)
    }
}


//  ASCENDING MAJOR VIEWS
//  ---------------------------------------------------------------------


impl < MatrixMajor, MatrixMinor >

    ViewRowAscend for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    IndicesAndCoefficients + ViewRowAscend,
        MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   MatrixMajor::ColIndex, 
                                ColIndex        =   MatrixMajor::RowIndex,
                                EntryMajor      =   MatrixMajor::EntryMinor,
                                EntryMinor      =   MatrixMajor::EntryMajor,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >,    
{
    type ViewMajorAscend            =   MatrixMajor::ViewMajorAscend;
    type ViewMajorAscendIntoIter    =   MatrixMajor::ViewMajorAscendIntoIter;    

    fn   view_major_ascend( &self, index: Self::RowIndex ) -> Self::ViewMajorAscend {
        self.matrix_major.view_major_ascend(index)
    }
}

//  DESCENDING MAJOR VIEWS
//  ---------------------------------------------------------------------


impl < MatrixMajor, MatrixMinor >

    ViewRowDescend for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    IndicesAndCoefficients + ViewRowDescend,
        MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   MatrixMajor::ColIndex, 
                                ColIndex        =   MatrixMajor::RowIndex,
                                EntryMajor      =   MatrixMajor::EntryMinor,
                                EntryMinor      =   MatrixMajor::EntryMajor,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >,    
{
    type ViewMajorDescend            =   MatrixMajor::ViewMajorDescend;
    type ViewMajorDescendIntoIter    =   MatrixMajor::ViewMajorDescendIntoIter;    

    fn   view_major_descend( &self, index: Self::RowIndex ) -> Self::ViewMajorDescend {
        self.matrix_major.view_major_descend(index)
    }
}


//  ASCENDING MINOR VIEWS
//  ---------------------------------------------------------------------


impl < MatrixMajor, MatrixMinor >

    ViewCol for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    IndicesAndCoefficients,
        MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   MatrixMajor::ColIndex, 
                                ColIndex        =   MatrixMajor::RowIndex,
                                EntryMajor      =   MatrixMajor::EntryMinor,
                                EntryMinor      =   MatrixMajor::EntryMajor,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >
                        + ViewRow,
{
    type ViewMinor            =   MatrixMinor::ViewMajor;
    type ViewMinorIntoIter    =   MatrixMinor::ViewMajorIntoIter;    

    fn   view_minor( &self, index: Self::ColIndex ) -> Self::ViewMinor {
        self.matrix_minor.view_major(index)
    }
}



//  ASCENDING MINOR VIEWS
//  ---------------------------------------------------------------------


impl < MatrixMajor, MatrixMinor >

    ViewColAscend for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    IndicesAndCoefficients,
        MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   MatrixMajor::ColIndex, 
                                ColIndex        =   MatrixMajor::RowIndex,
                                EntryMajor      =   MatrixMajor::EntryMinor,
                                EntryMinor      =   MatrixMajor::EntryMajor,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >
                        + ViewRowAscend,
{
    type ViewMinorAscend            =   MatrixMinor::ViewMajorAscend;
    type ViewMinorAscendIntoIter    =   MatrixMinor::ViewMajorAscendIntoIter;    

    fn   view_minor_ascend( &self, index: Self::ColIndex ) -> Self::ViewMinorAscend {
        self.matrix_minor.view_major_ascend(index)
    }
}

//  DESCENDING MINOR VIEWS
//  ---------------------------------------------------------------------


impl < MatrixMajor, MatrixMinor >

    ViewColDescend for 

    MatrixBimajor
        < MatrixMajor, MatrixMinor > 

    where 
        MatrixMajor:    IndicesAndCoefficients,
        MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   MatrixMajor::ColIndex, 
                                ColIndex        =   MatrixMajor::RowIndex,
                                EntryMajor      =   MatrixMajor::EntryMinor,
                                EntryMinor      =   MatrixMajor::EntryMajor,
                                Coefficient     =   MatrixMajor::Coefficient,
                            >
                        + ViewRowDescend,  
{
    type ViewMinorDescend            =   MatrixMinor::ViewMajorDescend;
    type ViewMinorDescendIntoIter    =   MatrixMinor::ViewMajorDescendIntoIter;    

    fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend {
        self.matrix_minor.view_major_descend(index)
    }
}





//  ==========================================================
//  BIMAJOR DATA (REQUIRES **DATA** FOR TWO ORACLES)
//  ==========================================================





//  ==========================================================
//  BIMAJOR (REQUIRES TWO ORACLES)
//  ==========================================================

/// Similar to [MatrixBimajor], but works with objects that only become oracles when referenced.
/// 
/// Stores two matrices, A and B; returns the major views of B as minor views.
/// For example, calling `view_minor_ascend` on this struct will return `B.view_major_ascend`.
/// 
/// The key distinction with [MatrixBimajor] is as follows:
/// - [MatrixBimajor] 
///   - is itself an oracle
///   - holds two oracles, A and B
/// - [MatrixBimajorData]
///   - is not an oracle, but [& MatrixBimajorData](MatrixBimajorData) is
///   - holds two objects A and B such that &A and &B *are* oracles (even if A and B are not)
/// 
/// # Motivation
/// 
/// We use this matrix type for several operations related to U-match factorization.
/// (Were it not for this important application, we would likely delete [MatrixBimajorData] from the library.)
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::matrices::query::{ViewRowAscend, ViewRowDescend, ViewColDescend};
/// use oat_rust::algebra::matrices::query::{ViewRow, ViewColAscend, ViewCol};        
/// use oat_rust::algebra::matrices::types::bimajor::MatrixBimajorData;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// 
/// use itertools::Itertools;   
/// 
/// 
/// // data for a matrix:
/// 
/// // |  1   0   0 |
/// // | -1   1   0 |
/// // |  0  -1   1 |
/// // |  0   0  -1 |           
/// let matrix_data = vec![ 
///         vec![ (0,1.),                       ], 
///         vec![ (0,-1.),  (1,1.),             ], 
///         vec![           (1,-1.), (2,1.),    ], 
///         vec![                    (2,-1.),   ],                 
///     ];
/// 
/// // define matrices
/// let matrix_data = VecOfVec::new(matrix_data);
/// let transpose_data   = matrix_data.transpose_deep(3).unwrap(); // a non-lazy transpose; copies the data
/// let bimajor_data     =   MatrixBimajorData{ matrix_major_data: matrix_data.clone(), matrix_minor_data: transpose_data };
/// let matrix  = & matrix_data;
/// 
/// // check that major views agree
/// for index in 0..4 {
///     assert_eq!(     matrix.view_major_ascend(index).collect_vec(), 
///                     ( & bimajor_data ).view_major_ascend(index).collect_vec(),     );
///     assert_eq!(     matrix.view_major_descend(index).collect_vec(), 
///                     ( & bimajor_data ).view_major_descend(index).collect_vec(),    );                            
///     assert_eq!(     matrix.view_major(index).collect_vec(), 
///                     ( & bimajor_data ).view_major(index).collect_vec(),            ); 
/// }
/// 
/// // check that minor  views agree
/// for index in 0..3 {
///     assert_eq!(     matrix.view_minor_ascend(index).collect_vec(), 
///                     ( & bimajor_data ).view_minor_ascend(index).collect_vec(),     );
///     assert_eq!(     matrix.view_minor_descend(index).collect_vec(), 
///                     ( & bimajor_data ).view_minor_descend(index).collect_vec(),    );                            
///     assert_eq!(     matrix.view_minor(index).collect_vec(), 
///                     ( & bimajor_data ).view_minor(index).collect_vec(),            );
/// }                         
/// ```        
#[derive(new,Clone,Copy,Debug,Getters,Dissolve,Eq,PartialEq)]
pub struct MatrixBimajorData< MatrixMajor, MatrixMinor >
{
    pub matrix_major_data:  MatrixMajor,
    pub matrix_minor_data:  MatrixMinor,
}


//  MATRIX ORACLE
//  ---------------------------------------------------------------------

impl < 'a, MatrixMajor, MatrixMinor >

    MatrixOracle for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    MatrixOracle,
        &'a MatrixMinor:    MatrixOracle< 
                                RowIndex        =   < &'a MatrixMajor as MatrixOracle >::ColumnIndex, 
                                ColumnIndex     =   < &'a MatrixMajor as MatrixOracle >::RowIndex,
                                Coefficient     =   < &'a MatrixMajor as MatrixOracle >::Coefficient,
                            >,    
{
    type Coefficient        =   <&'a MatrixMajor as MatrixOracle >::Coefficient;    
    type RowIndex           =   <&'a MatrixMajor as MatrixOracle >::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   <&'a MatrixMajor as MatrixOracle >::ColumnIndex;       // The type of column indices    
    type RowEntry           =   <&'a MatrixMajor as MatrixOracle >::RowEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   <&'a MatrixMinor as MatrixOracle >::RowEntry;          // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`                            
    type Row                =   <&'a MatrixMajor as MatrixOracle >::Row;       // What you get when you ask for a row.                            
    type RowIter            =   <&'a MatrixMajor as MatrixOracle >::RowIter;   // What you get when you call `row.into_iter()`, where `row` is a row                            
    type RowReverse         =   <&'a MatrixMajor as MatrixOracle >::RowReverse;// What you get when you ask for a row with the order of entries reversed                            
    type RowReverseIter     =   <&'a MatrixMajor as MatrixOracle >::RowReverseIter;// What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)                                
    type Column             =   <&'a MatrixMinor as MatrixOracle >::Row;// What you get when you ask for a column                            
    type ColumnIter         =   <&'a MatrixMinor as MatrixOracle >::RowIter;// What you get when you call `column.into_iter()`, where `column` is a column                            
    type ColumnReverse      =   <&'a MatrixMinor as MatrixOracle >::RowReverse;// What you get when you ask for a column with the order of entries reversed                                                         
    type ColumnReverseIter  =   <&'a MatrixMinor as MatrixOracle >::RowReverseIter;// What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)
                        
    fn entry(                   & self,  row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient >
        { ( & self.matrix_major_data ).entry( row, column ) }
    fn row(                     & self,  index: Self::RowIndex    )   -> Self::Row
        { ( & self.matrix_major_data ).row( index ) }
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>
        { ( & self.matrix_major_data ).row_opt( index ) }
    fn row_reverse(             & self,  index: Self::RowIndex    )   -> Self::RowReverse
        { ( & self.matrix_major_data ).row_reverse( index ) }
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>
        { ( & self.matrix_major_data ).row_reverse_opt( index ) }        
    fn column(                  & self,  index: Self::ColumnIndex )   -> Self::Column
        { ( & self.matrix_minor_data ).row(index) }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column>
        { ( & self.matrix_minor_data ).row_opt(index) }    
    fn column_reverse(          & self,  index: Self::ColumnIndex )   -> Self::ColumnReverse
        { ( & self.matrix_minor_data ).row_reverse(index) }
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse>
        { ( & self.matrix_minor_data ).row_reverse_opt(index) }    
        
}



//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------

impl < 'a, MatrixMajor, MatrixMinor >

    IndicesAndCoefficients for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    IndicesAndCoefficients,
        &'a MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex, 
                                ColIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex,
                                EntryMajor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor,
                                EntryMinor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor,
                                Coefficient     =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient,
                            >,    
{
    type Coefficient    =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient;
    type EntryMajor     =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor;
    type EntryMinor     =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor;
    type RowIndex       =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex;
    type ColIndex       =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex;
}


//  MAJOR VIEWS
//  ---------------------------------------------------------------------


impl < 'a, MatrixMajor, MatrixMinor >

    ViewRow for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    IndicesAndCoefficients + ViewRow,
        &'a MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex, 
                                ColIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex,
                                EntryMajor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor,
                                EntryMinor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor,
                                Coefficient     =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient,
                            >,   
{
    type ViewMajor            =   < &'a MatrixMajor as ViewRow>::ViewMajor;
    type ViewMajorIntoIter    =   < &'a MatrixMajor as ViewRow>::ViewMajorIntoIter;    

    fn   view_major( &self, index: Self::RowIndex ) -> Self::ViewMajor {
        return (&self.matrix_major_data).view_major(index)
    }
}


//  ASCENDING MAJOR VIEWS
//  ---------------------------------------------------------------------


impl < 'a, MatrixMajor, MatrixMinor >

    ViewRowAscend for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    IndicesAndCoefficients + ViewRowAscend,
        &'a MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex, 
                                ColIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex,
                                EntryMajor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor,
                                EntryMinor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor,
                                Coefficient     =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient,
                            >,  
{
    type ViewMajorAscend            =   < &'a MatrixMajor as ViewRowAscend>::ViewMajorAscend;
    type ViewMajorAscendIntoIter    =   < &'a MatrixMajor as ViewRowAscend>::ViewMajorAscendIntoIter;    

    fn   view_major_ascend( &self, index: Self::RowIndex ) -> Self::ViewMajorAscend {
        return (&self.matrix_major_data).view_major_ascend(index)
    }
}

//  DESCENDING MAJOR VIEWS
//  ---------------------------------------------------------------------


impl < 'a, MatrixMajor, MatrixMinor >

    ViewRowDescend for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    IndicesAndCoefficients + ViewRowDescend,
        &'a MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex, 
                                ColIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex,
                                EntryMajor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor,
                                EntryMinor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor,
                                Coefficient     =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient,
                            >,  
{
    type ViewMajorDescend            =   < &'a MatrixMajor as ViewRowDescend>::ViewMajorDescend;
    type ViewMajorDescendIntoIter    =   < &'a MatrixMajor as ViewRowDescend>::ViewMajorDescendIntoIter;      

    fn   view_major_descend( &self, index: Self::RowIndex ) -> Self::ViewMajorDescend {
        return (&self.matrix_major_data).view_major_descend(index)
    }
}


//  MINOR VIEWS
//  ---------------------------------------------------------------------


impl < 'a, MatrixMajor, MatrixMinor >

    ViewCol for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    IndicesAndCoefficients,
        &'a MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex, 
                                ColIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex,
                                EntryMajor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor,
                                EntryMinor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor,
                                Coefficient     =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient,
                            >
                            + ViewRow, 
{
    type ViewMinor            =   < &'a MatrixMinor as ViewRow>::ViewMajor;
    type ViewMinorIntoIter    =   < &'a MatrixMinor as ViewRow>::ViewMajorIntoIter;

    fn   view_minor( &self, index: Self::ColIndex ) -> Self::ViewMinor {
        return ( & self.matrix_minor_data ).view_major(index)
    }
}



//  ASCENDING MINOR VIEWS
//  ---------------------------------------------------------------------


impl < 'a, MatrixMajor, MatrixMinor >

    ViewColAscend for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    IndicesAndCoefficients,
        &'a MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex, 
                                ColIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex,
                                EntryMajor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor,
                                EntryMinor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor,
                                Coefficient     =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient,
                            >
                            + ViewRowAscend, 
{
    type ViewMinorAscend            =   < &'a MatrixMinor as ViewRowAscend>::ViewMajorAscend;
    type ViewMinorAscendIntoIter    =   < &'a MatrixMinor as ViewRowAscend>::ViewMajorAscendIntoIter;  

    fn   view_minor_ascend( &self, index: Self::ColIndex ) -> Self::ViewMinorAscend {
        return ( & self.matrix_minor_data ).view_major_ascend(index)
    }
}

//  DESCENDING MINOR VIEWS
//  ---------------------------------------------------------------------


impl < 'a, MatrixMajor, MatrixMinor >

    ViewColDescend for 

    &'a MatrixBimajorData
        < MatrixMajor, MatrixMinor > 

    where 
        &'a MatrixMajor:    IndicesAndCoefficients,
        &'a MatrixMinor:    IndicesAndCoefficients< 
                                RowIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::ColIndex, 
                                ColIndex        =   < &'a MatrixMajor as IndicesAndCoefficients >::RowIndex,
                                EntryMajor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMinor,
                                EntryMinor      =   < &'a MatrixMajor as IndicesAndCoefficients >::EntryMajor,
                                Coefficient     =   < &'a MatrixMajor as IndicesAndCoefficients >::Coefficient,
                            >
                            + ViewRowDescend, 
{
    type ViewMinorDescend            =   < &'a MatrixMinor as ViewRowDescend >::ViewMajorDescend;
    type ViewMinorDescendIntoIter    =   < &'a MatrixMinor as ViewRowDescend >::ViewMajorDescendIntoIter;   

    fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend {
        return ( & self.matrix_minor_data ).view_major_descend(index)
    }
}













//  ==========================================================
//  TESTS
//  ==========================================================







#[cfg(test)]
mod tests {
    use crate::algebra::matrices::types::bimajor::MatrixBimajorData;

    



    #[test]
    fn test_matrix_bimajor() {   

        use crate::algebra::matrices::operations::MatrixOperations;   
        use crate::algebra::matrices::query::{ViewRowAscend, ViewRowDescend, ViewColDescend};
        use crate::algebra::matrices::query::{ViewRow, ViewColAscend, ViewCol};        
        use crate::algebra::matrices::types::bimajor::MatrixBimajor;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;

        use itertools::Itertools;   


        // data for a matrix:

        // |  1   0   0 |
        // | -1   1   0 |
        // |  0  -1   1 |
        // |  0   0  -1 |           
        let matrix_data = vec![ 
                vec![ (0,1.),                       ], 
                vec![ (0,-1.),  (1,1.),             ], 
                vec![           (1,-1.), (2,1.),    ], 
                vec![                    (2,-1.),   ],                 
            ];
        
        // define matrices
        let matrix_data = VecOfVec::new(matrix_data);
        let matrix = & matrix_data;
        let transpose   = matrix.transpose();
        let bimajor     =   MatrixBimajor{ matrix_major: &matrix, matrix_minor: & transpose };
        
        // check that major views agree
        for index in 0..4 {
            assert_eq!(     matrix.view_major_ascend(index).collect_vec(), 
                            bimajor.view_major_ascend(index).collect_vec(),     );
            assert_eq!(     matrix.view_major_descend(index).collect_vec(), 
                            bimajor.view_major_descend(index).collect_vec(),    );                            
            assert_eq!(     matrix.view_major(index).collect_vec(), 
                            bimajor.view_major(index).collect_vec(),            ); 
        }

        // check that minor  views agree
        for index in 0..3 {
            assert_eq!(     matrix.view_minor_ascend(index).collect_vec(), 
                            bimajor.view_minor_ascend(index).collect_vec(),     );
            assert_eq!(     matrix.view_minor_descend(index).collect_vec(), 
                            bimajor.view_minor_descend(index).collect_vec(),    );                            
            assert_eq!(     matrix.view_minor(index).collect_vec(), 
                            bimajor.view_minor(index).collect_vec(),            );
        }                    
    }  








    #[test]
    fn test_matrix_bimajor_data() {   
 
        use crate::algebra::matrices::query::{ViewRowAscend, ViewRowDescend, ViewColDescend};
        use crate::algebra::matrices::query::{ViewRow, ViewColAscend, ViewCol};        
        use crate::algebra::matrices::types::bimajor::MatrixBimajorData;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;

        use itertools::Itertools;   


        // data for a matrix:

        // |  1   0   0 |
        // | -1   1   0 |
        // |  0  -1   1 |
        // |  0   0  -1 |           
        let matrix_data = vec![ 
                vec![ (0,1.),                       ], 
                vec![ (0,-1.),  (1,1.),             ], 
                vec![           (1,-1.), (2,1.),    ], 
                vec![                    (2,-1.),   ],                 
            ];
        
        // define matrices
        let matrix_data = VecOfVec::new(matrix_data);
        let transpose_data   = matrix_data.transpose_deep(3).unwrap(); // a non-lazy transpose; copies the data
        let bimajor_data     =   MatrixBimajorData{ matrix_major_data: matrix_data.clone(), matrix_minor_data: transpose_data };
        let matrix  = & matrix_data;

        // check that major views agree
        for index in 0..4 {
            assert_eq!(     matrix.view_major_ascend(index).collect_vec(), 
                            ( & bimajor_data ).view_major_ascend(index).collect_vec(),     );
            assert_eq!(     matrix.view_major_descend(index).collect_vec(), 
                            ( & bimajor_data ).view_major_descend(index).collect_vec(),    );                            
            assert_eq!(     matrix.view_major(index).collect_vec(), 
                            ( & bimajor_data ).view_major(index).collect_vec(),            ); 
        }

        // check that minor  views agree
        for index in 0..3 {
            assert_eq!(     matrix.view_minor_ascend(index).collect_vec(), 
                            ( & bimajor_data ).view_minor_ascend(index).collect_vec(),     );
            assert_eq!(     matrix.view_minor_descend(index).collect_vec(), 
                            ( & bimajor_data ).view_minor_descend(index).collect_vec(),    );                            
            assert_eq!(     matrix.view_minor(index).collect_vec(), 
                            ( & bimajor_data ).view_minor(index).collect_vec(),            );
        }                    
    }  





}

