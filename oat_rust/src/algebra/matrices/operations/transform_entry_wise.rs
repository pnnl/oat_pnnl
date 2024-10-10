//! Create a new matrix by applying a transformation to each entry of an existing matrix.

use std::{fmt::Debug, iter::Cloned, marker::PhantomData};

use crate::{algebra::matrices::query::{ViewRowAscend, IndicesAndCoefficients, ViewColDescend, MatrixOracle}, utilities::{functions::evaluate::EvaluateFunction, iterators::general::MapByTransform}};
use crate::algebra::vectors::entries::{KeyValGet, KeyValNew, ReindexEntry};


//  TRAIT: ENTRY-WISE TRANSFORMATIONS
//  ---------------------------------------------------------------------


/// Vector-wise transformations on matrix oracles.
///
pub trait TransformMajorAscend < ViewMajorAscend >

    where   Self:               ViewRowAscend,
            ViewMajorAscend:    IntoIterator

{
}


// struct MatrixPeekable {
//     matrix:             MatrixUnsimplified
// }



//  TRANSFORM ENTRY-WISE
//  ---------------------------------------------------------------------





//  CLONED ENTRIES
//  ---------------------------------------------------------------------

/// Wrapper for matrix oracles; applies `.simplify()` to each iterator returned by the matrix.
/// 
/// See the documentation for [IntoClonedRowEntries](crate::algebra::matrices::operations::transform_entry_wise::IntoClonedRowEntries), for an example.
pub struct ClonedRowEntries < MatrixUncloned >  
{
    matrix_uncloned:            MatrixUncloned,
}


impl < 'a, MatrixUncloned, RowEntryUnref, >

    MatrixOracle for

    ClonedRowEntries< MatrixUncloned >

    where
        MatrixUncloned:             MatrixOracle< RowEntry = &'a RowEntryUnref >,
        MatrixUncloned::RowIndex:   Clone,
        RowEntryUnref:              'a + Clone + Debug + KeyValGet< MatrixUncloned::ColumnIndex, MatrixUncloned::Coefficient >,
{

    type Coefficient        =   MatrixUncloned::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   MatrixUncloned::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   MatrixUncloned::ColumnIndex;       // The type of column indices
    
    type RowEntry           =   RowEntryUnref;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   MatrixUncloned::ColumnEntry;        // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Cloned< MatrixUncloned::RowIter >;               // What you get when you ask for a row.
    type RowIter            =   Cloned< MatrixUncloned::RowIter >;           // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse         =   Cloned< MatrixUncloned::RowReverseIter >;  // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter     =   Cloned< MatrixUncloned::RowReverseIter >;  // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    
    type Column             =   MatrixUncloned::Column;            // What you get when you ask for a column
    type ColumnIter         =   MatrixUncloned::ColumnIter;        // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse      =   MatrixUncloned::ColumnReverse;     // What you get when you ask for a column with the order of entries reversed 
    type ColumnReverseIter  =   MatrixUncloned::ColumnReverseIter; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    fn entry(                   & self,  row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient >
        { self.matrix_uncloned.entry( row, column ) }

    fn row(                     & self,  index: Self::RowIndex    )       -> Self::Row
        {   
            let a   =   (self.matrix_uncloned).row( index.clone() ).into_iter().cloned();
            let c: Vec<_> = a.collect();            
            println!("c = {:?}", c );
            let a   =   (self.matrix_uncloned).row( index.clone() ).into_iter().cloned();            
            return a
        }
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>
        { self.matrix_uncloned.row_opt( index ).map(|x| x.into_iter().cloned() ) }
    fn row_reverse(             & self,  index: Self::RowIndex    )       -> Self::RowReverse
        { self.matrix_uncloned.row_reverse( index ).into_iter().cloned() }
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>
        { self.matrix_uncloned.row_reverse_opt( index ).map(|x| x.into_iter().cloned() ) }    
    
    fn column(                  & self,  index: Self::ColumnIndex )       -> Self::Column
        { self.matrix_uncloned.column( index ) }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column>
        { self.matrix_uncloned.column_opt( index ) }
    fn column_reverse(          & self,  index: Self::ColumnIndex )       -> Self::ColumnReverse
        { self.matrix_uncloned.column_reverse( index ) }
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse> 
        { self.matrix_uncloned.column_reverse_opt( index ) }    
} 



// IndicesAndCoefficients
impl < 'a, T, MatrixUncloned >

    IndicesAndCoefficients for  

    ClonedRowEntries < MatrixUncloned >

    where
        MatrixUncloned:                         ViewRowAscend< EntryMajor = &'a T > + IndicesAndCoefficients,
        MatrixUncloned::ViewMajorAscend:        IntoIterator,
        T:                                      Clone,    // these requirements reflect the requirements for Cloned to implement Iterator (see the source code for Cloned) 
        T:                                      'a,       // these requirements reflect the requirements for Cloned to implement Iterator (see the source code for Cloned)           

{ 
    type EntryMajor = T;
    type EntryMinor = MatrixUncloned::EntryMinor;
    type ColIndex = MatrixUncloned::ColIndex; 
    type RowIndex = MatrixUncloned::RowIndex; 
    type Coefficient = MatrixUncloned::Coefficient; 
} 



//  ViewRowAscend
impl    < 'a, T, MatrixUncloned > 

        ViewRowAscend for 

        ClonedRowEntries
            < MatrixUncloned > 
            
        where
            MatrixUncloned:                         ViewRowAscend< EntryMajor = &'a T > + IndicesAndCoefficients,
            MatrixUncloned::ViewMajorAscend:        IntoIterator,
            T:                                      Clone,    // these requirements reflect the requirements for Cloned to implement Iterator (see the source code for Cloned) 
            T:                                      'a,       // these requirements reflect the requirements for Cloned to implement Iterator (see the source code for Cloned)   

{
    type ViewMajorAscend            =   Cloned < MatrixUncloned::ViewMajorAscendIntoIter >;
    type ViewMajorAscendIntoIter    =   Cloned < MatrixUncloned::ViewMajorAscendIntoIter >;

    fn view_major_ascend( & self, majkey: Self::RowIndex) -> Self::ViewMajorAscend { 
        self.matrix_uncloned.view_major_ascend( majkey ).into_iter().cloned()
    }
}

/// Trait to wrap a matrix oracle in a [`ClonedRowEntries`] struct, which wraps each view
/// of the matrix in a `Cloned` struct, thereby ensuring that each view iterator passes
/// a *clone* of its next item to the user, rather than the item itself.
pub trait IntoClonedRowEntries
        
    where 
        Self:                 Sized + ViewRowAscend,
        Self::EntryMajor:     Clone,
{

    /// Wraps a matrix oracle in a [`ClonedRowEntries`] struct, which wraps each view
    /// of the matrix in a `Cloned` struct, thereby ensuring that each view iterator passes
    /// a *clone* of its next item to the user, rather than the item itself.
    /// 
    /// # Examples
    /// 
    /// ```
    /// // import crates
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted_ref::VecOfVec;
    /// use oat_rust::algebra::matrices::operations::transform_entry_wise::IntoClonedRowEntries;
    /// use oat_rust::algebra::matrices::query::MatrixOracle;
    /// use oat_rust::utilities::order::OrderOperatorByKey;
    /// 
    /// // define the matrix and its entry-wise clone
    /// let matrix_origin   =   VecOfVec::new( 
    ///                                                             vec![ vec![ (0,0), (1,1), ] ],
    ///                                                         );
    /// let matrix_origin_ref   =   &matrix_origin;
    /// let matrix_cloned       =   ( matrix_origin_ref ).cloned_entries();
    /// 
    /// // check that iterators return the same sequence of entries
    /// assert!(    itertools::equal(   
    ///                     ( & matrix_origin ).row(0).cloned(),
    ///                     ( & matrix_cloned ).row(0),                            
    ///                 ));
    /// ```
    fn cloned_entries( self ) -> ClonedRowEntries < Self > {
        ClonedRowEntries::< Self > { matrix_uncloned: self }
    }
}

impl < MatrixUncloned, > 
    
    IntoClonedRowEntries for 

    MatrixUncloned where 

    MatrixUncloned:                         ViewRowAscend,
    MatrixUncloned::EntryMajor:   Clone

{ } // the function `cloned_entries` is auto-implemented






//  REINDEX
//  ---------------------------------------------------------------------

/// Reindexes a matrix given a pair of functions that map old indices to new indices, and vice versa.
/// 
/// Concretely, if the original matrix `M` satisfies 
/// ```text
/// M[i,j] = t
/// ```
/// then the new matrix `N` satisfies 
/// ```te
/// N[ f^{-1}(i), f^{-1}(j) ] = t
/// ```
pub struct ReindexSquareMatrix< MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
{
    matrix_unreindexed:             MatrixUnreindexed,
    reindex_function_old_to_new:    ReindexFunctionOldToNew,
    reindex_function_new_to_old:    ReindexFunctionNewToOld,    
    phantom_indexold:               PhantomData< IndexOld >,
    phantom_indexnew:               PhantomData< IndexNew >,
    phantom_entrynew:               PhantomData< EntryNew >,
}

//  Implement the struct
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
{
    pub fn new( matrix_unreindexed: MatrixUnreindexed, reindex_function_old_to_new: ReindexFunctionOldToNew, reindex_function_new_to_old: ReindexFunctionNewToOld )
        -> 
        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 
    {
        ReindexSquareMatrix{ matrix_unreindexed, reindex_function_old_to_new, reindex_function_new_to_old, phantom_indexnew: PhantomData, phantom_indexold: PhantomData, phantom_entrynew: PhantomData }
    }                    
}


//  IndicesAndCoefficients
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        IndicesAndCoefficients for 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        where
            MatrixUnreindexed:          IndicesAndCoefficients< ColIndex = IndexOld, RowIndex = IndexOld, >,
            ReindexFunctionOldToNew:    EvaluateFunction< MatrixUnreindexed::RowIndex, IndexNew >,
            ReindexFunctionNewToOld:    EvaluateFunction< IndexNew, MatrixUnreindexed::RowIndex >,            
{   
    type EntryMajor = EntryNew;    
    type EntryMinor = EntryNew;        
    type ColIndex = IndexNew;     
    type RowIndex = IndexNew;     
    type Coefficient = MatrixUnreindexed::Coefficient;    
}


//  ViewRowAscend
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        ViewRowAscend for 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        where
            MatrixUnreindexed:                          IndicesAndCoefficients< ColIndex = IndexOld, RowIndex = IndexOld, >,
            MatrixUnreindexed:                          ViewRowAscend,
            MatrixUnreindexed::ViewMajorAscend:         IntoIterator,
            MatrixUnreindexed::EntryMajor:    KeyValGet< MatrixUnreindexed::ColIndex, MatrixUnreindexed::Coefficient >,
            ReindexFunctionOldToNew:                    Clone + EvaluateFunction< MatrixUnreindexed::RowIndex, IndexNew >,
            ReindexFunctionNewToOld:                    EvaluateFunction< IndexNew, MatrixUnreindexed::RowIndex >,  
            EntryNew:                                   KeyValNew< IndexNew, MatrixUnreindexed::Coefficient >,                      
{
    type ViewMajorAscend           =   MapByTransform< 
                                                MatrixUnreindexed::ViewMajorAscendIntoIter, 
                                                EntryNew,
                                                ReindexEntry<   
                                                        MatrixUnreindexed::EntryMajor,
                                                        EntryNew,
                                                        MatrixUnreindexed::ColIndex,
                                                        IndexNew,
                                                        MatrixUnreindexed::Coefficient,
                                                        ReindexFunctionOldToNew,
                                                    >
                                            >;
    type ViewMajorAscendIntoIter   =   Self::ViewMajorAscend;

    fn view_major_ascend(&self, index: Self::ColIndex) -> Self::ViewMajorAscend {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index );
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.view_major_ascend( selecting_index_new ).into_iter(),
            entry_reindexer,
        )
    }
}

//  ViewColDescend
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        ViewColDescend for 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        where
            MatrixUnreindexed:                          IndicesAndCoefficients< ColIndex = IndexOld, RowIndex = IndexOld, >,
            MatrixUnreindexed:                          ViewColDescend,
            MatrixUnreindexed::ViewMinorDescend:        IntoIterator,
            MatrixUnreindexed::EntryMinor:   KeyValGet< MatrixUnreindexed::ColIndex, MatrixUnreindexed::Coefficient >,
            ReindexFunctionOldToNew:                    Clone + EvaluateFunction< MatrixUnreindexed::RowIndex, IndexNew >,
            ReindexFunctionNewToOld:                    EvaluateFunction< IndexNew, MatrixUnreindexed::RowIndex >,  
            EntryNew:                                   KeyValNew< IndexNew, MatrixUnreindexed::Coefficient >,          
{
    type ViewMinorDescend           =   MapByTransform< 
                                                MatrixUnreindexed::ViewMinorDescendIntoIter, 
                                                EntryNew,
                                                ReindexEntry<   
                                                        MatrixUnreindexed::EntryMinor,
                                                        EntryNew,
                                                        MatrixUnreindexed::ColIndex,
                                                        IndexNew,
                                                        MatrixUnreindexed::Coefficient,
                                                        ReindexFunctionOldToNew,
                                                    >
                                            >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;

    fn view_minor_descend(&self, index: Self::ColIndex) -> Self::ViewMinorDescend {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index );
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );

        MapByTransform::new(
            self.matrix_unreindexed.view_minor_descend( selecting_index_new ).into_iter(),
            entry_reindexer,
        )
    }
}       








//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {

}        
