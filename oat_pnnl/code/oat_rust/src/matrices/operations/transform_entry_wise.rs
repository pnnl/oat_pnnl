//! Create a new matrix by applying a transformation to each entry of an existing matrix.

use std::{iter::Cloned, marker::PhantomData};

use crate::{matrices::matrix_oracle_traits::{OracleMajorAscend, IndicesAndCoefficients, OracleMajor, OracleMinorDescend}, utilities::{functions::evaluate::EvaluateFunction, iterators::general::MapByTransform}, entries::{KeyValGet, KeyValNew, ReindexEntry}, vectors::operations::ChangeIndexSimple};


//  TRAIT: ENTRY-WISE TRANSFORMATIONS
//  ---------------------------------------------------------------------


/// Vector-wise transformations on matrix oracles.
///
pub trait TransformMajorAscend < ViewMajorAscend >

    where   Self:               OracleMajorAscend,
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
/// See the documentation for [IntoClonedEntries](crate::matrices::operations::transform_entry_wise::IntoClonedEntries), for an example.
pub struct ClonedEntries < MatrixUncloned >  
{
    matrix_uncloned:            MatrixUncloned,
}


// IndicesAndCoefficients
impl < MatrixUncloned: IndicesAndCoefficients >

    IndicesAndCoefficients for  

    ClonedEntries < MatrixUncloned >

    where
        MatrixUncloned:     IndicesAndCoefficients

{ type KeyMin = MatrixUncloned::KeyMin; type KeyMaj = MatrixUncloned::KeyMaj; type SnzVal = MatrixUncloned::SnzVal; } 



//  OracleMajorAscend
impl    < 'a, T, MatrixUncloned > 

        OracleMajorAscend for 

        ClonedEntries
            < MatrixUncloned > where

        MatrixUncloned:                         OracleMajorAscend< ViewMajorAscendEntry = &'a T > + IndicesAndCoefficients,
        MatrixUncloned::ViewMajorAscend:        IntoIterator,
        T:                                      Clone,    // these requirements reflect the requirements for Cloned to implement Iterator (see the source code for Cloned) 
        T:                                      'a,       // these requirements reflect the requirements for Cloned to implement Iterator (see the source code for Cloned)   

{
    type ViewMajorAscend            =   Cloned < MatrixUncloned::ViewMajorAscendIntoIter >;
    type ViewMajorAscendIntoIter    =   Cloned < MatrixUncloned::ViewMajorAscendIntoIter >;
    type ViewMajorAscendEntry       =   T;

    fn view_major_ascend( & self, majkey: Self::KeyMaj) -> Self::ViewMajorAscend { 
        self.matrix_uncloned.view_major_ascend( majkey ).into_iter().cloned()
    }
}

/// Trait to wrap a matrix oracle in a [`ClonedEntries`] struct, which wraps each view
/// of the matrix in a `Cloned` struct, thereby ensuring that each view iterator passes
/// a *clone* of its next item to the user, rather than the item itself.
pub trait IntoClonedEntries
        
    where 
        Self:                           Sized + OracleMajorAscend,
        Self::ViewMajorAscendEntry:     Clone,
{

    /// Wraps a matrix oracle in a [`ClonedEntries`] struct, which wraps each view
    /// of the matrix in a `Cloned` struct, thereby ensuring that each view iterator passes
    /// a *clone* of its next item to the user, rather than the item itself.
    /// 
    /// # Examples
    /// 
    /// ```
    /// // import crates
    /// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVec;
    /// use oat_rust::matrices::operations::transform_entry_wise::IntoClonedEntries;
    /// use oat_rust::matrices::matrix_oracle_traits::OracleMajorAscend;
    /// use oat_rust::utilities::partial_order::OrderComparatorAutoLtByKey;
    /// 
    /// // define the matrix and its entry-wise clone
    /// let matrix_origin   =   VecOfVec::new( 
    ///                                                             vec![ vec![ (0,0), (1,1), ] ],
    ///                                                             OrderComparatorAutoLtByKey::new(),
    ///                                                         );
    /// let matrix_origin_ref   =   &matrix_origin;
    /// let matrix_cloned   =   ( matrix_origin_ref ).cloned_entries();
    /// 
    /// // check that iterators return the same sequence of entries
    /// assert!(    itertools::equal(   
    ///                     ( & matrix_origin ).view_major_ascend(0).iter().cloned(),
    ///                     ( & matrix_cloned ).view_major_ascend(0),                            
    ///                 ));
    /// ```
    fn cloned_entries< 'a >( &'a self ) -> ClonedEntries < &'a Self > {
        ClonedEntries::< &'a Self > { matrix_uncloned: self }
    }
}

impl < 'a, MatrixUncloned, > 
    
    IntoClonedEntries for 

    MatrixUncloned where 

    MatrixUncloned:                         OracleMajorAscend,
    MatrixUncloned::ViewMajorAscendEntry:   Clone

{ } // the function `cloned_entries` is auto-implemented






//  REINDEX
//  ---------------------------------------------------------------------

/// Reindexes a matrix given a pair of functions that map old indices to new indices, and vice versa.
/// 
/// Concretely, if the original matrix `M` satisfies 
/// ```ignore
/// M[i,j] = t
/// ```
/// then the new matrix `N` satisfies 
/// ```ignore
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
            MatrixUnreindexed:          IndicesAndCoefficients< KeyMin = IndexOld, KeyMaj = IndexOld, >,
            ReindexFunctionOldToNew:    EvaluateFunction< MatrixUnreindexed::KeyMaj, IndexNew >,
            ReindexFunctionNewToOld:    EvaluateFunction< IndexNew, MatrixUnreindexed::KeyMaj >,            
{   type KeyMin = IndexNew;     type KeyMaj = IndexNew;     type SnzVal = MatrixUnreindexed::SnzVal;    }


//  OracleMajorAscend
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        OracleMajorAscend for 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        where
            MatrixUnreindexed:                          IndicesAndCoefficients< KeyMin = IndexOld, KeyMaj = IndexOld, >,
            MatrixUnreindexed:                          OracleMajorAscend,
            MatrixUnreindexed::ViewMajorAscend:         IntoIterator,
            MatrixUnreindexed::ViewMajorAscendEntry:    KeyValGet< MatrixUnreindexed::KeyMin, MatrixUnreindexed::SnzVal >,
            ReindexFunctionOldToNew:                    Clone + EvaluateFunction< MatrixUnreindexed::KeyMaj, IndexNew >,
            ReindexFunctionNewToOld:                    EvaluateFunction< IndexNew, MatrixUnreindexed::KeyMaj >,  
            EntryNew:                                   KeyValNew< IndexNew, MatrixUnreindexed::SnzVal >,                      
{
    type ViewMajorAscend           =   MapByTransform< 
                                                MatrixUnreindexed::ViewMajorAscendIntoIter, 
                                                EntryNew,
                                                ReindexEntry<   
                                                        MatrixUnreindexed::ViewMajorAscendEntry,
                                                        EntryNew,
                                                        MatrixUnreindexed::KeyMin,
                                                        IndexNew,
                                                        MatrixUnreindexed::SnzVal,
                                                        ReindexFunctionOldToNew,
                                                    >
                                            >;
    type ViewMajorAscendEntry      =   EntryNew;
    type ViewMajorAscendIntoIter   =   Self::ViewMajorAscend;

    fn view_major_ascend(&self, index: Self::KeyMin) -> Self::ViewMajorAscend {
        let selecting_index_new = (self.reindex_function_new_to_old).evaluate_function( index );
        let entry_reindexer = ReindexEntry::new( self.reindex_function_old_to_new.clone() );
        
        MapByTransform::new(
            self.matrix_unreindexed.view_major_ascend( selecting_index_new ).into_iter(),
            entry_reindexer,
        )
    }
}

//  OracleMinorDescend
impl    < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        OracleMinorDescend for 

        ReindexSquareMatrix
            < MatrixUnreindexed, ReindexFunctionOldToNew, ReindexFunctionNewToOld, IndexOld, IndexNew, EntryNew > 

        where
            MatrixUnreindexed:                          IndicesAndCoefficients< KeyMin = IndexOld, KeyMaj = IndexOld, >,
            MatrixUnreindexed:                          OracleMinorDescend,
            MatrixUnreindexed::ViewMinorDescend:        IntoIterator,
            MatrixUnreindexed::ViewMinorDescendEntry:   KeyValGet< MatrixUnreindexed::KeyMin, MatrixUnreindexed::SnzVal >,
            ReindexFunctionOldToNew:                    Clone + EvaluateFunction< MatrixUnreindexed::KeyMaj, IndexNew >,
            ReindexFunctionNewToOld:                    EvaluateFunction< IndexNew, MatrixUnreindexed::KeyMaj >,  
            EntryNew:                                   KeyValNew< IndexNew, MatrixUnreindexed::SnzVal >,          
{
    type ViewMinorDescend           =   MapByTransform< 
                                                MatrixUnreindexed::ViewMinorDescendIntoIter, 
                                                EntryNew,
                                                ReindexEntry<   
                                                        MatrixUnreindexed::ViewMinorDescendEntry,
                                                        EntryNew,
                                                        MatrixUnreindexed::KeyMin,
                                                        IndexNew,
                                                        MatrixUnreindexed::SnzVal,
                                                        ReindexFunctionOldToNew,
                                                    >
                                            >;
    type ViewMinorDescendEntry      =   EntryNew;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;

    fn view_minor_descend(&self, index: Self::KeyMin) -> Self::ViewMinorDescend {
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

    //  Entry-wise clone
    //  ---------------------------------------------------------------------
    
    #[test] 
    fn test_entrywise_clone() {
        // import crates
        use crate::matrices::matrix_types::vec_of_vec::VecOfVec;
        use crate::matrices::operations::transform_entry_wise::IntoClonedEntries;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use crate::utilities::partial_order::OrderComparatorAutoLtByKey;

        // define the matrix and its entry-wise clone
        let matrix_origin   =   VecOfVec::new( 
                                                                    vec![ vec![ (0,0), (1,1), ] ],
                                                                    OrderComparatorAutoLtByKey::new(),
                                                                );
        let matrix_origin_ref   =   &matrix_origin;
        let matrix_cloned   =   ( matrix_origin_ref ).cloned_entries();

        // check that iterators return the same sequence of entries
        assert!(    itertools::equal(   
                            ( & matrix_origin ).view_major_ascend(0).iter().cloned(),
                            ( matrix_cloned ).view_major_ascend(0),                            
                        ));
    }

}        
