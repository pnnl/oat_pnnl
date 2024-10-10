
//! Create a new matrix by applying a transformation to each row or column of an existing matrix.
//! 
//! See also [`transform_entry_wise`](oat_rust::algebra::matrices::operations::transform_entry_wise)



//  GENERAL VECTORWISE TRANSFORM
//  ---------------------------------------------------------------------

/// Wrapper for matrix oracles; applies `vector_transformer` to each sparse vector returned by the matrix.
/// 
/// # Examples
/// 
/// ```
/// // import crates
/// use oat_rust::algebra::matrices::operations::transform_vector_wise::VecWiseTransformed;        
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::ViewRowAscend;
/// use std::iter::Cloned;
/// use std::slice::Iter;
/// use itertools::Itertools;
/// 
/// // define a function that doubles the second entry in a tuple
/// fn double_coefficient( pair: (usize, usize) ) -> (usize, usize) { ( pair.0, 2 * pair.1) } 
/// fn vector_transformer( iter: Cloned< Iter< '_, (usize,usize) > > ) -> Vec< (usize, usize) > { 
///     iter
///     .map(|x| double_coefficient(x) )
///     .collect_vec() 
/// }
/// 
/// // define the untransformed matrix
/// let matrix_untransformed    =   VecOfVec::new( vec![ vec![ (0,1) ] ] );  
///     
/// // define the transformed matrix
/// let matrix_transformed      =   VecWiseTransformed::new( 
///                                         & matrix_untransformed,
///                                         vector_transformer,
///                                     );
/// let transformed_vector      =   ( matrix_transformed ).view_major_ascend( 0 ).into_iter().collect_vec();
/// assert_eq!( transformed_vector, vec![ (0,2) ] )
/// ```
pub struct VecWiseTransformed < MatrixUntransformed, VectorTransformer >  {
    untransformed_matrix:       MatrixUntransformed,
    vector_transformer:         VectorTransformer,
}

//  Implement the struct
impl < 'a, MatrixUntransformed, VectorTransformer >
        
    VecWiseTransformed 
    < MatrixUntransformed, VectorTransformer > 
{
    pub fn new( untransformed_matrix: MatrixUntransformed, vector_transformer: VectorTransformer )
        ->
        VecWiseTransformed < MatrixUntransformed, VectorTransformer > 
    {
        VecWiseTransformed{
            untransformed_matrix,
            vector_transformer,
        }        
    }
}

// IndicesAndCoefficients
impl < MatrixUntransformed, ViewMajorAscendTransformed, VectorTransformer, > 

    IndicesAndCoefficients for  

    VecWiseTransformed
        < MatrixUntransformed, VectorTransformer > where

    MatrixUntransformed:            IndicesAndCoefficients + ViewRowAscend,
    ViewMajorAscendTransformed:     IntoIterator,
    MatrixUntransformed:            ViewRowAscend,
    VectorTransformer:              Fn( MatrixUntransformed::ViewMajorAscend ) -> ViewMajorAscendTransformed,  

{
    type EntryMajor = ViewMajorAscendTransformed::Item;
    type EntryMinor = MatrixUntransformed::EntryMinor;
    type ColIndex = MatrixUntransformed::ColIndex; 
    type RowIndex = MatrixUntransformed::RowIndex; 
    type Coefficient = MatrixUntransformed::Coefficient; 
}  


//  ViewRowAscend
impl    < MatrixUntransformed, ViewMajorAscendTransformed, VectorTransformer, > 

        ViewRowAscend for 

        VecWiseTransformed
            < MatrixUntransformed, VectorTransformer > where

        MatrixUntransformed:            ViewRowAscend + IndicesAndCoefficients,
        ViewMajorAscendTransformed:     IntoIterator,
        MatrixUntransformed:            ViewRowAscend,
        VectorTransformer:              Fn( MatrixUntransformed::ViewMajorAscend ) -> ViewMajorAscendTransformed,                                    
{
    type ViewMajorAscend            =   ViewMajorAscendTransformed;
    type ViewMajorAscendIntoIter    =   ViewMajorAscendTransformed::IntoIter;

    fn view_major_ascend ( & self, majkey: Self::RowIndex ) -> Self::ViewMajorAscend { 
        (self.vector_transformer)(  self.untransformed_matrix.view_major_ascend( majkey )  )
    }
}




//  SIMPLIFY VECTORS
//  ---------------------------------------------------------------------




use crate::{algebra::matrices::query::{ViewRowAscend, ViewRow, IndicesAndCoefficients, ViewColDescend, ViewCol}, utilities::{iterators::general::PeekUnqualified, sets::MapHasKeyOrSequenceHasElement}, algebra::rings::operator_traits::Semiring, algebra::vectors::operations::{Simplify, VectorOperations}, };
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet};

/// Wrapper for matrix oracles; applies `.simplify()` to each iterator returned by the matrix.
struct MatrixSimplify 
        < MatrixUnsimplified, RingOperator >  
    where
        MatrixUnsimplified:     IndicesAndCoefficients
{
    ring_operator:              RingOperator,
    unsimplified_matrix:        MatrixUnsimplified,
}

// IndicesAndCoefficients

impl < MatrixUnsimplified: IndicesAndCoefficients, RingOperator >

    IndicesAndCoefficients for  

    MatrixSimplify
        < MatrixUnsimplified, RingOperator > where
                
    MatrixUnsimplified:                             ViewRowAscend + IndicesAndCoefficients,
    MatrixUnsimplified::ViewMajorAscend:            IntoIterator,
    // The next uncommented requirement is the same as the following commented one (Rust just has a hard time seeing that they are equivalent)
    // MatrixUnsimplified::EntryMajor:       KeyValGet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, > + 
    //                                                 KeyValSet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, >,        
    < MatrixUnsimplified::ViewMajorAscendIntoIter as IntoIterator >::Item:
                                                    KeyValGet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, > + 
                                                    KeyValSet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, >,        
    MatrixUnsimplified::ViewMajorAscendIntoIter:    PeekUnqualified,  
    MatrixUnsimplified::ColIndex:                     std::cmp::PartialEq,
    RingOperator:                                   Clone + Semiring< MatrixUnsimplified::Coefficient >,  

{ 
    type EntryMajor = MatrixUnsimplified::EntryMajor;
    type EntryMinor = MatrixUnsimplified::EntryMinor;    
    type RowIndex = MatrixUnsimplified::RowIndex; 
    type ColIndex = MatrixUnsimplified::ColIndex;     
    type Coefficient = MatrixUnsimplified::Coefficient; 
}  


// ViewRowAscend
impl    < MatrixUnsimplified, RingOperator > 

        ViewRowAscend for 

        MatrixSimplify
            < MatrixUnsimplified, RingOperator > where
                
        MatrixUnsimplified:                             ViewRowAscend + IndicesAndCoefficients,
        MatrixUnsimplified::ViewMajorAscend:            IntoIterator,
        // The next uncommented requirement is the same as the following commented one (Rust just has a hard time seeing that they are equivalent)
        // MatrixUnsimplified::EntryMajor:       KeyValGet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, > + 
        //                                                 KeyValSet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, >,        
        < MatrixUnsimplified::ViewMajorAscendIntoIter as IntoIterator >::Item:
                                                        KeyValGet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, > + 
                                                        KeyValSet < MatrixUnsimplified::ColIndex, MatrixUnsimplified::Coefficient, >,        
        MatrixUnsimplified::ViewMajorAscendIntoIter:    PeekUnqualified,  
        MatrixUnsimplified::ColIndex:                     std::cmp::PartialEq,
        RingOperator:                                   Clone + Semiring< MatrixUnsimplified::Coefficient >,
{
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;    
    type ViewMajorAscend            =   Simplify<
                                                MatrixUnsimplified::ViewMajorAscendIntoIter,
                                                MatrixUnsimplified::ColIndex,
                                                RingOperator,
                                                MatrixUnsimplified::Coefficient,
                                            > ;    

    fn view_major_ascend ( & self, majkey: Self::RowIndex) -> Self::ViewMajorAscend { 
        
        self.unsimplified_matrix
                .view_major_ascend( majkey )
                .into_iter()
                .simplify( self.ring_operator.clone() )
    }
}





//  MASK: ONLY INDICES OUTSIDE
//  -------------------------------------------------------------------------------------

//  MINOR KEYS
//  -----------------------------------------

use crate::algebra::vectors::operations::OnlyIndicesOutsideCollection;


/// A matrix oracle whose major views iterate only over entries that have indices **outside**
/// a given collection.
pub struct OnlyKeyMinOutsideCollection<
        Matrix, KeyMinToExclude,
    >
    where
        Matrix:                         IndicesAndCoefficients,    
        KeyMinToExclude:               MapHasKeyOrSequenceHasElement< Matrix::ColIndex >,                  
{
    matrix:                     Matrix,
    keymins_to_exclude:         KeyMinToExclude,
}

// implement the object
impl < Matrix, KeyMinToExclude, >

    OnlyKeyMinOutsideCollection
        < Matrix, KeyMinToExclude >

    where
        Matrix:                         IndicesAndCoefficients,    
        KeyMinToExclude:               MapHasKeyOrSequenceHasElement< Matrix::ColIndex >,                          

    {
        /// Create a new `OnlyKeyMinOutsideCollection` from a matrix and a collection
        /// of minor keys.
        /// 
        /// In order to operate properly, the collection of minor keys should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, keymins_to_exclude: KeyMinToExclude ) 
        -> 
        OnlyKeyMinOutsideCollection
            < Matrix, KeyMinToExclude >
        where
            Matrix:                         ViewRowAscend + IndicesAndCoefficients,    
            Matrix::ViewMajorAscend:        IntoIterator,
            KeyMinToExclude:               MapHasKeyOrSequenceHasElement< Matrix::ColIndex >,              
        { OnlyKeyMinOutsideCollection{ matrix, keymins_to_exclude, } }
    }    

// IndicesAndCoefficients

impl < Matrix: IndicesAndCoefficients, KeyMinToExclude >

    IndicesAndCoefficients for  

    OnlyKeyMinOutsideCollection
        < Matrix, KeyMinToExclude >

    where
        Matrix:                             IndicesAndCoefficients,
        KeyMinToExclude:                    Copy + MapHasKeyOrSequenceHasElement< Matrix::ColIndex >,            
{ 
    type EntryMajor = Matrix::EntryMajor;
    type EntryMinor = Matrix::EntryMinor;    
    type ColIndex = Matrix::ColIndex; 
    type RowIndex = Matrix::RowIndex; 
    type Coefficient = Matrix::Coefficient; 
}        


// ViewMajorAscend
impl < Matrix, KeyMinToExclude >

    ViewRowAscend for  

    OnlyKeyMinOutsideCollection
        < Matrix, KeyMinToExclude >

    where
        Matrix:                             ViewRowAscend + IndicesAndCoefficients,    
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::EntryMajor:                 KeyValGet< Matrix::ColIndex, Matrix::Coefficient >,
        KeyMinToExclude:                    Copy + MapHasKeyOrSequenceHasElement< Matrix::ColIndex >,  
        Matrix::ViewMajorAscendIntoIter:    Iterator,      
{
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscend            =   OnlyIndicesOutsideCollection< Matrix::ViewMajorAscendIntoIter, KeyMinToExclude, Self::ColIndex, Self::Coefficient >;  
    
    fn view_major_ascend( &self, keymaj: Self::RowIndex) 
        -> 
        OnlyIndicesOutsideCollection< Matrix::ViewMajorAscendIntoIter, KeyMinToExclude, Matrix::ColIndex, Matrix::Coefficient >  
    { 
        OnlyIndicesOutsideCollection::new( self.matrix.view_major_ascend(keymaj).into_iter(), self.keymins_to_exclude  ) 
    }
}

// ViewMajor
impl < Matrix, KeyMinToExclude >

    ViewRow for  

    OnlyKeyMinOutsideCollection
        < Matrix, KeyMinToExclude >

    where
        Matrix:                     ViewRow + IndicesAndCoefficients,    
        Matrix::ViewMajor:          IntoIterator,
        Matrix::EntryMajor:         KeyValGet< Matrix::ColIndex, Matrix::Coefficient >,
        KeyMinToExclude:            Copy + MapHasKeyOrSequenceHasElement< Matrix::ColIndex >,  
        Matrix::ViewMajorIntoIter:  Iterator,      
{
    type ViewMajorIntoIter    =   < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajor            =   OnlyIndicesOutsideCollection< Matrix::ViewMajorIntoIter, KeyMinToExclude, Self::ColIndex, Self::Coefficient >;  

    fn view_major( &self, keymaj: Self::RowIndex) -> Self::ViewMajor
    { 
        OnlyIndicesOutsideCollection::new( self.matrix.view_major(keymaj).into_iter(), self.keymins_to_exclude  ) 
    }
}



//  MAJOR KEYS
//  -----------------------------------------


/// A matrix oracle whose major views iterate only over entries that have indices **outside**
/// a given collection.
pub struct OnlyKeyMajOutsideCollection<
        Matrix, KeyMajToExclude,
    >
    where
        Matrix:                         IndicesAndCoefficients,    
        KeyMajToExclude:                MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                  
{
    matrix:                     Matrix,
    keymajs_to_exclude:         KeyMajToExclude,
}

// implement the object
impl < Matrix, KeyMajToExclude, >

    OnlyKeyMajOutsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                        IndicesAndCoefficients,    
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                          

    {
        /// Create a new `OnlyKeyMajOutsideCollection` from a matrix and a collection
        /// of minor keys.
        /// 
        /// In order to operate properly, the collection of minor keys should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, keymajs_to_exclude: KeyMajToExclude ) 
        -> 
        OnlyKeyMajOutsideCollection
            < Matrix, KeyMajToExclude >
        where
            Matrix:                         IndicesAndCoefficients,    
            KeyMajToExclude:                MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,              
        { OnlyKeyMajOutsideCollection{ matrix, keymajs_to_exclude, } }
    }    

// IndicesAndCoefficients

impl < Matrix: IndicesAndCoefficients, KeyMajToExclude >

    IndicesAndCoefficients for  

    OnlyKeyMajOutsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                        IndicesAndCoefficients,    
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,         

{ 
    type EntryMajor = Matrix::EntryMajor;
    type EntryMinor = Matrix::EntryMinor;
    type ColIndex = Matrix::ColIndex; 
    type RowIndex = Matrix::RowIndex; 
    type Coefficient = Matrix::Coefficient; 
}        


// ViewMinorDescend
impl < Matrix, KeyMajToExclude >

    ViewColDescend for  

    OnlyKeyMajOutsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                                 ViewColDescend + IndicesAndCoefficients,    
        Matrix::ViewMinorDescend:                                IntoIterator,
        Matrix::EntryMinor:                           KeyValGet< Matrix::RowIndex, Matrix::Coefficient >,
        KeyMajToExclude:                                       Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,  
        Matrix::ViewMinorDescendIntoIter:                        Iterator,      
{
    type ViewMinorDescendIntoIter    =   < Self::ViewMinorDescend as IntoIterator >::IntoIter;
    type ViewMinorDescend            =   OnlyIndicesOutsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Self::RowIndex, Self::Coefficient >;  
    
    fn view_minor_descend( &self, keymin: Self::ColIndex ) 
        -> 
        OnlyIndicesOutsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Matrix::RowIndex, Matrix::Coefficient >  
    { 
        OnlyIndicesOutsideCollection::new( self.matrix.view_minor_descend(keymin).into_iter(), self.keymajs_to_exclude  ) 
    }
}

// ViewMajor
impl < Matrix, KeyMajToExclude >

    ViewCol for  

    OnlyKeyMajOutsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                           ViewCol + IndicesAndCoefficients,    
        Matrix::ViewMinor:                                IntoIterator,
        Matrix::EntryMinor:                           KeyValGet< Matrix::RowIndex, Matrix::Coefficient >,
        KeyMajToExclude:                                 Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,  
        Matrix::ViewMinorIntoIter:                        Iterator,      
{
    type ViewMinorIntoIter    =   < Self::ViewMinor as IntoIterator >::IntoIter;
    type ViewMinor            =   OnlyIndicesOutsideCollection< Matrix::ViewMinorIntoIter, KeyMajToExclude, Self::RowIndex, Self::Coefficient >;  

    fn view_minor( &self, keymin: Self::ColIndex) -> Self::ViewMinor
    { 
        OnlyIndicesOutsideCollection::new( self.matrix.view_minor(keymin).into_iter(), self.keymajs_to_exclude  ) 
    }
}



//  MASK: ONLY INDICES INSIDE
//  -------------------------------------------------------------------------------------



//  MINOR KEYS
//  -----------------------------------------


use crate::algebra::vectors::operations::OnlyIndicesInsideCollection;


/// A matrix oracle whose major views iterate only over entries that have indices **outside**
/// a given collection.
#[derive(Copy, Clone)]
pub struct OnlyKeyMinInsideCollection
                < Matrix, KeyMinToInclude, >
{
    matrix:                     Matrix,
    keymins_to_include:         KeyMinToInclude,   
}

// implement the object
impl < Matrix, KeyMinToInclude >

    OnlyKeyMinInsideCollection
        < Matrix, KeyMinToInclude >

    {
        /// Create a new `OnlyKeyMinInsideColletion` from a matrix and a collection
        /// of minor keys.
        /// 
        /// In order to operate properly, the collection of minor keys should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, keymins_to_include: KeyMinToInclude ) 
        -> 
        OnlyKeyMinInsideCollection
            < Matrix, KeyMinToInclude, >
        where
            Matrix:                 IndicesAndCoefficients,    
            KeyMinToInclude:       MapHasKeyOrSequenceHasElement< Matrix::ColIndex >,              
        { OnlyKeyMinInsideCollection{ matrix, keymins_to_include } }
    }    

// IndicesAndCoefficients
impl < Matrix: IndicesAndCoefficients, KeyMinToInclude >

    IndicesAndCoefficients for  

    OnlyKeyMinInsideCollection
        < Matrix, KeyMinToInclude >

{ 
    type EntryMajor = Matrix::EntryMajor;
    type EntryMinor = Matrix::EntryMinor;
    type ColIndex = Matrix::ColIndex; 
    type RowIndex = Matrix::RowIndex; 
    type Coefficient = Matrix::Coefficient; 
} 


// ViewMajorAscend
impl < Matrix, KeyMinToInclude >

    ViewRowAscend for  

    OnlyKeyMinInsideCollection
        < Matrix, KeyMinToInclude >

    where
        Matrix:                             ViewRowAscend + IndicesAndCoefficients,    
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::ViewMajorAscendIntoIter:    Iterator,
        Matrix::EntryMajor:   KeyValGet< Matrix::ColIndex, Matrix::Coefficient >,        
        KeyMinToInclude:                   MapHasKeyOrSequenceHasElement< Self::ColIndex > + Copy,   
{
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscend            =   OnlyIndicesInsideCollection< Matrix::ViewMajorAscendIntoIter, KeyMinToInclude, Matrix::ColIndex, Matrix::Coefficient >;  

    fn view_major_ascend( &self, keymaj: Self::RowIndex ) 
        -> 
        Self::ViewMajorAscend
    { 
        OnlyIndicesInsideCollection::new( self.matrix.view_major_ascend(keymaj).into_iter(), self.keymins_to_include  ) 
    }
}


// ViewMajor
impl < Matrix, KeyMinToInclude, >

    ViewRow for  

    OnlyKeyMinInsideCollection
        < Matrix, KeyMinToInclude, >

    where
        Matrix:                     ViewRow + IndicesAndCoefficients,    
        Matrix::ViewMajor:          IntoIterator,
        Matrix::ViewMajorIntoIter:  Iterator,
        Matrix::EntryMajor:     KeyValGet< Matrix::ColIndex, Matrix::Coefficient >,        
        KeyMinToInclude:           MapHasKeyOrSequenceHasElement< Self::ColIndex > + Copy,     
{
    type ViewMajorIntoIter    =   < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajor            =   OnlyIndicesInsideCollection< Matrix::ViewMajorIntoIter, KeyMinToInclude, Matrix::ColIndex, Matrix::Coefficient >;  

    fn view_major( &self, keymaj: Self::RowIndex ) 
        -> 
        Self::ViewMajor
    { 
        OnlyIndicesInsideCollection::new( self.matrix.view_major(keymaj).into_iter(), self.keymins_to_include  ) 
    }
}




//  MAJOR KEYS
//  -----------------------------------------


/// A matrix oracle whose major views iterate only over entries that have indices **outside**
/// a given collection.
pub struct OnlyKeyMajInsideCollection<
        Matrix, KeyMajToExclude,
    >
    where
        Matrix:                        IndicesAndCoefficients,    
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                  
{
    matrix:                     Matrix,
    keymajs_to_exclude:         KeyMajToExclude,
}

// implement the object
impl < Matrix, KeyMajToExclude, >

    OnlyKeyMajInsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                        IndicesAndCoefficients,    
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,                          

    {
        /// Create a new `OnlyKeyMajInsideCollection` from a matrix and a collection
        /// of minor keys.
        /// 
        /// In order to operate properly, the collection of minor keys should implement
        /// the [`MapHasKeyOrSequenceHasElement`](crate::utilities::sets::MapHasKeyOrSequenceHasElement) 
        /// trait.
        pub fn new( matrix: Matrix, keymajs_to_exclude: KeyMajToExclude ) 
        -> 
        OnlyKeyMajInsideCollection
            < Matrix, KeyMajToExclude >
        where
            Matrix:                         IndicesAndCoefficients,    
            KeyMajToExclude:                MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,              
        { OnlyKeyMajInsideCollection{ matrix, keymajs_to_exclude, } }
    }    

// IndicesAndCoefficients

impl < Matrix: IndicesAndCoefficients, KeyMajToExclude >

    IndicesAndCoefficients for  

    OnlyKeyMajInsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                         IndicesAndCoefficients,    
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,         

{ 
    type EntryMinor = Matrix::EntryMinor; 
    type EntryMajor = Matrix::EntryMajor; 
    type ColIndex = Matrix::ColIndex; 
    type RowIndex = Matrix::RowIndex; 
    type Coefficient = Matrix::Coefficient; }        


// ViewMinorDescend
impl < Matrix, KeyMajToExclude >

    ViewColDescend for  

    OnlyKeyMajInsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                                 ViewColDescend + IndicesAndCoefficients,    
        Matrix::ViewMinorDescend:                                IntoIterator,
        Matrix::EntryMinor:                           KeyValGet< Matrix::RowIndex, Matrix::Coefficient >,
        KeyMajToExclude:                                       Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,  
        Matrix::ViewMinorDescendIntoIter:                        Iterator,      
{
    type ViewMinorDescendIntoIter    =   < Self::ViewMinorDescend as IntoIterator >::IntoIter;
    type ViewMinorDescend            =   OnlyIndicesInsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Self::RowIndex, Self::Coefficient >;  
    
    fn view_minor_descend( &self, keymin: Self::ColIndex) 
        -> 
        OnlyIndicesInsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Matrix::RowIndex, Matrix::Coefficient >  
    { 
        OnlyIndicesInsideCollection::new( self.matrix.view_minor_descend(keymin).into_iter(), self.keymajs_to_exclude  ) 
    }
}

// ViewMajor
impl < Matrix, KeyMajToExclude >

    ViewCol for  

    OnlyKeyMajInsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                           ViewCol + IndicesAndCoefficients,    
        Matrix::ViewMinor:                                IntoIterator,
        Matrix::EntryMinor:                           KeyValGet< Matrix::RowIndex, Matrix::Coefficient >,
        KeyMajToExclude:                                 Copy + MapHasKeyOrSequenceHasElement< Matrix::RowIndex >,  
        Matrix::ViewMinorIntoIter:                        Iterator,      
{
    type ViewMinorIntoIter    =   < Self::ViewMinor as IntoIterator >::IntoIter;
    type ViewMinor            =   OnlyIndicesInsideCollection< Matrix::ViewMinorIntoIter, KeyMajToExclude, Self::RowIndex, Self::Coefficient >;  

    fn view_minor( &self, keymin: Self::ColIndex) -> Self::ViewMinor
    { 
        OnlyIndicesInsideCollection::new( self.matrix.view_minor(keymin).into_iter(), self.keymajs_to_exclude  ) 
    }
}















//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {


    // Entry-wise clone
    // =====================================================================================================================
    
    #[test] 
    fn test_entrywise_clone() {

        // import crates
        use crate::algebra::matrices::operations::transform_vector_wise::VecWiseTransformed;        
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewRowAscend;
        use std::iter::Cloned;
        use std::slice::Iter;
        use itertools::Itertools;

        // define a function that doubles the second entry in a tuple
        fn double_coefficient( pair: (usize, usize) ) -> (usize, usize) { ( pair.0, 2 * pair.1) } 
        fn vector_transformer( iter: Cloned< Iter< '_, (usize,usize) > > ) -> Vec< (usize, usize) > { 
            iter
            .map(double_coefficient )
            .collect_vec() 
        }

        // define the untransformed matrix
        let matrix_untransformed    =   VecOfVec::new( vec![ vec![ (0,1) ] ] );  
        let matrix_untransformed_ref    =   & matrix_untransformed;   // a reference is required to implement the oracle trait   
            
        // define the transformed matrix
        let matrix_transformed      =   VecWiseTransformed::new( 
                                                matrix_untransformed_ref,
                                                vector_transformer,
                                            );
        let transformed_vector      =   ( matrix_transformed ).view_major_ascend( 0 ).into_iter().collect_vec();
        assert_eq!( transformed_vector, vec![ (0,2) ] )
    }

} 

