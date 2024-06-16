
//! Create a new matrix by applying a transformation to each row or column of an existing matrix.
//! 
//! See also [`transform_entry_wise`](oat_rust::matrices::operations::transform_entry_wise)



//  GENERAL VECTORWISE TRANSFORM
//  ---------------------------------------------------------------------

/// Wrapper for matrix oracles; applies `vector_transformer` to each sparse vector returned by the matrix.
/// 
/// # Examples
/// 
/// ```
/// // import crates
/// use oat_rust::matrices::operations::transform_vector_wise::VecWiseTransformed;        
/// use oat_rust::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
/// use oat_rust::matrices::matrix_oracle_traits::OracleMajorAscend;
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
/// let matrix_untransformed    =   VecOfVecSimple::new( vec![ vec![ (0,1) ] ] );  
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
impl < MatrixUntransformed: IndicesAndCoefficients, VectorTransformer >

    IndicesAndCoefficients for  

    VecWiseTransformed
        < MatrixUntransformed, VectorTransformer > 

{ type KeyMin = MatrixUntransformed::KeyMin; type KeyMaj = MatrixUntransformed::KeyMaj; type SnzVal = MatrixUntransformed::SnzVal; }  


//  OracleMajorAscend
impl    < MatrixUntransformed, ViewMajorAscendTransformed, VectorTransformer, > 

        OracleMajorAscend for 

        VecWiseTransformed
            < MatrixUntransformed, VectorTransformer > where

        MatrixUntransformed:            OracleMajorAscend + IndicesAndCoefficients,
        ViewMajorAscendTransformed:     IntoIterator,
        MatrixUntransformed:            OracleMajorAscend,
        VectorTransformer:              Fn( MatrixUntransformed::ViewMajorAscend ) -> ViewMajorAscendTransformed,                                    
{
    type ViewMajorAscend            =   ViewMajorAscendTransformed;
    type ViewMajorAscendEntry       =   ViewMajorAscendTransformed::Item;
    type ViewMajorAscendIntoIter    =   ViewMajorAscendTransformed::IntoIter;

    fn view_major_ascend ( & self, majkey: Self::KeyMaj ) -> Self::ViewMajorAscend { 
        return (self.vector_transformer)(  self.untransformed_matrix.view_major_ascend( majkey )  )
    }
}




//  SIMPLIFY VECTORS
//  ---------------------------------------------------------------------


use std::marker::PhantomData;

use crate::{matrices::matrix_oracle_traits::{OracleMajorAscend, OracleMajor, IndicesAndCoefficients, OracleMinorDescend, OracleMinor}, utilities::{iterators::general::PeekUnqualified, sets::MapHasKeyOrSequenceHasElement}, rings::operator_traits::Semiring, vectors::operations::{Simplify, Transforms}, entries::{KeyValGet, KeyValSet}};

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
        < MatrixUnsimplified, RingOperator >  

{ type KeyMin = MatrixUnsimplified::KeyMin; type KeyMaj = MatrixUnsimplified::KeyMaj; type SnzVal = MatrixUnsimplified::SnzVal; }  


// OracleMajorAscend
impl    < MatrixUnsimplified, RingOperator > 

        OracleMajorAscend for 

        MatrixSimplify
            < MatrixUnsimplified, RingOperator > 
        
        where
        MatrixUnsimplified:                             OracleMajorAscend + IndicesAndCoefficients,
        MatrixUnsimplified::ViewMajorAscend:            IntoIterator,
        // The next uncommented requirement is the same as the following commented one (Rust just has a hard time seeing that they are equivalent)
        // MatrixUnsimplified::ViewMajorAscendEntry:       KeyValGet < MatrixUnsimplified::KeyMin, MatrixUnsimplified::SnzVal, > + 
        //                                                 KeyValSet < MatrixUnsimplified::KeyMin, MatrixUnsimplified::SnzVal, >,        
        < MatrixUnsimplified::ViewMajorAscendIntoIter as IntoIterator >::Item:
                                                        KeyValGet < MatrixUnsimplified::KeyMin, MatrixUnsimplified::SnzVal, > + 
                                                        KeyValSet < MatrixUnsimplified::KeyMin, MatrixUnsimplified::SnzVal, >,        
        MatrixUnsimplified::ViewMajorAscendIntoIter:    PeekUnqualified,  
        MatrixUnsimplified::KeyMin:                     std::cmp::PartialEq,
        RingOperator:                                   Clone + Semiring< MatrixUnsimplified::SnzVal >,
{
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;    
    type ViewMajorAscendEntry       =   < Self::ViewMajorAscend as IntoIterator >::Item;
    type ViewMajorAscend            =   Simplify<
                                                MatrixUnsimplified::ViewMajorAscendIntoIter,
                                                MatrixUnsimplified::KeyMin,
                                                RingOperator,
                                                MatrixUnsimplified::SnzVal,
                                            > ;    

    fn view_major_ascend ( & self, majkey: Self::KeyMaj) -> Self::ViewMajorAscend { 
        let x = 
            self.unsimplified_matrix
                .view_major_ascend( majkey )
                .into_iter()
                .simplify( self.ring_operator.clone() );
        x
    }
}





//  MASK: ONLY INDICES OUTSIDE
//  -------------------------------------------------------------------------------------

//  MINOR KEYS
//  -----------------------------------------

use crate::vectors::operations::OnlyIndicesOutsideCollection;


/// A matrix oracle whose major views iterate only over entries that have indices **outside**
/// a given collection.
pub struct OnlyKeyMinOutsideCollection<
        Matrix, KeyMinToExclude,
    >
    where
        Matrix:                         IndicesAndCoefficients,    
        KeyMinToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMin >,                  
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
        KeyMinToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMin >,                          

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
            Matrix:                         OracleMajorAscend + IndicesAndCoefficients,    
            Matrix::ViewMajorAscend:        IntoIterator,
            KeyMinToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMin >,              
        { OnlyKeyMinOutsideCollection{ matrix: matrix, keymins_to_exclude: keymins_to_exclude, } }
    }    

// IndicesAndCoefficients

impl < Matrix: IndicesAndCoefficients, KeyMinToExclude >

    IndicesAndCoefficients for  

    OnlyKeyMinOutsideCollection
        < Matrix, KeyMinToExclude >

    where
        Matrix:                         IndicesAndCoefficients,    
        KeyMinToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMin >,         

{ type KeyMin = Matrix::KeyMin; type KeyMaj = Matrix::KeyMaj; type SnzVal = Matrix::SnzVal; }        


// ViewMajorAscend
impl < Matrix, KeyMinToExclude >

    OracleMajorAscend for  

    OnlyKeyMinOutsideCollection
        < Matrix, KeyMinToExclude >

    where
        Matrix:                                                 OracleMajorAscend + IndicesAndCoefficients,    
        Matrix::ViewMajorAscend:                                IntoIterator,
        Matrix::ViewMajorAscendEntry:                           KeyValGet< Matrix::KeyMin, Matrix::SnzVal >,
        KeyMinToExclude:                                       Copy + MapHasKeyOrSequenceHasElement< Matrix::KeyMin >,  
        Matrix::ViewMajorAscendIntoIter:                        Iterator,      
{
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscendEntry       =   < Self::ViewMajorAscend as IntoIterator >::Item;
    type ViewMajorAscend            =   OnlyIndicesOutsideCollection< Matrix::ViewMajorAscendIntoIter, KeyMinToExclude, Self::KeyMin, Self::SnzVal >;  
    
    fn view_major_ascend( &self, keymaj: Self::KeyMaj) 
        -> 
        OnlyIndicesOutsideCollection< Matrix::ViewMajorAscendIntoIter, KeyMinToExclude, Matrix::KeyMin, Matrix::SnzVal >  
    { 
        OnlyIndicesOutsideCollection::new( self.matrix.view_major_ascend(keymaj).into_iter(), self.keymins_to_exclude  ) 
    }
}

// ViewMajor
impl < Matrix, KeyMinToExclude >

    OracleMajor for  

    OnlyKeyMinOutsideCollection
        < Matrix, KeyMinToExclude >

    where
        Matrix:                                           OracleMajor + IndicesAndCoefficients,    
        Matrix::ViewMajor:                                IntoIterator,
        Matrix::ViewMajorEntry:                           KeyValGet< Matrix::KeyMin, Matrix::SnzVal >,
        KeyMinToExclude:                                 Copy + MapHasKeyOrSequenceHasElement< Matrix::KeyMin >,  
        Matrix::ViewMajorIntoIter:                        Iterator,      
{
    type ViewMajorIntoIter    =   < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajorEntry       =   < Self::ViewMajor as IntoIterator >::Item;
    type ViewMajor            =   OnlyIndicesOutsideCollection< Matrix::ViewMajorIntoIter, KeyMinToExclude, Self::KeyMin, Self::SnzVal >;  

    fn view_major( &self, keymaj: Self::KeyMaj) -> Self::ViewMajor
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
        KeyMajToExclude:                MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,                  
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
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,                          

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
            KeyMajToExclude:                MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,              
        { OnlyKeyMajOutsideCollection{ matrix: matrix, keymajs_to_exclude: keymajs_to_exclude, } }
    }    

// IndicesAndCoefficients

impl < Matrix: IndicesAndCoefficients, KeyMajToExclude >

    IndicesAndCoefficients for  

    OnlyKeyMajOutsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                        IndicesAndCoefficients,    
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,         

{ type KeyMin = Matrix::KeyMin; type KeyMaj = Matrix::KeyMaj; type SnzVal = Matrix::SnzVal; }        


// ViewMinorDescend
impl < Matrix, KeyMajToExclude >

    OracleMinorDescend for  

    OnlyKeyMajOutsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                                 OracleMinorDescend + IndicesAndCoefficients,    
        Matrix::ViewMinorDescend:                                IntoIterator,
        Matrix::ViewMinorDescendEntry:                           KeyValGet< Matrix::KeyMaj, Matrix::SnzVal >,
        KeyMajToExclude:                                       Copy + MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,  
        Matrix::ViewMinorDescendIntoIter:                        Iterator,      
{
    type ViewMinorDescendIntoIter    =   < Self::ViewMinorDescend as IntoIterator >::IntoIter;
    type ViewMinorDescendEntry       =   < Self::ViewMinorDescend as IntoIterator >::Item;
    type ViewMinorDescend            =   OnlyIndicesOutsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Self::KeyMaj, Self::SnzVal >;  
    
    fn view_minor_descend( &self, keymin: Self::KeyMin ) 
        -> 
        OnlyIndicesOutsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Matrix::KeyMaj, Matrix::SnzVal >  
    { 
        OnlyIndicesOutsideCollection::new( self.matrix.view_minor_descend(keymin).into_iter(), self.keymajs_to_exclude  ) 
    }
}

// ViewMajor
impl < Matrix, KeyMajToExclude >

    OracleMinor for  

    OnlyKeyMajOutsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                           OracleMinor + IndicesAndCoefficients,    
        Matrix::ViewMinor:                                IntoIterator,
        Matrix::ViewMinorEntry:                           KeyValGet< Matrix::KeyMaj, Matrix::SnzVal >,
        KeyMajToExclude:                                 Copy + MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,  
        Matrix::ViewMinorIntoIter:                        Iterator,      
{
    type ViewMinorIntoIter    =   < Self::ViewMinor as IntoIterator >::IntoIter;
    type ViewMinorEntry       =   < Self::ViewMinor as IntoIterator >::Item;
    type ViewMinor            =   OnlyIndicesOutsideCollection< Matrix::ViewMinorIntoIter, KeyMajToExclude, Self::KeyMaj, Self::SnzVal >;  

    fn view_minor( &self, keymin: Self::KeyMin) -> Self::ViewMinor
    { 
        OnlyIndicesOutsideCollection::new( self.matrix.view_minor(keymin).into_iter(), self.keymajs_to_exclude  ) 
    }
}



//  MASK: ONLY INDICES INSIDE
//  -------------------------------------------------------------------------------------



//  MINOR KEYS
//  -----------------------------------------


use crate::vectors::operations::OnlyIndicesInsideCollection;


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
            KeyMinToInclude:       MapHasKeyOrSequenceHasElement< Matrix::KeyMin >,              
        { OnlyKeyMinInsideCollection{ matrix, keymins_to_include } }
    }    

// IndicesAndCoefficients
impl < Matrix: IndicesAndCoefficients, KeyMinToInclude >

    IndicesAndCoefficients for  

    OnlyKeyMinInsideCollection
        < Matrix, KeyMinToInclude >

{ type KeyMin = Matrix::KeyMin; type KeyMaj = Matrix::KeyMaj; type SnzVal = Matrix::SnzVal; } 


// ViewMajorAscend
impl < Matrix, KeyMinToInclude >

    OracleMajorAscend for  

    OnlyKeyMinInsideCollection
        < Matrix, KeyMinToInclude >

    where
        Matrix:                             OracleMajorAscend + IndicesAndCoefficients,    
        Matrix::ViewMajorAscend:            IntoIterator,
        Matrix::ViewMajorAscendIntoIter:    Iterator,
        Matrix::ViewMajorAscendEntry:   KeyValGet< Matrix::KeyMin, Matrix::SnzVal >,        
        KeyMinToInclude:                   MapHasKeyOrSequenceHasElement< Self::KeyMin > + Copy,   
{
    type ViewMajorAscendIntoIter    =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscendEntry       =   < Self::ViewMajorAscend as IntoIterator >::Item;
    type ViewMajorAscend            =   OnlyIndicesInsideCollection< Matrix::ViewMajorAscendIntoIter, KeyMinToInclude, Matrix::KeyMin, Matrix::SnzVal >;  

    fn view_major_ascend( &self, keymaj: Self::KeyMaj ) 
        -> 
        Self::ViewMajorAscend
    { 
        OnlyIndicesInsideCollection::new( self.matrix.view_major_ascend(keymaj).into_iter(), self.keymins_to_include  ) 
    }
}


// ViewMajor
impl < Matrix, KeyMinToInclude, >

    OracleMajor for  

    OnlyKeyMinInsideCollection
        < Matrix, KeyMinToInclude, >

    where
        Matrix:                     OracleMajor + IndicesAndCoefficients,    
        Matrix::ViewMajor:          IntoIterator,
        Matrix::ViewMajorIntoIter:  Iterator,
        Matrix::ViewMajorEntry:     KeyValGet< Matrix::KeyMin, Matrix::SnzVal >,        
        KeyMinToInclude:           MapHasKeyOrSequenceHasElement< Self::KeyMin > + Copy,     
{
    type ViewMajorIntoIter    =   < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajorEntry       =   < Self::ViewMajor as IntoIterator >::Item;
    type ViewMajor            =   OnlyIndicesInsideCollection< Matrix::ViewMajorIntoIter, KeyMinToInclude, Matrix::KeyMin, Matrix::SnzVal >;  

    fn view_major( &self, keymaj: Self::KeyMaj ) 
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
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,                  
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
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,                          

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
            KeyMajToExclude:                MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,              
        { OnlyKeyMajInsideCollection{ matrix, keymajs_to_exclude, } }
    }    

// IndicesAndCoefficients

impl < Matrix: IndicesAndCoefficients, KeyMajToExclude >

    IndicesAndCoefficients for  

    OnlyKeyMajInsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                         IndicesAndCoefficients,    
        KeyMajToExclude:               MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,         

{ type KeyMin = Matrix::KeyMin; type KeyMaj = Matrix::KeyMaj; type SnzVal = Matrix::SnzVal; }        


// ViewMinorDescend
impl < Matrix, KeyMajToExclude >

    OracleMinorDescend for  

    OnlyKeyMajInsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                                 OracleMinorDescend + IndicesAndCoefficients,    
        Matrix::ViewMinorDescend:                                IntoIterator,
        Matrix::ViewMinorDescendEntry:                           KeyValGet< Matrix::KeyMaj, Matrix::SnzVal >,
        KeyMajToExclude:                                       Copy + MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,  
        Matrix::ViewMinorDescendIntoIter:                        Iterator,      
{
    type ViewMinorDescendIntoIter    =   < Self::ViewMinorDescend as IntoIterator >::IntoIter;
    type ViewMinorDescendEntry       =   < Self::ViewMinorDescend as IntoIterator >::Item;
    type ViewMinorDescend            =   OnlyIndicesInsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Self::KeyMaj, Self::SnzVal >;  
    
    fn view_minor_descend( &self, keymin: Self::KeyMin) 
        -> 
        OnlyIndicesInsideCollection< Matrix::ViewMinorDescendIntoIter, KeyMajToExclude, Matrix::KeyMaj, Matrix::SnzVal >  
    { 
        OnlyIndicesInsideCollection::new( self.matrix.view_minor_descend(keymin).into_iter(), self.keymajs_to_exclude  ) 
    }
}

// ViewMajor
impl < Matrix, KeyMajToExclude >

    OracleMinor for  

    OnlyKeyMajInsideCollection
        < Matrix, KeyMajToExclude >

    where
        Matrix:                                           OracleMinor + IndicesAndCoefficients,    
        Matrix::ViewMinor:                                IntoIterator,
        Matrix::ViewMinorEntry:                           KeyValGet< Matrix::KeyMaj, Matrix::SnzVal >,
        KeyMajToExclude:                                 Copy + MapHasKeyOrSequenceHasElement< Matrix::KeyMaj >,  
        Matrix::ViewMinorIntoIter:                        Iterator,      
{
    type ViewMinorIntoIter    =   < Self::ViewMinor as IntoIterator >::IntoIter;
    type ViewMinorEntry       =   < Self::ViewMinor as IntoIterator >::Item;
    type ViewMinor            =   OnlyIndicesInsideCollection< Matrix::ViewMinorIntoIter, KeyMajToExclude, Self::KeyMaj, Self::SnzVal >;  

    fn view_minor( &self, keymin: Self::KeyMin) -> Self::ViewMinor
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
        use crate::matrices::operations::transform_vector_wise::VecWiseTransformed;        
        use crate::matrices::matrix_types::vec_of_vec::VecOfVecSimple;
        use crate::matrices::matrix_oracle_traits::OracleMajorAscend;
        use std::iter::Cloned;
        use std::slice::Iter;
        use itertools::Itertools;

        // define a function that doubles the second entry in a tuple
        fn double_coefficient( pair: (usize, usize) ) -> (usize, usize) { ( pair.0, 2 * pair.1) } 
        fn vector_transformer( iter: Cloned< Iter< '_, (usize,usize) > > ) -> Vec< (usize, usize) > { 
            iter
            .map(|x| double_coefficient(x) )
            .collect_vec() 
        }

        // define the untransformed matrix
        let matrix_untransformed    =   VecOfVecSimple::new( vec![ vec![ (0,1) ] ] );  
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

