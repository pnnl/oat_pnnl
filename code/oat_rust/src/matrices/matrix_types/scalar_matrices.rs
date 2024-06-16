//! Memory-efficient scalar matrices (ie diagonal matrices where all entries on the diagonal are equal)

use crate::matrices::matrix_oracle_traits::{   OracleMajor,
                                        OracleMajorAscend,
                                        OracleMajorDescend,
                                        OracleMinor, 
                                        OracleMinorAscend,
                                        OracleMinorDescend, IndicesAndCoefficients,
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
/// use oat_rust::matrices::matrix_types::scalar_matrices::ScalarMatrix;
/// use oat_rust::matrices::matrix_oracle_traits::{OracleMajor};   
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
                scalar:             scalar,
                phantom_key:        PhantomData,
            }
    }
}


//  ---------------------
//  TRAIT IMPLEMENTATIONS
//  ---------------------

//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------

impl    < Key, Val >

    IndicesAndCoefficients      for 
    ScalarMatrix                < Key, Val >
{ type KeyMin = Key; type KeyMaj = Key; type SnzVal = Val; }


//  MAJOR VIEWS
//  ---------------------------------------------------------------------------


//  OracleMajor
//  
impl     < Key, Val >

        OracleMajor       for        
        ScalarMatrix      < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMajor = Once < (Key, Val) >;
    type ViewMajorIntoIter = < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajorEntry = (Key, Val);

    fn view_major( & self, index: Key ) -> Once < (Key, Val) > 
    { 
        once( ( index, self.scalar.clone() ) )
    }
}

//  OracleMajorAscend
//  
impl     < Key, Val >

        OracleMajorAscend   for
        ScalarMatrix        < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMajorAscend = Once < (Key, Val) >;
    type ViewMajorAscendIntoIter = < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscendEntry = (Key, Val);

    fn view_major_ascend( & self, index: Key ) -> Once < (Key, Val) >
    { 
        once( ( index, self.scalar.clone() ) )
    }
}


//  OracleMajorDescend
//  
impl     < Key, Val >

        OracleMajorDescend  for
        ScalarMatrix  < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMajorDescend           =   Once < (Key, Val) >;
    type ViewMajorDescendIntoIter   =   < Self::ViewMajorDescend as IntoIterator >::IntoIter;
    type ViewMajorDescendEntry      =   (Key, Val);    

    fn view_major_descend( & self, index: Key ) -> Once < (Key, Val) >
    { 
        once( ( index, self.scalar.clone() ) )
    }
}


//  MINOR VIEWS
//  ---------------------------------------------------------------------------


//  OracleMinor
//  
impl     < 'a, Key, Val >

        OracleMinor       for
        ScalarMatrix      < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMinor          =   Once < (Key, Val) >;
    type ViewMinorIntoIter  =   Once < (Key, Val) >;
    type ViewMinorEntry     =   (Key, Val);    

    fn view_minor( & self, index: Key ) -> Once < (Key, Val) > 
    { 
        once( ( index, self.scalar.clone() ) )
    }
}

//  OracleMinorAscend
//  
impl     < 'a, Key, Val >

        OracleMinorAscend   for
        ScalarMatrix        < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMinorAscend           =   Once < (Key, Val) >;
    type ViewMinorAscendIntoIter   =   < Self::ViewMinorAscend as IntoIterator >::IntoIter;
    type ViewMinorAscendEntry      =   (Key, Val);    

    fn view_minor_ascend( & self, index: Key ) -> Once < (Key, Val) >
    { 
        once( ( index, self.scalar.clone() ) )
    }
}


//  OracleMinorDescend
//  
impl     < 'a, Key, Val >

        OracleMinorDescend      for
        ScalarMatrix            < Key, Val > 
        
        where   Val: Clone, // hard to drop this requirement (tuples give move errors if no clone) 
                Key: Clone  // hard to drop this requirement (tuples give move errors if no clone) 
{
    type ViewMinorDescend           =   Once < (Key, Val) >;
    type ViewMinorDescendIntoIter   =   Once < (Key, Val) >;
    type ViewMinorDescendEntry      =   (Key, Val);        

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
        use crate::matrices::matrix_types::scalar_matrices::ScalarMatrix;
        use crate::matrices::matrix_oracle_traits::{OracleMajor};   

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
// /// use oat_rust::matrices::matrix_types::scalar_matrices::ScalarMatrixOracleUsize;
// /// use oat_rust::matrices::matrix_oracle_traits::{OracleMajor};   
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
// //  OracleMajor
// //  
// impl     < 'a, Val >
// 
//         OracleMajor < 'a, usize, Once< (usize, Val) > > for 
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
// //  OracleMajorAscend
// //  
// impl    < 'a, Val >
//         
//         OracleMajorAscend < 'a, usize, Once< (usize, Val) > > for 
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
// //  OracleMajorDescend
// //  
// impl     < 'a, Val >
// 
//         OracleMajorDescend          < 'a, usize, Once < (usize, Val) > > for         
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
// //  OracleMinor
// //  
// impl     < 'a, Val >
//         
//         OracleMinor                 < 'a, usize, Once < (usize, Val) > > for
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
// //  OracleMinorAscend
// //  
// impl     < 'a, Val >
//         
//         OracleMinorAscend           < 'a, usize, Once < (usize, Val) > > for
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
// //  OracleMinorDescend
// //  
// impl     < 'a, Val >
//         
//         OracleMinorDescend          < 'a, usize, Once < (usize, Val) > > for
//         ScalarMatrixOracleUsize     < Key, Val > 
//         
//         where   Val:    Clone, // hard to drop this requirement (tuples give move errors if no clone) 
// {
//     fn view_minor_descend<'b: 'a>( &'b self, index: usize ) -> Once < (usize, Val) >
//     { 
//         once( ( index, self.scalar.clone() ) )
//     }
// }











