
// OPTIONS

//  - macro to impl on refs (unexplained errors)
//  - use "impl Oracle for &T where T impl Oracle" (hard to explain errors)
//  - use "impl Oracle for &T where T impl Oracle + Inherit" + OracleRef
//  - 




//! Look up the (structural) nonzero entries in a row or column.
//! 
//! A matrix oracle is an object that allows one to ``read`` the rows and/or columns of a
//! matrix.
//! 
//! 
//! 
//! 
//! # Views and major dimensions
//! 
//! To understand matrix oracles in oat_rust, you need just two basic concepts: major dimension
//! and view.  These are easiest to understand in terms of an example:
//! 
//! **Example: Vector-of-Vectors**
//! 
//! Consider the following 2x2 matrix:
//! 
//! ```ignore
//! 5 6
//! 7 0
//! ```
//! 
//! We can represent this matrix as a vector-of-vectors, where each internal 
//! vector represents an ordered list of the nonzero entries in a given row
//! (we call this type of representation *row-major*):
//! 
//! ```ignore
//! vec![
//!     vec![ (0, 5), (1, 6) ],
//!     vec![ (0, 7) ],
//! ]
//! ```
//! 
//! Alternatively, we can represent the matrix as a vector-of-vectors where each internal 
//! vector represents an ordered list of the nonzero entries in a given column (we call
//! this type of repreentation *column-major*):
//! 
//! ```ignore
//! vec![
//!     vec![ (0, 5), (1, 7) ],
//!     vec![ (0, 6) ],
//! ]
//! ```
//! 
//! If you want to reconstruct a full row of the matrix from one of these vector-of-vector representations,
//! then you will probably use a different procedure than you would use to
//! reconstruct a full column of the matrix.
//! For example, given a *row-major* representation of a matrix,  
//! - you can reconstruct row `i` by iterating over the entries in the `i`th internal vector
//! - you can reconstruct column `i`, by searching through each and every internal vector,
//! to check if that vector contains an entry of form `(i, *)`.
//! 
//! Many sparse matrix data structures share this characteristic: they use different procedures to reconstruct
//! rows versus columns.  Typically, one of the two procedures will be significantly faster / more 
//! efficient than the other; in the case of vector-of-vectors, for example, it's much faster to iterate
//! over a single internal vector than it is to search across all of them, so reconstructing rows is
//! easier in row-major vector-of-vector representations, while reconstructing columns is easier in 
//! column-major vector-of-vector representations.
//! 
//! For the purpose of writing code, it's often irrelevant whether the internal vectors of a vector-of-vectors
//! represent rows versus columns.  All that matters is that 
//! - you have a data structure that can produce sparse vectors in one of two ways, and 
//! - one way is much faster / more efficient than the other.
//! 
//! A sparse vector constructed via the efficient method is called a **major view**.  A sparse vector 
//! constructed via the less efficient method is called a **minor view**.  A data structure that
//! represents a matrix is **row-major** if major views represent rows, and **column-major** if 
//! major-views represent columns.
//! 
//! That's it!  Now you understand about major views and row-major versus column-major reprsentations!
//! 
//! 
//! # Examples
//! 
//! 
//! The methods used to access major/minor views of a matrix (with entries listed in ascending or descending order)
//! are organized into traits.  A list of these traits can be found below.
//! See the [matrices](crate::matrices) module to find examples of how these methods are used, in practice.




use std::fmt::Debug;
use std::iter::IntoIterator;
use auto_impl::auto_impl; // auto-implement a trait on references to objects that implement the trait



// /// An enum with two values: `Row` and `Col`.
// #[derive(Clone, Debug)]
// pub enum MajorDimension{
//     Row,
//     Col
// }


//  ---------------------------------------------------------------------------
//  MAJOR DIMENSION 
//  ---------------------------------------------------------------------------

// /// Specifies a major dimension (row or column).
// pub trait WhichMajor{ fn major_dimension( &self ) -> MajorDimension; }


//  ---------------------------------------------------------------------------
//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------

/// Specifies the major indices, minor indices, and coefficients of a sparse matrix.
pub trait IndicesAndCoefficients {
    type KeyMin;
    type KeyMaj;
    type SnzVal;    
}

impl < 'a, T: IndicesAndCoefficients >
    
    IndicesAndCoefficients for 
    &'a T
{   type KeyMin = T::KeyMin; type KeyMaj = T::KeyMaj; type SnzVal = T::SnzVal;    }    


// pub trait MatrixTypeParameters {
//     type KeyMin;
//     type KeyMaj;
//     type SnzVal;    
//     type KeyValPairMin;
//     type KeyValPairMaj;
// }

// impl < 'a, T: MatrixTypeParameters >
    
//     MatrixTypeParameters for 
//     &'a T
// {   type KeyMin = T::KeyMin; type KeyMaj = T::KeyMaj; type SnzVal = T::SnzVal; type KeyValPairMaj = T::KeyValPairMaj; type KeyValPairMin = T::KeyValPairMin;   }    




//  ---------------------------------------------------------------------------
//  ORACLE REFERENCE INHERIT
//  ---------------------------------------------------------------------------

// pub trait OracleRefInherit {}

// /// If `T` implements `OracleMinorDesecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: IndicesAndCoefficients > IndicesAndCoefficients for &'a T {
//     type KeyMin     =   T::KeyMin;
//     type KeyMaj     =   T::KeyMaj;        
//     type SnzVal     =   T::SnzVal;                
// }

// /// If `T` implements `OracleInherit`, then so does `&'a T`.
// impl < 'a, T: OracleRefInherit > OracleRefInherit for & 'a T {}

// /// If `T` implements `OracleMajor` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + OracleMajor< KeyMaj, ViewMajor >, KeyMaj, ViewMajor:IntoIterator > OracleMajor< KeyMaj, ViewMajor > for &'a T {
//     fn view_major(&self, index: Self::KeyMaj) -> Self::ViewMajor {
//         (*self).view_major(index)
//     }
// }

// /// If `T` implements `OracleMajorAsecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + OracleMajorAscend< KeyMaj, ViewMajorAscend >, KeyMaj, ViewMajorAscend:IntoIterator > OracleMajorAscend< KeyMaj, ViewMajorAscend > for &'a T {
//     fn view_major_ascend(&self, index: Self::KeyMaj) -> Self::ViewMajorAscend {
//         (*self).view_major_ascend(index)
//     }
// }

// /// If `T` implements `OracleMajorDesecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + OracleMajorDescend< KeyMaj, ViewMajorDescend >, KeyMaj, ViewMajorDescend:IntoIterator > OracleMajorDescend< KeyMaj, ViewMajorDescend > for &'a T {
//     fn view_major_descend(&self, index: Self::KeyMaj) -> Self::ViewMajorDescend {
//         (*self).view_major_descend(index)
//     }
// }

// /// If `T` implements `OracleMinor` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + OracleMinor< KeyMaj, ViewMinor >, KeyMaj, ViewMinor > OracleMinor< KeyMaj, ViewMinor > for &'a T {
//     fn view_minor(&self, index: Self::KeyMaj) -> Self::ViewMinor {
//         (*self).view_minor(index)
//     }
// }

// /// If `T` implements `OracleMinorAsecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + OracleMinorAscend< KeyMaj, ViewMinorAscend >, KeyMaj, ViewMinorAscend:IntoIterator > OracleMinorAscend< KeyMaj, ViewMinorAscend > for &'a T {
//     fn view_minor_ascend(&self, index: Self::KeyMaj) -> Self::ViewMinorAscend {
//         (*self).view_minor_ascend(index)
//     }
// }

// impl < 'a, T: IndicesAndCoefficients + OracleMinorDescend > OracleMinorDescend for &'a T {
//     type ViewMinorDescend       =   T::ViewMinorDescend;
//     type ViewMinorDescendIntoIter   =   T::ViewMinorDescendIntoIter;        
//     type ViewMinorDescendEntry  =   T::ViewMinorDescendEntry;                
//     fn view_minor_descend(&self, index: T::KeyMin) -> T::ViewMinorDescend {
//         (*self).view_minor_descend(index)
//     }
// }



//  ---------------------------------------------------------------------------
//  ORACLE MAJOR
//  ---------------------------------------------------------------------------

//  OracleMajorInherit IS DEPRECATED IN FAVOR OF OracleRefInherit
//  
// /// A trait with no methods, used as a criterion for auto-implementation of the OracleMajor trait on references.
// pub trait OracleMajorInherit {}

// /// If `T` implements `OracleMajorInherit`, the so does `&'a T`.
// impl < 'a, T: OracleMajorInherit > 

//     OracleMajorInherit for 
    
//     & 'a T 
//     {}

/// Entries may not appear in sorted order.
// #[auto_impl(&)] 
pub trait OracleMajor: IndicesAndCoefficients
{
    type ViewMajor: IntoIterator< IntoIter = Self::ViewMajorIntoIter, Item = Self::ViewMajorEntry >;
    type ViewMajorIntoIter: Iterator< Item = Self::ViewMajorEntry >;
    type ViewMajorEntry;

    /// Get a major view.
    ///
    /// The order in which terms appear should be the same every time the
    /// function is called; however, the order need not be sorted.
    fn   view_major( &self, index: Self::KeyMaj ) -> Self::ViewMajor;
}



/// Entries appear in **strictly** ascending order, according to index.  Consecutive entries must have distinct indices.
// #[auto_impl(&)] 
pub trait OracleMajorAscend:    IndicesAndCoefficients
{
    type ViewMajorAscend: IntoIterator< IntoIter = Self::ViewMajorAscendIntoIter, Item = Self::ViewMajorAscendEntry >;
    type ViewMajorAscendIntoIter: Iterator< Item = Self::ViewMajorAscendEntry >;
    type ViewMajorAscendEntry;

    /// Get a major view with entries sorted in ascending order of index.
    fn   view_major_ascend( &self, index: Self::KeyMaj ) -> Self::ViewMajorAscend;
}

impl < 'a, T > 

    OracleMajorAscend for 
    
    &'a T

    where
        T:                  OracleMajorAscend + IndicesAndCoefficients,
        &'a T:              IndicesAndCoefficients< KeyMin=T::KeyMin, KeyMaj=T::KeyMaj, SnzVal=T::SnzVal >,        

    {
        type ViewMajorAscend            =   T::ViewMajorAscend;
        type ViewMajorAscendIntoIter    =   T::ViewMajorAscendIntoIter;
        type ViewMajorAscendEntry       =   T::ViewMajorAscendEntry;                

        fn   view_major_ascend( &self, index: T::KeyMaj ) -> T::ViewMajorAscend { (*self).view_major_ascend( index ) }
    }



/// Entries appear in **strictly** descending order, according to index.  Consecutive entries must have distinct indices.
pub trait OracleMajorDescend: IndicesAndCoefficients
{
    type ViewMajorDescend: IntoIterator< IntoIter = Self::ViewMajorDescendIntoIter, Item = Self::ViewMajorDescendEntry >;
    type ViewMajorDescendIntoIter: Iterator< Item = Self::ViewMajorDescendEntry >;
    type ViewMajorDescendEntry;
        
    /// Get a major view with entries sorted in descending order of index.
    fn   view_major_descend( &self, index: Self::KeyMaj ) -> Self::ViewMajorDescend;
}

/// Entries appear in **non-strictly** ascending order, according to index.  Consecutive entries may have identical indices.
// #[auto_impl(&)] 
pub trait OracleMajorAscendLax: IndicesAndCoefficients
{
    type ViewMajorAscendLax: IntoIterator< IntoIter = Self::ViewMajorAscendLaxIntoIter, Item = Self::ViewMajorAscendLaxEntry >;
    type ViewMajorAscendLaxIntoIter: Iterator< Item = Self::ViewMajorAscendLaxEntry >;
    type ViewMajorAscendLaxEntry;

    /// Get a major view with entries sorted in ascending order of index.
    fn   view_major_ascend( &self, index: Self::KeyMaj ) -> Self::ViewMajorAscendLax;
}

/// Entries appear in **non-strictly** descending order, according to index.   Consecutive entries may have identical indices.
// #[auto_impl(&)] 
pub trait OracleMajorDescendLax: IndicesAndCoefficients
{
    type ViewMajorDescendLax: IntoIterator< IntoIter = Self::ViewMajorDescendLaxIntoIter, Item = Self::ViewMajorDescendLaxEntry >;
    type ViewMajorDescendLaxIntoIter: Iterator< Item = Self::ViewMajorDescendLaxEntry >;
    type ViewMajorDescendLaxEntry;

    /// Get a major view with entries sorted in descending order of index.
    fn   view_major_descend( &self, index: Self::KeyMaj ) -> Self::ViewMajorDescendLax;
}

// FOR FUTURE CONSIDERATION
// pub trait OracleMajorAscendScoped < KeyMin, KeyMaj, SnzVal>
// {
//     type PairMajorAscendScoped: KeyValGet < KeyMin, SnzVal, >;
//     type ViewMajorAscendScoped: IntoIterator < Item = PairMajorAscendScoped >;
//     /// Get a major view with entries sorted in ascending order of index, clipped to range [min,
//     /// max).
//     fn   view_major_ascend_scoped( &self, index: Self::KeyMaj, min: KeyMin, max: KeyMin ) -> Self::ViewMajorAscendScoped;
// }

//  ---------------------------------------------------------------------------
//  ORACLE MINOR
//  ---------------------------------------------------------------------------

/// Entries may not appear in sorted order.
// #[auto_impl(&)] 
pub trait OracleMinor: IndicesAndCoefficients
{
    type ViewMinor: IntoIterator< IntoIter = Self::ViewMinorIntoIter, Item = Self::ViewMinorEntry >;
    type ViewMinorIntoIter: Iterator< Item = Self::ViewMinorEntry >;
    type ViewMinorEntry;

    /// Get a minor view.
    ///
    /// The order in which terms appear should be the same every time the
    /// function is called; however, the order need not be sorted.
    fn   view_minor( &self, index: Self::KeyMin ) -> Self::ViewMinor;
}

/// Entries appear in **strictly** ascending order, according to index.  Consecutive entries have distinct indices.
// #[auto_impl(&)] 
pub trait OracleMinorAscend: IndicesAndCoefficients
{
    type ViewMinorAscend: IntoIterator< IntoIter = Self::ViewMinorAscendIntoIter, Item = Self::ViewMinorAscendEntry >;
    type ViewMinorAscendIntoIter: Iterator< Item = Self::ViewMinorAscendEntry >;
    type ViewMinorAscendEntry;

    /// Get a minor view with entries sorted in ascending order of index.
    fn   view_minor_ascend( &self, index: Self::KeyMin ) -> Self::ViewMinorAscend;
}

/// Entries appear in **strictly** descending order, according to index.  Consecutive entries have distinct indices.
// #[auto_impl(&)] 
pub trait OracleMinorDescend: IndicesAndCoefficients
{
    type ViewMinorDescend: IntoIterator< IntoIter = Self::ViewMinorDescendIntoIter, Item = Self::ViewMinorDescendEntry >;
    type ViewMinorDescendIntoIter: Iterator< Item = Self::ViewMinorDescendEntry >;
    type ViewMinorDescendEntry;

    /// Get a minor view with entries sorted in descending order of index.
    fn   view_minor_descend( &self, index: Self::KeyMin ) -> Self::ViewMinorDescend;
}

/// Entries appear in **non-strictly** ascending order, according to index.  Consecutive entries may have identical indices.
// #[auto_impl(&)] 
pub trait OracleMinorAscendLax: IndicesAndCoefficients
{
    type ViewMinorAscendLax: IntoIterator< IntoIter = Self::ViewMinorAscendLaxIntoIter, Item = Self::ViewMinorAscendLaxEntry >;
    type ViewMinorAscendLaxIntoIter: Iterator< Item = Self::ViewMinorAscendLaxEntry >;
    type ViewMinorAscendLaxEntry;

    /// Get a major view with entries sorted in ascending order of index.
    fn   view_minor_ascend( &self, index: Self::KeyMin ) -> Self::ViewMinorAscendLax;
}

/// Entries appear in **non-strictly** descending order, according to index.   Consecutive entries may have identical indices.
// #[auto_impl(&)] 
pub trait OracleMinorDescendLax: IndicesAndCoefficients
{
    type ViewMinorDescendLax: IntoIterator< IntoIter = Self::ViewMinorDescendLaxIntoIter, Item = Self::ViewMinorDescendLaxEntry >;
    type ViewMinorDescendLaxIntoIter: Iterator< Item = Self::ViewMinorDescendLaxEntry >;
    type ViewMinorDescendLaxEntry;

    /// Get a major view with entries sorted in descending order of index.
    fn   view_minor_descend( &self, index: Self::KeyMin ) -> Self::ViewMinorDescendLax;
}




















//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    use crate::matrices::matrix_oracle_traits::IndicesAndCoefficients;


    
    #[test] 
    fn test_trait_implementation_demos() {
        // import crates
        use crate::matrices::matrix_oracle_traits::OracleMajor;          

        
        //  -----------------------------------------------------------------------------------------
        //  implement OracleMajor on a WRAPPER FOR &'a Vec< Vec< usize > > 
        //  -----------------------------------------------------------------------------------------

        // NB:  (1) the salient feature of this segment are that 
        //          (a) it uses an immutable reference to a vector-of-vectors.  This is because the
        //              `lifetime` of the reference can help us to overcome some issues concerning
        //              lifetimes and major views.
        //          (b) it introduces and uses a *wrapper* struct that contains a `&'a Vec< Vec< usize > >` 
        //              This is beecause Rust's "orphan rule" prevents the implementation of a
        //              foreign trait on a foreign type (see the following segment for details)

        // This struct is just a wrapper around `&'a Vec< Vec< (usize,usize)`
        struct VecOfVecReferenceWrapper< 'a > { vec_of_vec: &'a Vec< Vec< (usize,usize) > > }

        impl < 'a > IndicesAndCoefficients for VecOfVecReferenceWrapper< 'a > { 
            type KeyMin = usize;  type KeyMaj = usize;  type SnzVal = usize;
        }

        impl < 'a >

                OracleMajor for 

                VecOfVecReferenceWrapper< 'a >

        {
            type ViewMajor          = &'a [(usize,usize)];
            type ViewMajorIntoIter  = std::slice::Iter< 'a, (usize,usize) >;
            type ViewMajorEntry     = &'a (usize,usize);

            fn view_major( & self, index: Self::KeyMaj ) -> Self::ViewMajor {
                return self.vec_of_vec[index].as_slice()
            } 
        }

        // Get a major view
        let matrix = VecOfVecReferenceWrapper{ vec_of_vec: &vec![ vec![ (1,1), (2,2) ]  ] };
        let row = matrix.view_major( 0 );
        itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );

        // Get a major view from a *nested* immutable reference.
        let _row = (&(& matrix)).view_major( 0 );      


        //  -----------------------------------------------------------------------------------------
        //  implement OracleMajor on &'a Vec< Vec< usize > > 
        //  -----------------------------------------------------------------------------------------


        // NB:  SO FAR AS WE ARE CURRENTLY AWARE it is not possible to implement OracleMajor on 
        //      &'a Vec< Vec< usize > > **except when the code resides in the same file where this
        //      trait is defined.**  The reason for this is Rust's orphan rule 
        //      (c.f. https://doc.rust-lang.org/error-index.html#E0117), which prevents the
        //      implementation of a foreign trait on a foreign type (where "foreign" means a trait or
        //      type that is defined outside the file where you want to implement the trait).
        //      Therfore the following code compiles in the unit tests for the file
        //      src/oat_rust/matrices/entry_lookup.rs, but it does not compile in doc string tests or
        //      unit tests defined in other files.

        impl < 'a > IndicesAndCoefficients for &'a Vec < Vec < (usize,usize) > > { 
            type KeyMin = usize;  type KeyMaj = usize;  type SnzVal = usize;
        }

        impl < 'a >

                OracleMajor for 

                &'a Vec < Vec < (usize,usize) > > 

        {
            type ViewMajor          =   &'a [(usize,usize)];
            type ViewMajorIntoIter  =   std::slice::Iter< 'a, (usize,usize) >;
            type ViewMajorEntry     =   &'a (usize,usize);

            fn view_major( & self, index: Self::KeyMaj ) -> Self::ViewMajor {
                return self[index].as_slice()
            } 
        }

        // Get a major view from an immutable reference.     
        let mut matrix = vec![ vec![ (1,1), (2,2) ]  ];
        let row = (& matrix).view_major( 0 );
        itertools::assert_equal( row, vec![ (1,1), (2,2) ].iter() );

        // Get a major view from a *nested* immutable reference.
        let _row = (&(& matrix)).view_major( 0 );      
        
        // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        let _row_a = (& matrix).view_major( 0 );
        matrix.push( vec![ (4,4), (5,5)] );
        let row_b = (& matrix).view_major( 1 );   

        itertools::assert_equal(row_b, vec![ (4,4), (5,5)].iter() );
        // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );        
        
        
        //  -----------------------------------------------------------------------------------------
        //  implement OracleMajor on &'a mut Vec< Vec< usize > > 
        //  -----------------------------------------------------------------------------------------

        // NB:  SO FAR AS WE ARE CURRENTLY AWARE it is not possible to implement OracleMajor on &'a mut Vec< Vec< usize > > 
        //      in a manner similar to the example above.  Indeed, the following code is a facsimile of the 
        //      section above, with `&'a mut` replaceing `&'a`.  It does not seem to compile, due to 
        //      a lifetime conflict.

        // impl < 'a, 'b, (usize,usize) >

        //         OracleMajor
        //         <  usize,  &'a [(usize,usize)]  > for 

        //         &'a mut Vec < Vec < (usize,usize) > > 

        //         where   (usize,usize): Clone,
        //                 'b: 'a,
        // {
        //     fn view_major( & self, index: usize ) -> &'a [(usize,usize)] {
        //         return self[index].as_slice()
        //     } 
        // }

        // // Get a major view from an immutable reference.     
        // let mut matrix = vec![ vec![ 1, 2, 3 ]  ];
        // let row = (& matrix).view_major( 0 );
        // itertools::assert_equal( row, vec![ vec![ (1,1), (2,2) ]  ].iter() );

        // // Get a major view from a *nested* immutable reference.
        // let _row = (&(& matrix)).view_major( 0 );      
        
        // // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        // let _row_a = (& matrix).view_major( 0 );
        // matrix.push( vec![ (4,4), (5,5) ] );
        // let row_b = (& matrix).view_major( 0 );   

        // itertools::assert_equal(row_b.iter().cloned() , vec![ 1,2,3 ] );
        // // itertools::assert_equal( _row_a.iter().cloned() , vec![ 1,2,3] );           


        //  -----------------------------------------------------------------------------------------
        //  implement OracleMajor on a struct with a lifetime generic parameter
        //  -----------------------------------------------------------------------------------------

        //  NB: SO FAR AS WE ARE AWARE it is not straightforward to implement `OracleMajor` in a way
        //      that leverages a generic lifetime parameter associated with a struct.  This is in 
        //      contrast to some of the example shown below, which *do* leverage the lifetime
        //      parameter of a struct (it is key to note that different traits are implemented, and
        //      in the later examples the trait explicitly involves a lifetime parameter).  
        //      
        //  The following is a record of our "best attempt" at an implementation of `OracleMajor` that
        //  does make use of a lifetime parameter.  The attempt failed, but others may have more
        //  success in the future.

        // use std::marker::PhantomData;   

        // // Define the struct
        // pub struct VecOfVecWithLifetime
        // < 'a, IndexCoeffPair >
        // {
        //     pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
        //     phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
        // }

        // // Implement the trait
        // impl < 'a >

        //         OracleMajor
        //         <  usize, Cloned<Rev<std::slice::Iter<'a, (usize,usize)>>>  > for 

        //         VecOfVecWithLifetime< 'a, (usize,usize) >

        //         where   (usize,usize): Clone,
        // {
        //     fn view_major( & self, index: usize ) -> Cloned<Rev<std::slice::Iter<'a, (usize,usize)>>> {
        //         return self.vec_of_vec[index].iter().rev().cloned()
        //     } 
        // }


        //  -----------------------------------------------------------------------------------------
        //  -----------------------------------------------------------------------------------------
        //  Documentation of alternate oracle trait: with lifetimes
        //  -----------------------------------------------------------------------------------------
        //  -----------------------------------------------------------------------------------------        
        use std::marker::PhantomData;

        /// Entries may not appear in sorted order.
        pub trait OracleMajorWithLifetime< 'a, MajKey, ViewMajor >
        {
            /// Get a major view.
            ///
            /// The order in which terms appear should be the same every time the
            /// function is called; however, the order need not be sorted.
            fn   view_major_with_life<'b: 'a>( &'b self, index: MajKey ) -> ViewMajor;
        }

        //  -----------------------------------------------------------------------------------------
        //  implement OracleMajorWithLifetime on &'a mut Vec< Vec< (usize,usize) > >
        //  -----------------------------------------------------------------------------------------

        impl < 'a >

                OracleMajorWithLifetime
                <  'a, usize,  &'a [(usize,usize)]  > for 

                &'a mut Vec < Vec < (usize,usize) > > 

                where   (usize,usize): Clone,
        {
            fn view_major_with_life<'b: 'a>( &'b self, index: usize ) -> &'a [(usize,usize)] {
                return self[index].as_slice()
            } 
        }

        // Define a mutable reference to a matrix
        let mut matrix_core = vec![ vec![ (1,1), (2,2) ]  ];
        let matrix_ref = &mut matrix_core;

        // Get a major view from the mutable reference.     
        let row = matrix_ref.view_major_with_life( 0 );
        itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );

        // Get a major view from a *nested* immutable reference.
        let _row = (& matrix_ref).view_major_with_life( 0 );      
        
        // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` 
        // below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        let _row_a = matrix_ref.view_major_with_life( 0 );
        matrix_ref.push( vec![ (4,4), (5,5) ] );
        let row_b = matrix_ref.view_major_with_life( 1 );   

        itertools::assert_equal(row_b.iter().cloned() , vec![ (4,4), (5,5) ] );
        // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] ); 
        
        
        //  -----------------------------------------------------------------------------------------
        //  implement OracleMajorWithLifetime on VecOfVecWithLifetime< 'a, (usize,usize) >
        //  -----------------------------------------------------------------------------------------

        // DEFINE THE STRUCT
        pub struct VecOfVecWithLifetime
                    < 'a, IndexCoeffPair >
        
        {
            // pub major_dimension: MajorDimension, 
            pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
            pub phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
        }

        // DEFINE ITS METHODS
        impl    < 'a, IndexCoeffPair >
                
                VecOfVecWithLifetime 
                < 'a, IndexCoeffPair > 

        {
            // Make new (empty) VecOfVec. 
            pub fn new( vecvec: Vec < Vec < IndexCoeffPair > > ) -> Self  
            {
                VecOfVecWithLifetime{   
                            // major_dimension: major_dimension,
                            vec_of_vec:     vecvec,                    
                            phantom_kvpair: PhantomData,                   
                        }
            }
        }

        // IMPLEMENT OracleMajorWithLifetime
        impl < 'a >

                OracleMajorWithLifetime
                <  'a, usize,  &'a Vec< (usize,usize) >  > for 

                VecOfVecWithLifetime< 'a, (usize,usize) >

                where   (usize,usize): Clone,
        {
            fn view_major_with_life<'b: 'a>( &'b self, index: usize ) -> &'a Vec< (usize,usize) > {
                return &self.vec_of_vec[index]
            } 
        }        
        
        // Define a mutable reference to a matrix
        let mut matrix = VecOfVecWithLifetime::new( vec![ vec![ (1,1), (2,2) ]  ] );

        // Get a major view from the mutable reference.     
        let row = matrix.view_major_with_life( 0 );
        itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );

        // Get a major view from a *nested* immutable reference.
        let _row = (& matrix).view_major_with_life( 0 );      
        
        // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        let _row_a = matrix.view_major_with_life( 0 );
        matrix.vec_of_vec.push( vec![ (4,4), (5,5) ] );
        let row_b = matrix.view_major_with_life( 1 );   

        itertools::assert_equal(row_b.iter().cloned() , vec![ (4,4), (5,5) ] );
        // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );         


    }

}        
