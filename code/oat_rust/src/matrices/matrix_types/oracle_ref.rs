//! Wrapper for an immutable reference to a matrix oracle; implements every oracle trait that the wrapped object implements; also implements `Copy` and `Clone`.

use crate::matrices::matrix_oracle_traits::{OracleMajor, OracleMajorAscend, OracleMajorDescend, OracleMinor, OracleMinorAscend, OracleMinorDescend, IndicesAndCoefficients};



/// Wrapper for an immutable reference to a matrix oracle; implements every oracle trait that the wrapped object implements; also implements `Copy` and `Clone`.
/// 
/// If `T` implements a matrix oracle trait, it does not necessarily follow that `& T` implements the trait also.  However, `OracleRef< 'a, T >` does.
pub struct OracleRef< 'a, T > { oracle: & 'a T }

impl < 'a, T > OracleRef< 'a, T > {
    /// Wrap a reference to an oracle inside an `OracleRef`.
    pub fn new( oracle: & 'a T ) -> OracleRef< 'a, T > { OracleRef{ oracle } }
}

//  Clone and copy
//  NOTE: We implement these manually because the derive method of implementation imposes an unnecessary constraint that every type parameter be clone and copy, c.f. https://doc.rust-lang.org/std/marker/trait.Copy.html
impl < 'a, T > Clone for OracleRef< 'a, T > { fn clone(&self) -> Self { OracleRef::new( self.oracle ) } }
impl < 'a, T > Copy for OracleRef< 'a, T > {}


//  IndicesAndCoefficients
impl < 'a, T: IndicesAndCoefficients > IndicesAndCoefficients for OracleRef< 'a, T >
    {   type KeyMin = T::KeyMin; type KeyMaj = T::KeyMaj; type SnzVal = T::SnzVal;  }


//  Oracle major
impl < 'a, T >

    OracleMajor for

    OracleRef< 'a, T >

    where  
            T:              OracleMajor + IndicesAndCoefficients,
            T::ViewMajor:   IntoIterator,

{   
    type ViewMajor          =   T::ViewMajor;
    type ViewMajorEntry     =   < Self::ViewMajor as IntoIterator >::Item;
    type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;     
    
    fn   view_major( &self, index: Self::KeyMaj ) -> Self::ViewMajor { self.oracle.view_major( index ) }    
}


//  Oracle major ascend
impl < 'a, T >

    OracleMajorAscend for

    OracleRef< 'a, T >

    where  
            T:                  OracleMajorAscend + IndicesAndCoefficients,
            T::ViewMajorAscend: IntoIterator,            

{   
    type ViewMajorAscend          =   T::ViewMajorAscend;
    type ViewMajorAscendEntry     =   < Self::ViewMajorAscend as IntoIterator >::Item;
    type ViewMajorAscendIntoIter  =   < Self::ViewMajorAscend as IntoIterator >::IntoIter;     
    
    fn   view_major_ascend( &self, index: Self::KeyMaj ) -> Self::ViewMajorAscend { self.oracle.view_major_ascend( index ) }    
}


//  Oracle major descend
impl < 'a, T >

    OracleMajorDescend for

    OracleRef< 'a, T >

    where  
            T:                      OracleMajorDescend,
            T::ViewMajorDescend:    IntoIterator,

{   
    type ViewMajorDescend          =   T::ViewMajorDescend;
    type ViewMajorDescendEntry     =   < Self::ViewMajorDescend as IntoIterator >::Item;
    type ViewMajorDescendIntoIter  =   < Self::ViewMajorDescend as IntoIterator >::IntoIter; 

    fn   view_major_descend( &self, index: Self::KeyMaj ) -> Self::ViewMajorDescend { self.oracle.view_major_descend( index ) }    
}


//  Oracle minor
impl < 'a, T >

    OracleMinor for

    OracleRef< 'a, T >

    where  
            T:              OracleMinor,
            T::ViewMinor:   IntoIterator,  

{   
    type ViewMinor          =   T::ViewMinor;
    type ViewMinorEntry     =   < Self::ViewMinor as IntoIterator >::Item;
    type ViewMinorIntoIter  =   < Self::ViewMinor as IntoIterator >::IntoIter;     
    
    fn   view_minor( &self, index: Self::KeyMin ) -> Self::ViewMinor { self.oracle.view_minor( index ) }    
}


//  Oracle minor ascend
impl < 'a, T >

    OracleMinorAscend for

    OracleRef< 'a, T >

    where  
            T:                      OracleMinorAscend,
            T::ViewMinorAscend:     IntoIterator,

{   
    type ViewMinorAscend          =   T::ViewMinorAscend;
    type ViewMinorAscendEntry     =   < Self::ViewMinorAscend as IntoIterator >::Item;
    type ViewMinorAscendIntoIter  =   < Self::ViewMinorAscend as IntoIterator >::IntoIter;    

    fn   view_minor_ascend( &self, index: Self::KeyMin ) -> Self::ViewMinorAscend { self.oracle.view_minor_ascend( index ) }    
}


//  Oracle minor descend
impl < 'a, T >

    OracleMinorDescend for

    OracleRef< 'a, T >

    where  
            T:                      OracleMinorDescend,
            T::ViewMinorDescend:    IntoIterator,

{   
    type ViewMinorDescend          =   T::ViewMinorDescend;
    type ViewMinorDescendEntry     =   < Self::ViewMinorDescend as IntoIterator >::Item;
    type ViewMinorDescendIntoIter  =   < Self::ViewMinorDescend as IntoIterator >::IntoIter;     
    
    fn   view_minor_descend( &self, index: Self::KeyMin ) -> Self::ViewMinorDescend { self.oracle.view_minor_descend( index ) }    
}