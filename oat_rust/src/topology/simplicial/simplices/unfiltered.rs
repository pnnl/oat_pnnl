//!  Data structures for unfiltered simplices
//! 
//! # Advisory: consider ordinary Rust vectors
//! 





use crate::utilities::sequences_and_ordinals::{SortedVec};

use std::cmp::Ordering;

use std::hash::Hash;



use std::fmt::Debug;



//  ===========================================================================
//  COMBINATORIAL SIMPLEX (UNWEIGHTED) -- DEFINITION
//  ===========================================================================


/// An unweighted simplex
/// 
/// The vertices should sorted in ascending order, however this condition is not
/// gauranteed.  Compare with [SimplexSafe], where strict ordering *is* gauranteed.
/// 
/// # Advisory
/// 
/// Plain Rust vectors can serve quite well to represent simplices, OAT provides
/// a substantial tool set to work with vector representations.  Using this more
/// primitive data format may simplifiy your code, and provide access to a 
/// broader tool set.
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Simplex< Vertex > 
{ 
    pub vertices: Vec< Vertex >     //  vertices should be sorted in ascending order
} 

impl    < Vertex > 
        Simplex
        < Vertex >   
        {
    
    pub fn num_vertices( &self ) -> usize { self.vertices.len() }
    pub fn dimension( &self ) -> usize { self.vertices.len() - 1 }
}        


impl    < Vertex >           
        PartialOrd for Simplex
        < Vertex >

    where   Vertex: Ord     {

    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl    < Vertex >   
        Ord for Simplex
        < Vertex >

    where Vertex: Ord   {

    fn cmp(&self, other: &Self) -> Ordering {

        // next compare simplex dimensions
        let comp = self.num_vertices().cmp( & other.vertices.len() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        self.vertices.cmp( & other.vertices )
    }
}

impl    < Vertex >   
        IntoIterator for Simplex
        < Vertex >      {

    type Item = Vertex;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter { self.vertices.into_iter() }
}



//  ================================================================================
//  COMBINATORIAL SIMPLEX (UNWEIGHTED) -- DEFINITION
//  ================================================================================


/// An unweighted simplex; the vertices should sorted in ascending order.
/// 
/// This is simply a wrapper around a `SortedVec`.  However, this wrapper has a special
/// order (first by length, second by lexicographic order of elements) which is useful
/// in homology computations.
/// 
/// # Advisory
/// 
/// Plain Rust vectors can serve quite well to represent simplices, OAT provides
/// a substantial tool set to work with vector representations.  Using this more
/// primitive data format may simplifiy your code, and provide access to a 
/// broader tool set.  Moreover, the safety gaurantees associated with this data
/// struture can often be replicated by simpler means.
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct SimplexSafe< Vertex: Ord + PartialEq + Eq > 
{ 
    vertices: SortedVec< Vertex >     //  vertices should be sorted in ascending order
} 

impl    < Vertex: Ord > 
        SimplexSafe
        < Vertex >   
        {
    
    /// Construct a new simplex
    pub fn new( vertices: SortedVec< Vertex > ) -> Self { SimplexSafe { vertices } }
    /// Construct from a vector; throws an error if the vector is not sorted in striclty ascending order
    pub fn from_vec( vertices: Vec< Vertex > ) -> Result< Self, Vec<Vertex > > where Vertex: Ord {
        SortedVec::new( vertices ).map(|x| SimplexSafe{ vertices: x } )        
    }        
    /// Number of vertices
    pub fn num_vertices( &self ) -> usize { self.vertices.len() }
    /// Dimension of the simplex
    pub fn dimension( &self ) -> usize { self.vertices.len() - 1 }
    /// Immutable reference to the internally stored sorted vector
    pub fn sorted_vec( &self ) -> & SortedVec< Vertex > { & self.vertices }
    /// Immutable reference to the internally stored vector 
    pub fn vertices( &self ) -> & Vec< Vertex > { self.vertices.vec() }    
    /// Consumes the simplex, returning the sorted vector of vertices
    pub fn into_sorted_vec( self ) -> SortedVec< Vertex > { self.vertices }    
    /// Consumes the simplex, returning the vector of vertices
    pub fn into_vec( self ) -> Vec< Vertex > { self.vertices.into_vec() }        
}        


impl    < Vertex >           
        PartialOrd for SimplexSafe
        < Vertex >

    where   Vertex: Ord     {

    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl    < Vertex >   
        Ord for SimplexSafe
        < Vertex >

    where Vertex: Ord   {

    fn cmp(&self, other: &Self) -> Ordering {

        // next compare simplex dimensions
        let comp = self.num_vertices().cmp( & other.num_vertices() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        return self.vertices().cmp( other.vertices() )
    }
}

impl    < Vertex: Ord >   
        IntoIterator for SimplexSafe
        < Vertex >      {

    type Item = Vertex;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter { self.into_vec().into_iter() }
}