//! Tools to define and enumerate filtered and unfiltered simplices
//! 
//! In this library,
//! - filtered simplices (meaning simplices with assigned birth times) are represented by the [`SimplexFiltered`] struct
//! - unfiltered simplices (meaning simplices without birth times) are represented by vectors with *strictly ascending entries*

use derive_getters::Dissolve;
use derive_new::new;



use crate::utilities::order::{JudgeOrder, JudgePartialOrder};


use std::cmp::Ordering;

use std::hash::Hash;



use std::fmt::Debug;


type Vertex = u16;




//  ================================================================================
//  TWIST ORDER OPERATOR - FILTERED SIMPLEX
//  ================================================================================


/// Orders simplices first by dimension (descending) then by filtration (ascending)
/// then by lexicographic order (ascending)
/// 
/// This order is specially designed to facilitate the `twist` optimization in 
/// persistent homology calculations.
/// 
/// This object implements the `JudgeOrder` and `JudgePartialOrder` traits.
#[derive(Clone,Debug,new)]
pub struct OrderOperatorTwistSimplexFiltered;

impl < FilVal: Ord + Clone + Debug >

    JudgePartialOrder
        < SimplexFiltered< FilVal > > for 

    OrderOperatorTwistSimplexFiltered 
{
    /// Orders simplices first by dimension (descending) then by filtration (ascending)
    /// then by lexicographic order (ascending)
    /// 
    /// This order is specially designed to facilitate the `twist` optimization in 
    /// persistent homology calculations.    
    fn judge_partial_cmp( & self, lhs: &SimplexFiltered< FilVal >, rhs: &SimplexFiltered< FilVal > ) -> Option<Ordering> {
        Some(self.judge_cmp(lhs,rhs))
    }
}

impl < FilVal: Ord + Clone + Debug >

    JudgeOrder
        < SimplexFiltered< FilVal > > for 

    OrderOperatorTwistSimplexFiltered 
{

    /// Orders simplices first by dimension (descending) then by filtration (ascending)
    /// then by lexicographic order (ascending)
    /// 
    /// This order is specially designed to facilitate the `twist` optimization in 
    /// persistent homology calculations.    
    fn judge_cmp( & self, lhs: &SimplexFiltered< FilVal >, rhs: &SimplexFiltered< FilVal > ) -> Ordering {
        // compare simplex dimensions
        let mut comp = rhs.vertices.len().cmp( & lhs.vertices.len() ); // note we put RHS on the LEFT here
        if comp != Ordering::Equal { return comp }

        // first compare filtration values
        comp = lhs.filtration.cmp( & rhs.filtration );
        if comp != Ordering::Equal { return comp }        

        // finally, compare simplices lexicographically
        lhs.vertices.cmp( & rhs.vertices )
    }
}







//  ================================================================================
//  COMBINATORIAL SIMPLEX (WEIGHTED) -- DEFINITION
//  ================================================================================



/// A simplex associated with a filtration value
#[derive(Debug, PartialEq, Clone, Eq, Hash, Dissolve)]
pub struct SimplexFiltered<FilVal: Clone + Debug>{
    pub filtration: FilVal,     // the filtration value of simplex
    pub vertices:   Vec<Vertex> // a sorted vector in strictly descending order recording vertices
}

impl <FilVal: Clone + Debug>
    SimplexFiltered<FilVal> {

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize { self.vertices.len() }

    /// Dimension of the simplex.
    pub fn dimension(&self) -> usize { self.vertices.len() - 1 }

    /// Filtraiton value of the simplex.
    pub fn filtration(&self) -> FilVal { self.filtration.clone() }

    /// The vertices of the simplex.
    pub fn vertices(&self) -> &Vec<Vertex> { &self.vertices }
}

impl <FilVal>
    PartialOrd for SimplexFiltered<FilVal> 
    where FilVal: Eq + PartialEq + Ord + PartialOrd + Clone + Debug

{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl <FilVal>
    Ord for SimplexFiltered<FilVal> 
    where FilVal: Eq + PartialEq + Ord + PartialOrd + Clone + Debug
{
    /// Compares simplices first by filtration value, and second by lexicographic order.
    fn cmp(&self, other: &Self) -> Ordering {

        // first compare filtration values
        let comp = self.filtration.cmp( & other.filtration );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        self.vertices.cmp( & other.vertices )
    }
}




