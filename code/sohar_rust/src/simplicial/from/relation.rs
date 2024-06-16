//! Simplicial complex by a superset of the maximal simplices (c.f. Dowker complex)
//! 
//! This modules includes data structures for
//! - the boundary matrix in row-major form
//! - the major views of this matrix

use oat_rust::matrices::matrix_oracle_traits::{IndicesAndCoefficients, OracleMajorAscend, OracleMajorDescend, OracleMinorAscend, OracleMinorDescend};
use oat_rust::matrices::operations::umatch::row_major::ParetoShortCircuit;
use oat_rust::rings::operator_traits::{Semiring, Ring};

use std::collections::HashSet;
use std::hash::Hash;
use std::iter::FromIterator;
use std::marker::PhantomData;

// use oat_rust::boundary_matrices::{SimplexBoundaryAscend, SimplexBoundaryDescend};


//  COUBOUNDARY VECTOR DESCENDING
//  -----------------------------------------------------------------------------------



/// Iterates over the terms of the coboundary of a simplex in descending lexicographic order.
/// 
/// The coefficient is calculated as `(-1)^k`, where `xk` is the vertex which appears twice in the template.
/// 
/// # Examples
/// 
/// ```
/// use sohar_rust::simplicial::from::relation::CoboundaryDowkerDescend;
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use std::collections::HashSet;
/// use std::iter::FromIterator;
/// 
/// let ring_operator = PrimeOrderFieldOperator::new(3);
/// let simplex = vec![1,3];
/// let dowker_simplices = vec![ HashSet::from_iter( vec![0,1,2,3,4] ) ];
/// let coboundary = CoboundaryDowkerDescend::from_vec_of_dowker_simplices( simplex, &dowker_simplices, ring_operator );
/// 
/// itertools::assert_equal( coboundary, vec![ (vec![1,3,4], 1), (vec![1,2,3], 2), (vec![0,1,3], 1) ]);
/// ```
#[derive(Clone, Debug)]
pub struct CoboundaryDowkerDescend< Vertex, RingOperator, RingElement >
    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,
{
    next_cofacet_opt:           Option< Vec< Vertex > >,
    next_coefficient:           RingElement,
    vertices_to_insert:         Vec< Vertex >,              // should be sorted in ascending order
    retrieval_locus:            usize,                      // the place where the vertex inserted into `next_cofacet` was drawn from 
    insertion_locus:            usize,                      // this points to the place in `next_cofacet` where we inserted a vertex    
    ring_operator:              RingOperator,
}   


impl < Vertex, RingOperator, RingElement >

CoboundaryDowkerDescend
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Clone + Hash + Ord,
{
    /// Generates a [CoboundaryDowkerAscend] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_vec_of_dowker_simplices( facet: Vec< Vertex >, dowker_simplices: &Vec< HashSet< Vertex > >, ring_operator: RingOperator ) -> Self {

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return CoboundaryDowkerDescend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator:              ring_operator,         // these are arbitrary values                    
            }                      
        }

        //  obtain a list of vertices to insert
        //  ---------------------------------------------        
        let facet_vertex_set = HashSet::from_iter( facet.iter().cloned() );
        let mut vertices_to_insert = HashSet::new();

        // take the union of all dowker super-simplices
        for dowker_simplex in dowker_simplices
                                        .iter() // iterate over dowker siplices
                                        .filter( |x| x.is_superset( &facet_vertex_set) ) // only keep those that contain the facet
            {
                vertices_to_insert.extend( dowker_simplex.iter().cloned() )
            }

        // take a set difference: vertices_to_insert \ facet
        for vertex in facet.iter() { vertices_to_insert.remove( vertex ); }

        // convert to a vector
        let mut vertices_to_insert: Vec< Vertex > = vertices_to_insert.into_iter().collect();

        // sort the vector
        vertices_to_insert.sort();
        

        //  if there are no vertices to insert, then the coboundary is zero
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `next_cofacet_opt` is None
            return CoboundaryDowkerDescend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),     // these are abitrary values            
                vertices_to_insert:         vec![],                  // these are abitrary values  
                retrieval_locus:            0,                       // these are abitrary values  
                insertion_locus:            0,                       // these are abitrary values  
                ring_operator:              ring_operator,           // these are abitrary values                    
            }                      
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = ring_operator.minus_one_to_power( facet.len() );
        let mut next_cofacet = facet.clone();
        let mut insertion_locus = facet.len();
        let retrieval_locus = vertices_to_insert.len() - 1;        
        let inserted_vertex = vertices_to_insert[retrieval_locus].clone();
        while   ( insertion_locus > 0 ) 
                && 
                ( inserted_vertex < next_cofacet[ insertion_locus - 1 ] ) {
            insertion_locus -= 1;
            coefficient = ring_operator.negate( coefficient );
        }
        next_cofacet.insert( insertion_locus, inserted_vertex );

        CoboundaryDowkerDescend{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert:         vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus:            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus:            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator:              ring_operator,             
        }

    }

    // /// Generates a [CoboundaryDowkerAscend] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< Vertex >, dowker_matrix_csr: Vec< HashSet< Vertex > >, dowker_matrix_csc: Vec< HashSet< Vertex > > ) {        
    // }
}


impl < Vertex, RingOperator, RingElement >

    Iterator for

    CoboundaryDowkerDescend
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,   
{
    type Item = (Vec<Vertex>, RingElement);

    fn next( &mut self ) -> Option< Self::Item >{

        // println!("{:?} -- DELETE THIS AND ALL DEBUG REQUIREMENTS ON THIS IMPLEMENTATION OF ITER", &self);

        match self.next_cofacet_opt {
            None => { return None }
            Some( ref mut next_cofacet ) => {

                // make a copy of the return value, before we change the internal state of the iterator
                let return_value    =   ( next_cofacet.clone(), self.next_coefficient.clone() );

                // if there are more vertices to insert, update the internal state with data for the next entry
                if self.retrieval_locus > 0 {
                    
                    // grab the next vertex
                    self.retrieval_locus -= 1; 
                    let inserted_vertex     =   self.vertices_to_insert[ self.retrieval_locus ].clone();
                    
                    // update pointers and coefficients
                    while   ( self.insertion_locus > 0 ) 
                            && 
                            ( inserted_vertex < next_cofacet[ self.insertion_locus - 1 ] ) {
                        next_cofacet[ self.insertion_locus ] = next_cofacet[ self.insertion_locus - 1 ].clone();
                        self.insertion_locus -= 1;
                        self.next_coefficient = self.ring_operator.negate( self.next_coefficient.clone() );
                    }
                    // update the cofacet
                    next_cofacet[ self.insertion_locus ] = inserted_vertex;
                } 
                // otherwise flag the fact that there will be no more entries, after we return the one that we have already extracted
                else {
                    self.next_cofacet_opt = None;
                }

                return Some( return_value )
            }
        }

    }
}       



//  COUBOUNDARY VECTOR ASCENDING
//  -----------------------------------------------------------------------------------

/// Iterates over the terms of the coboundary of a simplex in ascending lexicographic order.
/// 
/// The coefficient is calculated as `(-1)^k`, where `xk` is the vertex which appears twice in the template.
/// 
/// # Examples
/// 
/// ```
/// use sohar_rust::simplicial::from::relation::CoboundaryDowkerAscend;
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use std::collections::HashSet;
/// use std::iter::FromIterator;
/// 
/// let ring_operator = PrimeOrderFieldOperator::new(3);
/// let simplex = vec![1,3];
/// let dowker_simplices = vec![ HashSet::from_iter( vec![0,1,2,3,4] ) ];
/// let coboundary = CoboundaryDowkerAscend::from_vec_of_dowker_simplices( simplex, &dowker_simplices, ring_operator );
/// 
/// itertools::assert_equal( coboundary, vec![ (vec![0,1,3], 1), (vec![1,2,3], 2), (vec![1,3,4], 1) ]);
/// ```
#[derive(Clone, Debug)]
pub struct CoboundaryDowkerAscend< Vertex, RingOperator, RingElement >
    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,
{
    next_cofacet_opt:           Option< Vec< Vertex > >,
    next_coefficient:           RingElement,
    vertices_to_insert:         Vec< Vertex >,              // should be sorted in ascending order
    retrieval_locus:            usize,                      // the place where the vertex inserted into `next_cofacet` was drawn from 
    insertion_locus:            usize,                      // this points to the place in `next_cofacet` where we inserted a vertex    
    ring_operator:              RingOperator,
}   


impl < Vertex, RingOperator, RingElement >

    CoboundaryDowkerAscend
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Clone + Hash + Ord,
{
    /// Generates a [CoboundaryDowkerAscend] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_vec_of_dowker_simplices( facet: Vec< Vertex >, dowker_simplices: &Vec< HashSet< Vertex > >, ring_operator: RingOperator ) -> Self {

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return CoboundaryDowkerAscend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator:              ring_operator,         // these are arbitrary values                    
            }                      
        }

        //  obtain a list of vertices to insert
        //  ---------------------------------------------        
        let facet_vertex_set = HashSet::from_iter( facet.iter().cloned() );
        let mut vertices_to_insert = HashSet::new();

        // take the union of all dowker super-simplices
        for dowker_simplex in dowker_simplices
                                        .iter() // iterate over dowker siplices
                                        .filter( |x| x.is_superset( &facet_vertex_set) ) // only keep those that contain the facet
            {
                vertices_to_insert.extend( dowker_simplex.iter().cloned() )
            }

        // take a set difference: vertices_to_insert \ facet
        for vertex in facet.iter() { vertices_to_insert.remove( vertex ); }

        // convert to a vector
        let mut vertices_to_insert: Vec< Vertex > = vertices_to_insert.into_iter().collect();

        // sort the vector
        vertices_to_insert.sort();
        

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `vertices_to_insert` is empty
            return CoboundaryDowkerAscend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),
                vertices_to_insert:         vec![],     
                retrieval_locus:            0,          
                insertion_locus:            0,          
                ring_operator:              ring_operator,                
            }                      
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = RingOperator::one();
        let mut next_cofacet = facet.clone();
        let mut insertion_locus = 0;
        let retrieval_locus = 0;        
        let inserted_vertex = vertices_to_insert[retrieval_locus].clone();
        while   ( insertion_locus < next_cofacet.len() ) 
                && 
                ( inserted_vertex > next_cofacet[ insertion_locus ] ) {
            insertion_locus += 1;
            coefficient = ring_operator.negate( coefficient );
        }
        next_cofacet.insert( insertion_locus, inserted_vertex );

        CoboundaryDowkerAscend{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert:         vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus:            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus:            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator:              ring_operator,             
        }

    }

    // /// Generates a [CoboundaryDowkerAscend] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< Vertex >, dowker_matrix_csr: Vec< HashSet< Vertex > >, dowker_matrix_csc: Vec< HashSet< Vertex > > ) {        
    // }
}


impl < Vertex, RingOperator, RingElement >

    Iterator for

    CoboundaryDowkerAscend
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,   
{
    type Item = (Vec<Vertex>, RingElement);

    fn next( &mut self ) -> Option< Self::Item >{

        println!("THIS COULD BE IMPROVED WITH SHORT-CIRCUITING");

        // println!("{:?} -- DELETE THIS AND ALL DEBUG REQUIREMENTS ON THIS IMPLEMENTATION OF ITER", &self);

        match self.next_cofacet_opt {
            None => { return None }
            Some( ref mut next_cofacet ) => {

                // make a copy of the return value, before we change the internal state of the iterator
                let return_value    =   ( next_cofacet.clone(), self.next_coefficient.clone() );

                // if there are more vertices to insert, update the internal state with data for the next entry
                if self.retrieval_locus + 1 < self.vertices_to_insert.len() {
                    
                    // grab the next vertex
                    self.retrieval_locus += 1; 
                    let inserted_vertex     =   self.vertices_to_insert[ self.retrieval_locus ].clone();
                    
                    // update pointers and coefficients
                    while   ( self.insertion_locus + 1 < next_cofacet.len() ) 
                            && 
                            ( inserted_vertex > next_cofacet[ self.insertion_locus + 1 ] ) {
                        next_cofacet[ self.insertion_locus ] = next_cofacet[ self.insertion_locus +1 ].clone();
                        self.insertion_locus += 1;
                        self.next_coefficient = self.ring_operator.negate( self.next_coefficient.clone() );
                    }
                    // update the cofacet
                    next_cofacet[ self.insertion_locus ] = inserted_vertex;
                } 
                // otherwise flag the fact that there will be no more entries, after we return the one that we have already extracted
                else {
                    self.next_cofacet_opt = None;
                }

                return Some( return_value )
            }
        }

    }
}    



impl < Vertex, RingOperator, RingElement >

    ParetoShortCircuit< (Vec<usize>, RingElement) > for

    CoboundaryDowkerAscend
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,   
{
    fn pareto_short_circuit(& self) -> Option< (Vec<usize>, RingElement) > {
        return None
    }
}








//  ---------------------------------------------------------------------------
//  DOWKER BOUNDARY MATRIX
//  ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct DowkerComplexBoundaryMatrixRowMajor
            < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    dowker_simplices:       Vec< HashSet< Vertex > >,
    ring_operator:          RingOperator,
    phantom_ringelement:    PhantomData< RingElement >,
}


impl < Vertex, RingOperator, RingElement >
    
    DowkerComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    pub fn new( dowker_simplices: Vec< HashSet< Vertex > >, ring_operator: RingOperator ) -> Self {
        DowkerComplexBoundaryMatrixRowMajor{ dowker_simplices, ring_operator, phantom_ringelement: PhantomData }
    }
}


//  INDICES AND COEFFICIENTS
//  ------------------------------------------


impl < Vertex, RingOperator, RingElement >
    
    IndicesAndCoefficients for

    DowkerComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type KeyMaj = Vec< Vertex >; type KeyMin = Vec< Vertex >; type SnzVal = RingElement;
}


//  ORACLE MAJOR ASCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    OracleMajorAscend for

    DowkerComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMajorAscend            =   CoboundaryDowkerAscend< Vertex, RingOperator, RingElement >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;
    type ViewMajorAscendEntry       =   ( Self::KeyMaj, Self::SnzVal );

    fn view_major_ascend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMajorAscend {
        CoboundaryDowkerAscend::from_vec_of_dowker_simplices( keymaj, & self.dowker_simplices, self.ring_operator.clone() )
    }
}

//  ORACLE MAJOR DESCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
OracleMajorDescend for

    DowkerComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMajorDescend           =   CoboundaryDowkerDescend< Vertex, RingOperator, RingElement >;
    type ViewMajorDescendIntoIter   =   Self::ViewMajorDescend;
    type ViewMajorDescendEntry      =   ( Self::KeyMaj, Self::SnzVal );

    fn view_major_descend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMajorDescend {

        println!("GOT THROUGH THIS FILE ANE REMOVE THE DEBUG REQUIREMENTS");
        CoboundaryDowkerDescend::from_vec_of_dowker_simplices( keymaj, & self.dowker_simplices, self.ring_operator.clone() )
    }
}



//  ORACLE MINOR ASCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    OracleMinorAscend for

    DowkerComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMinorAscend            =   SimplexBoundaryAscend< Vertex, RingOperator, RingElement >;
    type ViewMinorAscendIntoIter    =   Self::ViewMinorAscend;
    type ViewMinorAscendEntry       =   ( Self::KeyMaj, Self::SnzVal );

    fn view_minor_ascend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMinorAscend {
        SimplexBoundaryAscend::new( keymaj, self.ring_operator.clone() )
    }
}

//  ORACLE MINOR DESCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    OracleMinorDescend for

    DowkerComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMinorDescend           =   SimplexBoundaryDescend< Vertex, RingOperator, RingElement >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;
    type ViewMinorDescendEntry      =   ( Self::KeyMaj, Self::SnzVal );

    fn view_minor_descend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMinorDescend {
        SimplexBoundaryDescend::new( keymaj, self.ring_operator.clone() )
    }
}













#[cfg(test)]
mod tests {
    

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use std::collections::HashSet;
    use std::iter::FromIterator;
    use oat_rust::matrices::matrix_types::oracle_ref::OracleRef;    
    use oat_rust::matrices::matrix_oracle_traits::OracleMinorDescend;
    use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor};
    use crate::simplicial::from::relation::DowkerComplexBoundaryMatrixRowMajor;
    use crate::simplicial::simplices::{subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend};
    use oat_rust::utilities::partial_order::OrderComparatorAutoLt;        
    use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    use oat_rust::matrices::debug::{verify_viewmajorascend_compatible_with_viewminordescend, verify_viewmajorascend_compatible_with_viewmajordescend, verify_viewminorascend_compatible_with_viewminordescend};
    use oat_rust::utilities::random::rand_sequences;

    fn test_dowker_boundary_helper( dowker_simplices_vec: Vec<Vec<usize>> ){

        let ring_operator                =   PrimeOrderFieldOperator::new(47);
        let dowker_simplices: Vec<_>    =   dowker_simplices_vec.iter().cloned().map(|x| HashSet::from_iter(x) ).collect();    

        // define an iterator to run over all simplices; simplices are ordered first by dimension (ascending), then by lexicographic order (descending)
        let iter_keymaj = subsimplices_dim_0_thru_d_iter_descend( &dowker_simplices_vec, 2 );        
        
        // define the boundary matrix
        let boundary_matrix = DowkerComplexBoundaryMatrixRowMajor::new( dowker_simplices, ring_operator.clone() );

        // check that ascending major views agree with descending MINOR views
        verify_viewmajorascend_compatible_with_viewminordescend(
                OracleRef::new( & boundary_matrix ),
                iter_keymaj.clone(),
                iter_keymaj.clone(),
            );

        // check that ascending major views agree with descending MAJOR views            
        verify_viewmajorascend_compatible_with_viewmajordescend(
                OracleRef::new( & boundary_matrix ),
                iter_keymaj.clone(),
            );

        // check that ascending MINOR views agree with descending MINOR views            
        verify_viewminorascend_compatible_with_viewminordescend(
            OracleRef::new( & boundary_matrix ),
            iter_keymaj.clone(),            
        )        
    }


    #[test]
    fn test_dowker_boundary_small(){
        
        // used for enumerating simplices
        let dowker_simplices_vec =  
                vec![ 
                        vec![0,1,2],
                        vec![0,3],                      
                        vec![1,3], 
                        vec![2,3]                                           
                    ];      

        test_dowker_boundary_helper(dowker_simplices_vec)
    }



    #[test]
    fn test_dowker_boundary_big(){
        
        for _ in 0..10 {
            
            let dowker_simplices_vec    =   rand_sequences(10, (0..7).step_by(2), 0.1 );

            test_dowker_boundary_helper(dowker_simplices_vec)
        }

    }    



}    


























//  =================================================================================
//  MATRIX FACOTRIZATION AND HOMOLOGY COMPUTATION
//  =================================================================================



// use std::collections::HashSet;
// use std::iter::FromIterator;
use itertools::Itertools;
use num::rational::Ratio;


use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor, new_umatchrowmajor_with_clearing, UmatchRowMajor};
use crate::simplicial::boundary::{SimplexBoundaryAscend, SimplexBoundaryDescend};
use crate::simplicial::simplices::{subsimplices_dim_0_thru_d_iter_ascend, subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend, subsimplices_dim_d_iter_descend};
use oat_rust::utilities::partial_order::{OrderComparatorAutoLt, OrderComparatorLtByKey};        
use oat_rust::rings::operator_structs::field_prime_order::{PrimeOrderFieldOperator};
use oat_rust::rings::operator_structs::ring_native::{DivisionRingNative, field_rational_i64};
use oat_rust::utilities::iterators::is_sorted::IsSortedBy;


pub enum UseClearing{ Yes, No }

type RingOperator = DivisionRingNative< Ratio<i64> >;
type RingElement = Ratio<i64>;

/// Compute a U-match factorization of the (all-dimensions) boundary matrix of a Dowker complex, over the rational numbers.
/// 
/// Input is a relation formatted as a `m x n` vec-of-rowvec matrix `M`.  The simplicial complex consists of
/// every subset `S` of `{0,..,n-1}` such that `M[r,S]= [1,..,1]` for some row index `r`
pub fn umatch_from_dowker( 
        dowker_simplices: & Vec<Vec<usize>>, 
        maxdim: usize,
        clearing: UseClearing,
        ) 
        -> 
        // (
        //     DowkerComplexBoundaryMatrixRowMajor<usize, RingOperator, RingElement>,
            UmatchRowMajor<
                    DowkerComplexBoundaryMatrixRowMajor<
                            usize, 
                            RingOperator, 
                            RingElement
                        >, 
                        RingOperator, 
                    OrderComparatorLtByKey<Vec<usize>, RingElement, (Vec<usize>, RingElement), OrderComparatorAutoLt<Vec<usize>>>, 
                    OrderComparatorLtByKey<Vec<usize>, RingElement, (Vec<usize>, RingElement), OrderComparatorAutoLt<Vec<usize>>>
                >
        // )
                
    {

    // Verify that inputs are sorted
    for vec in dowker_simplices {
        if ! vec.windows(2).all(|w| w[0] <= w[1]) { panic!("The vertices of each input simplex must be sorted.") }
    }

    // Define the ring operator for the rational numbers.
    // let ring_operator = RingOperator::new();
    let ring_operator = field_rational_i64();

    // We will build a dowker complex.
    // A dowker complex is defined by a vertex set V and a family S
    // of subsets of V.  A subset of V forms a simplex iff it is 
    // a subset of some element of S.  We refer to the elements 
    // of S as "dowker simplices".     

    // Here is another format; it replaces inner vectors with sets; each format has its own use.
    let dowker_simplices_hashset_format  =  dowker_simplices
                                                .iter()
                                                .cloned()
                                                .map( |x| HashSet::from_iter( x ) )
                                                .collect_vec();

    // Build the boundary matrix.
    // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
    let boundary_matrix = DowkerComplexBoundaryMatrixRowMajor::new( dowker_simplices_hashset_format.clone(), ring_operator.clone() );

    // This iterates over simplices in ascending order of 
    // dimension (first) and descending lexicographic order (second).
    // When computing homology of dimension d, we only need to iterate
    // over simplices of dimension d and below.   

    let iter_keymaj = 
        (0..maxdim+1)
            .map(|x| subsimplices_dim_d_iter_descend(&dowker_simplices, x) )
            .flatten();

    let fun = match clearing { 
        UseClearing::No   => new_umatchrowmajor, 
        UseClearing::Yes  => new_umatchrowmajor_with_clearing  
    };    

    // Compute a umatch factorization of the boundary matrix.
    // For details on what this factorization entails, see the paper 
    // "U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co) homology"
    // by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.  You can also check out the
    // oat_rust documentation for `umatch`.
    let umatch = fun(
            boundary_matrix, 
            iter_keymaj, 
            ring_operator.clone(), 
            OrderComparatorAutoLt::new(), 
            OrderComparatorAutoLt::new(), 
        );   

    return umatch

}


/// Compute betti numbers and cycle representatives for a , over the rational numbers.
/// 
/// Input is a relation formatted as a `m x n` vec-of-rowvec matrix `M`.  The simplicial complex consists of
/// every subset `S` of `{0,..,n-1}` such that `M[r,S]= [1,..,1]` for some row index `r`
pub fn homology_basis_from_dowker( 
    dowker_simplices: & Vec<Vec<usize>>, 
    maxdim: usize,
    clearing: UseClearing,
    ) 
    -> 
    Vec<Vec<Vec<(Vec<usize>,RingElement)>>>
    {

    // Get the U-match factorization
    // let (boundary_matrix, umatch) = umatch_from_dowker( dowker_simplices, maxdim, clearing );
    let umatch = umatch_from_dowker( dowker_simplices, maxdim, clearing );    
    
    // Get the domain COMB (cf the  paper on umatch factorization)
    let array_comb_domain = umatch.array_comb_domain();
    
    // The set {columns of the domain comb that are not matched upward or downward} 
    // forms a basis for homology
    let mut cycles = vec![Vec::new(); maxdim+1];

    let mut dim;
    for simplex in subsimplices_dim_0_thru_d_iter_ascend( &dowker_simplices, maxdim ) 
                        .filter(|x| umatch.array_matching_ref().lacks_key(x) )
    {        
        dim = simplex.len()-1;        
        cycles[dim].push( array_comb_domain.view_minor_descend( simplex ).collect_vec() )
    }

    for vec in cycles.iter_mut() { vec.shrink_to_fit() }

    return cycles  
}




#[cfg(test)]
mod tests_for_factorization_and_homology {
    use super::*;


    #[test]
    fn test_homology_basis_from_dowker() {
        let dowker_simplices = vec![ vec![0,1], vec![1,2], vec![0,2] ];
        let basis = homology_basis_from_dowker( &dowker_simplices, 1, UseClearing::No );
        let basis = homology_basis_from_dowker( &dowker_simplices, 1, UseClearing::Yes );        
    }     
    
    
    fn test_umatch() {

    }
}    