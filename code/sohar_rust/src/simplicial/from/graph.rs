//! The Clique complex (aka Flag complex) of a simple combinatorial graph

use itertools::Itertools;

use oat_rust::matrices::matrix_oracle_traits::{IndicesAndCoefficients, OracleMajorAscend, OracleMajorDescend, OracleMinorAscend, OracleMinorDescend};
use oat_rust::rings::operator_traits::{Semiring, Ring};
use oat_rust::utilities::iterators::general::IntersectOrderedIteratorsUnsafe;

use std::collections::HashSet;
use std::hash::Hash;
use std::iter::FromIterator;
use std::marker::PhantomData;

use crate::simplicial_complex::boundary::{SimplexBoundaryAscend, SimplexBoundaryDescend};

// use super::boundary_matrices::{SimplexBoundaryAscend, SimplexBoundaryDescend};




// pub struct FlagComplexSimplexIterAscend{ 
//     proximity_matrix:   & Vec<HashSet<usize>>,
//     cardinality:        usize,
//     current_face:       Vec<usize>,
// }

// impl    Iterator for 
//         FlagComplexSimplexIterAscend {

//     type Item = Vec< usize >;

//     fn next(&mut self) -> Option<Self::Item> {
//         let last_index = self.cardinality - 1;
//         let last_vertex     =   self.current_face[ last_index ];
//         for next_option in last_vertex + 1 .. self.proximity_matrix.len() {
//             if self.current_face.iter().map( |x| self.proximity_matrix.contains(x) ).take( last_index ).all() {
//                 self.current_face[ last_index ] = next_option;
//                 return Some(self.current_face.clone());
//             }
//         }
        
//     }
// }





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
/// let coboundary = CoboundaryDowkerDescend::from_proximity_matrix( simplex, &dowker_simplices, ring_operator );
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


impl < RingOperator, RingElement >

CoboundaryDowkerDescend
        < usize, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        // Vertex:             Clone + Hash + Ord,
{
    /// Generates a [CoboundaryDowkerAscend] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_proximity_matrix( facet: Vec< usize >, proximity_matrix: &Vec< Vec< usize > >, ring_operator: RingOperator ) -> Self {

        println!("NOTE: we assume that the proximity graph is just a sparsity pattern and DOES NOT CONTAIN the diagonal");
        println!("NOTE: THIS IS INEFFICIENT; CAN REPLACE vertices-to-insert with an iterator");
        println!("NOTE: consider merging the constructors for the ascending and descending ... or maybe not that; consider merging the constructor with the actual Dowker version ");

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
        
        let iters = facet.iter().map(|x| proximity_matrix[*x].iter().cloned() );
        let vertices_to_insert = IntersectOrderedIteratorsUnsafe::new( iters ).collect_vec();
        

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
    // fn from_csr_and_csc_matrices( facet: Vec< usize >, dowker_matrix_csr: Vec< HashSet< usize > >, dowker_matrix_csc: Vec< HashSet< usize > > ) {        
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
/// let coboundary = CoboundaryDowkerAscend::from_proximity_matrix( simplex, &dowker_simplices, ring_operator );
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


impl < RingOperator, RingElement >

    CoboundaryDowkerAscend
        < usize, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    /// Generates a [CoboundaryDowkerAscend] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_proximity_matrix( facet: Vec< usize >, proximity_matrix: &Vec< Vec< usize > >, ring_operator: RingOperator ) -> Self {

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
        let iters = facet.iter().map(|x| proximity_matrix[*x].iter().cloned() );
        let vertices_to_insert = IntersectOrderedIteratorsUnsafe::new( iters ).collect_vec();
        

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
    // fn from_csr_and_csc_matrices( facet: Vec< usize >, dowker_matrix_csr: Vec< HashSet< usize > >, dowker_matrix_csc: Vec< HashSet< usize > > ) {        
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












//  ---------------------------------------------------------------------------
//  DOWKER BOUNDARY MATRIX
//  ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct FlagComplexBoundaryMatrixRowMajor
            < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    proximity_matrix:       Vec< Vec< Vertex > >,
    ring_operator:          RingOperator,
    phantom_ringelement:    PhantomData< RingElement >,
}


impl < Vertex, RingOperator, RingElement >
    
    FlagComplexBoundaryMatrixRowMajor
        < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    pub fn new( proximity_matrix: Vec< Vec< Vertex > >, ring_operator: RingOperator ) -> Self {
        FlagComplexBoundaryMatrixRowMajor{ proximity_matrix, ring_operator, phantom_ringelement: PhantomData }
    }
}


//  INDICES AND COEFFICIENTS
//  ------------------------------------------


impl < Vertex, RingOperator, RingElement >
    
    IndicesAndCoefficients for

    FlagComplexBoundaryMatrixRowMajor
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

impl < RingOperator, RingElement >
    
    OracleMajorAscend for

    FlagComplexBoundaryMatrixRowMajor
        < usize, RingOperator, RingElement >

    where
        // Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMajorAscend            =   CoboundaryDowkerAscend< usize, RingOperator, RingElement >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;
    type ViewMajorAscendEntry       =   ( Self::KeyMaj, Self::SnzVal );

    fn view_major_ascend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMajorAscend {
        CoboundaryDowkerAscend::from_proximity_matrix( keymaj, & self.proximity_matrix, self.ring_operator.clone() )
    }
}

//  ORACLE MAJOR DESCEND
//  ------------------------------------------

impl < RingOperator, RingElement >
    
OracleMajorDescend for

    FlagComplexBoundaryMatrixRowMajor
        < usize, RingOperator, RingElement >

    where
        // Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMajorDescend           =   CoboundaryDowkerDescend< usize, RingOperator, RingElement >;
    type ViewMajorDescendIntoIter   =   Self::ViewMajorDescend;
    type ViewMajorDescendEntry      =   ( Self::KeyMaj, Self::SnzVal );

    fn view_major_descend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMajorDescend {

        println!("GOT THROUGH THIS FILE ANE REMOVE THE DEBUG REQUIREMENTS");
        CoboundaryDowkerDescend::from_proximity_matrix( keymaj, & self.proximity_matrix, self.ring_operator.clone() )
    }
}



//  ORACLE MINOR ASCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    OracleMinorAscend for

    FlagComplexBoundaryMatrixRowMajor
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

    FlagComplexBoundaryMatrixRowMajor
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
    use crate::simplicial_complex::from::graph::FlagComplexBoundaryMatrixRowMajor;
    use crate::simplicial_complex::simplices::{subsimplices_dim_0_thru_d_iter_descend, subsimplices_dim_d_iter_ascend};    
    use oat_rust::matrices::matrix_types::oracle_ref::OracleRef;
    use oat_rust::matrices::matrix_oracle_traits::OracleMinorDescend;
    use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor};
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
        let boundary_matrix = FlagComplexBoundaryMatrixRowMajor::new( dowker_simplices, ring_operator.clone() );

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