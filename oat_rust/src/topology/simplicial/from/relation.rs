//! The Dowker complex of a set relation.
//! 
//! Equivalently, a simplicial complex `K` represented by a 
//! subset of simplices `Q` containing every maximal simplex of `K`.
//! 
//! This modules includes data structures for
//! - the boundary matrix in row-major form
//! - the major views of this matrix
//! 
//! 
//! # Example
//! 
//! Here we compute a basis for homology in dimension 1,
//! for the simplicial complex on vertex set `{0,1,2,3}`
//! whose maximal faces are `{0,1,2}, {0,3}, {1,3}, {2,3}`.
//! It "looks" like a tetrahedron where the 3-face and 
//! all but one of the 2-faces have been removed.
//! The coefficient field is the finite field of order 3.
//! 
//! **To run this example on your desktop computer** first checkout the
//! [quick start tutorial in OAT]() for instructions
//! on installing Rust and running a program.  As part of this process,
//! you'll create a new folder that contains a file called `main.rs`.  Inside
//! `main.rs` is some text that reads `fn main{ .. }`.  Delete everything
//! between `{` and `}`, and paste in the following:
//! 
//! ```
//! use oat_rust::algebra::chains::factored::factor_boundary_matrix;    
//! use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
//!             
//! use oat_rust::topology::simplicial::from::relation::BoundaryMatrixDowker;
//! use oat_rust::topology::simplicial::simplices::vector::{subsimplices_dim_d_iter_ascend, subsimplices_dim_d_iter_descend}; 
//! use oat_rust::topology::simplicial::simplices::vector::{subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex}; 
//!     
//! use oat_rust::utilities::order::OrderOperatorAuto;        
//! use oat_rust::utilities::sequences_and_ordinals::SortedVec;
//! 
//! use itertools::Itertools;
//!     
//! // Parameters
//! // ----------
//!     
//! // Define the maximum homology dimensino we want to compute.
//! let max_homology_dimension                  =   2;
//!     
//! // Define the ring operator for the finite field of order 3.
//! // You can use this object to perform arithmetic operations, e.g., 
//! // to add 1 and 1, you can run `ring_operator.add(1,1)`.
//! let ring_operator          =   PrimeOrderFieldOperator::new(3);
//!     
//! // We will build a dowker complex.
//! // A dowker complex is defined by a vertex set V and a family S
//! // of subsets of V.  A subset of V forms a simplex iff it is 
//! // a subset of some element of S.  We refer to the elements 
//! // of S as "dowker simplices".
//!     
//! // Each dowker simplex is represented by a SortedVec of vertices.
//! // We store the list of all such simplices inside a larger vector.
//! let dowker_simplices 
//!     =   vec![    
//!                 vec![0,1,2], 
//!                 vec![0,3], 
//!                 vec![1,3], 
//!                 vec![2,3]  
//!             ]
//!             .into_iter()
//!             .map( |x| SortedVec::new(x).unwrap() )  // we unwrap because `new` can return an error
//!             .collect_vec();
//!     
//! //  Boundary matrix
//! //  ---------------
//!     
//! // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
//! let boundary_matrix = BoundaryMatrixDowker::new( dowker_simplices.clone(), ring_operator.clone() );
//!     
//! //  Simplex iterators
//! //  -----------------
//!     
//! // An iterator that runs over all triangles in the complex, in ascending 
//! // lexicographic order
//! let triangles_ascending_order = subsimplices_dim_d_iter_ascend( &dowker_simplices, 2);
//!     
//! // An iterator that runs over all edges in the complex, in descending 
//! // lexicographic order
//! let triangles_descending_order = subsimplices_dim_d_iter_descend( &dowker_simplices, 2);   
//!     
//! // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
//! // ordered first by dimension (ascending) and second by lexicographic order (descending)
//! let row_indices = boundary_matrix.row_indices_in_descending_order( max_homology_dimension );
//! 
//! // row_indices contains a reference to boundary_matrix, which will cause problems. Instead,
//! // we can construct an iterator that runs over the same sequence of simplices, like so:
//! let row_indices     =   subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex(&dowker_simplices, max_homology_dimension);
//!     
//! //  Homology computation (by matrix factorization)
//! //  ----------------------------------------------
//!     
//! // Factor the boundary matrix
//! let factored    =   factor_boundary_matrix( 
//!                         boundary_matrix, 
//!                         ring_operator, 
//!                         OrderOperatorAuto, 
//!                         row_indices,
//!                     );
//!                 
//! //  Printing results
//! //  ----------------
//!                 
//! // Betti numbers.  For this computation we have to provide a
//! // function that assigns a dimension to each index (i.e. to each simplex)
//! let betti_numbers   =   factored.betti_numbers(|x| x.len() as isize -1 ); 
//! for dim in 0 .. 2 {
//!     println!(
//!             // we'll insert two values into this string
//!             "The betti number in dimension {:?} is {:?}.",
//!             dim,                  // the dimension
//!             betti_numbers         // and the betti number
//!                 .get( & dim )     // this looks up a value in the hashmap
//!                 .unwrap_or( & 0)  // if the hashmap doesn't have the value, then use 0
//!         );            
//! }
//! println!(""); // insert line break
//! 
//! // Cycle representatives for homology
//! println!(
//!         "The following are basis vectors for homology in dimensions 0 through {:?}",
//!         max_homology_dimension,
//!     );
//! for (cycle_number, cycle) in factored.basis_harmonic().enumerate() {
//!     // `cycle` is an iterator.  For convenience, collect the elements of the
//!     // iterator into a Rust vector.
//!     let cycle: Vec<_> = cycle.collect();
//!     println!("Basis vector {:?} = {:?}", cycle_number, cycle);
//! }
//! ```
//! 
//! ### Try changing the coefficient ring
//! 
//! OAT has a number of different [predefined coefficient rings](oat_rust::rings), which you can substitute into
//! the example above in order to calculate homology with different coefficients.  Simply replace the
//! line 
//! ```
//! # use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;       
//! let ring_operator   =   PrimeOrderFieldOperator::new(3);
//! ``` 
//! with one of the `let ring_operator = ...` lines listed under *Predefined Rings*, [here](oat_rust::rings).
//! 
//! ### Try changing the complex
//! 
//! These are the maximal faces in a (minimal) triangulation of the real projective plane:
//! 
//! ```
//! # use oat_rust::utilities::sequences_and_ordinals::SortedVec;
//! # use std::iter::FromIterator;
//! let dowker_simplices    =   Vec::from_iter(
//!     vec![    
//!             vec![1,2,6], vec![1,3,6], vec![1,4,5], 
//!             vec![1,2,5], vec![2,3,5], vec![2,4,6], 
//!             vec![1,3,4], vec![2,3,4], vec![3,5,6], 
//!             vec![4,5,6],  
//!         ]
//!         .into_iter()
//!         .map( |x| SortedVec::new(x) ) 
//! );
//! ```
//! 
//! # About this module
//! 
//! This module uses a data structure [`oat_rust::utilities::sequences_and_ordinals::SortedVec`]
//! for a variety of internal opertaions, e.g. binary search, set union, set difference, evaluation of
//! set containment, etc.

use crate::algebra::matrices::query::{IndicesAndCoefficients, ViewRowAscend, ViewRowDescend, ViewColAscend, ViewColDescend};
use crate::algebra::matrices::operations::umatch::row_major::ParetoShortCircuit;
use crate::algebra::rings::operator_traits::{Semiring, Ring,};
use crate::utilities::iterators::merge::hit::HitMerge;
use crate::utilities::sequences_and_ordinals::{SortedVec, CombinationsReverse};

use std::fmt::Debug;
use std::hash::Hash;
use std::iter::{ Flatten};
use std::marker::PhantomData;

use itertools::{Itertools, Dedup};


use crate::topology::simplicial::boundary::{SimplexBoundaryAscend, SimplexBoundaryDescend};
use crate::topology::simplicial::simplices::vector::{subsimplices_dim_d_iter_descend, subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex};
use crate::utilities::order::{ OrderOperatorAutoReverse, OrderOperatorAuto};        
use crate::algebra::rings::operator_structs::field_prime_order::{PrimeOrderFieldOperator};

// use oat_rust::boundary_matrices::{SimplexBoundaryAscend, SimplexBoundaryDescend};








//  ===================================================================================
//  DOWKER COMPLEX
//  ===================================================================================


/// Represents the Dowker complex of a family of subsets.
/// 
/// We require the elements of the sets to implement `Ord`, as this facilitates numerous
/// counting and look-up operations.
#[derive(Clone)]
pub struct BoundaryMatrixDowker
            < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    dowker_sets:            Vec< SortedVec< Vertex > >,
    ring_operator:          RingOperator,
    phantom_ringelement:    PhantomData< RingElement >,
}


impl < Vertex, RingOperator, RingElement >
    
    BoundaryMatrixDowker
        < Vertex, RingOperator, RingElement >
    where
        Vertex:             Clone + Debug + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    /// Construct a new Dowker complex
    pub fn new( dowker_sets: Vec< SortedVec< Vertex > >, ring_operator: RingOperator ) -> Self {
        BoundaryMatrixDowker{ dowker_sets, ring_operator, phantom_ringelement: PhantomData }
    }

    /// Construct a new Dowker complex from vectors; panics if one of the vectors is not sorted in strictly ascending order.
    pub fn from_vectors( dowker_simplices: Vec< Vec< Vertex > >, ring_operator: RingOperator ) -> Result< Self, Vec< Vertex > > {
        let mut dowker_sets = Vec::with_capacity( dowker_simplices.len() );
        for ( counter, simplex ) in dowker_simplices.into_iter().enumerate() {
            match SortedVec::new(simplex) {
                Err( vec ) => { 
                    println!("Error: attempted to create a Dowker boundary matrix from a sequence of simplices, but simplex number {:?}, which is {:?}, is not sorted.", counter, &vec); 
                    return Err( vec );
                } Ok(set) => { 
                    dowker_sets.push( set ); 
                }
            }            
        }
        // let dowker_sets = dowker_simplices.into_iter().map(|sequence| SortedVec::new(sequence)? ).collect();
        Ok( BoundaryMatrixDowker{ dowker_sets, ring_operator, phantom_ringelement: PhantomData } )
    }    
    
    /// A superset of the family of maximal simplices (represented by hash sets)
    pub fn dowker_sets( &self ) -> & Vec< SortedVec< Vertex > > { & self.dowker_sets }
    
    /// A superset of the family of maximal simplices (represented by vectors whose vertices
    /// are sorted in strictly ascending order)
    pub fn dowker_simplices( &self ) -> & Vec< SortedVec < Vertex > > { 
        & self.dowker_sets
        // self.dowker_sets
        //     .iter()
        //     .map(
        //             |x|
        //             {
        //                 let mut a: Vec<_> = x.iter().cloned().collect();
        //                 a.sort();
        //                 a
        //             }
        //         )
        //     .collect()
    }   

    /// Returns the sequence of row indices of the boundary matrix sorted (first) 
    /// in ascending order of dimension, and (second) in descending lexicographic
    /// order (excluding simplices of dimension `> max_simplex_dimension`)
    pub fn row_indices_in_descending_order( 
                    &self, 
                    max_simplex_dimension: isize 
                ) -> 
            // Vec< Vec< Vertex > > 
            Flatten< std::vec::IntoIter<
                    Dedup< 
                            HitMerge< 
                                    CombinationsReverse< Vertex, & Vec< Vertex > >,
                                    OrderOperatorAutoReverse,
                                >,  
                        >     
                >>                
    {
        return subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex( self.dowker_simplices(), max_simplex_dimension )
        // let iter_keymaj = 
        //     (0..max_simplex_dimension+1)
        //         .map(|x| subsimplices_dim_d_iter_descend(&dowker_simplices, x).unwrap() )
        //         .flatten();
    } 




}


//  INDICES AND COEFFICIENTS
//  ------------------------------------------


impl < Vertex, RingOperator, RingElement >
    
    IndicesAndCoefficients for

    BoundaryMatrixDowker
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type EntryMinor       =   ( Self::RowIndex, Self::Coefficient );    
    type EntryMajor =   ( Self::RowIndex, Self::Coefficient );
    type RowIndex = Vec< Vertex >; 
    type ColIndex = Vec< Vertex >; 
    type Coefficient = RingElement;
}


//  ORACLE MAJOR ASCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    ViewRowAscend for

    BoundaryMatrixDowker
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Debug + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMajorAscend            =   CoboundaryDowkerAscend< Vertex, RingOperator, RingElement >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymaj: Self::RowIndex ) -> Self::ViewMajorAscend {
        CoboundaryDowkerAscend::from_vec_of_dowker_sets( keymaj, & self.dowker_sets, self.ring_operator.clone() ).unwrap()
    }
}

//  ORACLE MAJOR DESCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
ViewRowDescend for

    BoundaryMatrixDowker
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Debug + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMajorDescend           =   CoboundaryDowkerDescend< Vertex, RingOperator, RingElement >;
    type ViewMajorDescendIntoIter   =   Self::ViewMajorDescend;

    fn view_major_descend( &self, keymaj: Self::RowIndex ) -> Self::ViewMajorDescend {

        println!("GOT THROUGH THIS FILE ANE REMOVE THE DEBUG REQUIREMENTS");
        CoboundaryDowkerDescend::from_vec_of_dowker_sets( keymaj, & self.dowker_sets, self.ring_operator.clone() ).unwrap()
    }
}



//  ORACLE MINOR ASCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    ViewColAscend for

    BoundaryMatrixDowker
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMinorAscend            =   SimplexBoundaryAscend< Vertex, RingOperator, RingElement >;
    type ViewMinorAscendIntoIter    =   Self::ViewMinorAscend;

    fn view_minor_ascend( &self, keymaj: Self::RowIndex ) -> Self::ViewMinorAscend {
        SimplexBoundaryAscend::new( keymaj, self.ring_operator.clone() )
    }
}

//  ORACLE MINOR DESCEND
//  ------------------------------------------

impl < Vertex, RingOperator, RingElement >
    
    ViewColDescend for

    BoundaryMatrixDowker
        < Vertex, RingOperator, RingElement >

    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
{
    type ViewMinorDescend           =   SimplexBoundaryDescend< Vertex, RingOperator, RingElement >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;

    fn view_minor_descend( &self, keymaj: Self::RowIndex ) -> Self::ViewMinorDescend {
        SimplexBoundaryDescend::new( keymaj, self.ring_operator.clone() )
    }
}









//  ===================================================================================
//  MAJOR AND MINOR VIEWS
//  ===================================================================================


//  COUBOUNDARY VECTOR DESCENDING
//  -----------------------------------------------------------------------------------


/// Iterates over the terms of the coboundary of a simplex in descending lexicographic order.
/// 
/// The coefficient is calculated as `(-1)^k`, where `xk` is the vertex which appears twice in the template.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::from::relation::CoboundaryDowkerDescend;
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let ring_operator = PrimeOrderFieldOperator::new(3);
/// let simplex = vec![1,3];
/// let dowker_sets = vec![ SortedVec::from_iter( vec![0,1,2,3,4] ) ];
/// // note we have to unwrap, because the constructor returns a Result
/// let coboundary = CoboundaryDowkerDescend::from_vec_of_dowker_sets( simplex, &dowker_sets, ring_operator ).unwrap();
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
        Vertex:             Clone + Debug + Hash + Ord,
{
    /// Generates a [CoboundaryDowkerAscend] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_vec_of_dowker_sets( facet: Vec< Vertex >, dowker_sets: &Vec< SortedVec< Vertex > >, ring_operator: RingOperator ) -> Result< Self, Vec< Vertex > > {

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return Ok( CoboundaryDowkerDescend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator,         // these are arbitrary values                    
            } )                   
        }

        let facet = SortedVec::new( facet ); // place the facet in a safe wrapper; this will panic if the facet isn't strictly sorted
        if let Err(vec) = facet { 
            println!("Error: attempted to compute the coboundary of simplex {:?}, but the vertices of the simplex are not sorted in strictly ascending order", &vec );
            return Err( vec );
        }
        let facet = facet.unwrap(); // we have handled the error, so it's safe to unwrap
        let vertices_to_insert: Vec< Vertex >   =   dowker_sets
                                                        .iter()
                                                        .filter( |x|  x.contains_subset( &facet ) )
                                                        .map(|x| x.vec().iter())
                                                        .kmerge()
                                                        .dedup()
                                                        .filter(|x| ! facet.contains(x) )
                                                        .cloned()
                                                        .collect();


        // //  obtain a list of vertices to insert
        // //  ---------------------------------------------        
        // let facet_vertex_set = SortedVec::from_iter( facet.iter().cloned() );
        // let mut vertices_to_insert = SortedVec::new();

        // // take the union of all dowker super-simplices
        // for dowker_set in dowker_sets
        //                                 .iter() // iterate over dowker siplices
        //                                 .filter( |x| x.is_superset( &facet_vertex_set) ) // only keep those that contain the facet
        //     {
        //         vertices_to_insert.extend( dowker_set.iter().cloned() )
        //     }

        // // take a set difference: vertices_to_insert \ facet
        // for vertex in facet.iter() { vertices_to_insert.remove( vertex ); }

        // // convert to a vector
        // let mut vertices_to_insert: Vec< Vertex > = vertices_to_insert.into_iter().collect();

        // // sort the vector
        // vertices_to_insert.sort();
        

        //  if there are no vertices to insert, then the coboundary is zero
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `next_cofacet_opt` is None
            return Ok( CoboundaryDowkerDescend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),     // these are abitrary values            
                vertices_to_insert:         vec![],                  // these are abitrary values  
                retrieval_locus:            0,                       // these are abitrary values  
                insertion_locus:            0,                       // these are abitrary values  
                ring_operator,           // these are abitrary values                    
            } )                  
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = ring_operator.minus_one_to_power( facet.len() );
        let mut next_cofacet = facet.into_vec();
        let mut insertion_locus = next_cofacet.len();
        let retrieval_locus = vertices_to_insert.len() - 1;        
        let inserted_vertex = vertices_to_insert[retrieval_locus].clone();
        while   ( insertion_locus > 0 ) 
                && 
                ( inserted_vertex < next_cofacet[ insertion_locus - 1 ] ) {
            insertion_locus -= 1;
            coefficient = ring_operator.negate( coefficient );
        }
        next_cofacet.insert( insertion_locus, inserted_vertex );

        Ok( CoboundaryDowkerDescend{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator,             
        } )

    }

    // /// Generates a [CoboundaryDowkerAscend] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< Vertex >, dowker_matrix_csr: Vec< SortedVec< Vertex > >, dowker_matrix_csc: Vec< SortedVec< Vertex > > ) {        
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

        // println!("{:?} -- !!!! DELETE THIS AND ALL DEBUG REQUIREMENTS ON THIS IMPLEMENTATION OF ITER", &self);

        match self.next_cofacet_opt {
            None => { None }
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

                Some( return_value )
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
/// use oat_rust::topology::simplicial::from::relation::CoboundaryDowkerAscend;
/// use oat_rust::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let ring_operator = PrimeOrderFieldOperator::new(3);
/// let simplex = vec![1,3];
/// let dowker_sets = vec![ SortedVec::from_iter( vec![0,1,2,3,4] ) ];
/// 
/// // note we have to unwrap, because the constructor returns a Result
/// let coboundary = CoboundaryDowkerAscend::from_vec_of_dowker_sets( simplex, &dowker_sets, ring_operator ).unwrap();
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
        Vertex:             Clone + Debug + Hash + Ord,
{
    /// Generates a [CoboundaryDowkerAscend] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_vec_of_dowker_sets( facet: Vec< Vertex >, dowker_sets: &Vec< SortedVec< Vertex > >, ring_operator: RingOperator ) -> Result< Self, Vec< Vertex > > {

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return Ok( CoboundaryDowkerAscend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator,         // these are arbitrary values                    
            } )                
        }


        let facet = SortedVec::new( facet ); // place the facet in a safe wrapper; this will panic if the facet isn't strictly sorted
        if let Err(vec) = facet { 
            println!("Error: attempted to compute the coboundary of simplex {:?}, but the vertices of the simplex are not sorted in strictly ascending order", &vec );
            return Err( vec );
        }
        let facet = facet.unwrap(); // we have handled the error, so it's safe to unwrap
        let vertices_to_insert: Vec< Vertex >   =   dowker_sets
                                                        .iter()
                                                        .filter( |x|  x.contains_subset( &facet ) )
                                                        .map(|x| x.vec().iter())
                                                        .kmerge()
                                                        .dedup()
                                                        .filter(|x| ! facet.contains(x) )
                                                        .cloned()
                                                        .collect();
        // //  obtain a list of vertices to insert
        // //  ---------------------------------------------        
        // let facet_vertex_set = SortedVec::from_iter( facet.iter().cloned() );
        // let mut vertices_to_insert = SortedVec::new();

        // // take the union of all dowker super-simplices
        // for dowker_simplex in dowker_sets
        //                                 .iter() // iterate over dowker siplices
        //                                 .filter( |x| x.is_superset( &facet_vertex_set) ) // only keep those that contain the facet
        //     {
        //         vertices_to_insert.extend( dowker_simplex.iter().cloned() )
        //     }

        // // take a set difference: vertices_to_insert \ facet
        // for vertex in facet.iter() { vertices_to_insert.remove( vertex ); }

        // // convert to a vector
        // let mut vertices_to_insert: Vec< Vertex > = vertices_to_insert.into_iter().collect();

        // // sort the vector
        // vertices_to_insert.sort();
        

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `vertices_to_insert` is empty
            return Ok( CoboundaryDowkerAscend{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),
                vertices_to_insert:         vec![],     
                retrieval_locus:            0,          
                insertion_locus:            0,          
                ring_operator,                
            } )                   
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = RingOperator::one();
        let mut next_cofacet = facet.into_vec();
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

        Ok( CoboundaryDowkerAscend{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator,             
        } )

    }

    // /// Generates a [CoboundaryDowkerAscend] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< Vertex >, dowker_matrix_csr: Vec< SortedVec< Vertex > >, dowker_matrix_csc: Vec< SortedVec< Vertex > > ) {        
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
            None => { None }
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

                Some( return_value )
            }
        }

    }
}    



impl < Vertex, RingOperator, RingElement >

    ParetoShortCircuit< (Vec<Vertex>, RingElement) > for

    CoboundaryDowkerAscend
        < Vertex, RingOperator, RingElement >

    where 
        RingOperator:       Semiring< RingElement > + Ring< RingElement >,
        RingElement:        Clone,
        Vertex:             Ord + Clone,   
{
    fn pareto_short_circuit(& self) -> Option< (Vec<Vertex>, RingElement) > {
        None
    }
}






//  ===================================================================================
//  OPTIMIZATION
//  ===================================================================================



pub fn optimize() { println!("NOTE THAT WE REALLY NEED A TOOL TO IMINIMIZE A CYCLE WITHIN A FIXED HOMOLOGY CLASS"); }


//  ===================================================================================
//  TESTING
//  ===================================================================================




use crate::algebra::matrices::debug::{verify_view_minor_descend_is_sorted_strictly, verify_view_major_ascend_is_sorted_strictly, verify_viewmajorascend_compatible_with_viewminordescend, verify_viewmajorascend_compatible_with_viewmajordescend, verify_viewminorascend_compatible_with_viewminordescend};


/// A tool to check that 
pub fn dowker_boundary_diagnostic<T: Clone + Hash + Debug + Ord>( dowker_simplices_vec: Vec<Vec<T>>, maxdim: isize, ) -> Result<(), Vec< T > > {

    // first ensure simplices are sorted
    let dowker_sets: Result< Vec<SortedVec<T>>, Vec<T> > = dowker_simplices_vec.into_iter().map(|x| SortedVec::new(x) ).collect();
    let dowker_sets 
        =   dowker_sets.map_err(
                |vector|
                {
                    println!("Error: attempted to convert a Vec< Vec<T> > to a Vec< SortedVec< T > >, but one of the inner vectors wasn't sorted: {:?}", & vector );
                    vector         
                }
            )?;

    let ring_operator                =   PrimeOrderFieldOperator::new(47);
    // define the boundary matrix
    let boundary_matrix = BoundaryMatrixDowker::new( dowker_sets.clone(), ring_operator );        

    // define an iterator to run over all simplices; simplices are ordered first by dimension (ascending), then by lexicographic order (descending)
    let iter_keymaj = boundary_matrix.row_indices_in_descending_order( maxdim );     


    let keys = iter_keymaj.clone().collect_vec();
    println!("keys: {:?}", keys );

    for dim in 0 .. maxdim+1 {
        let v: Vec<_> = subsimplices_dim_d_iter_descend( &dowker_sets, dim ).unwrap().collect();
        println!("KEYS OF DIMENSION {:?} === {:?}", dim, v );
    }

    // check that the input is strictly sorted, first by dimension 
    let mut old_opt: Option< Vec< T > > = None;
    for new in iter_keymaj.clone() {
        if let Some( old ) = old_opt {
            if ( old.len() == new.len() ) && ( new >= old ) {
                panic!("The keymaj iterator is not strictly sorted, first by dimension (ascending) then lexicographically (descending)")
            }            
        }
        old_opt = Some( new );
    }
    


    // if verbose {
    //     println!("major views");
    //     for keymaj in iter_keymaj.clone() {
    //         println!("row {:?} = ", boundary_matrix.view_major_ascend(keymaj).collect_vec() );
    //     }
    //     println!("minor views");
    //     for keymaj in iter_keymaj.clone() {
    //         println!("column {:?} = ", boundary_matrix.view_minor_descend(keymaj).collect_vec() );
    //     }        
    // }    

    let print_fn = |x: Vec< Vec<T> >| println!("{:?}", x );

    // check that major ascending views are sorted in strictly ascending order / minor descending views are sorted in strictly descending order
    verify_view_major_ascend_is_sorted_strictly(  & boundary_matrix, iter_keymaj.clone(), OrderOperatorAuto, print_fn );
    verify_view_minor_descend_is_sorted_strictly( & boundary_matrix, iter_keymaj.clone(), OrderOperatorAutoReverse::new(), print_fn );    // NOTE we use GT for this one, because order is reversed

    // check that ascending major views agree with descending MINOR views
    verify_viewmajorascend_compatible_with_viewminordescend(
            & boundary_matrix,
            iter_keymaj.clone(),
            iter_keymaj.clone(),
        );

    // check that ascending major views agree with descending MAJOR views            
    verify_viewmajorascend_compatible_with_viewmajordescend(
            & boundary_matrix,
            iter_keymaj.clone(),
        );

    // check that ascending MINOR views agree with descending MINOR views            
    verify_viewminorascend_compatible_with_viewminordescend(
        & boundary_matrix,
        iter_keymaj.clone(),            
    );

    Ok(())      
}




#[cfg(test)]
mod tests {
    

use crate::utilities::random::rand_sequences;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

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

        let _ = dowker_boundary_diagnostic( dowker_simplices_vec, 2 );
    }



    #[test]
    fn test_dowker_boundary_big(){
        
        for _ in 0..10 {
            
            let dowker_simplices_vec    =   rand_sequences(10, (0..7).step_by(2), 0.1 );
            let _verbose = false;
            let _ = dowker_boundary_diagnostic( dowker_simplices_vec, 2 );
        }

    }    



}   



#[cfg(test)]
mod docstring_tests {
    use crate::{topology::simplicial::simplices::vector::{subsimplices_dim_d_iter_ascend, subsimplices_dim_d_iter_descend, subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex}, utilities::order::OrderOperatorAuto};

    



    #[test]
    fn docstring_test_dowker_homology() {
        use itertools::Itertools;
            
        use crate::topology::simplicial::from::relation::BoundaryMatrixDowker;
            
        use crate::algebra::chains::factored::factor_boundary_matrix;       
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;        
      
        use crate::utilities::sequences_and_ordinals::SortedVec;
            
        // Parameters
        // ----------
            
        // Define the maximum homology dimensino we want to compute.
        let max_homology_dimension                  =   2;
            
        // Define the ring operator for the finite field of order 3.
        // You can use this object to perform arithmetic operations, e.g., 
        // to add 1 and 1, you can run `ring_operator.add(1,1)`.
        let ring_operator          =   PrimeOrderFieldOperator::new(3);
            
        // We will build a dowker complex.
        // A dowker complex is defined by a vertex set V and a family S
        // of subsets of V.  A subset of V forms a simplex iff it is 
        // a subset of some element of S.  We refer to the elements 
        // of S as "dowker simplices".
            
        // Each dowker simplex is represented by a SortedVec of vertices.
        // We store the list of all such simplices inside a larger vector.
        let dowker_simplices 
            =   vec![    
                        vec![0,1,2], 
                        vec![0,3], 
                        vec![1,3], 
                        vec![2,3]  
                    ]
                    .into_iter()
                    .map( |x| SortedVec::new(x).unwrap() ) 
                    .collect_vec();
            
        //  Boundary matrix
        //  ---------------
            
        // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
        let boundary_matrix = BoundaryMatrixDowker::new( dowker_simplices.clone(), ring_operator );
            
        //  Simplex iterators
        //  -----------------
            
        // An iterator that runs over all triangles in the complex, in ascending 
        // lexicographic order
        let _triangles_ascending_order = subsimplices_dim_d_iter_ascend( &dowker_simplices, 2);
            
        // An iterator that runs over all edges in the complex, in descending 
        // lexicographic order
        let _triangles_descending_order = subsimplices_dim_d_iter_descend( &dowker_simplices, 2);   
            
        // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
        // ordered first by dimension (ascending) and second by lexicographic order (descending)
        let _row_indices = boundary_matrix.row_indices_in_descending_order( max_homology_dimension );

        // row_indices contains a reference to boundary_matrix, which will cause problems. Instead,
        // we can construct an iterator that runs over the same sequence of simplices, like so:
        let row_indices     =   subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex(&dowker_simplices, max_homology_dimension);
            
        //  Homology computation (by matrix factorization)
        //  ----------------------------------------------
            
        // Factor the boundary matrix
        let factored    =   factor_boundary_matrix( 
                                boundary_matrix, 
                                ring_operator, 
                                OrderOperatorAuto, 
                                row_indices,
                            );
                        
        //  Printing results
        //  ----------------
                        
        // Betti numbers.  For this computation we have to provide a
        // function that assigns a dimension to each index (i.e. to each simplex)
        let betti_numbers   =   factored.betti_numbers(|x| x.len() as isize -1 ); 
        for dim in 0 .. 2 {
            println!(
                    // we'll insert two values into this string
                    "The betti number in dimension {:?} is {:?}.",
                    dim,                  // the dimension
                    betti_numbers         // and the betti number
                        .get( & dim )     // this looks up a value in the hashmap
                        .unwrap_or( & 0)  // if the hashmap doesn't have the value, then use 0
                );            
        }
        println!(); // insert line break
        
        // Cycle representatives for homology
        println!(
                "The following are basis vectors for homology in dimensions 0 through {:?}",
                max_homology_dimension,
            );
        for (cycle_number, cycle) in factored.basis_harmonic().enumerate() {
            // `cycle` is an iterator.  For convenience, collect the elements of the
            // iterator into a Rust vector.
            let cycle: Vec<_> = cycle.collect();
            println!("Basis vector {:?} = {:?}", cycle_number, cycle);
        }
    }  


    fn test() {
        use crate::topology::simplicial::from::relation::CoboundaryDowkerAscend;
        use crate::algebra::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
        use crate::utilities::sequences_and_ordinals::SortedVec;

        let ring_operator = PrimeOrderFieldOperator::new(3);
        let simplex = vec![1,3];
        let dowker_sets = vec![ SortedVec::from_iter( vec![0,1,2,3,4] ) ];
        let coboundary = CoboundaryDowkerAscend::from_vec_of_dowker_sets( simplex, &dowker_sets, ring_operator ).unwrap();

        itertools::assert_equal( coboundary, vec![ (vec![0,1,3], 1), (vec![1,2,3], 2), (vec![1,3,4], 1) ]);        
    }    

}





