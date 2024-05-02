/*!
Data structures for filtered clique complexes.

# Clique chain complex oracles
Mathematically, a filtered clique complex is determined by two pieces of data:
* a **dissimilarity matrix**
    * we represent this by a symmetric matrix *S* that has been flattened into a vector *v*
* a **threshold**  that determines where we stop growing the filtration
    * we represent this as a real number *t*
The boundary matrices of this chain complex oracle are indexed by `SimplexFiltered` objects which are
(essentially) pairs of form `(simplex, filtration_value)`.
```
/*
Example:
    - construct a dissimilarity matrix
    - construct the associated clique complex
    - access a row + column of the boundary matrix
    - compute the dimension 1 barcode
///////////////////////////////////////////////////////////////////////////////////////////////////
// Setup ring data
	let ring_operator = RingMetadata{
		ringspec: RingSpec::Modulus(field),
		identity_additive: 0,
		identity_multiplicative: 1,
	};
	//println!("Ring: {:?}", ring_operator.identity_additive);
///////////////////////////////////////////////////////////////////////////////////////////////////
// Compute persistence of clique complex
/*
	if maxdis == OrderedFloat(0.0) { maxdis = radius; }
	let chx = CliqueComplex {
		dissimilarity_matrix: dismat_sample,
		dissimilarity_value_max: maxdis,
		safe_homology_degrees_to_build_boundaries: (1..(dim+1)).collect(),
		major_dimension: MajorDimension::Row,
		ring_operator: ringdata,
		simplex_count: Vec::new()
	};
*/
///////////////////////////////////////////////////////////////////////////////////////////////////
*/
```
*/



use std::ops::Bound;
use std::{ops::Neg, marker::PhantomData};
use std::convert::TryInto;
use std::hash::Hash;
use std::fmt::Debug;
use itertools::Itertools;
use ordered_float::OrderedFloat;

// use crate::matrix::{SmOracle, RingMetadata, MajorDimension};
// use crate::chx::{ChainComplex, ChxTransformKind};

use std::cmp::Ordering;

use oat_rust::matrices::operations::umatch::row_major::ParetoShortCircuit;
use oat_rust::rings::operator_traits::DivisionRing;
use oat_rust::utilities::iterators::general::{PeekUnqualified, find_min};
use oat_rust::utilities::partial_order::StrictlyLess;
use oat_rust::{rings::operator_traits::{Semiring, Ring}, matrices::matrix_oracle_traits::{IndicesAndCoefficients, OracleMajorAscend, OracleMinorDescend}};

use crate::simplicial::simplices::SimplexFiltered;

use super::user_interface::CliqueParams;

// The specific type we use to store vertices
type Vertex = u16;




/// Boundary matrix represented as an [`SmOracle`] trait object.
pub struct CliqueBoundaryMatrix<FilVal, SnzVal, RingOperator>
    where 
        SnzVal:         Clone,
        RingOperator:   Semiring< SnzVal > + Ring< SnzVal >
{
    pub ring_operator: RingOperator, // ring meta data
    /// The "maxdis" value represents the maximum of filtration value
    pub maxdis: FilVal,
    /// A vector representing the dissimilarity matrix by listing all its rows
	pub dismat: Vec<Vec<FilVal>>,
    /// A vector representing the neighborhood within "maxdis" of each vertex
    pub cutoff_matrix: Vec<Vec<Vertex>>,
    pub phantomsnzval: PhantomData< SnzVal >,
}

/// A function that produces "cutoff_matrix"
///
/// # Parameters
/// - `dismat`: A dissimilarity/distance matrix
/// - `maxdis`: The radius of neighborhood
pub fn get_cutoff_matrix<FilVal: PartialOrd + Debug>(dismat: &Vec<Vec<FilVal>>, maxdis: FilVal) -> Vec<Vec<Vertex>> {
    let mut output: Vec<Vec<Vertex>> = Vec::new();
    for row in dismat.iter() {
        let mut index_vec = Vec::new();
        for jj in 0..row.len() {
            if row[jj] <= maxdis {
                index_vec.push(jj.try_into().unwrap());
            }
        }
        output.push(index_vec);
    }
    return output;
}

/// Methods of CliqueBoundaryMatrix struct
impl <FilVal: Clone + Debug + PartialOrd + Ord, SnzVal: Clone, RingOperator: Semiring< SnzVal > + Ring< SnzVal > > 

    CliqueBoundaryMatrix
        < FilVal, SnzVal, RingOperator > {

    /// Output the diameter of a simplex
    ///
    /// # Parameters
    /// -`vertices`: vertices of the simplex
    /// # Returns
    /// Return None if some edge of simplex is out of cutoff value; return Some(diam) otherwise
    fn diam(&self, vertices: &Vec<Vertex>) -> Option<FilVal> {
        if vertices.is_empty() { return None };

        let mut a = usize::from( vertices[0] );     // first vertex
        let mut b = a.clone();                      // second vertex

        let mut diam_bound = &self.dismat[ a ][ b ];   // lower bound on diameter
        let mut diam = diam_bound.clone(); 
        for ii in 0..vertices.len() {
            a = usize::from( vertices[ii] ); 
            for jj in 0..ii {
                b = usize::from( vertices[jj] ); 
                diam_bound = &self.dismat[ a ][ b ];
                if diam_bound > &diam { diam = diam_bound.clone(); }
            }
        }
        return Some(diam);
    }    

}

/// Methods of CliqueBoundaryMatrix struct
impl < SnzVal: Clone, RingOperator: Semiring< SnzVal > + Ring< SnzVal > > 

    CliqueBoundaryMatrix
        < OrderedFloat<f64>, SnzVal, RingOperator > {


    pub fn new( 
                params: & CliqueParams<RingOperator, SnzVal> 
            ) 
            -> 
            Self
        where 
            RingOperator:   Clone + Semiring< SnzVal > + Ring< SnzVal > + DivisionRing< SnzVal >,
            SnzVal:         Clone + Debug,        
        {
            let cutoff_matrix = get_cutoff_matrix( params.dismat_ref(), params.maxdis() );
            CliqueBoundaryMatrix
                {
                    ring_operator: params.ring_operator(),
                    /// The "maxdis" value represents the maximum of filtration value
                    maxdis: params.maxdis(),
                    /// A vector representing the dissimilarity matrix by listing all its rows
                    dismat: params.dismat_ref().clone(),
                    /// A vector representing the neighborhood within "maxdis" of each vertex
                    cutoff_matrix,
                    phantomsnzval: PhantomData,
                }
    }
}


/// An iterator that iterates all cofacets of a given simplex
struct CofacetIter
        <'a, FilVal: Clone, SnzVal: Clone, RingOperator: Semiring< SnzVal > + Ring< SnzVal > >  {
	pub clique: &'a CliqueBoundaryMatrix<FilVal, SnzVal, RingOperator >,
	pub next_cofacet_vertices: Vec<Vertex>,
    pub simplex_filvalue: FilVal,
	pub insertion_location: usize,
    pub candidate_location: usize,
    pub first_vert: Vertex,
	pub coeff: SnzVal,
    pub ring_operator: RingOperator,
}

/// implement standard methods of Iterator for CofacetIter struct
impl<'a, FilVal, SnzVal, RingOperator: Semiring< SnzVal > + Ring< SnzVal >, > 

    Iterator for 
    
    CofacetIter
        <'a, FilVal, SnzVal, RingOperator, > 
        
    where
        FilVal: Clone + Debug + PartialOrd,
        SnzVal: Clone
{
	type Item = (SimplexFiltered<FilVal>, SnzVal);

    fn next( &mut self ) -> Option<Self::Item> {

        let candidacies = &self.clique.cutoff_matrix[usize::from(self.first_vert)];
        let vertices = &mut self.next_cofacet_vertices;

        loop {
            if self.candidate_location >= candidacies.len(){
                return None;
            }
            let new_vertex = candidacies[self.candidate_location];
            let mut flag = false;
            vertices[self.insertion_location] = new_vertex;
            let mut new_filval = self.simplex_filvalue.clone();
            for vert in vertices.iter() {
                if self.clique.dismat[usize::from(new_vertex)][usize::from(*vert)] > self.clique.maxdis {
                    flag = true;
                    break;
                } else if self.clique.dismat[usize::from(new_vertex)][usize::from(*vert)] > new_filval {
                    new_filval = self.clique.dismat[usize::from(new_vertex)][usize::from(*vert)].clone();
                }
            }
            if flag { self.candidate_location += 1; continue; }

            while self.insertion_location < vertices.len()-1 &&
                                new_vertex >= vertices[self.insertion_location+1] {
                if new_vertex == vertices[self.insertion_location+1] {
                    flag = true;
                    break;
                }
                vertices[self.insertion_location] = vertices[self.insertion_location+1];
                self.insertion_location += 1;
                self.coeff = self.ring_operator.negate(self.coeff.clone());
            }
            if flag { self.candidate_location += 1; continue; }

            if self.insertion_location >= vertices.len() {
                return None;
            }
            vertices[self.insertion_location] = new_vertex;
        //println!("new_filval={:?}, simplex_filvalue={:?}", new_filval, self.simplex_filvalue);
            let simp = SimplexFiltered{ vertices: vertices.clone(), filvalue: new_filval };
            self.candidate_location += 1;
            return Some((simp, self.coeff.clone()));
        }
        //println!("candidate_location={:?}", self.candidate_location);
	}
}


/// A iterator that iterates all facets of a given simplex
struct FacetIter
        <'a, FilVal: Clone + Debug, SnzVal: Clone, RingOperator: Semiring< SnzVal > + Ring< SnzVal > >  {
    pub clique: &'a CliqueBoundaryMatrix<FilVal, SnzVal, RingOperator>,
	pub simp: SimplexFiltered<FilVal>,
	pub removal_location: usize,
    pub ring_operator: RingOperator,
}

/// implement standard methods of Iterator for FacetIter struct
impl<FilVal, SnzVal, RingOperator> 

    Iterator for 
    
    FacetIter<'_, FilVal, SnzVal, RingOperator> 
    
    where
        FilVal: Clone + PartialOrd + Debug + Ord,
        SnzVal: Clone,
        RingOperator: Semiring< SnzVal > + Ring< SnzVal >,
{
	type Item = (SimplexFiltered<FilVal>, SnzVal);

	fn next( &mut self ) -> Option<Self::Item> {
        if self.simp.vertices.len() == 1 { return None } // the boundary of a 0-simplex is empty
        if self.removal_location == self.simp.vertices.len() { return None; }
        
        let mut simplex = self.simp.clone();
        simplex.vertices.remove(self.removal_location);
        simplex.vertices.shrink_to_fit();
        simplex.filvalue = self.clique.diam(&simplex.vertices).unwrap();
        let coeff = self.ring_operator.minus_one_to_power( self.removal_location );
        self.removal_location += 1;
        return Some((simplex, coeff));

        // if let Some(diam) = self.clique.diam(&simplex.vertices) {
        //     simplex.filvalue = diam;
        //     self.coeff = self.ring_operator.minus_one_to_power( self.removal_location );
        //     self.removal_location += 1;
        //     return Some((simplex, self.coeff.clone()));
        // } else {
        //     return None;
        // }
    }
}


fn get_coboundary_as_vec< FilVal, SnzVal, RingOperator >( 
            matrix: & CliqueBoundaryMatrix< FilVal, SnzVal, RingOperator >, 
            keymaj: SimplexFiltered<FilVal> 
        ) 
        -> Vec< ( SimplexFiltered<FilVal>, SnzVal ) >
        
    where
        FilVal:             Clone + Debug + PartialOrd + Ord, 
        SnzVal:             Clone, 
        RingOperator:       Clone + Semiring< SnzVal > + Ring< SnzVal >   {

    let mut vertices = keymaj.vertices.clone();
    let first_vert = vertices[0];
    vertices.insert(0,0);
    vertices.shrink_to_fit();
    let iter = 
    CofacetIter{
        clique: &matrix,
        next_cofacet_vertices: vertices,
        simplex_filvalue: keymaj.filvalue.clone(),
        insertion_location: 0,
        candidate_location: 0,
        coeff: RingOperator::one(),
        first_vert,
        ring_operator: matrix.ring_operator.clone(),
    };   
    let mut vec = iter.collect_vec();
    vec.shrink_to_fit();
    vec.sort_by(|x,y| x.0.cmp(&y.0) );
    return vec
}


//  ---------------------------------------------------------------------------
//  CLIQUE BOUNDARY MATRIX
//  ---------------------------------------------------------------------------


//  INDICES AND COEFFICIENTS
//  ------------------------------------------


impl < FilVal, SnzVal, RingOperator >
    
    IndicesAndCoefficients for    

    CliqueBoundaryMatrix
        < FilVal, SnzVal, RingOperator >

    where
        FilVal:             Clone + Debug + PartialOrd + Ord, 
        SnzVal:             Clone, 
        RingOperator:       Semiring< SnzVal > + Ring< SnzVal >      
{
    type KeyMaj = SimplexFiltered<FilVal>; type KeyMin = SimplexFiltered<FilVal>; type SnzVal = SnzVal;
}


//  ORACLE MAJOR ASCEND
//  ------------------------------------------


impl < 'a, FilVal, SnzVal, RingOperator >
    
OracleMajorAscend for 

    &'a CliqueBoundaryMatrix
        < FilVal, SnzVal, RingOperator >

    where
        FilVal:             Clone + Debug + PartialOrd + Ord, 
        SnzVal:             Clone, 
        RingOperator:       Clone + Semiring< SnzVal > + Ring< SnzVal >,
{
    type ViewMajorAscend            =   LazyOrderedCoboundary<  'a, FilVal, SnzVal, RingOperator >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;
    type ViewMajorAscendEntry       =   ( Self::KeyMaj, Self::SnzVal );

    fn view_major_ascend( &self, keymaj: Self::KeyMaj ) -> Self::ViewMajorAscend {

        let first_v = keymaj.vertices[0];
        let mut shortcircuit_vertex = None;        
        for outer_v in self.cutoff_matrix[first_v as usize ].iter().cloned() {
            if outer_v >= first_v { break }
            let disvec = & self.dismat[ outer_v as usize ];
            let diam = keymaj.fil();
            let mut diam_ok = true;
            for inner_v in & keymaj.vertices {
                if disvec[ inner_v.clone() as usize ] > diam { diam_ok = false; break }                                
            }
            if diam_ok { shortcircuit_vertex=Some( outer_v ); break }
        }

        if let Some( outer_v ) = shortcircuit_vertex {
            // IN THIS CASE WE ONLY BUILD THE FIRST ENTRY
            let mut first_cofacet = Vec::with_capacity( keymaj.vertices.len()+1 );
            first_cofacet.push( outer_v );
            for inner_v in keymaj.vertices.iter().cloned() { first_cofacet.push(inner_v) }
            let cofacet = SimplexFiltered{ vertices: first_cofacet, filvalue: keymaj.fil() };
            let first_entry = ( cofacet, RingOperator::one() );
            LazyOrderedCoboundary{
                facet:              keymaj,
                boundary_matrix:    self,
                coboundary:         vec![ first_entry ],
                next_index:         0,
                built_all:          false,
            }
        } else {
            let coboundary = get_coboundary_as_vec(self, keymaj.clone() );
            LazyOrderedCoboundary{
                facet:              keymaj,
                boundary_matrix:    self,
                coboundary:         coboundary,
                next_index:         0,
                built_all:          true,
            }            
        }

    }
}


#[derive(Clone)]
pub struct LazyOrderedCoboundary
        <'a, FilVal, SnzVal, RingOperator >
    where 
        FilVal:         Clone + Debug,
        SnzVal:         Clone,
        RingOperator:   Semiring< SnzVal > + Ring< SnzVal >,
{
    facet:              SimplexFiltered<FilVal>,
    boundary_matrix:    &'a CliqueBoundaryMatrix < FilVal, SnzVal, RingOperator >,
    coboundary:         Vec< (SimplexFiltered<FilVal>, SnzVal) >,
    next_index:         usize,
    built_all:          bool
}

impl < 'a, FilVal, SnzVal, RingOperator >
    
    Iterator for
    
    LazyOrderedCoboundary
        < 'a, FilVal, SnzVal, RingOperator >

    where 
        FilVal:         Clone + Debug + Ord,
        SnzVal:         Clone,
        RingOperator:   Clone + Semiring< SnzVal > + Ring< SnzVal >,      
{
    type Item = (SimplexFiltered<FilVal>, SnzVal);

    fn next(&mut self) -> Option<Self::Item> {

        while self.next_index >= self.coboundary.len() {
            if self.built_all { return None }
            else{
                let mut v = get_coboundary_as_vec( self.boundary_matrix, self.facet.clone() );
                v.sort_by(|x,y| x.0.cmp(&y.0)); 
                self.coboundary = v;
                self.built_all = true;                
            }
        }

        let return_value = self.coboundary[ self.next_index ].clone();
        self.next_index += 1;

        return Some( return_value );
    }
}   


impl < 'a, FilVal, SnzVal, RingOperator >
    
    ParetoShortCircuit
        < (SimplexFiltered<FilVal>, SnzVal) > for
    
    LazyOrderedCoboundary
        < 'a, FilVal, SnzVal, RingOperator >

    where 
        FilVal:         Clone + Debug + Ord,
        SnzVal:         Clone,
        RingOperator:   Clone + Semiring< SnzVal > + Ring< SnzVal >,      
{
    fn pareto_short_circuit(& self) -> Option< (SimplexFiltered<FilVal>, SnzVal) > {
        if ! self.built_all { return Some( self.coboundary[0].clone() ) }
        return None
    }
} 


//  ORACLE MINOR ASCEND
//  ------------------------------------------



impl < 'a, FilVal, SnzVal, RingOperator >
    
OracleMinorDescend for 

    &'a CliqueBoundaryMatrix
        < FilVal, SnzVal, RingOperator >

    where
        FilVal:             Clone + Debug + PartialOrd + Ord, 
        SnzVal:             Clone, 
        RingOperator:       Clone + Semiring< SnzVal > + Ring< SnzVal >  
{
    type ViewMinorDescend            =   Vec< ( SimplexFiltered<FilVal>, SnzVal ) >;
    type ViewMinorDescendIntoIter    =   std::vec::IntoIter<(SimplexFiltered<FilVal>, SnzVal)>;
    type ViewMinorDescendEntry       =   ( Self::KeyMaj, Self::SnzVal );

    fn view_minor_descend( &self, keymin: Self::KeyMin ) -> Self::ViewMinorDescend {

        let iter = FacetIter {
            clique: &self,
            simp: keymin.clone(),
            removal_location: 0,
            ring_operator: self.ring_operator.clone(),
        };
        let mut vec = iter.collect_vec();     
        vec.shrink_to_fit();
        vec.sort_by(|x,y| y.0.cmp(&x.0));  // NOTE WE GIVE REVERSE ORDER

        return vec
    }
}


// See matrix.rs file for specific definition of SmOracle trait
// impl<FilVal, SnzVal> 

//     SmOracle
//         <SimplexFiltered<FilVal>, SimplexFiltered<FilVal>, SnzVal> for 
        
//     CliqueBoundaryMatrix
//         <FilVal, SnzVal> 
        
//     where
//         FilVal: PartialOrd + Clone + Debug + Eq + Hash + Ord,
//         SnzVal: Neg<Output = SnzVal> + Clone
// {

//     fn ring( &self ) -> &RingMetadata<SnzVal> {
//         &self.ring_operator
//     }

// 	fn maj_dim( &self ) -> MajorDimension { MajorDimension::Row }

// 	fn maj_itr(&self, majkey: &SimplexFiltered<FilVal>) -> Box<dyn Iterator<Item=(SimplexFiltered<FilVal>, SnzVal)> + '_> {
//         let mut vertices = majkey.vertices.clone();
//         let first_vert = vertices[0];
//         vertices.insert(0,0);
//         vertices.shrink_to_fit();
//         Box::new(CofacetIter{
//         	clique: &self,
//         	next_cofacet_vertices: vertices,
//             simplex_filvalue: majkey.filvalue.clone(),
//         	insertion_location: 0,
//             candidate_location: 0,
//         	coeff: self.ring_operator.identity_multiplicative.clone(),
//             first_vert,
//         })
// 	}

// 	fn min_itr( &self, minkey: &SimplexFiltered<FilVal> ) -> Box<dyn Iterator<Item=(SimplexFiltered<FilVal>, SnzVal)> + '_> {
//         Box::new( FacetIter {
//             clique: &self,
//             simp: minkey.clone(),
//             removal_location: 0,
//             coeff: -self.ring_operator.identity_multiplicative.clone()
//         })
// 	}

// 	// Total # of nonzeros in the matrix (if this number can actually be calculated)
//     // Fill later
// 	fn countsnz( &self ) -> Option<usize> {
//         //for
//         //for iter in self.maj_itr()
// 		None
// 	}

// 	fn finiteminors( &self ) -> Option<bool> {
// 		Some(true)
// 	}

// 	fn finitemajors( &self ) -> Option<bool> {
// 		Some(true)
// 	}

//     fn is_pivot(&self, majkey: &SimplexFiltered<FilVal>, minkey: &SimplexFiltered<FilVal>) -> Option<bool> {
//         if majkey.filvalue == minkey.filvalue && majkey.vertices[0] < minkey.vertices[0]{
//             return Some(true);
//         }
//         else {
//             return Some(false);
//         }
//     }

// }


// /// Based, filtered chain complex implementing the [`ChainComplex`](ChainComplex) trait
// pub struct CliqueComplex<SnzVal: Clone, FilVal> {
//     pub ring_operator: RingMetadata<SnzVal>,
//     pub dissimilarity_matrix: Vec<Vec<FilVal>>,
//     pub dissimilarity_value_max: FilVal,
//     pub safe_homology_degrees_to_build_boundaries: Vec<usize>,
//     pub major_dimension: MajorDimension,
//     pub simplex_count: Vec<(usize,usize)>
// }

/// An iterator that iterates all simplices of given dimension
pub struct SimplexIter<'a, FilVal: Clone + Debug, > {
    pub dissimilarity_matrix: &'a Vec< Vec< FilVal > >,
    pub dissimilarity_value_max: FilVal,
    pub filvec: Vec<FilVal>,
    pub vec: Vec<Vertex>,
    pub val: Vertex,
    pub loc: usize,
}

impl <'a, FilVal: Clone + Debug >

    SimplexIter
        <'a, FilVal > 
    {

    pub fn new( 
            dim:usize, 
            dissimilarity_matrix: &'a Vec< Vec< FilVal > >,
            dissimilarity_value_max: FilVal,
        ) 
        -> Self {
        SimplexIter {
            dissimilarity_matrix,
            dissimilarity_value_max,
            filvec: vec![dissimilarity_matrix[0][0].clone(); dim+2],
            vec: vec![0; dim+1],
            val: 0,
            loc: 0,
        }
    }
}


/// implement standard methods of Iterator for SimplexIter struct
impl<'a, FilVal > Iterator for SimplexIter<'a, FilVal > where
FilVal: Clone + PartialOrd + Debug,
{
    type Item = SimplexFiltered<FilVal>;

    fn next(&mut self) -> Option<Self::Item> {
        let size = self.dissimilarity_matrix.len();
        if self.vec.len() > size { return None; }
        loop {
            while usize::from(self.val) <= size-self.vec.len()+self.loc {
                self.vec[self.loc] = self.val;
                self.filvec[self.loc+1] = self.filvec[self.loc].clone();
                for ii in 0..self.loc {
                    if self.filvec[self.loc+1] < self.dissimilarity_matrix[usize::from(self.val)][usize::from(self.vec[ii])] {
                        self.filvec[self.loc+1] = self.dissimilarity_matrix[usize::from(self.val)][usize::from(self.vec[ii])].clone();
                    }
                }
                if self.filvec[self.loc+1] <= self.dissimilarity_value_max {
                    if self.loc == self.vec.len()-1 {
                        self.val += 1;
                        return Some(SimplexFiltered{
                                vertices: self.vec.clone(),
                                filvalue: self.filvec[self.loc+1].clone()
                            });
                    } else {
                        self.loc += 1;
                    }
                }
                self.val += 1;
            }
            if self.loc == 0 { return None; }
            self.loc = 	self.loc - 1;
            self.val =	self.vec[self.loc] + 1;

        }
    }
}

// // Implimentation of ChainComplex trait. See chx.rs for definition of ChainComplex trait.
// impl<SnzVal, FilVal> ChainComplex<SimplexFiltered<FilVal>, SnzVal, FilVal> for CliqueComplex<SnzVal, FilVal> where
// SnzVal: Clone + Neg<Output = SnzVal> + Debug,
// FilVal: PartialOrd + Copy + Debug + Eq + Hash + Ord
// {
//     type Matrix = CliqueBoundaryMatrix<FilVal, SnzVal>;

//     fn get_smoracle(
//         &self,
//         major_dim       :   MajorDimension, // row or column
//         transform       :   ChxTransformKind
//     ) -> CliqueBoundaryMatrix<FilVal, SnzVal> {
//         CliqueBoundaryMatrix {
//             ring_operator: self.ring_operator.clone(),
//             maxdis: self.dissimilarity_value_max,
//         	dismat: self.dissimilarity_matrix.clone(),
//             cutoff_matrix: cutoff_matrix(&self.dissimilarity_matrix, self.dissimilarity_value_max)
//         }
//     }

//     /// Return an iterator of all simplices for a given dimension
//     ///
//     /// # Parameters
//     /// -`dim`: The given dimension
//     fn keys_unordered_itr(&self, dim:usize) -> Box<dyn Iterator<Item=SimplexFiltered<FilVal>> + '_> {
//         Box::new( SimplexIter {
//             chx: &self,
//             filvec: vec![self.dissimilarity_matrix[0][0].clone(); dim+2],
//             vec: vec![0; dim+1],
//             val: 0,
//             loc: 0
//         })
//     }

//     fn keys_ordered(&self, dim: usize) -> Vec<SimplexFiltered<FilVal>> {
//         let mut vector: Vec<SimplexFiltered<FilVal>> = Vec::new();

//         for simp in self.keys_unordered_itr(dim) {
//             vector.push(simp.clone());
//         }

//         vector.sort();
//         return vector;
//     }

//     fn key_2_filtration(&self, key: &SimplexFiltered<FilVal>) -> FilVal {
//         return key.filvalue;
//     }

//     fn max_filtration(&self) -> FilVal {
//         return self.dissimilarity_value_max;
//     }
// }






//  ==============================================================================================
//  LAZY SORTED ITERATOR
//  ==============================================================================================




// pub struct LazySortedCofacets< 'a, FilVal, OrderComparator >{
//     dismat:             &'a Vec< Vec< FilVal > >,
//     facet_ver:           Vec< Vertex >,
//     facet_fil:          FilVal,
//     nbr_iters_ordfil:   &'a Vec< Vec< (Vertex,FilVal) > >, 
//     nbr_iters_ordlex:   &'a Vec< Vec< (Vertex,FilVal) > >,     
//     maxdis:             FilVal,
//     order_comparator:   OrderComparator,
// }

// impl < 'a, FilVal, OrderComparator >

//     Iterator for 

//     LazySortedCofacets
//         < 'a, FilVal, OrderComparator >

//     where
//         FilVal:             Clone + Ord,
//         OrderComparator:    StrictlyLess <  usize >,
//         Iter:               Iterator< Item = usize > + PeekUnqualified,  !!!! LOOK AT RECENT CRATE ON MEMORY USAGE
// {    
//     type Item = usize;

//     fn next(&mut self) -> Option<Self::Item> {

//         println!("move this function to another file");        

//         println!("we require neighbors to be sorted lexicographically, after sorting by diameter; this ensures correct behavior");

//         let numiter = self.nbr_iters_ordfil.len();
//         let numiter_minus1 = numiter-1;
//         // abreviations:
//         //  ii = iterator index
//         let mut maxii=0; // iterator index of max length edge
//         let mut maxfil; // max filtration value    
//         let mut maxnbr;                
//         // let mut tryiio; // offset between maxii and tryii
//         let mut tryii; // iterator index of an iterator we will "try out"
//         let mut tryfil;
//         let mut trynbr;
//         let mut iter;
//         let mut restart_for_loop = false;
//         let mut obtained_next_output = false;
//         match self.nbr_iters_ordfil[0].next() {
//             None => return None,
//             Some( x ) => { maxnbr = x; maxfil = self.dismat[0][maxnbr].clone() }
//         }

//         println!("FIRST ITERATE IN LEXICOGRAPHIC ORDER OVER NEIGHBORS OF FIRST VERTEX");
//         let disvec;
//         while maxfil < self.facet_fil {
//             disvec = self.dismat[maxnbr];
//             if self.facet_ver.iter().map(|x| disvec[*x as usize] ).max().unwrap() < self.facet_fil {
                
//             }
//         }


//         let most_remote = find_min(  // find the vertex with the fewest neighbors
//                     self.facet_ver.iter().map(
//                             |x| 
//                             self.nbr_iters_ordfil[*x as usize].len() 
//                         ) 
//                 ).unwrap();

//         for (outer_v, outer_dis) in self.nbr_iters_ordlex[most_remote as usize ].iter().cloned() {
//             if outer_v >= most_remote { break }
//             let disvec = & self.dismat[ outer_v as usize ];
//             let diam = self.facet_fil();
//             let mut diam_ok = true;
//             for inner_v in & self.facet_ver {
//                 if disvec[ inner_v as usize ] > diam { diam_ok = false; break }                                
//             }
//             if diam_ok { shortcircuit_vertex=Some( outer_v ); break }
//         }        

//         println!("THEN SKIP FORWARD THE OTHER ITERATOR UNTIL YOU GET TO GREATER DIAMETERS");        

//         loop { // loop A
//             // let maxii=0;
//             // let maxnbr be the closest remaining neighbor to maxii
//             // let maxfil = dismat[ self.facet_ver[maxii], maxnbr ]
//             // for every other vetex v_i in self.facet_ver = [v_0, .., v_n]
//             // loop {
//             //   - let n_i be the closest remaining neihbor to v_i, and d_i be its distance
//             //     - if n_i does not exist then terminate/return none: there are no more cofacets
//             //     - if d_i > maxfil
//             //       - reset maxii = i; maxfil = d_i; maxnbr = n_i; [ note that n_i couled equal maxnbr]
//             //       - restart the for-loop
//             //     - if d_i <= maxfil && n_i != maxnbr
//             //       - remove n_i and restart the inner loop
//             // }
//             // if you make it here, then maxnbr can be added to create the next facet! return this value!
//             //
//             // proof of correctness: 
//             // - proceed by induction
//             // - follows that on step k, the smallest simplex we can add has filtration value maxfil
//             // - each simplex can be associated with its longest edge(s) enter;            
//             // - any subsequent simplices added in a true/correct enumeration will have equal or greater diam; the associated edges that bear witness to those creations will of course not be lost if we throw out edges of length less than maxfil
//             // - so we are gauranteed that throwing out edges will leave all the relevant witness edges
//             // - clearly none of the edges we throw out are witness edges
//             // - therefore the edges that we don't throw out are exactly the witness edges
//             // - and clearly those edges are identified in order of length
//             // - there are some minor details conserning edges of equal length; this should be handled with lexicographic order; rather we order neighbors first by edge length then by vertex num; this effectively gives every edge a unique length
//             restart_for_loop = false;
//             for tryiio in 0..numiter_minus1 { // loop B
//                 tryii  = maxii + tryiio % numiter;
//                 iter = &mut self.nbr_iters_ordfil[ tryii ];
//                 loop { // loop C
//                     match iter.peek_unqualified() {
//                         None => return None, // there exist no more cofacets
//                         Some ( x ) => {
//                             trynbr = x.clone();
//                             tryfil = self.dismat[ self.facet_ver[ tryii ] ][ trynbr ];
//                             if tryfil > maxfil { // we'll have to restart the forloop
//                                 maxii = tryii; maxfil = tryfil; maxnbr = trynbr;
//                                 restart_for_loop = true;
//                                 break; // now loop B will restart
//                             } else if tryfil <= maxfil && trynbr != maxnbr {
//                                 let _ = iter.next(); // remove the nearest remaining neighbor (of v_tryii) from iter
//                                 // now the loop C restarts
//                             }
//                         }
//                     }
//                 }
//                 if maxfil <= simplex_diameter { restart_for_loop = true } // skip over cofacets with filtration value equal to the facet
//                 if restart_for_loop { break }  // breaking the for-loop essentially restarts it              
//             }
//             if ! restart_for_loop { 
//                 // we have found a valid result! 
//                 // the next element of every iterator is maxnbr; clear these elements before proceeding
//                 for iter in self.nbr_iters_ordfil.iter_mut() { _ = iter.next(); }
//                 return Some( maxnbr ) 
//             }
//         }
//     }
// }

// /// Given a sorted vector v = [v0,..,vm] and an inteter u, find the
// /// minimum k such that [v0,..,u, vk,..,vm] is sorted
// pub fn get_insertion_locus( vertices: & Vec<usize>, vertex: usize ) -> usize {
//     println!("move this function to another file");
//     let mut insertion_locus = 0;
//     let nvertices = vertices.len();
//     while vertices[insertion_locus] < vertex && insertion_locus < nvertices {
//         insertion_locus += 1;
//     }
//     return insertion_locus
// }

// /// Inserts a vertex at the given location, returning a new vec and leaving the input vec unchanged.
// pub fn insert_out_of_place( vertices: & Vec<usize>, vertex: usize, insertion_locus: usize ) 
//     ->
//     Vec<usize>
//     {
//     println!("move this function to another file");        
//     let nvertices = vertices.len();
//     let mut cofacet = Vec::with_capacity( nvertices + 1 );
//     for v in vertices[0..insertion_locus].iter().cloned() { 
//         cofacet.push(v);
//     }
//     cofacet.push(vertex);
//     for v in vertices[insertion_locus..nvertices-1].iter().cloned() {
//         cofacet.push(v);
//     }
//     return cofacet
// }


// /// Returns the entry of the coboundary matrix corresponding to the cofacet
// /// obtained by adding vertex `v` to `facet`.
// pub fn get_coboundary_entry< SnzVal, RingOperator >( 
//             facet: Vec<usize>, 
//             vertex: usize, 
//             ring_operator: RingOperator,
//         )
//         -> ( Vec<usize>, SnzVal )
//     where
//         RingOperator:       Semiring< SnzVal > + Ring< SnzVal >,
// {
//     let insertion_locus     =   get_insertion_locus( &facet, vertex );
//     let cofacet = insert_out_of_place( &facet, vertex, insertion_locus);
//     let coeff = ring_operator.minus_one_to_power( insertion_locus );
//     return (cofacet, coeff)
// }        

























//  =========================================================================
//  UNIT TESTS
//  =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    use itertools::Itertools;
    use num::Float;
    use num::rational::Ratio;
    use crate::simplicial::from::graph_weighted::user_interface::{CliqueParams, cliques_in_order};
    // use oat_rust::utilities::homology::clique::{CliqueParams, cliques_in_order};
    use crate::simplicial::from::relation::{homology_basis_from_dowker, UseClearing};
    // use oat_rust::utilities::cw_complexes::simplices_weighted::clique::{CliqueBoundaryMatrix, get_cutoff_matrix, SimplexIter, Simplex};
    use crate::point_cloud::unit_circle;    
    use oat_rust::rings::operator_structs::field_prime_order::BooleanFieldOperator;
    use oat_rust::rings::operator_structs::ring_native::{DivisionRingNative, field_rational_i64};    
    use oat_rust::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
    use oat_rust::matrices::matrix_oracle_traits::{OracleMinorDescend, OracleMajorAscend};
    use oat_rust::matrices::matrix_types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
    use oat_rust::matrices::matrix_types::oracle_ref::OracleRef;
    use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor_with_clearing, ParetoShortCircuit};    
    use oat_rust::utilities::partial_order::{OrderComparatorAutoLt, is_sorted_strictly, OrderComparatorLtByKey, OrderComparatorAutoGt};
    use oat_rust::utilities::random::random_matrix;
    use oat_rust::utilities::distances::{rowwise_distances, minmax};
    use std::f64::INFINITY;
    use std::f64::consts::PI;
    use std::io;
    use std::marker::PhantomData;
    use std::time::{Duration, Instant};
    use ordered_float::OrderedFloat;    
    
    #[test]
    fn check_that_some_basic_functions_run_without_error() {

        let m = 20;
        let maxdim = 1;
        let maxdis = None;

        let pcloud = unit_circle( m, Some(-1.0 .. 1.0));

        let dissimilarity_matrix = rowwise_distances(pcloud);        

        for i in 0 .. dissimilarity_matrix.len() {
            for j in i .. dissimilarity_matrix.len() {
                assert_eq!( dissimilarity_matrix[i][j], dissimilarity_matrix[j][i] );
            }
        }
     
    
        let ring_operator = field_rational_i64();
        let params = CliqueParams::new(
                    & dissimilarity_matrix, 
                    maxdim, 
                    maxdis, 
                    ring_operator.clone(),
                );
        let boundary_matrix = CliqueBoundaryMatrix::new( & params );
        let boundary_matrix_ref = & boundary_matrix;    
        let keymaj_vec = cliques_in_order( & params );
        let keymin_vec = cliques_in_order( & CliqueParams::new(
                    & dissimilarity_matrix, maxdim+1, maxdis, ring_operator.clone()
                ) ) ;
    
    
        verify_viewmajorascend_compatible_with_viewminordescend(
                boundary_matrix_ref,
                keymin_vec.iter().cloned(),
                keymaj_vec.iter().cloned(),
            );       
    
        // println!("(oat_rust_python) keymaj_vec = "); 
        // for entry in keymaj_vec.iter() { println!("{:?}", entry) }
    
    
        // println!("press any key to continue");
        // let mut guess = String::new();
        // io::stdin()
        //     .read_line(&mut guess)
        //     .expect("Failed to read line");    
    
        let iter_keymaj = keymaj_vec.iter().cloned();    
    
        println!("check that oracle has strictly sorted rows");
        // print_indexed_major_views( & boundary_matrix_ref, iter_keymaj.clone() );  // print the major views       
        for keymaj in iter_keymaj.clone() {        
            assert!( is_sorted_strictly( 
                                            &boundary_matrix_ref.view_major_ascend(keymaj.clone()).collect_vec() , 
                                            &OrderComparatorLtByKey::new( OrderComparatorAutoLt::new() ) 
                                        ) );
        }
    
        // println!("press enter to continue");
        // let mut guess = String::new();
        // io::stdin()
        //     .read_line(&mut guess)
        //     .expect("Failed to read line");        
    
        println!("check that oracle has strictly sorted columns");
        // print_indexed_minor_views( & boundary_matrix_ref, iter_keymaj.clone() );  // print the major views        
        for keymaj in iter_keymaj.clone() {
            assert!( is_sorted_strictly(    &boundary_matrix_ref.view_minor_descend(keymaj).iter().cloned().collect_vec() , 
                                            &OrderComparatorLtByKey::new( OrderComparatorAutoGt::new() )  // NOTE THAT HERE WE USE GT 
                                        ) );
        }    
    
        // println!("press enter to continue");
        // let mut guess = String::new();
        // io::stdin()
        //     .read_line(&mut guess)
        //     .expect("Failed to read line");       
    
        println!("starting umatch");
        let umatch = new_umatchrowmajor_with_clearing(
                boundary_matrix_ref, 
                iter_keymaj.clone(), 
                ring_operator.clone(), 
                OrderComparatorAutoLt::new(), 
                OrderComparatorAutoLt::new(), 
            );      
    
        // println!("start build bd matrix");
        // let boundary_matrix = oat_rust::utilities::homology::clique::get_clique_boundary_matrix(
        //     dissimilarity_matrix,
        //     maxdis,
        //     field_rational_i64()        
        // );
    
        // println!("start umatch");    
        // let umatch = oat_rust::utilities::homology::clique::umatch_from_clique(
        //     & boundary_matrix,
        //     maxdim,
        // );
    
           
        // let iter_keymaj = 
        //     (0..maxdim+1).map(
        //         |dim|
        //         {
        //             let mut vec = SimplexIter::new( 
        //                         dim,
        //                         & umatch.array_mapping_ref().dismat,
        //                         umatch.array_mapping_ref().maxdis,                
        //                     )
        //                     .collect_vec();
        //             vec.sort_by(|x,y| y.filvalue.cmp(&x.filvalue));
        //             vec
        //         }
        //     )
        //     .flatten();    
    
    
        println!("setting up to unpack");  
        let dim_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.dim();
        let fil_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.fil().into_inner();    
        let barcode = umatch.barcode( iter_keymaj, dim_fn, fil_fn, true , true);
    
        println!("getting intervals, reps, bounding chains");
        let (intervals, representatives, bounding_chains ) = barcode.unwrap();
    
        let convert_coefficients = |x: Ratio<i64> | (x.numer().clone(), x.denom().clone());
        let convert_simplex = |x: SimplexFiltered< OrderedFloat<f64> > | x.vertices.iter().map(|y| y.clone() as usize ).collect_vec();
        
        println!("start reformat reps");     
        let represntatives_new = representatives.unwrap().iter().map(
                    |x| // each x is a list of cycle reps
                    x.iter().map(
                        |y| // each y is a cycle rep
                        y.iter().map(
                            |z| // each z is an entry in a cycle rep
                            ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
                        ).collect_vec()
                    ).collect_vec()
                ).collect_vec();
    
        println!("start reformat bounding chains");                 
        let bounding_chains_new = bounding_chains.unwrap().iter().map(
                    |x| // each x is a list of cycle reps
                    x.iter().map(
                        |y| // each y is a cycle rep
                        y.iter().map(
                            |z| // each z is an entry in a cycle rep
                            ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
                        ).collect_vec()
                    ).collect_vec()
                ).collect_vec();              
    }
}    