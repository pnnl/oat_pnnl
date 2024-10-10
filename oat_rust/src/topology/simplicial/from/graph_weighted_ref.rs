//! The filtered Clique complex (aka filtered Flag complex, filtered Vietoris-Rips complex); often constructed from a metric space or point cloud


//! Data structures for filtered clique complexes.
//! 
//! # Clique chain complex oracles
//! Mathematically, a filtered clique complex is determined by two pieces of data:
//! * a **dissimilarity matrix**
//!     * we represent this by a symmetric matrix *S* that has been flattened into a vector *v*
//! * a **threshold**  that determines where we stop growing the filtration
//!     * we represent this as a real number *t*
//! The boundary matrices of this chain complex oracle are indexed by `SimplexFiltered` objects which are
//! (essentially) pairs of form `(simplex, filtration_value)`.
//! 
//! Example:
//!     - construct a dissimilarity matrix
//!     - construct the associated clique complex
//!     - access a row + column of the boundary matrix
//!     - compute the dimension 1 barcode




use std::ops::Bound;
use std::sync::Arc;
use std::{ops::Neg, marker::PhantomData};
use std::convert::TryInto;
use std::hash::Hash;
use std::fmt::Debug;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use crate::utilities::iterators::merge::hit::HitMerge;
use crate::algebra::vectors::entries::KeyValGet;


// use crate::matrix::{SmOracle, RingMetadata, MajorDimension};
// use crate::chx::{ChainComplex, ChxTransformKind};

use std::cmp::Ordering;

use crate::algebra::chains::factored::{factor_boundary_matrix, FactoredBoundaryMatrix};
use crate::algebra::matrices::operations::umatch::row_major::{ParetoShortCircuit, Umatch, Umatch::factor_with_clearing};
use crate::algebra::matrices::query::{IndicesAndCoefficients, ViewRowAscend, ViewColDescend, MatrixEntry};
use crate::algebra::rings::operator_traits::DivisionRing;
use crate::algebra::rings::operator_traits::{Semiring, Ring};
use crate::utilities::iterators::general::{PeekUnqualified, find_min};
use crate::utilities::order::{ JudgePartialOrder, OrderOperatorByKey, OrderOperatorByKeyCutsom, JudgeOrder };
use crate::utilities::iterators::general::{minmax, IntersectOrderedIterators};
// use oat_rust::utilities::indexing_and_bijection::sparsevec_lookup;

use crate::topology::simplicial::misc::simplex_permutation::Simplex;
use crate::topology::simplicial::simplices::{SimplexFiltered, OrderOperatorTwistSimplexFiltered};

// The specific type we use to store vertices
type Vertex = u16;
// type OrderOperatorSimplexFiltered =   OrderOperatorByKey< 
//                                                 SimplexFiltered< OrderedFloat< f64 > >, 
//                                                 SnzValRational,
//                                                 (
//                                                     SimplexFiltered< OrderedFloat< f64 > >, 
//                                                     SnzValRational,
//                                                 )
//                                             >;
type OrderOperatorSimplexFiltered< FilVal, Coefficient > =   OrderOperatorByKeyCutsom< 
                                                                SimplexFiltered< FilVal >, 
                                                                Coefficient,
                                                                (
                                                                    SimplexFiltered< FilVal>, 
                                                                    Coefficient,
                                                                ),
                                                                OrderOperatorAuto,
                                                            >;                                            


//  ===========================================================
//  SORTED ENUMERATION OF CLIQUES
//  ===========================================================


// pub trait RowIter
//             < 'a, T, I: Iterator< Item=T > > {
//     fn row_iter( &'a self) -> I;
// }

// impl < 'a, T >

//     RowIter
//         < 'a, &'a T, std::slice::Iter< 'a, T > > for 
//     Vec< Vec< T > >
// {
//     fn iter( &'a self ) -> std::slice::Iter< 'a, T > { self.iter() }
// }

// pub trait RowIter
//             < T, I: Iterator< Item=T > > {
//     fn row_iter( & self) -> I;
// }

// impl < T >

//     RowIter
//         < & T, std::slice::Iter< T > > for 
//     Vec< Vec< T > >
// {
//     fn iter( & self ) -> std::slice::Iter< T > { self.iter() }
// }

pub fn test() {
    let v = Arc::new( vec![ vec![0,1], vec![0,1] ] );
    let u = v.iter();
    let a = minmax(u);

    let w = vec![ vec![0,1], vec![0,1] ];
    let z = w.iter();
    let a = minmax(z);
}


/// A vector of simplices sorted first by dimension (ascending `0..dimension_max`) 
/// then by diameter (descending) then by lexicographic order (descending)
pub fn filtered_cliques_in_order< FilVal, DisMat >( 
            dimension_max:                  isize,
            dissimilarity_value_max:        FilVal,
            dissimilarity_matrix:           DisMat,
            dissimilarity_matrix_size:      usize,            
        ) 
        -> 
        Vec< SimplexFiltered< FilVal > >  
    where
        FilVal:                         Clone + Debug + Ord,
        DisMat:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > + ViewRowAscend + MatrixEntry,
        DisMat::EntryMajor:   KeyValGet< usize, FilVal >,
    {
        let order_operator = OrderOperatorTwistSimplexFiltered{};
        let dissimilarity_matrix_ref = & dissimilarity_matrix;
        let vec: Vec<_> = (0..dimension_max+1).map(
            |dim|
            {
                let mut vec =   SimplexIter::new(                          
                                        dissimilarity_matrix_ref, // a simple trick to make a copy of the dissimilarity matrix that implements the Copy trait
                                        dissimilarity_matrix_size,                
                                        dissimilarity_value_max.clone(),                            
                                        dim,                               
                                    )
                                    .collect_vec();
                vec.sort_by(|x,y| order_operator.judge_cmp(y,x) );
                vec
            }
        )
        .flatten()
        .collect();
        return vec
}





//  ===========================================================
//  CUTOFF AND DISMAT
//  ===========================================================


/// The "cutoff" of a dissimilarity matrix
///
/// # Parameters
/// - `dissimilarity_matrix`: Iterable that runs over rows of the dissimilarity matrix
/// - `dissimilarity_value_max`: The radius of neighborhood
pub fn get_cutoff_matrix
            < 'a, FilVal, DisMat, DisMatEntry >
            (dissimilarity_matrix: DisMat, dissimilarity_value_max: FilVal) 
            -> 
            Vec<Vec<Vertex>> 
    where
        FilVal:             'a + PartialOrd + Debug,
        DisMat:             IntoIterator,
        DisMat::Item:       IntoIterator< Item=DisMatEntry >,
        DisMatEntry:        KeyValGet< usize, FilVal >,

{
    let mut output: Vec<Vec<Vertex>> = Vec::new();
    for row in dissimilarity_matrix {
        let mut index_vec = Vec::new();
        for entry in row {
            let k = entry.key(); let v = entry.val();
            if v <= dissimilarity_value_max { index_vec.push( k as Vertex ) }
        }
        output.push(index_vec);
    }
    return output;
}


// UNDER CONSTRUCTION -- WE PROBABLY HAVE PLANS TO USE THIS
// /// A sparse symmetric matrix with efficient look-up of individual entries.
// /// 
// /// The underlying data strcutre is a row-major vec-of-vec matrix for the 
// /// upper-triangular-including-diagonal part of the matrix.
// pub struct SymmetricVecVecIntegerIndexMatrix< T > {
//     vecvec:    Vec< Vec< (usize, T) > >
// }

// impl < T: Clone > 

//     SymmetricVecVecIntegerIndexMatrix< T > {  

//     /// Create a new `SymmetricVecVecIntegerIndexMatrix` whose underlying data structure (a row-major vec-of-vec matrix
//     /// representing the upper-triangular-including-diagonal part of the matrix ) is `vecvec`.
//     /// 
//     /// Throws an error the entries of any vector in `vecvec` either (i) contain an entry that lies strictly to the left of
//     /// the diagonal, or (ii) is not soted in strictly ascending order by index
//     pub fn new( vecvec: Vec< Vec< (usize,T) > >) -> Self {            
//         for (vec_num, vec) in vecvec.iter().enumerate() {
//             // check that the entries in this row lie above the diagonal
//             if ! vec.is_empty() && vec[0].0 < vec_num { panic!("Construction of SymmetricVecVecIntegerIndexMatrix failed because row {:?} of the input data contains an entry with index < {:?}.", vec_num, vec_num ); }            
//             // check that input is sorted
//             let is_sorted_strictly = (0..vec.len()-1).all(|i| vec[i].0 < vec[i+1].0 );
//             if ! is_sorted_strictly { panic!("Construction of SymmetricVecVecIntegerIndexMatrix failed because row {:?} of the input data is not sorted.", vec_num); }            
//         }
//         return SymmetricVecVecIntegerIndexMatrix{ vecvec }
//     }

//     // Return entry `Some( (i,j) )`, if there is a value stored for `(i,j)`; otherwise return `None`.
//     pub fn entry(&self, i: usize, j: usize) -> Option< T > {
//         let pq = match i < j {
//             true => (i,j), false => (j,i)
//         };
//         return sparsevec_lookup( &self.vecvec[pq.0], pq.1 ).map(|x| x.1.clone() )
//     }  

// }



//  ===========================================================
//  WEIGHTED CLIQUE BOUNDARY MATRIX 
//  ===========================================================

// pub struct ChainComplexVrFilteredUmatch<FilVal, Coefficient, RingOperator>
//     where 
//         Coefficient:         Clone,
//         RingOperator:   Semiring< Coefficient > + Ring< Coefficient >,
//         FilVal:         Hash,
// {
//     umatch:     Umatch<  
//                         ChainComplexVrFiltered<FilVal, Coefficient, RingOperator>,
//                         RingOperator, 
//                         OrderOperatorSimplexFiltered,
//                         OrderOperatorSimplexFiltered,
//                     >,
//     dimension_max:     isize,

// }        

//  ===========================================================
//  FILTERED VIETORIS RIPS BOUNDARY MATRIX 
//  ===========================================================

/// Boundary matrix represented by a matrix oracle.
#[derive(Clone,Copy,Debug,new)]
pub struct ChainComplexVrFiltered< DisMat, FilVal, Coefficient, RingOperator >
    where 
        Coefficient:         Clone,
        RingOperator:   Semiring< Coefficient > + Ring< Coefficient >,
        DisMat:         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                        ViewRowAscend + 
                        MatrixEntry,
{
    pub ring_operator: RingOperator, // ring meta data
    /// The "dissimilarity_value_max" value represents the maximum of filtration value
    pub dissimilarity_value_max: FilVal,
    /// The dissimilarity matrix
    pub dismin: FilVal,
    /// The dissimilarity score (or "diameter") of the emtpy simplex
	pub dissimilarity_matrix: DisMat,
    /// The number of rows (respectively, columns) of the dissimilarity matrix
    pub dissimilarity_matrix_size: usize,
    /// A vector representing the neighborhood within "dissimilarity_value_max" of each vertex
    pub cutoff_matrix: Vec<Vec<Vertex>>,
    pub phantomsnzval: PhantomData< Coefficient >,
}


/// Methods of ChainComplexVrFiltered struct
impl < DisMat, FilVal, Coefficient, RingOperator  > 

    ChainComplexVrFiltered
        < DisMat, FilVal, Coefficient, RingOperator > 

    where
        FilVal:                         Clone + Debug + PartialOrd + Ord,
        Coefficient:                         Clone,
        RingOperator:                   Semiring< Coefficient > + Ring< Coefficient >,
        DisMat:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                                        ViewRowAscend + 
                                        MatrixEntry, 
        DisMat::EntryMajor:   KeyValGet< usize, FilVal >,       
{

    /// Construct a (lazy) boundary matrix for the filtered VR complex of a symmetric matrix.
    /// - `dissimilarity_matrix` a symmetric matrix representing pairwise dissimilarity values
    /// - `dissimilarity_matrix_size` an integer representing the number of rows/columns of the dissimilarity matrix; this parameter is necessary because the matrix oracle traits are passive, and do not give a means to directly count or generate the number of row/column indices
    /// - `dissimilarity_value_max` maximum dissimilarity threshold
    /// - `dismin` diameter of the empty simplex; this can be any number that is equal or less than every structural nonzero entry in the dissimilarity matrix
    /// - `ring_operator`: ring operator for the coefficient ring
    /// The filtration parameter of
    /// - the empty simplex `[]` is `dismin`
    /// - a singletime `[i]` is `dissimilarity_matrix[i,i]`
    /// - an edge `[i,j]` is `dissimilarity_matrix[i,j]`
    /// - an arbitrary simplex `s` is the maximum of the filtration parameters of cardinality 0, 1, and 2 subsets of `s`
    /// subject to the condition that the filtration parameter of
    /// - a singleton `[i]` is undefined if `dissimilarity_matrix[i,i]` is structurally zero
    /// - an edge `[i,j]` is undefined if `dissimilarity_matrix[i,j]` is structurally zero
    /// - an arbitrary simplex `s` is undefined if the filtration parameter of any cardinality 0, 1, or 2 subset of `s` is undefined
    /// This means that some vertices may never enter the filtration.
    /// 
    /// 
    /// The constructor checks that 
    /// - the input matrix is symmetric
    /// - the filtration parameter of every edge equals or exceeds the filtration parameter of its vertices.
    pub fn new(
            dissimilarity_matrix:                         DisMat,
            dissimilarity_matrix_size:                    usize,
            dissimilarity_value_max:        Option< FilVal >,
            dismin:                         FilVal,
            ring_operator:                  RingOperator,
        ) -> Self
    {

        println!("In future the ChainComplexVrFiltered struct will be updated to use the SymmetricVecVecIntegerIndexMatrix in the dissimilarity_matrix field.");        

        let m = dissimilarity_matrix_size;
        let disemp_nonegreater = NoneGreater::from_val( dismin.clone() );

        println!("checking that (1) dissimilarity matrix is symmetric, (2) all nonzero entries in a given row, i, match or exceed entry (i,i), and (3) and all nonzero entries match or exceed the diameter of the empty simplex");
        for i in 0 .. m {
            for entry in dissimilarity_matrix.view_major_ascend(i) {   
                let j = entry.key(); 
                let v = entry.val();             
                if Some( v ) != dissimilarity_matrix.entry( j, i ) {
                    panic!("Input matrix is not symmetric: entry {:?} doesn't equal entry {:?}", (i,j), (j,i));                    
                }
            }

            // check that the filtration function respects face relations; this means that diagonal elements
            // (1) are minimal in their respective rows
            // (2) equal or exceed the diameter of the empty simplex
            let min = dissimilarity_matrix.view_major_ascend(i).into_iter().map(|x| x.val() ).min();    // min entry this row
            let dia = dissimilarity_matrix.entry(i,i);  
            if min != dia {
                panic!("Entry ({:?},{:?}) of the dissimilarity matrix is {:?} but the minimum structural nonzero entry in row {:?} is {:?}.  In this case `None` represents a value greater strictly greater than `Some(x)`, for every filtration value `x`.", i, i, &dia, i, &min );
            }
            if NoneGreater::new(min.clone() ) < disemp_nonegreater { // we wrap filtration values in `NoneGreater` so that None values will be regarded as maximal instead of minimal
                panic!("The dissimilarity matrix assigns  ({:?},{:?}) a value of {:?}, which is strictly less than the value {:?} assigned to the empty simplex.  In this case `None` represents a value greater strictly greater than `x`, for every filtration value `x`.", i, i, min.clone(), dismin.clone() );
            }            
        }        
    
        // if no dissimilarity_value_max is provided then use the enclosing radius
        let dissimilarity_value_max = dissimilarity_value_max.unwrap_or( 
                    minmax( 
                            (0..dissimilarity_matrix_size).map(
                                    |x| 
                                    dissimilarity_matrix.view_major_ascend(x).into_iter().map(
                                            |x| 
                                            x.val()
                                        ) 
                                ) 
                        )
                        .unwrap_or( dismin.clone() )
                );
        let cutoff_matrix = get_cutoff_matrix( dissimilarity_matrix.views_major_ascend(0..dissimilarity_matrix_size), dissimilarity_value_max.clone() );
        
        ChainComplexVrFiltered { ring_operator, dissimilarity_value_max, dissimilarity_matrix, dissimilarity_matrix_size, cutoff_matrix, dismin, phantomsnzval: PhantomData }
    }            

    /// Output the diameter of a simplex
    ///
    /// # Parameters
    /// -`vertices`: vertices of the simplex
    /// # Returns
    /// Return 
    /// - `None` if any pairwise distance is not recorded in the
    /// - `self.dismin` if the `vertices` is empty
    /// - the maximum of the pairwise dissimilarities, otherwise
    fn diameter( & self, vertices: & Vec< Vertex > ) -> Option< FilVal > {
        let mut diam = self.dismin.clone(); 
        let mut a;
        let mut b;
        for ii in 0..vertices.len() {
            a = usize::from( vertices[ii] ); 
            for jj in ii .. vertices.len() {
                b = usize::from( vertices[jj] ); 
                match self.dissimilarity_matrix.entry( a, b ) {
                    None => { return None } // the simplex never enters the filtration
                    Some( diam_bound ) => {
                        if diam_bound > diam { diam = diam_bound.clone(); }
                    }
                }
            }
        }
        return Some(diam);
    }    

    /// Returns a reference to the internally stored symmetric dissimilarity matrix, wrapped in a struct that allows query operations to be performed.
    pub fn dissimilarity_matrix_ref(&self) -> DisMat { &self.dissimilarity_matrix }

    /// Returns the dissimilarity value of two vertices, wrapped in a `NoneGreater` struct
    /// 
    /// The purpose of the `NoneGreater` struct is to "switch" the order of None values,
    /// so that they are regarded as strictly greater (rather than strictly less) than all `Some(x)` values.
    pub fn dissimilarity_as_nonegreater( &self, a: Vertex, b: Vertex ) -> NoneGreater<FilVal> 
        { NoneGreater { opt: self.dissimilarity_matrix.entry(a as usize, b as usize) } }
    
    /// The maximum dissimilarity threshold
    pub fn dissimilarity_value_max( &self ) -> FilVal { self.dissimilarity_value_max.clone() }

    /// The ring operator for the coefficient ring.
    pub fn ring_operator( &self ) -> RingOperator where RingOperator: Clone { self.ring_operator.clone() }

    /// A vector of simplices sorted first by dimension (ascending `0..dimension_max`,
    /// including `dimension_max`) then by diameter (descending) then by lexicographic order
    pub fn cliques_in_order( &self, dimension_max: isize ) -> Vec< SimplexFiltered< FilVal > > { 
        return filtered_cliques_in_order( dimension_max, self.dissimilarity_value_max(), self.dissimilarity_matrix_ref(), self.dissimilarity_matrix_size )
    }    

    /// Returns a factorization of the boundary matrix restricted to simplices
    /// of dimension 0 .. dimension_max+1 (inclusive)
    /// 
    /// This factorization can be used to compute betti numbers, cycle representatives,
    /// and more.  See the documenation for `FactoredChainComplex`.
    pub fn factor_from_ref( &self, max_homology_dim: isize ) 
            ->
            FactoredBoundaryMatrix< 
                    &Self, 
                    RingOperator, 
                    OrderOperatorSimplexFiltered< FilVal, Coefficient >, 
                    SimplexFiltered<FilVal>,
                    (SimplexFiltered<FilVal>, Coefficient),                    
                    Vec< SimplexFiltered<FilVal> >,                     
                >    
                // Umatch<  
                //         &Self,                        
                //         RingOperator, 
                //         OrderOperatorSimplexFiltered< FilVal, Coefficient >,
                //         OrderOperatorSimplexFiltered< FilVal, Coefficient >,
                //     >
        // where
        //     FilVal:         Hash,
        //     RingOperator:   DivisionRing< Coefficient >,
        where
            FilVal:             Clone + Debug + PartialOrd + Ord + Hash, 
            Coefficient:             Clone + Debug + PartialEq, 
            RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient > + DivisionRing< Coefficient >,            
    {
        let iter_keymaj     =   self.cliques_in_order( max_homology_dim );
        // Umatch::factor_with_clearing(
        //     self, 
        //     iter_keymaj, 
        //     self.ring_operator.clone(), 
        //     OrderOperatorAutoLt::new(),  // this will be wrapped in an AutoComparatroLtByKey
        //     OrderOperatorAutoLt::new(),  // this will be wrapped in an AutoComparatroLtByKey
        // )
        factor_boundary_matrix(
            & self,
            self.ring_operator.clone(),
            OrderOperatorAuto,
            iter_keymaj,
        )
    }

}


//  ===========================================================
//  COFACET ITERATOR
//  ===========================================================


/// An iterator that runs over all cofacets of a given simplex
struct IterCoboundary
        < 'a, DisMat, FilVal, Coefficient, RingOperator >
    where
        FilVal:             Clone + Debug,
        Coefficient:             Clone,
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,     
{
	pub clique: &'a ChainComplexVrFiltered< DisMat, FilVal, Coefficient, RingOperator >,
	pub next_cofacet_vertices: Vec<Vertex>,
    pub simplex_filtration: FilVal,
	pub insertion_location: usize,
    pub candidate_location: usize,
    pub first_vert: Vertex,
	pub coeff: Coefficient,
    pub ring_operator: RingOperator,
}

/// implement standard methods of Iterator for IterCoboundary struct
impl< 'a, DisMat, FilVal, Coefficient, RingOperator >

    Iterator for 
    
    IterCoboundary
        < 'a, DisMat, FilVal, Coefficient, RingOperator >
        
    where
        FilVal:             Clone + Debug + Ord,
        Coefficient:             Clone,
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,    
{
	type Item = (SimplexFiltered<FilVal>, Coefficient);

    fn next( &mut self ) -> Option<Self::Item> {

        let candidacies = &self.clique.cutoff_matrix[usize::from(self.first_vert)]; // vertices we will try to insert
        let vertices = &mut self.next_cofacet_vertices;

        loop {
            if self.candidate_location >= candidacies.len(){
                return None;
            }
            let new_vertex = candidacies[self.candidate_location];
            let mut flag = false; // flag intidcates that we need to increment the insertion location
            vertices[self.insertion_location] = new_vertex;
            let mut new_filval = self.simplex_filtration.clone();
            for vert in vertices.iter() {
                // get the dissimilarity between the new vertex and the next-old-vertex
                let newdis_opt  =   self.clique.dissimilarity_matrix.entry( 
                                        new_vertex as usize, 
                                        *vert as usize
                                    );
                match newdis_opt {
                    None            =>  {
                                            flag = true; // this pair of vertices is non-adjacent, so skip to next insertion locus
                                            break;                        
                                        },
                    Some(newdis)    => {
                        if newdis > self.clique.dissimilarity_value_max {
                            flag = true; // simplex would be too big, so skip to next insertion locus
                            break;
                        } else if newdis > new_filval {
                            new_filval = newdis;
                        }
                    }
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
        //println!("new_filval={:?}, simplex_filtration={:?}", new_filval, self.simplex_filtration);
            let simp = SimplexFiltered{ vertices: vertices.clone(), filtration: new_filval };
            self.candidate_location += 1;
            return Some((simp, self.coeff.clone()));
        }
        //println!("candidate_location={:?}", self.candidate_location);
	}
}



/// Returns the coboundary of a simplex `keymaj` as a vector **with entries sorted in ascending order**
/// 
/// The order is the usual lexicographic order: first by birth time, second by lexicographic order.
fn get_coboundary_as_vec< DisMat, FilVal, Coefficient, RingOperator >( 
            matrix: & ChainComplexVrFiltered< DisMat, FilVal, Coefficient, RingOperator >, 
            keymaj: SimplexFiltered<FilVal> 
        ) 
        -> Vec< ( SimplexFiltered<FilVal>, Coefficient ) >
        
    where
        FilVal:             Clone + Debug + PartialOrd + Ord, 
        Coefficient:             Clone, 
        RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient >,
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,        
{

    let mut vertices = keymaj.vertices.clone();
    let first_vert = vertices[0];
    vertices.insert(0,0);
    vertices.shrink_to_fit();
    let iter = 
    IterCoboundary{
        clique: &matrix,
        next_cofacet_vertices: vertices,
        simplex_filtration: keymaj.filtration.clone(),
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



//  ===========================================================
//  FACET ITERATORS
//  ===========================================================


/// A iterator that runs over all entries in the boundary of a simplex, in **descending lexicographic order**
/// 
/// If the original simplex is `[0,1,2]`, it return entries in the following order: `([1,2],1), ([0,2],-1), ([0,1],1)`.
struct IterBoundary
        < 'a, DisMat, FilVal, Coefficient, RingOperator >
    where
        FilVal:                         Clone + Debug,
        Coefficient:                         Clone,
        DisMat:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                                        ViewRowAscend + 
                                        MatrixEntry,
        DisMat::EntryMajor:   KeyValGet< usize, FilVal >,
        RingOperator:                   Semiring< Coefficient > + Ring< Coefficient >, 
{
    pub clique: &'a ChainComplexVrFiltered< DisMat, FilVal, Coefficient, RingOperator >,
	pub simp: SimplexFiltered<FilVal>,
	pub removal_location: usize,
    pub ring_operator: RingOperator,
}

/// implement standard methods of Iterator for IterBoundary struct
impl< DisMat, FilVal, Coefficient, RingOperator >

    Iterator for 
    
    IterBoundary
        < '_, DisMat, FilVal, Coefficient, RingOperator >
    where
        FilVal:                         Clone + Debug + Ord + PartialOrd,
        Coefficient:                         Clone,
        DisMat:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                                        ViewRowAscend + 
                                        MatrixEntry,
        DisMat::EntryMajor:   KeyValGet< usize, FilVal >,
        RingOperator:                   Semiring< Coefficient > + Ring< Coefficient >, 
{
	type Item = (SimplexFiltered<FilVal>, Coefficient);

	fn next( &mut self ) -> Option<Self::Item> {
        if self.simp.vertices.len() == 1 { return None } // the boundary of a 0-simplex is empty
        if self.removal_location == self.simp.vertices.len() { return None; }
        
        let mut simplex = self.simp.clone();
        simplex.vertices.remove(self.removal_location);
        simplex.vertices.shrink_to_fit();
        simplex.filtration = self.clique.diameter(&simplex.vertices).unwrap();
        let coeff = self.ring_operator.minus_one_to_power( self.removal_location );
        self.removal_location += 1;
        return Some((simplex, coeff));

        // if let Some(diam) = self.clique.diameter(&simplex.vertices) {
        //     simplex.filtration = diam;
        //     self.coeff = self.ring_operator.minus_one_to_power( self.removal_location );
        //     self.removal_location += 1;
        //     return Some((simplex, self.coeff.clone()));
        // } else {
        //     return None;
        // }
    }
}




//  ---------------------------------------------------------------------------
//  IMPLEMENT MATRIX ORACLE CLIQUE BOUNDARY MATRIX + DEFINE SOME SORTED (CO)BOUNDARY ITERATORS
//  ---------------------------------------------------------------------------


//  INDICES AND COEFFICIENTS
//  ------------------------------------------


impl < DisMat, FilVal, Coefficient, RingOperator >
    
    IndicesAndCoefficients for    

    ChainComplexVrFiltered
        < DisMat, FilVal, Coefficient, RingOperator >

    where
        FilVal:             Clone + Debug + PartialOrd + Ord, 
        Coefficient:             Clone, 
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,         
{
    type RowIndex = SimplexFiltered<FilVal>; type ColIndex = SimplexFiltered<FilVal>; type Coefficient = Coefficient;
}


//  ORACLE MAJOR ASCEND
//  ------------------------------------------


impl < 'a, DisMat, FilVal, Coefficient, RingOperator >
    
ViewRowAscend for 

    &'a ChainComplexVrFiltered
        < DisMat, FilVal, Coefficient, RingOperator >

    where
        FilVal:                         Clone + Debug + Ord + PartialOrd,
        Coefficient:                         Clone,
        DisMat:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                                        ViewRowAscend + 
                                        MatrixEntry,
        DisMat::EntryMajor:   KeyValGet< usize, FilVal >,
        RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >,       
{
    type ViewMajorAscend            =   LazyOrderedCoboundary< 'a, DisMat, FilVal, Coefficient, RingOperator >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;
    type EntryMajor =   ( Self::RowIndex, Self::Coefficient );

    fn view_major_ascend( &self, keymaj: Self::RowIndex ) -> Self::ViewMajorAscend {

        let first_v = keymaj.vertices[0];
        let mut shortcircuit_vertex = None;        


        for outer_v in self.cutoff_matrix[first_v as usize ].iter().cloned() {
            // search for a vertex, `outer_v` whose id number is strictly below that of all 
            // vertices inside the simplex, such that `outer_v` could be added to the simplex 
            // without increasing its diameter
            // NB: we specifically need the vertex with the *smallest possible id number*
            //     that satisfies this property

            if outer_v >= first_v { break }
            let diam = NoneGreater::from_val( keymaj.filtration() );
            let mut diam_ok = true;
            // exclude outer_v if simplex [outer_v] is born too late
            if self.dissimilarity_as_nonegreater( outer_v, outer_v ) > diam 
                { diam_ok = false; break }
            // exclude outer_v if it is too far from any vertex in the simplex
            for inner_v in keymaj.vertices.iter().cloned() {
                if self.dissimilarity_as_nonegreater( outer_v, inner_v ) > diam 
                    { diam_ok = false; break }
            }
            if diam_ok { shortcircuit_vertex=Some( outer_v ); break }
        }

        if let Some( outer_v ) = shortcircuit_vertex {
            // if the preceding for-loop returns a valid value for outer_v, then construct
            // only the first entry

            let mut first_cofacet = Vec::with_capacity( keymaj.vertices.len()+1 );
            first_cofacet.push( outer_v );
            for inner_v in keymaj.vertices.iter().cloned() { first_cofacet.push(inner_v) }
            let cofacet = SimplexFiltered{ vertices: first_cofacet, filtration: keymaj.filtration() };
            let first_entry = ( cofacet, RingOperator::one() );
            LazyOrderedCoboundary{
                facet:              keymaj,
                boundary_matrix:    self,
                coboundary:         vec![ first_entry ],
                next_index:         0,
                built_all:          false,
            }
        } else {
            // otherwise construct precompute *all* entries, and store them in a Vec in sorted order

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


/// Iterator returning entries of a coboundary in order
/// 
/// This iterator is "lazy" because it will try to stop building the coboundary after generating,
/// the first element, using the Morse optimization.
/// 
/// After building the first element it will generate all other elements all at once, place them in a sorted vector,
/// and return elements of the vector one at a time.
#[derive(Clone,Copy,Debug,new)]
pub struct LazyOrderedCoboundary
        <'a, DisMat, FilVal, Coefficient, RingOperator >
    where 
        FilVal:             Clone + Debug,
        Coefficient:             Clone,
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,        
{
    facet:              SimplexFiltered<FilVal>,
    boundary_matrix:    &'a ChainComplexVrFiltered < DisMat, FilVal, Coefficient, RingOperator >,
    coboundary:         Vec< (SimplexFiltered<FilVal>, Coefficient) >,
    next_index:         usize,
    built_all:          bool
}

impl < 'a, DisMat, FilVal, Coefficient, RingOperator >
    
    Iterator for
    
    LazyOrderedCoboundary
        < 'a, DisMat, FilVal, Coefficient, RingOperator >

    where 
        FilVal:             Clone + Debug + Ord,
        Coefficient:             Clone,
        RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient >,     
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,         
{
    type Item = (SimplexFiltered<FilVal>, Coefficient);

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


impl < 'a, DisMat, FilVal, Coefficient, RingOperator >
    
    ParetoShortCircuit
        < (SimplexFiltered<FilVal>, Coefficient) > for
    
    LazyOrderedCoboundary
        < 'a, DisMat, FilVal, Coefficient, RingOperator >

    where 
        FilVal:             Clone + Debug + Ord,
        Coefficient:             Clone,
        RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient >,   
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,                    
{
    fn pareto_short_circuit(& self) -> Option< (SimplexFiltered<FilVal>, Coefficient) > {
        if ! self.built_all { return Some( self.coboundary[0].clone() ) }
        return None
    }
} 


//  ORACLE MINOR DESCEND
//  ------------------------------------------


impl < 'a, DisMat, FilVal, Coefficient, RingOperator >
    
ViewColDescend for 

    &'a ChainComplexVrFiltered
        < DisMat, FilVal, Coefficient, RingOperator >

    where
        FilVal:             Clone + Debug + PartialOrd + Ord, 
        Coefficient:             Clone, 
        RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient >,
        DisMat:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                            ViewRowAscend + 
                            MatrixEntry,
        DisMat::EntryMajor:   KeyValGet< usize, FilVal >,                                                                 
{
    type ViewMinorDescend            =   Vec< ( SimplexFiltered<FilVal>, Coefficient ) >;
    type ViewMinorDescendIntoIter    =   std::vec::IntoIter<(SimplexFiltered<FilVal>, Coefficient)>;
    type EntryMinor       =   ( Self::RowIndex, Self::Coefficient );

    fn view_minor_descend( &self, keymin: Self::ColIndex ) -> Self::ViewMinorDescend {

        let iter = IterBoundary {
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


// /// Based, filtered chain complex implementing the [`ChainComplex`](ChainComplex) trait
// pub struct CliqueComplex<Coefficient: Clone, FilVal> {
//     pub ring_operator: RingMetadata<Coefficient>,
//     pub dissimilarity_matrix: Vec<Vec<FilVal>>,
//     pub dissimilarity_value_max: FilVal,
//     pub safe_homology_degrees_to_build_boundaries: Vec<usize>,
//     pub major_dimension: MajorDimension,
//     pub simplex_count: Vec<(usize,usize)>
// }


//  SIMPLEX ITERATOR
//  ------------------------------------------

/// An iterator that iterates all simplices of given dimension
#[derive(Clone,Copy,Debug,new)]
pub struct SimplexIter< FilVal, DissimilarityMatrix > {
    pub dissimilarity_matrix: DissimilarityMatrix,
    pub dissimilarity_matrix_size: usize,    
    pub dissimilarity_value_max: FilVal,
    pub filvec: Vec<FilVal>,
    pub vec: Vec<Vertex>,
    pub val: Vertex,
    pub loc: usize,
}

impl < FilVal, DissimilarityMatrix >

    SimplexIter
        < FilVal, DissimilarityMatrix > 

    where
        FilVal:                 Clone + Debug + Ord,
    {

    pub fn new( 
            dissimilarity_matrix:           DissimilarityMatrix,
            dissimilarity_matrix_size:      usize,
            dissimilarity_value_max:        FilVal,
            dim:                            isize,             
        ) 
        -> 
        Self 
        where
            DissimilarityMatrix:                            IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                                                            ViewRowAscend,
            DissimilarityMatrix::EntryMajor:      KeyValGet< usize, FilVal >,
    {

        let dim = dim as usize;
        let minval 
            = (0..dissimilarity_matrix_size)
                .map(|x| dissimilarity_matrix.view_major_ascend(x) )
                .flatten().map(|x| x.val() ).min().unwrap();
        SimplexIter {
            dissimilarity_matrix,
            dissimilarity_matrix_size,
            dissimilarity_value_max,
            filvec: vec![minval; dim + 2 ], // we will maintain an internal state of the iterator such that filvec[p] = diameter of the simplex spanned by the first (p-1) simplices; this is useful as it allows us to store a value for the -1 simplex
            vec: vec![0; dim +1 ], // the vertices of the simplex we're building
            val: 0, // the vertex we insert / have inserted
            loc: 0, // where in the simplex we insert / have inserted
        }
    }
}


/// implement standard methods of Iterator for SimplexIter struct
impl< FilVal, DissimilarityMatrix > 

    Iterator for 
    
    SimplexIter
        < FilVal, DissimilarityMatrix > 
        
    where
        FilVal: Clone + PartialOrd + Ord + Debug,
        DissimilarityMatrix:                            IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=FilVal > +
                                                        ViewRowAscend +
                                                        MatrixEntry,
        DissimilarityMatrix::EntryMajor:      KeyValGet< usize, FilVal >,
{
    type Item = SimplexFiltered<FilVal>;

    /// Get the next simplex
    /// 
    /// This iterator works on the following principles:
    /// - we iterate over simplices in lexicographic order
    /// - we add vertices in a loop
    ///   - suppose we are at a point in the loop where we have filled in the first 3 slots: [1 2 3 * * *]
    ///   - the next simplex we add will have vertices inserted in strictly ascending order, so we can essentially
    ///     do a for-loop over vertices 4... inserting the first sequence of vertices possible
    ///   - after we do that we can try to increment the last vertex
    ///   - then if the last vertex is at ceiling we can try to incremeent the second-to-last vertex
    ///   - etc.
    fn next(&mut self) -> Option<Self::Item> {
        let size = self.dissimilarity_matrix_size;
        if self.vec.len() > size { return None; }
        loop {
            while usize::from(self.val) <= size - (self.vec.len() - self.loc) {
                
                // write in a new vertex value
                self.vec[self.loc] = self.val; 

                // compute the diameter of the simplex created by adding a vertex just now
                self.filvec[self.loc+1] = self.filvec[self.loc].clone(); // initialize the diameter as the diameter of the previous simplex; our constructor function inializes the diameter of the emepty simplex as the minimum entry of the dissimilarity_matrix
                
                
                for ii in 0..self.loc+1 {
                    match self.dissimilarity_matrix.entry( self.val as usize, self.vec[ii] as usize ) {
                        Some(filval) => {
                            if self.filvec[self.loc+1]  <  filval {
                                self.filvec[self.loc+1] = filval; // increase the recorded diameter if necessary
                            }
                        }
                        None => { continue } // if no dissimilarity is recorded, then we assume this pair of vertices to be NON-adjacent
                    }
                }

                // if the simplex has the correct dimension and a sufficiently small diameter, then return it
                if self.filvec[self.loc+1] <= self.dissimilarity_value_max {
                    if self.loc == self.vec.len()-1 {
                        self.val += 1;
                        return Some(SimplexFiltered{
                                vertices: self.vec.clone(),
                                filtration: self.filvec[self.loc+1].clone()
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






//  ================================================================================
//  COBOUNDARY / COFACET ITERATORS
//  ================================================================================





// pub struct ShortcircuitVertices
//         < FilVal, I >
//     where
//         FilVal:             Clone + Debug + Ord,
//         I:                  Iterator,
//         I::Item:            KeyValGet< Vertex, FilVal >,
// {
//     iter:                   IntersectOrderedIterators< I, OrderOperatorByKey< Vertex, FilVal, I::Item > >,
//     filtration_max:         FilVal,
//     facet_vertices:         Vec< Vertex >,
//     facet_cardinality:      usize,
//     duplication_pointer:    usize,
// }        


// impl < FilVal, I >

//     Iterator for

//     ShortcircuitVertices
//             < FilVal, I >
//         where
//             FilVal:         Clone + Debug + Ord,
//             I:              Iterator,
//             I::Item:        KeyValGet< Vertex, FilVal >,
//     {
//         type Item = Vertex;

//         fn next(&mut self) -> Option<Self::Item> {

//             println!("modify the HitMerge operator to ");
            
//             let mut multiplicity = 0;
//             let mut vertex_new;
//             let mut vertex_new_opt = None;
//             let mut vertex_old_opt = None;            

//             loop {
//                 match self.iter.next() {
//                     None => { 
//                         return None; 
//                     } Some( edge ) => {
//                         match edge.val() < self.filtration_max {
//                             false => { 
//                                 continue;
//                             } true => {
//                                 vertex_new      = edge.key();
//                                 vertex_new_opt  = Some( vertex_new );
                                
//                                 match vertex_new_opt == vertex_old_opt {
//                                     // check if the new vertex matches the old vertex
//                                     true => { 
//                                         // if it does then increment our mulitplicity counter

//                                         multiplicity += 1;
//                                         if multiplicity == self.facet_cardinality {
//                                             // if we've checked every vertex in the facet, then return the vertex
//                                             return Some( vertex_new )
//                                         }


//                                     } false => {

//                                         // update the duplication pointer, so that it points at the vertex that 
//                                         // vertex_new might duplicate
//                                         while ( self.duplication_pointer < self.facet_cardinality )
//                                             &&
//                                             (   vertex_new < self.facet_vertices[ self.duplication_pointer ] ) {
//                                             self.duplication_pointer += 1;
//                                         }
//                                         if vertex_new == self.facet_vertices[ self.duplication_pointer ] {
//                                             // if vertex_new duplicates a vertex in the facet, then skip it
//                                             // and move on
//                                             vertex_old_opt = None;
//                                             continue;
//                                         } else {
//                                             // otherwise use it to initialize a new search
//                                             multiplicity = 1; 
//                                             vertex_old_opt = vertex_new_opt;
//                                         }
//                                     }
//                                 }                                
//                             }                            
//                         }
//                     }
//                 }
//             }
//         }
// }






// /// Iterates over pairs (v,f) such that v can be added to a simplex s, producing a cofacet of diameter f.
// /// 
// /// Pairs are returned in strictly sorted order, first by filtration (ascending) and second by lexicographic order of vertices in v.
// /// 
// /// 
// /// How it works
// /// 
// /// - first run through neighbors in lexicographic order to get equal-diameter cofacets; this must run through the end of the nbr_iters_ordlex in order to respect order
// /// - second run through neighbors in filtration order
// pub struct CofacetEdgeIter
//                 < FilVal, I, J >
//     where
//         I:      Iterator< Item=(Vertex,FilVal) > + PeekUnqualified,
//         J:      Iterator< Item= Vertex >,
// {
//     facet_fil:                              FilVal,
//     neighbor_iters_filtration_order:        Vec< I >,
//     neighbors_within_diameter_lex_order:    J,
// }



// impl < FilVal, I, J >

//     Iterator for 

//     CofacetEdgeIter
//         < FilVal, I, J >
    
//     where
//         I:                  Iterator< Item=(Vertex,FilVal) > + PeekUnqualified,
//         J:                  Iterator< Item=Vertex >,
//         FilVal:             Clone + Ord,
// {    
//     type Item = ( Vertex, FilVal );

//     fn next(&mut self) -> Option<Self::Item> {

//         println!("J is a merge of several iterators; design J so that it stops when one of its internal iterators is empty, and returns only vertices which meet diameter bounds for all vertices of the facet, and excludes vertices in the facet");

//         // first check if we can generate a cofacet of equal diameter
//         // ----------------------------------------------------------
//         let a = self.neighbors_within_diameter_lex_order.next();
//         if a.is_some() { return Some( ( a.unwrap(), self.facet_fil.clone() ) ) }

//         // if not, then start iterating over edges incident to the vertices in the facet
//         // -----------------------------------------------------------------------------
        
//         // initialize some variables
//         let facet_len = self.neighbor_iters_filtration_order.len();
//         let mut iter;
//         let mut ii_old = 0;
//         let mut ii_new = 0;
//         let e_old_opt = self.neighbor_iters_filtration_order[ii_old].next();
//         match e_old_opt {
//             None => { 
//                 // in this case we've exhausted all neighbors the ii_old'th vertex in our facet, so no more cofacets exist
//                 return None 
//             } Some( mut e_old ) => {
//                 // in this case we have an edge e_old of length greater than the diameter of the facet
//                 // the edge connects the ii_old'th vertex in the facet to another vertex, v_other.                  
                
//                 loop {
//                     // in this loop we check to see if v_other is (1) incident to every other vertex in 
//                     // the facet, and (2) no farther than than length-of-e_old; if not, then pick a
//                     // new edge to try to add as the longest edge of a new cofacet, and start again

//                     // increment to the next vertex in the facet
//                     ii_new = ( ii_new + 1 ) % facet_len;

//                     // if we've successfully checked every vertex of the facet and returned to the first one 
//                     // that we originally checked, then return the new edge!
//                     if ii_new == ii_old { return Some(e_old) }

//                     // otherwise we have to do some work
//                     iter = &mut self.neighbor_iters_filtration_order[ ii_new ];                    
//                     loop {
//                         // check that v_other lies within distance length-of-e_old of the ii_new'th 
//                         // vertex in the facet; if not, then either (i) identify a different edge
//                         // to add, or (ii) delare that we have finished generating all cofacets
//                         match iter.next() {
//                             None => { return None } // iterator is emtpy
//                             Some( e_new ) => {
//                                 if e_new.1 > e_old.1 {
//                                     // in this case e_old can't be added as the longest edge in a new facet
//                                     if e_new.0 != e_old.0 {
//                                         // if e_new has a different "outer" vertex than e_old, then we have
//                                         // to re-start the process of checking that it's a viable choice
//                                         // for every vertex in the facet
//                                         ii_old = ii_new % facet_len;
//                                         ii_new = ii_old;
//                                     }
//                                     // record that we are now trying to add e_new as the longest edge
//                                     e_old = e_new;     
//                                     // now move on to the next neighbor iterator                               
//                                     break;
//                                 } else if e_new.0 == e_old.0 {
//                                     // in this case e_old still has a chance of being added as (one of the)
//                                     // longest edges in a new cofacet     
//                                     break;
//                                 }
//                             }
//                         }

//                     }
//                 }
//             }
//         }
//     }
// }   




// /// Iterates over cofacets in strictly ascending order of filtration value, with ties broken by lexicographic order.
// pub struct CofacetIter
//                 < FilVal, I, J >
//     where
//         I:      Iterator< Item=(Vertex,FilVal) > + PeekUnqualified,
//         J:      Iterator< Item= Vertex >,
// {
//     cofacet_edge_iter:      CofacetEdgeIter< FilVal, I, J >,
//     facet_vertices:         Vec< Vertex >,
// }

// impl < FilVal, I, J >

//     Iterator for 

//     CofacetIter
//         < FilVal, I, J >
    
//     where
//         I:                  Iterator< Item=(Vertex,FilVal) > + PeekUnqualified,
//         J:                  Iterator< Item=Vertex >,
//         FilVal:             Clone + Debug + Ord,
// {    
//     type Item = SimplexFiltered< FilVal >;

//     fn next(&mut self) -> Option<Self::Item> {

//         match self.cofacet_edge_iter.next() {
//             None => {
//                  return None 
//             } Some( (vertex, filtration_value) ) =>  {
//                 let facet_cardinality = self.facet_vertices.len();
//                 let mut insertion_locus = facet_cardinality;
//                 for p in 0 .. facet_cardinality {
//                     if vertex < self.facet_vertices[p] {
//                         insertion_locus = p;
//                         break;
//                     }
//                 }
//                 let mut cofacet_vertices = self.cofacet_edge_iter.facet_vertices.clone();
//                 cofacet_vertices.insert( insertion_locus, vertex );
//                 return Some( SimplexFiltered{ vertices: cofacet_vertices, filtration: filtration_value } );
//             }
//         }        
//     }
// }

// /// Iterates over entries of the coboundary in strictly ascending order of filtration value.
// pub struct CoboundaryIter
//                 < RingOperator, Coefficient, FilVal, I, J, >
//     where
//         I:                  Iterator< Item=(Vertex,FilVal) > + PeekUnqualified,
//         J:                  Iterator< Item= Vertex >,
//         RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
// {
//     cofacet_edge_iter:      CofacetEdgeIter< FilVal, I, J >,
//     facet_vertices:         Vec< Vertex >,
//     ring_operator:          RingOperator,
//     phantom_snzval:         PhantomData<Coefficient>,
// }

// impl < RingOperator, Coefficient, FilVal, I, J, >

//     Iterator for 

//     CoboundaryIter
//         < RingOperator, Coefficient, FilVal, I, J, >
    
//     where
//         I:                  Iterator< Item=(Vertex,FilVal) > + PeekUnqualified,
//         J:                  Iterator< Item=Vertex >,
//         FilVal:             Clone + Debug + Ord,
//         RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
// {    
//     type Item = ( SimplexFiltered< FilVal >, Coefficient );

//     fn next(&mut self) -> Option<Self::Item> {

//         match self.cofacet_edge_iter.next() {
//             None => {
//                  return None 
//             } Some( (vertex, filtration) ) =>  {
//                 let facet_cardinality = self.facet_vertices.len();
//                 let mut insertion_locus = facet_cardinality;
//                 for p in 0 .. facet_cardinality {
//                     if vertex < self.facet_vertices[p] {
//                         insertion_locus = p;
//                         break;
//                     }
//                 }
//                 let mut vertices = self.facet_vertices.clone();
//                 vertices.insert( insertion_locus, vertex );
//                 let cofacet = SimplexFiltered{ vertices: vertices, filtration };
//                 let coefficient = self.ring_operator.minus_one_to_power( insertion_locus );
//                 return Some( ( cofacet, coefficient ) );
//             }
//         }        
//     }
// }



























//  =========================================================================
//  UNIT TESTS
//  =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    use itertools::Itertools;
    use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    use crate::topology::point_cloud::unit_circle;    
use crate::algebra::rings::operator_structs::ring_native::{DivisionRingNative, FieldRational64};    
    use crate::algebra::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
    use crate::algebra::matrices::query::{ViewColDescend, ViewRowAscend};
    use crate::algebra::matrices::operations::umatch::row_major::{Umatch::factor_with_clearing, ParetoShortCircuit};    
    use crate::algebra::matrices::types::third_party::IntoCSR;
use crate::utilities::order::{ is_sorted_strictly, OrderOperatorByKeyCutsom, OrderOperatorAutoReverse};
use crate::utilities::random::random_matrix;
use crate::utilities::distances::{rowwise_distances};
  
    
    #[test]
    fn check_that_some_basic_functions_run_without_error() {

        let m = 20;
        let dimension_max = 1;
        let dissimilarity_value_max = None;
        let dissimilarity_value_empty = OrderedFloat(0.0);

        let pcloud = unit_circle( m, Some(-1.0 .. 1.0));

        let dissimilarity_matrix_data
            = rowwise_distances(pcloud)
                .into_iter()
                .map(|x| x.into_iter().enumerate().collect_vec() )
                .collect_vec()
                .into_csr( m, m );
        let dissimilarity_matrix = & dissimilarity_matrix_data;
    
        let ring_operator = FieldRational64::new();
        let chain_complex = ChainComplexVrFiltered::new( & dissimilarity_matrix, m, dissimilarity_value_max, dissimilarity_value_empty, ring_operator.clone() );
        let chain_complex_ref = & chain_complex;    
        let keymaj_vec = chain_complex.cliques_in_order(dimension_max);
        let keymin_vec = chain_complex.cliques_in_order(dimension_max+1);
    
    
        verify_viewmajorascend_compatible_with_viewminordescend(
                chain_complex_ref,
                keymin_vec.iter().cloned(),
                keymaj_vec.iter().cloned(),
            );       
    
        let iter_keymaj = keymaj_vec.iter().cloned();    
    
        println!("check that oracle has strictly sorted rows");
        // print_indexed_major_views( & chain_complex_ref, iter_keymaj.clone() );  // print the major views       
        for keymaj in iter_keymaj.clone() {        
            assert!( is_sorted_strictly( 
                                            &chain_complex_ref.view_major_ascend(keymaj.clone()).collect_vec() , 
                                            &OrderOperatorByKeyCutsom::new( OrderOperatorAuto ) 
                                        ) );
        }
    
        // println!("press enter to continue");
        // let mut guess = String::new();
        // io::stdin()
        //     .read_line(&mut guess)
        //     .expect("Failed to read line");        
    
        println!("check that oracle has strictly sorted columns");
        // print_indexed_minor_views( & chain_complex_ref, iter_keymaj.clone() );  // print the major views        
        for keymaj in iter_keymaj.clone() {
            assert!( is_sorted_strictly(    &chain_complex_ref.view_minor_descend(keymaj).iter().cloned().collect_vec() , 
                                            &OrderOperatorByKeyCutsom::new( OrderOperatorAutoReverse::new() )  // NOTE THAT HERE WE USE GT 
                                        ) );
        }    
    
        // println!("press enter to continue");
        // let mut guess = String::new();
        // io::stdin()
        //     .read_line(&mut guess)
        //     .expect("Failed to read line");       
    
        println!("starting umatch");
        let umatch = Umatch::factor_with_clearing(
                chain_complex_ref, 
                iter_keymaj.clone(), 
                ring_operator.clone(), 
                OrderOperatorAuto, 
                OrderOperatorAuto, 
            );      
    
        // println!("start build bd matrix");
        // let chain_complex = oat_rust::utilities::homology::clique::get_clique_chain_complex(
        //     dissimilarity_matrix,
        //     dissimilarity_value_max,
        //     FieldRational64::new()        
        // );
    
        // println!("start umatch");    
        // let umatch = oat_rust::utilities::homology::clique::umatch_from_clique(
        //     & chain_complex,
        //     dimension_max,
        // );
    
           
        // let iter_keymaj = 
        //     (0..dimension_max+1).map(
        //         |dim|
        //         {
        //             let mut vec = SimplexIter::new( 
        //                         dim,
        //                         & umatch.mapping_ref().dissimilarity_matrix,
        //                         umatch.mapping_ref().dissimilarity_value_max,                
        //                     )
        //                     .collect_vec();
        //             vec.sort_by(|x,y| y.filtration.cmp(&x.filtration));
        //             vec
        //         }
        //     )
        //     .flatten();    
    
    
        println!("setting up to unpack");  
        let dim_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.dimension() as isize;
        let fil_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.filtration();    
        let barcode = crate::barcode::barcode( &umatch, iter_keymaj, dim_fn, fil_fn, true , true);
    
        println!("getting intervals, reps, bounding chains");
        // let (intervals, representatives, bounding_chains ) = barcode.unwrap();
    
        // let convert_coefficients = |x: Ratio<i64> | (x.numer().clone(), x.denom().clone());
        // let convert_simplex = |x: SimplexFiltered< OrderedFloat<f64> > | x.vertices.iter().map(|y| y.clone() as usize ).collect_vec();
        
        // println!("start reformat reps");     
        // let represntatives_new = representatives.unwrap().iter().map(
        //             |x| // each x is a list of cycle reps
        //             x.iter().map(
        //                 |y| // each y is a cycle rep
        //                 y.iter().map(
        //                     |z| // each z is an entry in a cycle rep
        //                     ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
        //                 ).collect_vec()
        //             ).collect_vec()
        //         ).collect_vec();
    
        // println!("start reformat bounding chains");                 
        // let bounding_chains_new = bounding_chains.unwrap().iter().map(
        //             |x| // each x is a list of cycle reps
        //             x.iter().map(
        //                 |y| // each y is a cycle rep
        //                 y.iter().map(
        //                     |z| // each z is an entry in a cycle rep
        //                     ( convert_simplex(z.0.clone()), convert_coefficients(z.1.clone()) )
        //                 ).collect_vec()
        //             ).collect_vec()
        //         ).collect_vec();              
    }
}    