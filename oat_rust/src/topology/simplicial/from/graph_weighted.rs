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





use std::sync::Arc;
use std::{marker::PhantomData};

use std::hash::Hash;
use std::fmt::Debug;

use derive_new::new;
use itertools::Itertools;

use crate::topology::simplicial::simplices::filtered::{SimplexFiltered, OrderOperatorTwistSimplexFiltered};

use crate::algebra::vectors::entries::KeyValGet;
use crate::topology::simplicial::simplices::unfiltered::Simplex;


// use crate::matrix::{SmOracle, RingMetadata, MajorDimension};
// use crate::chx::{ChainComplex, ChxTransformKind};

use std::cmp::Ordering;

use crate::algebra::chains::factored::{factor_boundary_matrix, FactoredBoundaryMatrix};
use crate::algebra::matrices::operations::umatch::row_major::{ParetoShortCircuit,};
use crate::algebra::matrices::query::{IndicesAndCoefficients, ViewRowAscend, ViewColDescend, MatrixEntry};
use crate::algebra::rings::operator_traits::DivisionRing;
use crate::algebra::rings::operator_traits::{Semiring, Ring};

use crate::utilities::order::{ OrderOperatorByKeyCutsom, JudgeOrder, NoneGreater, OrderOperatorAuto};
use crate::utilities::iterators::general::{minmax, IntersectOrderedIterators, PeekUnqualified};
// use oat_rust::utilities::indexing_and_bijection::sparsevec_lookup;



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
type OrderOperatorSimplexFiltered< Filtration, Coefficient > =   OrderOperatorByKeyCutsom< 
                                                                SimplexFiltered< Filtration >, 
                                                                Coefficient,
                                                                (
                                                                    SimplexFiltered< Filtration>, 
                                                                    Coefficient,
                                                                ),
                                                                OrderOperatorAuto
                                                            >;                                            



//  ===========================================================
//  TRAIT FOR DISSIMILARITY MATRIX
//  ===========================================================






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
    let _a = minmax(u);

    let w = vec![ vec![0,1], vec![0,1] ];
    let z = w.iter();
    let _a = minmax(z);
}


/// A vector of simplices sorted first by dimension (ascending `0..dimension_max`, including `dimension_max`) 
/// then by diameter (descending) then by lexicographic order (descending)
pub fn filtered_cliques_in_order< Filtration, DissimilarityMatrix >( 
            dimension_max:                  isize,
            dissimilarity_value_max:        Filtration,
            dissimilarity_matrix:           DissimilarityMatrix,
            dissimilarity_matrix_size:      usize,            
        ) 
        -> 
        Vec< SimplexFiltered< Filtration > >  
    where
        Filtration:                         Clone + Debug + Ord,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > + ViewRowAscend + MatrixEntry,
        DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
    {
        let order_operator = OrderOperatorTwistSimplexFiltered{};
        let dissimilarity_matrix_ref = & dissimilarity_matrix;
        let vec: Vec<_> = (0..dimension_max+1).flat_map(|dim|
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
            })
        .collect();
        vec
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
            < 'a, Filtration, DissimilarityMatrix, DisMatEntry >
            (dissimilarity_matrix: DissimilarityMatrix, dissimilarity_value_max: Filtration) 
            -> 
            Vec<Vec<Vertex>> 
    where
        Filtration:             'a + PartialOrd + Debug,
        DissimilarityMatrix:             IntoIterator,
        DissimilarityMatrix::Item:       IntoIterator< Item=DisMatEntry >,
        DisMatEntry:        KeyValGet< usize, Filtration >,

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
    output
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

// pub struct ChainComplexVrFilteredUmatch<Filtration, Coefficient, RingOperator>
//     where 
//         Coefficient:         Clone,
//         RingOperator:   Semiring< Coefficient > + Ring< Coefficient >,
//         Filtration:         Hash,
// {
//     umatch:     Umatch<  
//                         ChainComplexVrFiltered<Filtration, Coefficient, RingOperator>,
//                         RingOperator, 
//                         OrderOperatorSimplexFiltered,
//                         OrderOperatorSimplexFiltered,
//                     >,
//     dimension_max:     isize,

// }        

//  ===========================================================
//  FILTERED VIETORIS RIPS BOUNDARY MATRIX 
//  ===========================================================

/// # Examples 
/// Boundary matrix represented by a matrix oracle.
/// 
/// 
/// ```
/// // Import modules
/// 
/// use oat_rust::topology::simplicial::{simplices::filtered::SimplexFiltered, from::graph_weighted::{ChainComplexVrFiltered}};        
/// use oat_rust::topology::point_cloud::unit_circle;    
/// use oat_rust::algebra::vectors::entries::KeyValGet;
/// use oat_rust::algebra::rings::operator_structs::ring_native::{FieldRational64};    
/// use oat_rust::algebra::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
/// use oat_rust::algebra::matrices::query::{ViewColDescend, ViewRowAscend};
/// use oat_rust::algebra::matrices::operations::umatch::row_major::{Umatch};    
/// use oat_rust::algebra::matrices::types::third_party::IntoCSR;
/// use oat_rust::utilities::order::{ is_sorted_strictly, OrderOperatorByKeyCutsom, OrderOperatorAuto, OrderOperatorAutoReverse};
/// use oat_rust::utilities::iterators::general::minmax;    
/// use oat_rust::utilities::distances::{rowwise_distances};
/// 
/// use std::sync::Arc;
/// use ordered_float::OrderedFloat;  
/// use itertools::Itertools;      
/// 
/// 
/// // Define a point cloud (points sampled from a unit circle)
/// 
/// let npoints = 20;
/// let pcloud = unit_circle( npoints, Some(-1.0 .. 1.0));
/// 
/// // Get the distance matrix
/// 
/// let dissimilarity_matrix_data
///     = rowwise_distances(pcloud)
///         .into_iter()
///         .map(|x| x.into_iter().enumerate().collect_vec() )
///         .collect_vec()
///         .into_csr( npoints, npoints );
/// let dissimilarity_matrix = & dissimilarity_matrix_data;
/// 
/// // Determine the upper and lower bounds for filtration values we will consider
/// 
/// let dissimilarity_value_min = OrderedFloat(0.0);        
/// let dissimilarity_value_max = 
///     minmax( 
///             (0..npoints).map(
///                 |x| 
///                 dissimilarity_matrix.view_major_ascend(x).into_iter().map(
///                     |x| 
///                     x.val()
///                 ) 
///             ) 
///     ).unwrap_or( dissimilarity_value_min.clone() 
/// );         
/// 
/// // Define the coefficient ring
/// 
/// let ring_operator = FieldRational64::new();
/// 
/// // Define the chain complex
/// 
/// let chain_complex_data = ChainComplexVrFiltered::new( & dissimilarity_matrix, npoints, dissimilarity_value_max, dissimilarity_value_min, ring_operator );
/// let chain_complex = Arc::new(chain_complex_data);    
/// 
/// let dimension_max = 2;       
/// let keymaj_vec = chain_complex.cliques_in_order(dimension_max); // runs over the row indices we want to use when computing the factorization
/// let keymin_vec = chain_complex.cliques_in_order(dimension_max+1);
/// 
/// 
/// 
/// let iter_keymaj = keymaj_vec.iter().cloned();    
/// 
/// // Check that oracle has strictly sorted rows
/// 
/// for keymaj in iter_keymaj.clone() {        
///     assert!( is_sorted_strictly( 
///                                     & chain_complex.view_major_ascend(keymaj.clone()).collect_vec() , // this line accesses a row 
///                                     & OrderOperatorByKeyCutsom::new( OrderOperatorAuto ) 
///                                 ) );
/// }      
/// 
/// // Check that oracle has strictly sorted columns
/// 
/// for keymaj in iter_keymaj.clone() {
///     assert!( is_sorted_strictly(    & chain_complex.view_minor_descend(keymaj).iter().cloned().collect_vec() ,  // this line accesses a column
///                                     & OrderOperatorByKeyCutsom::new( OrderOperatorAutoReverse::new() )  
///                                 ) );
/// }    
///       
/// 
/// // Get a Umatch factorization
/// 
/// let umatch = Umatch::factor_with_clearing(
///         chain_complex.clone(), 
///         iter_keymaj.clone(), 
///         ring_operator, 
///         OrderOperatorAuto, 
///         OrderOperatorAuto, 
///     );         
/// 
/// // Get the barcode
/// 
/// let dim_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.dimension() as isize;
/// let fil_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.filtration();    
/// let _barcode = oat_rust::algebra::chains::barcode::barcode( &umatch, iter_keymaj, dim_fn, fil_fn, true , true);
/// ```

#[derive(Clone,Debug,)]
pub struct ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    where 
        Coefficient:         Clone,
        RingOperator:   Semiring< Coefficient > + Ring< Coefficient >,
        DissimilarityMatrix:         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                        ViewRowAscend + 
                        MatrixEntry,
{
    pub ring_operator: RingOperator, // ring meta data
    /// The "dissimilarity_value_max" value represents the maximum of filtration value
    pub dissimilarity_value_max: Filtration,
    /// The dissimilarity matrix
    pub dissimilarity_value_min: Filtration,
    /// The dissimilarity score (or "diameter") of the emtpy simplex
	pub dissimilarity_matrix: DissimilarityMatrix,
    /// The number of rows (respectively, columns) of the dissimilarity matrix
    pub dissimilarity_matrix_size: usize,
    /// A vector representing the neighborhood within "dissimilarity_value_max" of each vertex
    pub cutoff_matrix: Vec<Vec<Vertex>>,
    pub phantomsnzval: PhantomData< Coefficient >,
}


/// Methods of ChainComplexVrFiltered struct
impl < DissimilarityMatrix, Filtration, Coefficient, RingOperator  > 

    ChainComplexVrFiltered
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator > 

    where
        Filtration:                         Clone + Debug + PartialOrd + Ord,
        Coefficient:                    Clone,
        RingOperator:                   Semiring< Coefficient > + Ring< Coefficient >,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                        ViewRowAscend + 
                                        MatrixEntry, 
        DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,       
{

    /// Construct a (lazy) boundary matrix for the filtered VR complex of a symmetric matrix.
    /// - `dissimilarity_matrix` a symmetric matrix representing pairwise dissimilarity values
    /// - `dissimilarity_matrix_size` an integer representing the number of rows/columns of the dissimilarity matrix; this parameter is necessary because the matrix oracle traits are passive, and do not give a means to directly count or generate the number of row/column indices
    /// - `dissimilarity_value_max` maximum dissimilarity threshold
    /// - `dissimilarity_value_min` diameter of the empty simplex; this can be any number that is equal or less than every structural nonzero entry in the dissimilarity matrix
    /// - `ring_operator`: ring operator for the coefficient ring
    /// The filtration parameter of
    /// - the empty simplex `[]` is `dissimilarity_value_min`
    /// - a singletime `[i]` is `dissimilarity_matrix[i,i]`
    /// - an edge `[i,j]` is `dissimilarity_matrix[i,j]`
    /// - an arbitrary simplex `s` is the maximum of the filtration parameters of cardinality 0, 1, and 2 subsets of `s`
    /// subject to the condition that the filtration parameter of
    /// - a singleton `[i]` is undefined if `dissimilarity_matrix[i,i]` is structurally zero
    /// - an edge `[i,j]` is undefined if `dissimilarity_matrix[i,j]` is structurally zero
    /// - an arbitrary simplex `s` is undefined if the filtration parameter of any cardinality 0, 1, or 2 subset of `s` is undefined
    /// This means that some vertices may never enter the filtration.
    /// 
    /// # Safety checks
    /// 
    /// The constructor checks that 
    /// - the input matrix is symmetric
    /// - the filtration parameter of every edge equals or exceeds the filtration parameter of its vertices.
    pub fn new(
            dissimilarity_matrix:                         DissimilarityMatrix,
            dissimilarity_matrix_size:                    usize,
            dissimilarity_value_max:        Filtration,
            dissimilarity_value_min:        Filtration,
            ring_operator:                  RingOperator,
        ) -> Self
    {

        let npoints = dissimilarity_matrix_size;
        let disemp_nonegreater = NoneGreater::from_val( dissimilarity_value_min.clone() );

        for i in 0 .. npoints {
            for entry in dissimilarity_matrix.view_major_ascend(i) {   
                let j = entry.key(); 
                let v = entry.val();             
                if Some( v ) != dissimilarity_matrix.entry_major_at_minor( j, i ) {
                    panic!("\n\nError: The dissimilarity matrix passed to the `ChainComplexVrFiltered` constructor is not symmetric: entry {:?} doesn't equal entry {:?}.\nThis message is generated by OAT.\n\n", (i,j), (j,i));                    
                }
            }

            // check that the filtration function respects face relations; this means that diagonal elements
            // (1) are minimal in their respective rows
            // (2) equal or exceed the diameter of the empty simplex
            let min = dissimilarity_matrix.view_major_ascend(i).into_iter().map(|x| x.val() ).min();    // min entry this row
            let dia = dissimilarity_matrix.entry_major_at_minor(i,i);  
            if min != dia {
                panic!("\n\nError: Entry ({:?},{:?}) of the dissimilarity matrix passed to the `ChainComplexVrFiltered` constructor is {:?} but the minimum structural nonzero entry in row {:?} is {:?}.  In this case `None` represents a value greater strictly greater than `Some(x)`, for every filtration value `x`.\nThis message is generated by OAT.\n\n", i, i, &dia, i, &min );
            }
            if NoneGreater::new(min.clone() ) < disemp_nonegreater { // we wrap filtration values in `NoneGreater` so that None values will be regarded as maximal instead of minimal
                panic!("\n\nError: The dissimilarity matrix passed to the `ChainComplexVrFiltered` constructor assigns  ({:?},{:?}) a value of {:?}, which is strictly less than the value {:?} assigned to the empty simplex.  In this case `None` represents a value greater strictly greater than `x`, for every filtration value `x`.\nThis message is generated by OAT.\n\n", i, i, min, dissimilarity_value_min );
            }            
        }        
    
        // -----------------------------------------------------------------------------
        // THIS BLOCK IS NOW DEPRECATED
        // It made sense when we had dense dissimilarity matrices, but now that the input
        // is a sparse matrix, this block will return the correct enclosing radius ONLY
        // if there exists a vertex that is adjacent to every other vertex.  This is too
        // much of a complication; we'll just ask the user to pass the matrix that they
        // want.
        //
        // if no dissimilarity_value_max is provided then use the enclosing radius
        // let dissimilarity_value_max = dissimilarity_value_max.unwrap_or( 
        //             minmax( 
        //                     (0..dissimilarity_matrix_size).map(
        //                             |x| 
        //                             dissimilarity_matrix.view_major_ascend(x).into_iter().map(
        //                                     |x| 
        //                                     x.val()
        //                                 ) 
        //                         ) 
        //                 )
        //                 .unwrap_or( dissimilarity_value_min.clone() )
        //         );
        // -----------------------------------------------------------------------------
        let cutoff_matrix = get_cutoff_matrix( (& dissimilarity_matrix).views_major_ascend(0..dissimilarity_matrix_size), dissimilarity_value_max.clone() );
        
        ChainComplexVrFiltered { ring_operator, dissimilarity_value_max, dissimilarity_matrix, dissimilarity_matrix_size, cutoff_matrix, dissimilarity_value_min, phantomsnzval: PhantomData }
    }            

    /// Output the diameter of a simplex
    ///
    /// # Parameters
    /// -`vertices`: vertices of the simplex
    /// # Returns
    /// Return 
    /// - `None` if any pairwise distance is not recorded in the
    /// - `self.dissimilarity_value_min` if the `vertices` is empty
    /// - the maximum of the pairwise dissimilarities, otherwise
    pub fn diameter( & self, vertices: & Vec< Vertex > ) -> Option< Filtration > {
        let mut diam = self.dissimilarity_value_min.clone(); 
        let mut a;
        let mut b;
        for ii in 0..vertices.len() {
            a = usize::from( vertices[ii] ); 
            for jj in ii .. vertices.len() {
                b = usize::from( vertices[jj] ); 
                match self.dissimilarity_matrix.entry_major_at_minor( a, b ) {
                    None => { return None } // the simplex never enters the filtration
                    Some( diam_bound ) => {
                        if diam_bound > diam { diam = diam_bound.clone(); }
                    }
                }
            }
        }
        Some(diam)
    }    

    /// Returns a reference to the internally stored symmetric dissimilarity matrix, wrapped in a struct that allows query operations to be performed.
    pub fn dissimilarity_matrix_ref(&self) -> & DissimilarityMatrix { & self.dissimilarity_matrix }

    /// Returns a reference to the internally stored symmetric dissimilarity matrix, wrapped in a struct that allows query operations to be performed.
    pub fn dissimilarity_matrix_size(&self) -> usize { self.dissimilarity_matrix_size }    

    /// Returns the dissimilarity value of two vertices, wrapped in a `NoneGreater` struct
    /// 
    /// The purpose of the `NoneGreater` struct is to "switch" the order of None values,
    /// so that they are regarded as strictly greater (rather than strictly less) than all `Some(x)` values.
    pub fn dissimilarity_as_nonegreater( &self, a: Vertex, b: Vertex ) -> NoneGreater<Filtration> 
        { NoneGreater { opt: self.dissimilarity_matrix.entry_major_at_minor(a as usize, b as usize) } }
    
    /// The maximum dissimilarity threshold
    pub fn dissimilarity_value_max( &self ) -> Filtration { self.dissimilarity_value_max.clone() }

    /// The ring operator for the coefficient ring.
    pub fn ring_operator( &self ) -> RingOperator where RingOperator: Clone { self.ring_operator.clone() }

    /// A vector of simplices sorted first by dimension (ascending `0..dimension_max`, including `dimension_max`) 
    /// then by diameter (descending) then by lexicographic order (descending)
    pub fn cliques_in_order( &self, dimension_max: isize ) -> Vec< SimplexFiltered< Filtration > > { 
        return filtered_cliques_in_order( dimension_max, self.dissimilarity_value_max(), self.dissimilarity_matrix_ref(), self.dissimilarity_matrix_size )
    }    

    /// An iterator that runs over simplices in ascending lexicographic order (not by filtration, etc.).
    pub fn cliques_in_lexicographic_order_fixed_dimension( &self, dimension: isize ) -> SimplexIter< Filtration, & DissimilarityMatrix > { 
        SimplexIter::new(                          
            self.dissimilarity_matrix_ref(), // a simple trick to make a copy of the dissimilarity matrix that implements the Copy trait
            self.dissimilarity_matrix_size,                
            self.dissimilarity_value_max(),
            dimension,                               
        )
    }        

    /// Returns a actorization of the boundary matrix restricted to simplices
    /// of dimension 0 .. dimension_max+1 (inclusive)
    /// 
    /// This factorization can be used to compute betti numbers, cycle representatives,
    /// and more.  See the documenation for `FactoredChainComplex`.
    pub fn factor_from_arc( self, max_homology_dim: isize ) 
            ->
            FactoredBoundaryMatrix< 
                    Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >, 
                    RingOperator, 
                    OrderOperatorSimplexFiltered< Filtration, Coefficient >, 
                    SimplexFiltered<Filtration>,
                    (SimplexFiltered<Filtration>, Coefficient),                    
                    Vec< SimplexFiltered<Filtration> >,                     
                >    
                // Umatch<  
                //         &Self,                        
                //         RingOperator, 
                //         OrderOperatorSimplexFiltered< Filtration, Coefficient >,
                //         OrderOperatorSimplexFiltered< Filtration, Coefficient >,
                //     >
        // where
        //     Filtration:         Hash,
        //     RingOperator:   DivisionRing< Coefficient >,
        where
            Filtration:             Clone + Debug + PartialOrd + Ord + Hash, 
            Coefficient:             Clone + Debug + PartialEq, 
            RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient > + DivisionRing< Coefficient >,            
    {
        let iter_keymaj     =   self.cliques_in_order( max_homology_dim );
        let arc = Arc::new(self);
        let ring_operator = arc.ring_operator();
        factor_boundary_matrix(
            arc,
            ring_operator,
            OrderOperatorAuto,
            iter_keymaj,
        )
    }

}


//  ===========================================================
//  COFACET ITERATOR
//  ===========================================================


/// An iterator that runs over all cofacets of a given simplex, in ascending lexicographic order, regardless of filtration value
pub struct IterCoboundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    where
        Filtration:             Clone + Debug,
        Coefficient:             Clone,
        DissimilarityMatrix:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                            ViewRowAscend + 
                            MatrixEntry,
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,     
{
	pub clique: Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
	pub next_cofacet_vertices: Vec<Vertex>,
    pub simplex_filtration: Filtration,
	pub insertion_location: usize,
    pub candidate_location: usize,
    pub first_vert: Vertex,
	pub coeff: Coefficient,
    pub ring_operator: RingOperator,
}

/// implement standard methods of Iterator for IterCoboundary struct
impl< 'a, DissimilarityMatrix, Filtration, Coefficient, RingOperator >

    Iterator for 
    
    IterCoboundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
        
    where
        Filtration:             Clone + Debug + Ord,
        Coefficient:             Clone,
        DissimilarityMatrix:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                            ViewRowAscend + 
                            MatrixEntry,
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,    
{
	type Item = (SimplexFiltered<Filtration>, Coefficient);

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
                let newdis_opt  =   self.clique.dissimilarity_matrix.entry_major_at_minor( 
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


//  ===========================================================
//  CLIQUE COBOUNDARY TRAIT
//  ===========================================================


pub trait GetBoundaryIters< Filtration >{
    fn coboundary_lexicographic_ascend_iter
            < DissimilarityMatrix, Coefficient, RingOperator >
            ( 
                self, 
                matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
            ) 
    -> 
    IterCoboundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    where
        Filtration:                         Clone + Debug + PartialOrd + Ord, 
        Coefficient:                         Clone, 
        RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                        ViewRowAscend + 
                                        MatrixEntry,        
        ;   

    fn coboundary_filtration_ascend_vec
            < DissimilarityMatrix, Coefficient, RingOperator >
            ( 
                self, 
                matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
            ) 
    -> Vec< ( SimplexFiltered<Filtration>, Coefficient ) >
    where
        Filtration:                         Clone + Debug + PartialOrd + Ord, 
        Coefficient:                         Clone, 
        RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                        ViewRowAscend + 
                                        MatrixEntry,        
        ; 


    fn boundary_lexicographic_descend_iter
            < DissimilarityMatrix, Coefficient, RingOperator >
            ( 
                self, 
                matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
            ) 
    -> 
    IterBoundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    where
        Filtration:                         Clone + Debug + PartialEq + Eq + Ord + PartialOrd,
        Coefficient:                    Clone,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                        ViewRowAscend + 
                                        MatrixEntry,
        DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
        RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >;   


    fn boundary_filtration_descend_vec
                < DissimilarityMatrix, Coefficient, RingOperator >
                ( 
                    self, 
                    matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
                ) 
        -> 
        Vec< ( SimplexFiltered<Filtration>, Coefficient ) >
        where
            Filtration:                         Clone + Debug + PartialEq + Eq + Ord + PartialOrd,
            Coefficient:                    Clone,
            DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                            ViewRowAscend + 
                                            MatrixEntry,
            DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
            RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >;          
}

impl < Filtration >

    GetBoundaryIters< Filtration > for 

    SimplexFiltered< Filtration > 
    
    where 
        Filtration:     Clone + Debug,

{
    /// Returns entries of the boundary matrix column for this simplex, in descending lexicographic order (filtration values are not sorted)
    fn boundary_lexicographic_descend_iter
                < DissimilarityMatrix, Coefficient, RingOperator >
                ( 
                    self, 
                    matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
                ) 
        -> 
        IterBoundary< DissimilarityMatrix, Filtration, Coefficient, RingOperator >
        where
            Filtration:                         Clone + Debug + PartialEq + Eq + Ord + PartialOrd,
            Coefficient:                    Clone,
            DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                            ViewRowAscend + 
                                            MatrixEntry,
            DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
            RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient > 
    {
        IterBoundary {
            clique: matrix.clone(),
            simp: self,
            removal_location: 0,
            ring_operator: matrix.ring_operator(),
        }
    }

    /// Returns entries of the boundary matrix column for this simplex, in descending lexicographic order, first by filtration, and then by vertices
    fn boundary_filtration_descend_vec
                < DissimilarityMatrix, Coefficient, RingOperator >
                ( 
                    self, 
                    matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
                ) 
        -> 
        Vec< ( SimplexFiltered<Filtration>, Coefficient ) >
        where
            Filtration:                         Clone + Debug + PartialEq + Eq + Ord + PartialOrd,
            Coefficient:                    Clone,
            DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                            ViewRowAscend + 
                                            MatrixEntry,
            DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
            RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient > 
    {
        let mut boundary = self.boundary_lexicographic_descend_iter(matrix).collect_vec();
        boundary.shrink_to_fit();
        boundary.sort_by(|x,y| x.0.cmp(&y.0) );
        boundary
    }    

    /// Returns the coboundary of `self` as a vector **with entries sorted in ascending lexicograph order (filtration values are not sorted)**
    fn coboundary_lexicographic_ascend_iter
                < DissimilarityMatrix, Coefficient, RingOperator >
                ( 
                    self, 
                    matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
                ) 
        -> 
        IterCoboundary< DissimilarityMatrix, Filtration, Coefficient, RingOperator >
        where
            Filtration:                         Clone + Debug + PartialOrd + Ord, 
            Coefficient:                         Clone, 
            RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >,
            DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                            ViewRowAscend + 
                                            MatrixEntry, 
    {
        let ( filtration, mut vertices) = self.dissolve();
        let first_vert = vertices[0];
        vertices.insert(0,0);
        vertices.shrink_to_fit();
        let ring_operator = matrix.ring_operator.clone();
        IterCoboundary{
            clique: matrix,
            next_cofacet_vertices: vertices,
            simplex_filtration: filtration,
            insertion_location: 0,
            candidate_location: 0,
            coeff: RingOperator::one(),
            first_vert,
            ring_operator,
        }        
    }

    /// Returns the coboundary of `self` as a vector **with entries sorted in ascending lexicograph order, first by filtration, then by vertices.**
    fn coboundary_filtration_ascend_vec
                < DissimilarityMatrix, Coefficient, RingOperator >
                ( 
                    self, 
                    matrix:  Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
                ) 
        -> Vec< ( SimplexFiltered<Filtration>, Coefficient ) >
        where
            Filtration:                         Clone + Debug + PartialOrd + Ord, 
            Coefficient:                         Clone, 
            RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >,
            DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                            ViewRowAscend + 
                                            MatrixEntry, 
    {
        let mut coboundary = self.coboundary_lexicographic_ascend_iter(matrix).collect_vec();
        coboundary.shrink_to_fit();
        coboundary.sort_by(|x,y| x.0.cmp(&y.0) );
        coboundary
    }
}    




// /// Returns the coboundary of a simplex `keymaj` as a vector **with entries sorted in ascending order**
// /// 
// /// The order is the usual lexicographic order: first by birth time, second by lexicographic order.
// fn get_coboundary_as_vec< DissimilarityMatrix, Filtration, Coefficient, RingOperator >( 
//             matrix: Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >, 
//             keymaj: SimplexFiltered<Filtration> 
//         ) 
//         -> Vec< ( SimplexFiltered<Filtration>, Coefficient ) >
        
//     where
//         Filtration:                         Clone + Debug + PartialOrd + Ord, 
//         Coefficient:                         Clone, 
//         RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >,
//         DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
//                                         ViewRowAscend + 
//                                         MatrixEntry,
// {

//     let mut vertices = keymaj.vertices.clone();
//     let first_vert = vertices[0];
//     vertices.insert(0,0);
//     vertices.shrink_to_fit();
//     let ring_operator = matrix.ring_operator.clone();
//     let iter = 
//     IterCoboundary{
//         clique: matrix,
//         next_cofacet_vertices: vertices,
//         simplex_filtration: keymaj.filtration,
//         insertion_location: 0,
//         candidate_location: 0,
//         coeff: RingOperator::one(),
//         first_vert,
//         ring_operator,
//     };   
//     let mut vec = iter.collect_vec();
//     vec.shrink_to_fit();
//     vec.sort_by(|x,y| x.0.cmp(&y.0) );
//     vec
// }



//  ===========================================================
//  FACET ITERATORS
//  ===========================================================


/// A iterator that runs over all entries in the boundary of a simplex, in **descending lexicographic order**
/// 
/// If the original simplex is `[0,1,2]`, it return entries in the following order: `([1,2],1), ([0,2],-1), ([0,1],1)`.
pub struct IterBoundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    where
        Filtration:                         Clone + Debug,
        Coefficient:                         Clone,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                        ViewRowAscend + 
                                        MatrixEntry,
        DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
        RingOperator:                   Semiring< Coefficient > + Ring< Coefficient >, 
{
    pub clique: Arc< ChainComplexVrFiltered< DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
	pub simp: SimplexFiltered<Filtration>,
	pub removal_location: usize,
    pub ring_operator: RingOperator,
}

/// implement standard methods of Iterator for IterBoundary struct
impl< DissimilarityMatrix, Filtration, Coefficient, RingOperator >

    Iterator for 
    
    IterBoundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    where
        Filtration:                         Clone + Debug + Ord + PartialOrd,
        Coefficient:                         Clone,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                        ViewRowAscend + 
                                        MatrixEntry,
        DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
        RingOperator:                   Semiring< Coefficient > + Ring< Coefficient >, 
{
	type Item = (SimplexFiltered<Filtration>, Coefficient);

	fn next( &mut self ) -> Option<Self::Item> {
        if self.simp.vertices.len() == 1 { return None } // the boundary of a 0-simplex is empty
        if self.removal_location == self.simp.vertices.len() { return None; }
        
        let mut simplex = self.simp.clone();
        // println!("simplex we are bounding: {:?}", &simplex.vertices);
        simplex.vertices.remove(self.removal_location);
        simplex.vertices.shrink_to_fit();
        // println!("print out of the vertices: {:?}", &simplex.vertices);
        simplex.filtration = self.clique.diameter(&simplex.vertices).unwrap();
        let coeff = self.ring_operator.minus_one_to_power( self.removal_location );
        self.removal_location += 1;
        Some((simplex, coeff))

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


impl < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    
    IndicesAndCoefficients for    

    // ChainComplexVrFilteredArc
    //     < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    Arc< 
            ChainComplexVrFiltered
                < DissimilarityMatrix, Filtration, Coefficient, RingOperator > 
        >    

    where
        Filtration:             Clone + Debug + PartialOrd + Ord, 
        Coefficient:             Clone, 
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
        DissimilarityMatrix:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                            ViewRowAscend + 
                            MatrixEntry,         
{
    type EntryMajor =   ( Self::RowIndex, Self::Coefficient );
    type EntryMinor       =   ( Self::RowIndex, Self::Coefficient );    
    type RowIndex = SimplexFiltered<Filtration>; 
    type ColIndex = SimplexFiltered<Filtration>; 
    type Coefficient = Coefficient;
}


//  ORACLE MAJOR ASCEND
//  ------------------------------------------


impl < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    
    ViewRowAscend for 

    // ChainComplexVrFilteredArc
    //             < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    Arc< 
            ChainComplexVrFiltered
                < DissimilarityMatrix, Filtration, Coefficient, RingOperator > 
        >        

    where
        Filtration:                         Clone + Debug + Ord + PartialOrd,
        Coefficient:                         Clone,
        DissimilarityMatrix:                         IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                        ViewRowAscend + 
                                        MatrixEntry,
        DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,
        RingOperator:                   Clone + Semiring< Coefficient > + Ring< Coefficient >,       
{
    type ViewMajorAscend            =   LazyOrderedCoboundary< DissimilarityMatrix, Filtration, Coefficient, RingOperator >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymaj: Self::RowIndex ) -> Self::ViewMajorAscend {

        let first_v = keymaj.vertices[0];
        // let mut shortcircuit_vertex = None;

        'search_pareto: for term in keymaj.clone().coboundary_lexicographic_ascend_iter(self.clone()) {
            if term.0.filtration == keymaj.filtration { // if the cofacet has the same filtration as the facet
                if term.0.vertices[0] < first_v { // and it's obtained by adding a vertex with number lower than all vertices in the facet
                    return LazyOrderedCoboundary{ // then it's on the pareto frontier
                        facet:              keymaj,
                        boundary_matrix:    self.clone(),
                        coboundary:         vec![ term ], 
                        next_index:         0,
                        built_all:          false,
                    }                    
                } else {
                    let bottom_of_column 
                        =   term.0
                                .clone()
                                .boundary_lexicographic_descend_iter(self.clone())
                                .map(|x| x.0)
                                .max()
                                .unwrap();
                    if bottom_of_column == keymaj { // if the facet lives at the bottom of the column indexed by `term`
                        return  LazyOrderedCoboundary{ // then it's on the pareto frontier
                                    facet:              keymaj,
                                    boundary_matrix:    self.clone(),
                                    coboundary:         vec![ term ], 
                                    next_index:         0,
                                    built_all:          false,
                                }                           
                    } else {
                        break 'search_pareto // we have to break the search after finding the first `term` with equal filtration, if the term fails the test
                    }
                }
            }
        }
        LazyOrderedCoboundary{
            facet:              keymaj.clone(),
            boundary_matrix:    self.clone(),
            coboundary:         keymaj.coboundary_filtration_ascend_vec(self.clone()),
            next_index:         0,
            built_all:          true,
        }  


        // 'outer_vertex_loop: for outer_v in self.arc.cutoff_matrix[first_v as usize ].iter().cloned() {
        //     // search for a vertex, `outer_v` whose id number is strictly below that of all 
        //     // vertices inside the simplex, such that `outer_v` could be added to the simplex 
        //     // without increasing its diameter
        //     // NB: we specifically need the vertex with the *smallest possible id number*
        //     //     that satisfies this property

        //     // -- record the first vertex
        //     // -- loop through all cofacets
        //     // -- first with same diameter: 
        //     //    -- check if has same first vertex; if so we're good 
        //     //    -- otherwise check to see if is the last facet

        //     if outer_v >= first_v { break }
        //     let diam = NoneGreater::from_val( keymaj.filtration() );
        //     // exclude outer_v if simplex [outer_v] is born too late
        //     if self.arc.dissimilarity_as_nonegreater( outer_v, outer_v ) > diam 
        //         { continue }
        //     // exclude outer_v if it is too far from any vertex in the simplex
        //     for inner_v in keymaj.vertices.iter().cloned() {
        //         if self.arc.dissimilarity_as_nonegreater( outer_v, inner_v ) > diam 
        //             { continue 'outer_vertex_loop }
        //     }
        //     shortcircuit_vertex=Some( outer_v ); 
        //     break
        // }

        // if let Some( outer_v ) = shortcircuit_vertex {
        //     // if the preceding for-loop returns a valid value for outer_v, then construct
        //     // only the first entry

        //     let mut first_cofacet = Vec::with_capacity( keymaj.vertices.len()+1 );
        //     first_cofacet.push( outer_v );
        //     for inner_v in keymaj.vertices.iter().cloned() { first_cofacet.push(inner_v) }
        //     let cofacet = SimplexFiltered{ vertices: first_cofacet, filtration: keymaj.filtration() };
        //     let first_entry = ( cofacet, RingOperator::one() );
        //     LazyOrderedCoboundary{
        //         facet:              keymaj,
        //         boundary_matrix:    self.arc.clone(),
        //         coboundary:         vec![ first_entry ],
        //         next_index:         0,
        //         built_all:          false,
        //     }
        // } else {
        //     // otherwise construct precompute *all* entries, and store them in a Vec in sorted order

        //     let coboundary = keymaj.clone().coboundary_filtration_ascend_vec(self.matrix_arc().clone());//get_coboundary_as_vec(self.arc.clone(), keymaj.clone() );
        //     LazyOrderedCoboundary{
        //         facet:              keymaj,
        //         boundary_matrix:    self.arc.clone(),
        //         coboundary,
        //         next_index:         0,
        //         built_all:          true,
        //     }            
        // }

    }
}


/// Iterator returning entries of a coboundary in order
/// 
/// This iterator is "lazy" because it will try to stop building the coboundary after generating,
/// the first element, using the Morse optimization.
/// 
/// After building the first element it will generate all other elements all at once, place them in a sorted vector,
/// and return elements of the vector one at a time.
#[derive(Clone,Debug,new)]
pub struct LazyOrderedCoboundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    where 
        Filtration:             Clone + Debug,
        Coefficient:             Clone,
        RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
        DissimilarityMatrix:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                            ViewRowAscend + 
                            MatrixEntry,        
{
    facet:              SimplexFiltered<Filtration>,
    boundary_matrix:    Arc< ChainComplexVrFiltered < DissimilarityMatrix, Filtration, Coefficient, RingOperator > >,
    coboundary:         Vec< (SimplexFiltered<Filtration>, Coefficient) >,
    next_index:         usize,
    built_all:          bool
}

impl < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    
    Iterator for
    
    LazyOrderedCoboundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >

    where 
        Filtration:             Clone + Debug + Ord,
        Coefficient:             Clone,
        RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient >,     
        DissimilarityMatrix:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                            ViewRowAscend + 
                            MatrixEntry,         
{
    type Item = (SimplexFiltered<Filtration>, Coefficient);

    fn next(&mut self) -> Option<Self::Item> {

        while self.next_index >= self.coboundary.len() {
            if self.built_all { return None }
            else{
                let mut v = self.facet.clone().coboundary_filtration_ascend_vec(self.boundary_matrix.clone());//get_coboundary_as_vec( self.boundary_matrix.clone(), self.facet.clone() );
                v.sort_by(|x,y| x.0.cmp(&y.0)); 
                self.coboundary = v;
                self.built_all = true;                
            }
        }

        let return_value = self.coboundary[ self.next_index ].clone();
        self.next_index += 1;

        Some( return_value )
    }
}   


impl < 'a, DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    
    ParetoShortCircuit
        < (SimplexFiltered<Filtration>, Coefficient) > for
    
    LazyOrderedCoboundary
        < DissimilarityMatrix, Filtration, Coefficient, RingOperator >

    where 
        Filtration:             Clone + Debug + Ord,
        Coefficient:             Clone,
        RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient >,   
        DissimilarityMatrix:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                            ViewRowAscend + 
                            MatrixEntry,                    
{
    fn pareto_short_circuit(& self) -> Option< (SimplexFiltered<Filtration>, Coefficient) > {
        if ! self.built_all { return Some( self.coboundary[0].clone() ) }
        None
    }
} 


//  ORACLE MINOR DESCEND
//  ------------------------------------------


impl < 'a, DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    
ViewColDescend for 

    // ChainComplexVrFilteredArc
    //     < DissimilarityMatrix, Filtration, Coefficient, RingOperator >
    Arc< 
            ChainComplexVrFiltered
                < DissimilarityMatrix, Filtration, Coefficient, RingOperator > 
        >        

    where
        Filtration:             Clone + Debug + PartialOrd + Ord, 
        Coefficient:             Clone, 
        RingOperator:       Clone + Semiring< Coefficient > + Ring< Coefficient >,
        DissimilarityMatrix:             IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                            ViewRowAscend + 
                            MatrixEntry,
        DissimilarityMatrix::EntryMajor:   KeyValGet< usize, Filtration >,                                                                 
{
    type ViewMinorDescend            =   Vec< ( SimplexFiltered<Filtration>, Coefficient ) >;
    type ViewMinorDescendIntoIter    =   std::vec::IntoIter<(SimplexFiltered<Filtration>, Coefficient)>;

    fn view_minor_descend( &self, keymin: Self::ColIndex ) -> Self::ViewMinorDescend {

        let iter = IterBoundary {
            clique: self.clone(),
            simp: keymin,
            removal_location: 0,
            ring_operator: self.ring_operator.clone(),
        };
        let mut vec = iter.collect_vec();     
        vec.shrink_to_fit();
        vec.sort_by(|x,y| y.0.cmp(&x.0));  // NOTE WE GIVE REVERSE ORDER

        vec
    }
}


// /// Based, filtered chain complex implementing the [`ChainComplex`](ChainComplex) trait
// pub struct CliqueComplex<Coefficient: Clone, Filtration> {
//     pub ring_operator: RingMetadata<Coefficient>,
//     pub dissimilarity_matrix: Vec<Vec<Filtration>>,
//     pub dissimilarity_value_max: Filtration,
//     pub safe_homology_degrees_to_build_boundaries: Vec<usize>,
//     pub major_dimension: MajorDimension,
//     pub simplex_count: Vec<(usize,usize)>
// }


//  SIMPLEX ITERATOR
//  ------------------------------------------

/// An iterator that runs over all simplices of given dimension in a clique complex, in **ascending lexicographic order**.
/// 
/// # Examples
/// 
/// ```
/// use sprs::CsMatBase;
/// use ordered_float::OrderedFloat;
/// use oat_rust::topology::simplicial::from::graph_weighted::SimplexIter;
/// use oat_rust::topology::simplicial::simplices::filtered::SimplexFiltered;
/// 
/// 
/// let dissimilarity_matrix = CsMatBase::new( (2,2), vec![0,1,2], vec![0,1], vec![OrderedFloat(0.0),OrderedFloat(0.0)] );
///     
/// let simplices: Vec<_>    =   SimplexIter::new( 
///                     & dissimilarity_matrix, 
///                     2,          //  matrix size
///                     OrderedFloat(1.0),        //  dissimilarity value max
///                     1,          //  dimension
///                 ).collect();
/// let empty =  Vec::< SimplexFiltered<OrderedFloat<f64>> >::new();
///             
/// assert!( empty.eq( &simplices ) ) 
/// ```
#[derive(Clone,Debug)]
pub struct SimplexIter< Filtration, DissimilarityMatrix > {
    pub dissimilarity_matrix: DissimilarityMatrix,
    pub dissimilarity_matrix_size: usize,    
    pub dissimilarity_value_max: NoneGreater< Filtration >,
    pub filvec: Vec< NoneGreater< Filtration > >,
    pub vec: Vec<Vertex>,
    pub val: Vertex,
    pub loc: usize,
}

impl < Filtration, DissimilarityMatrix >

    SimplexIter
        < Filtration, DissimilarityMatrix > 

    where
        Filtration:                 Clone + Debug + Ord,
    {

    pub fn new( 
            dissimilarity_matrix:           DissimilarityMatrix,
            dissimilarity_matrix_size:      usize,
            dissimilarity_value_max:        Filtration,
            dim:                            isize,             
        ) 
        -> 
        Self 
        where
            DissimilarityMatrix:                            IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                                            ViewRowAscend,
            DissimilarityMatrix::EntryMajor:      KeyValGet< usize, Filtration >,
    {

        let dim = dim as usize;
        let minval 
            = (0..dissimilarity_matrix_size)
                .flat_map(|x| dissimilarity_matrix.view_major_ascend(x)).map(|x| x.val() ).min().unwrap();
        SimplexIter {
            dissimilarity_matrix,
            dissimilarity_matrix_size,
            dissimilarity_value_max: NoneGreater::from_val( dissimilarity_value_max ),
            filvec: vec![NoneGreater::from_val(minval); dim + 2 ], // we will maintain an internal state of the iterator such that filvec[p] = diameter of the simplex spanned by the first (p-1) simplices; this is useful as it allows us to store a value for the -1 simplex
            vec: vec![0; dim +1 ], // the vertices of the simplex we're building
            val: 0, // the vertex we insert / have inserted
            loc: 0, // where in the simplex we insert / have inserted
        }
    }
}


/// implement standard methods of Iterator for SimplexIter struct
impl< Filtration, DissimilarityMatrix > 

    Iterator for 
    
    SimplexIter
        < Filtration, DissimilarityMatrix > 
        
    where
        Filtration:                                 Clone + PartialOrd + Ord + Debug,
        DissimilarityMatrix:                            IndicesAndCoefficients< ColIndex=usize, RowIndex=usize, Coefficient=Filtration > +
                                                        ViewRowAscend +
                                                        MatrixEntry,
        DissimilarityMatrix::EntryMajor:        KeyValGet< usize, Filtration >,
{
    type Item = SimplexFiltered<Filtration>;

    /// Get the next simplex
    /// 
    /// This iterator works on the following principles:
    /// - we iterate over simplices in ascending lexicographic order
    /// - we add vertices in a loop
    ///   - suppose we are at a point in the loop where we have filled in the first 3 slots: [1 2 3 * * *]
    ///   - the next simplex we add will have vertices inserted in strictly ascending order, so we can essentially
    ///     do a for-loop over vertices 4... inserting the first sequence of vertices possible
    ///   - after we do that we can try to increment the last vertex
    ///   - then if the last vertex is at ceiling we can try to incremeent the second-to-last vertex
    ///   - etc.
    fn next(&mut self) -> Option<Self::Item> {
        let size = self.dissimilarity_matrix_size;
        let mut filval;
        if self.vec.len() > size { return None; }
        loop {
            'insert_larger_vertex_number: while usize::from(self.val) <= size - (self.vec.len() - self.loc) {
                
                // write in a new vertex value
                self.vec[self.loc] = self.val; 

                // compute the diameter of the simplex created by adding a vertex just now
                self.filvec[self.loc+1] = self.filvec[self.loc].clone(); // initialize the diameter as the diameter of the previous simplex; our constructor function inializes the diameter of the emepty simplex as the minimum entry of the dissimilarity_matrix
                
                
                for ii in 0..self.loc+1 {                    
                    filval = NoneGreater::from_opt( // wrap the dissimilarity score in a NoneGreater
                        self.dissimilarity_matrix.entry_major_at_minor( self.val as usize, self.vec[ii] as usize ) 
                    );
                    // println!("II {:?}", &ii );
                    // println!("self.loc {:?}", &self.loc);
                    // println!("self.val {:?}", &self.val);
                    // println!("self.vec[ii] {:?}", &self.vec[ii] );
                    // println!("filval {:?}", &filval);
                    // println!("self.dissimilarity_value_max {:?}", &self.dissimilarity_value_max);                    
                    if filval > self.dissimilarity_value_max {  // if the dissimilarity is too great
                        self.val +=1;   // increment the id number of the vertex we try to insert
                        // println!("restarting while loop");                   
                        continue 'insert_larger_vertex_number   // and start over
                    } else if filval > self.filvec[self.loc+1] {
                        self.filvec[self.loc+1] = filval; // increase the recorded diameter if necessary
                    }
                    // match self.dissimilarity_matrix.entry_major_at_minor( self.val as usize, self.vec[ii] as usize ) {
                    //     Some(filval) => {
                    //         if self.filvec[self.loc+1]  <  filval {
                    //             self.filvec[self.loc+1] = filval; // increase the recorded diameter if necessary
                    //         }
                    //     }
                    //     None => { continue } // if no dissimilarity is recorded, then we assume this pair of vertices to be NON-adjacent
                    // }
                }

                // if the simplex has the correct dimension and a sufficiently small diameter, then return it
                if self.filvec[self.loc+1] <= self.dissimilarity_value_max {
                    if self.loc == self.vec.len()-1 {
                        self.val += 1;
                        return Some(SimplexFiltered{
                                vertices: self.vec.clone(),
                                filtration: self.filvec[self.loc+1].val()
                            });
                    } else {
                        self.loc += 1;
                    }
                }
                self.val += 1;
            }
            if self.loc == 0 { return None; }
            self.loc -= 1;
            self.val =	self.vec[self.loc] + 1;

        }
    }
}






//  ================================================================================
//  COBOUNDARY / COFACET ITERATORS
//  ================================================================================




// /// I think this is supposed to run over all vertices that can be added to the facet without changing its diameter, in lexicographic order.
// pub struct ShortcircuitVertices
//         < Filtration, I >
//     where
//         Filtration:             Clone + Debug + Ord,
//         I:                  Iterator,
//         I::Item:            KeyValGet< Vertex, Filtration >,
// {
//     iter:                   IntersectOrderedIterators< I, OrderOperatorByKey< Vertex, Filtration, I::Item > >,
//     filtration_max:         Filtration,
//     facet_vertices:         Vec< Vertex >,
//     facet_cardinality:      usize,
//     duplication_pointer:    usize,
// }        


// impl < Filtration, I >

//     Iterator for

//     ShortcircuitVertices
//             < Filtration, I >
//         where
//             Filtration:         Clone + Debug + Ord,
//             I:              Iterator,
//             I::Item:        KeyValGet< Vertex, Filtration >,
//     {
//         type Item = Vertex;

//         fn next(&mut self) -> Option<Self::Item> {

//             println!("modify the HitMerge operator to ");
            
//             let mut multiplicity = 0;
//             let mut vertex_new;
//             let mut vertex_new_opt = None;
//             let mut vertex_old_opt = None;            

//             loop {
//                 match self.iter.next() { // self.iter is obtained by merging I_1,..,I_m, where m = number of vertices in the facet and I_i is the set of neighbors of facet[i], ordered lexicographically
//                     None => { 
//                         return None; 
//                     } Some( edge ) => {
//                         match edge.val() < self.filtration_max { // if the vertex is too far, then exclude it
//                             false => { 
//                                 continue;
//                             } true => {
//                                 vertex_new      = edge.key();
//                                 vertex_new_opt  = Some( vertex_new ); // otherwise we'll check to see if it's close enought to every other vertex in the facet, by counting multiplicity
                                
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


// /// A wrapper struct; used instead of a tuple of form (vertex_outer, filtration), because the fields are more readable
// /// 
// /// In practice this struct will represent an edge of form (vertex_inner, vertex_outer) with length = filtration. Thanks
// /// to the context in which this struct is used, we never have to store vertex_inner explicitly
// struct EdgeToAdd< Filtration, >{
//     outer_vertex:   Vertex,    
//     filtration:     Filtration,
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
//                 < Filtration, I, J >
//     where
//         I:      Iterator< Item=(Vertex,Filtration) > + PeekUnqualified,
//         J:      Iterator< Item= Vertex >,
// {
//     facet_fil:                              Filtration,
//     neighbor_iters_filtration_order:        Vec< I >,
//     neighbors_within_diameter_lex_order:    J,
// }


// impl < Filtration, I, J >

//     Iterator for 

//     CofacetEdgeIter
//         < Filtration, I, J >
    
//     where
//         I:                  Iterator< Item=EdgeToAdd< Filtration > + PeekUnqualified,
//         J:                  Iterator< Item=Vertex >,
//         Filtration:             Clone + Ord,
// {    
//     type Item = ( Vertex, Filtration );

//     fn next(&mut self) -> Option<Self::Item> {

//         println!("J is a merge of several iterators; design J so that it stops when one of its internal iterators is empty, and returns only vertices which meet diameter bounds for all vertices of the facet, and excludes vertices in the facet");

//         // first check if we can generate a cofacet of equal diameter
//         // ----------------------------------------------------------
//         let a = self.neighbors_within_diameter_lex_order.next();
//         if a.is_some() { return Some( ( a.unwrap(), self.facet_fil.clone() ) ) }

//         // if not, then start iterating over edges incident to the vertices in the facet
//         // -----------------------------------------------------------------------------
        
//         // initialize some variables
//         let facet_len = self.neighbor_iters_filtration_order.len(); // number of vertices in the facet
//         let mut iter;        // a dummy variable
//         let mut ii_old = 0;  // an integer in 0..len(facet).  it marks a vertex within the facet who has a neighber-that-we-want-to-add-to-form-a-cofacet
//         let mut ii_new = 0;  // an integer in 0..len(facet).  it marks a vertex within the facet that we want to compare to ii_old
//         let e_old_opt = self.neighbor_iters_filtration_order[ii_old].next(); // see explanation of e_old, below
        
//         match e_old_opt {
//             None => { 
//                 // in this case we've exhausted all neighbors the ii_old'th vertex in our facet, so no more cofacets exist
//                 return None 
//             } Some( mut e_old ) => { 
//                 // e_old is a pair of form 
//                 //     (vertex_outer_old, `filtration_old`); 
//                 // it represents a edge of form 
//                 //     (facet[ii_old], vertex_outer_old) with filtration value `filtration_old`
//                 // this edge has length greater than the diameter of the facet, because (assuming that we've
//                 // done our homework) we will have already pulled vertices closer than diameter-of-the-facet
//                 // out of their respective iterators
                
//                 'loop_label_outer: loop {
//                     // in this loop we check to see if e_old.outer_vertex is (1) incident to every other vertex in 
//                     // the facet, and (2) no farther from these vertices than e_old.filtration; 
//                     // if not, then pick a new edge to try to add as the longest edge of a new cofacet, and start again

//                     // increment to the next vertex in the facet
//                     ii_new = ( ii_new + 1 ) % facet_len;

//                     // if we've successfully checked every vertex of the facet and circled back to the first one 
//                     // that we originally checked, then return the new edge!
//                     if ii_new == ii_old { return Some(e_old) }

//                     // otherwise we have to do some work
//                     iter = &mut self.neighbor_iters_filtration_order[ ii_new ];                    
//                     'loop_label_inner: loop {
//                         // check that e_old.vertex_outer lies within distance e_old.filtration of the ii_new'th 
//                         // vertex in the facet, that is facet[ii_new]; if not, then either 
//                         // (i) identify a different edge to add, or 
//                         // (ii) delare that we have finished generating all cofacets
//                         match iter.next() {
//                             None => { 
//                                 // in this case facet[ii_new] has no more neighbors; therefore there exist no more cofacets!  we're done!    
//                                 return None 
//                             } Some( e_new ) => { 
//                                 // e_new is a pair of form 
//                                 //    (vertex_outer_new, `filtration_new`).  
//                                 // it represents an edge of form 
//                                 //     (facet[ii_new], e_new.vertex_outer) whose filtration value is `filtration_new`
//                                 if e_new.filtration > e_old.filtration {  
//                                     // in this case e_old would not be the longest edge in a new facet, because we've just found a longer one
//                                     if e_new.outer_vertex != e_old.outer_vertex { 
//                                         // if e_new.outer_vertex != e_old.outer_vertex, then we have to check that
//                                         // every vertex in facet lies within distance e_new.filtration of e_new.outer_vertex;
//                                         // that is, we need to check vertices 
//                                         // facet[ii_new+1 .. ii_new]  <-- cyclice indexing
//                                         // we can ensure that all these vertices will be checked, by resetting
//                                         ii_old = ii_new;

//                                         // however, if e_new.outer_vertex == e_old.outer_vertex, then we only have to check
//                                         // facet[ii_new+1 .. ii_old], so we don't have to re-assign the value of ii_old

//                                         // NOTA BENE: I don't think we have to place this re-assignment inside an if-clause, but it
//                                         // benefits us to do so, because we will do fewer operations / less checking if we
//                                         // don't reassign ii_old to ii_new.  The point is that if e_new.outer_vertex = e_old.outer_vertex,
//                                         // then we already have a certificate that this vertex lies within distance e_new.filtration
//                                         // of the vertices in facet[ii_old .. ii_new]

//                                     }
//                                     // record that we are now trying to add e_new as the longest edge
//                                     e_old = e_new;     
//                                     // now move on to the next neighbor iterator                               
//                                     break;
//                                 } else if e_new.vertex_outer == e_old.vertex_outer {
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
//                 < Filtration, I, J >
//     where
//         I:      Iterator< Item=(Vertex,Filtration) > + PeekUnqualified,
//         J:      Iterator< Item= Vertex >,
// {
//     cofacet_edge_iter:      CofacetEdgeIter< Filtration, I, J >,
//     facet_vertices:         Vec< Vertex >,
// }

// impl < Filtration, I, J >

//     Iterator for 

//     CofacetIter
//         < Filtration, I, J >
    
//     where
//         I:                  Iterator< Item=(Vertex,Filtration) > + PeekUnqualified,
//         J:                  Iterator< Item=Vertex >,
//         Filtration:             Clone + Debug + Ord,
// {    
//     type Item = SimplexFiltered< Filtration >;

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
//                 < RingOperator, Coefficient, Filtration, I, J, >
//     where
//         I:                  Iterator< Item=(Vertex,Filtration) > + PeekUnqualified,
//         J:                  Iterator< Item= Vertex >,
//         RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
// {
//     cofacet_edge_iter:      CofacetEdgeIter< Filtration, I, J >,
//     facet_vertices:         Vec< Vertex >,
//     ring_operator:          RingOperator,
//     phantom_snzval:         PhantomData<Coefficient>,
// }

// impl < RingOperator, Coefficient, Filtration, I, J, >

//     Iterator for 

//     CoboundaryIter
//         < RingOperator, Coefficient, Filtration, I, J, >
    
//     where
//         I:                  Iterator< Item=(Vertex,Filtration) > + PeekUnqualified,
//         J:                  Iterator< Item=Vertex >,
//         Filtration:             Clone + Debug + Ord,
//         RingOperator:       Semiring< Coefficient > + Ring< Coefficient >,
// {    
//     type Item = ( SimplexFiltered< Filtration >, Coefficient );

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
    use crate::{topology::simplicial::simplices::filtered::SimplexFiltered, algebra::matrices::query::MatrixEntry, utilities::order::OrderOperatorAuto};

    


  
    
    #[test]
    fn check_that_some_basic_functions_run_without_error() {


        use crate::topology::simplicial::{simplices::filtered::SimplexFiltered, from::graph_weighted::{ChainComplexVrFiltered}};        
        use crate::topology::point_cloud::unit_circle;    

        use crate::algebra::vectors::entries::KeyValGet;
        use crate::algebra::rings::operator_structs::ring_native::{FieldRational64};    
        use crate::algebra::matrices::debug::verify_viewmajorascend_compatible_with_viewminordescend;
        use crate::algebra::matrices::query::{ViewColDescend, ViewRowAscend};
        use crate::algebra::matrices::operations::umatch::row_major::{Umatch};    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        use crate::utilities::order::{ is_sorted_strictly, OrderOperatorByKeyCutsom, OrderOperatorAutoReverse};
        use crate::utilities::iterators::general::minmax;    
        use crate::utilities::distances::{rowwise_distances};

        use std::sync::Arc;
        use ordered_float::OrderedFloat;  
        use itertools::Itertools;      

        let npoints = 20;
        let dimension_max = 2;


        let pcloud = unit_circle( npoints, Some(-1.0 .. 1.0));

        let dissimilarity_matrix_data
            = rowwise_distances(pcloud)
                .into_iter()
                .map(|x| x.into_iter().enumerate().collect_vec() )
                .collect_vec()
                .into_csr( npoints, npoints );
        let dissimilarity_matrix = & dissimilarity_matrix_data;

        let dissimilarity_value_min = OrderedFloat(0.0);        
        let dissimilarity_value_max = 
        minmax( 
                (0..npoints).map(
                        |x| 
                        dissimilarity_matrix.view_major_ascend(x).into_iter().map(
                                |x| 
                                x.val()
                            ) 
                    ) 
            ).unwrap_or( dissimilarity_value_min.clone() );         

    
        let ring_operator = FieldRational64::new();
        let chain_complex_data = ChainComplexVrFiltered::new( & dissimilarity_matrix, npoints, dissimilarity_value_max, dissimilarity_value_min, ring_operator );
        // let chain_complex_ref = & chain_complex; 
        // let chain_complex = ChainComplexVrFilteredArc{ arc: Arc::new(chain_complex_data) };   
        let chain_complex = Arc::new(chain_complex_data);           
        let keymaj_vec = chain_complex.cliques_in_order(dimension_max);
        let keymin_vec = chain_complex.cliques_in_order(dimension_max+1);
    
    
        verify_viewmajorascend_compatible_with_viewminordescend(
                chain_complex.clone(),
                keymin_vec.iter().cloned(),
                keymaj_vec.iter().cloned(),
            );        
    
        let iter_keymaj = keymaj_vec.iter().cloned();    
    
        println!("check that oracle has strictly sorted rows");
        // print_indexed_major_views( & chain_complex_ref, iter_keymaj.clone() );  // print the major views       
        for keymaj in iter_keymaj.clone() {        
            assert!( is_sorted_strictly( 
                                            & chain_complex.view_major_ascend(keymaj.clone()).collect_vec() , 
                                            & OrderOperatorByKeyCutsom::new( OrderOperatorAuto ) 
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
            assert!( is_sorted_strictly(    & chain_complex.view_minor_descend(keymaj).iter().cloned().collect_vec() , 
                                            & OrderOperatorByKeyCutsom::new( OrderOperatorAutoReverse::new() )  // NOTE THAT HERE WE USE GT 
                                        ) );
        }    
    
        // println!("press enter to continue");
        // let mut guess = String::new();
        // io::stdin()
        //     .read_line(&mut guess)
        //     .expect("Failed to read line");       
    
        println!("starting umatch");
        let umatch = Umatch::factor_with_clearing(
                chain_complex.clone(), 
                iter_keymaj.clone(), 
                ring_operator, 
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
        let _barcode = crate::algebra::chains::barcode::barcode( &umatch, iter_keymaj, dim_fn, fil_fn, true , true);
    
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




    #[test]
    fn test_empty_simplex_iter() {
        use sprs::CsMatBase;
        use ordered_float::OrderedFloat;
        use crate::topology::simplicial::from::graph_weighted::SimplexIter;

        
        let dissimilarity_matrix = CsMatBase::new( (2,2), vec![0,1,2], vec![0,1], vec![OrderedFloat(0.0),OrderedFloat(0.0)] );

        println!("entry (0,1) = {:?}", (& dissimilarity_matrix).entry_major_at_minor(0, 1));
            
        let simplices: Vec<_>    =   SimplexIter::new( 
                            & dissimilarity_matrix, 
                            2,          //  matrix size
                            OrderedFloat(1.0),        //  dissimilarity value max
                            1,          //  dimension
                        ).collect();
        println!("generated simplices = {:?}", &simplices );
        let empty =  Vec::< SimplexFiltered<OrderedFloat<f64>> >::new();
                    
        assert!( empty.eq( &simplices ) )        
    }

}    