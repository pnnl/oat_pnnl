//! Define and enumerate filtered and unfiltered simplices

use oat_rust::utilities::iterators::general::{MapByTransform, MapByTransformTrait};
use oat_rust::utilities::functions::misc_functions::ReverseVector;
use oat_rust::utilities::iterators::merge::heap_of_iterators::{HitMergeByPredicateTrait, HitMerge};
use oat_rust::utilities::partial_order::OrderComparatorAutoGt;
use itertools::{Dedup, KMerge, Itertools, Combinations};
use std::cmp::Ordering;
use std::collections::{HashSet};
use std::hash::Hash;
use std::iter::{FromIterator, Flatten, Cloned, Rev};
use std::slice::Iter;
use std::vec::IntoIter;

use std::fmt::Debug;


type Vertex = u16;

//  ================================================================================
//  COMBINATORIAL SIMPLEX (UNWEIGHTED) -- DEFINITION
//  ================================================================================


/// An unweighted simplex; the vertices should sorted in ascending order.
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
    pub fn dim( &self ) -> usize { self.vertices.len() - 1 }
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
        return self.vertices.cmp( & other.vertices )
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
//  COMBINATORIAL SIMPLEX (WEIGHTED) -- DEFINITION
//  ================================================================================



/// A simplex associated with a filtration value
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct SimplexFiltered<FilVal: Clone + Debug>{
    pub filvalue: FilVal,     // the filtration value of simplex
    pub vertices: Vec<Vertex> // a sorted vector in strictly descending order recording vertices
}

impl <FilVal: Clone + Debug>
    SimplexFiltered<FilVal> {
    
    /// Dimension of the simplex.
    pub fn dim(&self) -> usize { self.vertices.len() - 1 }

    /// Filtraiton value of the simplex.
    pub fn fil(&self) -> FilVal { self.filvalue.clone() }
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
    fn cmp(&self, other: &Self) -> Ordering {

        // first compare filtration values
        let mut comp = self.filvalue.cmp( & other.filvalue );
        if comp != Ordering::Equal { return comp }

        // next compare simplex dimensions
        comp = self.vertices.len().cmp( & other.vertices.len() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        return self.vertices.cmp( & other.vertices )
    }
}




//  ================================================================================
//  FACETS OF SIMPLICES
//  ================================================================================





//  ---------------------------------------------------------------------------
//  FACES FROM FACETS-OF-THE-COMPLEX (NEW)
//  ---------------------------------------------------------------------------


/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// an iterator that runs over `dim`-dimensional subsimplices in **descending** lexicographic 
/// order.
/// 
/// # Examples
/// 
/// ```
/// use sohar_rust::simplicial::simplices::subsimplices_dim_d_iter_descend;
/// 
/// let dowker_simplices    =   vec![ vec![0,1,2], vec![0,3] ];
/// let one_simplices       =   subsimplices_dim_d_iter_descend( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![1,2], vec![0,3], vec![0,2], vec![0,1] ] );
/// ```
pub fn subsimplices_dim_d_iter_descend< Vertex >(
    complex_facets: & Vec< Vec< Vertex >>, 
    dim: usize 
) 
    -> 
    Dedup< 
            HitMerge< 
                    MapByTransform< 
                            Combinations< Cloned< Rev< Iter<'_, Vertex >>>>, Vec< Vertex >, ReverseVector,
                        >,
                    OrderComparatorAutoGt,
                >,         
        > 
    
    where Vertex: Ord + Clone
{
    complex_facets
        .iter()
        .map( |x| x.iter().rev().cloned().combinations( dim + 1 ).map_by_transform( ReverseVector::new() ) )
        .hit_merge_by_predicate( OrderComparatorAutoGt )
        .dedup()
}

/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// an iterator that runs over `dim`-dimensional subsimplices in **ascending** lexicographic 
/// order.
/// 
/// # Examples
/// 
/// ```
/// use sohar_rust::simplicial::simplices::subsimplices_dim_d_iter_ascend;
/// 
/// let dowker_simplices    =   vec![ vec![0,1,2], vec![0,3] ];
/// let one_simplices       =   subsimplices_dim_d_iter_ascend( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![0,1], vec![0,2], vec![0,3], vec![1,2], ] );
/// ```
pub fn subsimplices_dim_d_iter_ascend< Vertex >(
    complex_facets: & Vec< Vec< Vertex >>, 
    dim: usize 
) 
    -> 
    Dedup< KMerge<  Combinations<Cloned<Iter<Vertex>>> > >
    where Vertex: Ord + Clone
{
    // let x = vec![0,1,2].iter().rev().cloned().combinations(2);
    complex_facets
        .iter()
        .map( |x| x.iter().cloned().combinations( dim + 1 )  )
        .kmerge()
        .dedup()
}



/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension, then lexicographically, in ascending order.
/// 
/// # Examples
/// 
/// ```
/// use sohar_rust::simplicial::simplices::subsimplices_dim_0_thru_d_iter_ascend;
/// 
/// let dowker_simplices    =   vec![ vec![0,1,2], vec![0,3] ];
/// let one_simplices       =   subsimplices_dim_0_thru_d_iter_ascend( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![0], vec![1], vec![2], vec![3], vec![0,1], vec![0,2], vec![0,3], vec![1,2], ] );
/// ```
pub fn  subsimplices_dim_0_thru_d_iter_ascend< Vertex >( 
    complex_facets: & Vec< Vec< Vertex >>, 
    max_dim: usize 
) 
    -> 
    Flatten<IntoIter<Dedup< KMerge<  itertools::Combinations< Cloned< Iter<'_, Vertex>>  > > >>>
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .map(   |dim|
                subsimplices_dim_d_iter_ascend(complex_facets, dim)
            )
        .collect_vec()
        .into_iter()
        .flatten()
}

/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension, then lexicographically, in descending order.
/// 
/// # Examples
/// 
/// ```
/// use sohar_rust::simplicial::simplices::subsimplices_dim_0_thru_d_iter_descend;
/// 
/// let dowker_simplices    =   vec![ vec![0,1,2], vec![0,3] ];
/// let one_simplices       =   subsimplices_dim_0_thru_d_iter_descend( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![1,2], vec![0,3], vec![0,2], vec![0,1], vec![3], vec![2], vec![1], vec![0], ] );
/// ```
pub fn  subsimplices_dim_0_thru_d_iter_descend< Vertex >( 
    complex_facets: & Vec< Vec< Vertex >>, 
    max_dim: usize 
) 
    -> 
    Flatten<IntoIter<
            Dedup< 
                    HitMerge< 
                            MapByTransform< 
                                    Combinations< Cloned< Rev< Iter<'_, Vertex >>>>, Vec< Vertex >, ReverseVector,
                                >,
                            OrderComparatorAutoGt,
                        >,         
                >     
        >>
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .rev()
        .map(   |dim|
                subsimplices_dim_d_iter_descend( complex_facets, dim )
            )
        .collect_vec()
        .into_iter()
        .flatten()
}



/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// a vector of vectors whose kth entry contains all dimension-k subsimplices,
/// sorted in in lexicographic order.
pub fn  subsimplices_dim_0_thru_d_vecvec< Vertex >( 
    complex_facets: & Vec< Vec< Vertex >>, 
    max_dim: usize 
) 
-> 
Vec< Vec< Vec< Vertex >>> 
where Vertex: Ord + Clone
{
    let mut seq             =   Vec::with_capacity( max_dim );
    for dim in 0 .. max_dim + 1  {
        let vec: Vec<_>     =   subsimplices_dim_d_iter_ascend(
                                    complex_facets,
                                    dim
                                )
                                .collect();
        seq.push( vec );
    }
    seq
}

/// Similar to `subsimplices_dim_0_thru_d_vecvec`, but stores all simplices in a single vector, sorted first by dimension and second by lexicographic order.
pub fn subsimplices_dim_0_thru_d_concatenated< Vertex >( 
    complex_facets: & Vec< Vec< Vertex >>, 
    max_dim: usize 
) 
-> 
Vec< Vec< Vertex >>
    where Vertex: Ord + Clone
{
    let mut a = subsimplices_dim_0_thru_d_vecvec( complex_facets, max_dim );
    let mut b = Vec::new();
    for i in 0 .. a.len() {
        b.append( &mut a[ i ] );
    }
    b
}




//  ---------------------------------------------------------------------------
//  FACES FROM FACETS-OF-THE-COMPLEX ( OLD )
//  ---------------------------------------------------------------------------

/// Given something that iterates over vectors (each of which represents a strictly 
/// ascending sequence of vertices), return a HashSet containing all nonempty subsequences.
pub fn  set_of_subsequences< IterFacet, Vertex >( facets: IterFacet ) -> HashSet< Vec< Vertex > > 
    where   IterFacet:      IntoIterator< Item = Vec< Vertex > >,
            Vertex:    Ord + Hash + Clone
{
    println!("THIS FUNCTION COULD PROBABLY BE MADE MUCH MORE EFFICIENT");    
    let mut faces       =   HashSet::new();
    for facet in facets {
        for seq_length in 1 .. facet.len() {
            for comb in facet.iter().cloned().combinations( seq_length ) {
                faces.insert( comb );
            }
        }
    }
    faces
}

/// Given something that iterates over vectors (each of which represents a strictly 
/// ascending sequence of vertices), return a vector V containing all nonempty ordered
/// subsequences; V is strictly ascending under the order that first compares length of 
/// a sequence, then compares equal-length sequences lexicographically.
/// 
//  NB: THE USE OF SIMPLICES RATHER THAN VECTORS IS IMPORTANT HERE, BECAUSE THE TWO STRUCTS HAVE
//      **DIFFERENT** TOTAL ORDERS
pub fn  ordered_sequence_of_faces< IterFacet, Vertex >( facets: IterFacet ) -> Vec< Simplex< Vertex > > 
    where   IterFacet:  IntoIterator< Item = Vec< Vertex > >,
            Vertex:     Ord + Hash + Clone
{
    println!("THIS FUNCTION COULD PROBABLY BE MADE MUCH MORE EFFICIENT");
    let mut faces   =   set_of_subsequences(facets);
    let mut faces   =   Vec::from_iter( faces.drain().map(|x| Simplex{vertices: x}) );
    faces.sort();
    faces
}   

//  ---------------------------------------------------------------------------
//  FACETS-OF-A-SIMPLEX
//  ---------------------------------------------------------------------------

/// Maintains an "inner state" that steps through the facets of a simplex in 
/// ascending lexicographic order; only returns `Some(())` or `None`.
/// 
/// # Examples
/// 
/// ```

/// use sohar_rust::simplicial::simplices::{Simplex, FacetIteratorNoReturnAscending};

/// // Create the iterator
/// let mut facet_iterator_noreturn     =   FacetIteratorNoReturnAscending::new(
///                                             Simplex{ vertices: vec![0, 1, 2] },
///                                             None
///                                         );
///
/// // Test it                                                
/// let mut answers = vec![
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 2] }, deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![1, 2] }, deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![]     }, deleted_vertex_index: None    },            
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },                                                                        
/// ];
///
/// for i in 0..5 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }      
//
/// // Re-initialize with a new simplex
///
/// facet_iterator_noreturn.reinitialize_with_simplex( Simplex{ vertices: vec![0 ,3]} );
///
/// // Test again        
///
/// answers = vec![
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![3] }, deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![]  }, deleted_vertex_index: None    },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },                                                                                                      
/// ];    
///
/// for i in 0..4 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }   
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct FacetIteratorNoReturnAscending< Vertex >
{
    pub simplex: Simplex< Vertex > ,
    pub facet: Simplex< Vertex >,
    pub deleted_vertex_index: Option<usize>
}

impl < Vertex > FacetIteratorNoReturnAscending < Vertex > 
    where Vertex: Clone
{

    /// Initialize a no-return facet iterator.  
    /// 
    /// When it is first constructed, the internal state of the iterator does not represent a facet.
    /// The internal state will be updated to represent the first facet when `next()` is called for the first time.
    pub fn new( simplex: Simplex< Vertex >, buffer: Option< Simplex<Vertex> > ) -> Self {
        let buff = 
            if let Some( vec ) = buffer { vec } 
            else { 
                Simplex{ 
                    vertices: Vec::with_capacity( simplex.dim() ) // = 1 less than num_vertices = NUM VERTICES OF A FACET
                }             
            };
        
        FacetIteratorNoReturnAscending {
            simplex: simplex,
            facet: buff,
            deleted_vertex_index: None
        }
    }

    /// Reinitialize a no-return facet iterator with a new simplex.
    /// 
    /// The result is essentially the same as replacing our iterator with `FacetIteratorNoReturnAscending::new( simplex )`.
    pub fn reinitialize_with_simplex( &mut self, simplex: Simplex< Vertex > ) {
        // if necessary, expand the capacity of the facet vector
        if simplex.dim() > self.facet.vertices.capacity() { 
            self.facet.vertices.reserve_exact(
                simplex.dim() - self.facet.vertices.capacity()
            ) 
        }
        // replace the old simplex with the new
        self.simplex    =   simplex;
        // update the state to indicate that it does not represent a facet
        self.facet.vertices.clear();
        self.deleted_vertex_index = None;
    }
}


impl < Vertex >
    Iterator for 
    FacetIteratorNoReturnAscending < Vertex >     
    where Vertex : Clone
{
    type Item    =   ();

    fn next( &mut self ) -> Option<()> {

        if let Some( deleted_index ) = self.deleted_vertex_index {

            if deleted_index == 0 {
                // if we start from the facet obtained by deleting vertex 0, then the 
                // next state should **not** represent a facet
                self.deleted_vertex_index   =   None;
                self.facet.vertices.clear();
                return None
                
            } else {
                // if we start from the facet obtained by deleting vertex k > 0, then 
                // the next state should represent the facet obtained by deleting vertex k-1
                let next_deleted_index  =   deleted_index - 1;
                self.facet.vertices[ next_deleted_index ] = self.simplex.vertices[ deleted_index ].clone(); // replace the deleted index and remove the next one
                self.deleted_vertex_index = Some( next_deleted_index );
                return Some( () )
            }
        
        } else {

            self.facet.vertices.clear();
            for i in 0..self.simplex.dim() {   // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
                self.facet.vertices.push( self.simplex.vertices[ i ].clone() ) 
            }      
            // set deleted vertex equal to last
            self.deleted_vertex_index = Some( self.simplex.dim() );  // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
            return Some(())            
        }
    }
}




//  ===========================================================================
//  ===========================================================================
//  TESTS
//  ===========================================================================
//  ===========================================================================





#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_ordered_subsimplices_up_thru_dim() {

        let complex_facets          =   vec![ vec![0, 1, 2] ];

        assert_eq!(         subsimplices_dim_0_thru_d_vecvec( & complex_facets, 2),
                            vec![
                                vec![   vec![0],     vec![1],    vec![2]         ],                                
                                vec![   vec![0,1],   vec![0,2],  vec![1,2]       ],
                                vec![   vec![0,1,2]                              ]
                            ]
        );

        assert_eq!(         subsimplices_dim_0_thru_d_concatenated( & complex_facets, 2),
                            vec![
                                        vec![0],     vec![1],    vec![2],                                         
                                        vec![0,1],   vec![0,2],  vec![1,2],       
                                        vec![0,1,2]                              
                            ]
        ) ;       


    }


    #[test]
    fn test_ascending_facet_iterator_no_return()
    {

        // Create the iterator
        let mut facet_iterator_noreturn     =   FacetIteratorNoReturnAscending::new(
                                                    Simplex{ vertices: vec![0, 1, 2] },
                                                    None
                                                );

        // Test it                                                
        let mut answers = vec![
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 2] }, deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![1, 2] }, deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![]     }, deleted_vertex_index: None    },            
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },                                                                        
        ];

        for i in 0..5 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }      

        // Re-initialize with a new simplex

        facet_iterator_noreturn.reinitialize_with_simplex( Simplex{ vertices: vec![0 ,3]} );

        // Test again        

        answers = vec![
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![3] }, deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![]  }, deleted_vertex_index: None    },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },                                                                                                      
        ];    

        for i in 0..4 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }     
                
    }       

}