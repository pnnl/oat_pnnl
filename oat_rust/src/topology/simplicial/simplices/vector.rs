//! Simplices represented by Rust vectors

use derive_new::new;


use crate::utilities::iterators::merge::hit::{HitMergeByPredicateTrait, HitMerge};
use crate::utilities::order::{OrderOperatorAutoReverse, JudgeOrder, JudgePartialOrder};
use crate::utilities::sequences_and_ordinals::{CombinationsReverse, SortedVec};
use itertools::{Dedup, KMerge, Itertools, Combinations};
use std::cmp::Ordering;
use std::collections::{HashSet};
use std::hash::Hash;
use std::iter::{Flatten, Cloned};
use std::slice::Iter;

use std::fmt::Debug;



//  ================================================================================
//  TWIST ORDER OPERATOR - UNFILTERED SIMPLEX
//  ================================================================================

/// Orders simplices first by dimension (descending) then by lexicographic order (ascending)
/// 
/// This order is specially designed to facilitate the `twist` optimization in 
/// persistent homology calculations.
/// 
/// This object implements the `JudgeOrder` and `JudgePartialOrder` traits.
#[derive(Clone,Debug,new)]
pub struct OrderOperatorTwistSimplex;

impl < Vertex: Ord >

    JudgePartialOrder
        < Vec< Vertex > > for 

    OrderOperatorTwistSimplex
{
    fn judge_partial_cmp
        ( & self, lhs: & Vec< Vertex >, rhs: & Vec< Vertex > ) -> Option<Ordering>
    {   Some( self.judge_cmp(lhs,rhs) )      }
}


impl < Vertex: Ord >

    JudgeOrder
        < Vec< Vertex > > for 

    OrderOperatorTwistSimplex
{
    fn judge_cmp
        ( & self, lhs: & Vec< Vertex >, rhs: & Vec< Vertex > ) -> Ordering
    {
        // compare simplex dimensions
        let comp = rhs.len().cmp( & lhs.len() ); // note we put RHS on the LEFT here
        if comp != Ordering::Equal { return comp }
        // then compare simplices lexicographically
        lhs.cmp( rhs )        
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
/// use oat_rust::topology::simplicial::simplices::vector::subsimplices_dim_d_iter_descend;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let one_simplices       =   subsimplices_dim_d_iter_descend( &dowker_simplices, 1 ).unwrap();
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![1,2], vec![0,3], vec![0,2], vec![0,1] ] );
/// ```
pub fn subsimplices_dim_d_iter_descend< Vertex >(
    complex_facets: &Vec< SortedVec< Vertex >>, 
    dim: isize 
) 
    -> 
    Option< 
        Dedup< 
                HitMerge< 
                        CombinationsReverse< Vertex, &Vec< Vertex > >,
                        OrderOperatorAutoReverse,
                    >,         
            > 
        >
    
    where Vertex: Ord + Clone
{
    if dim < -1 { return None }
    Some( 
            complex_facets
                .iter()
                .map( |vertices| CombinationsReverse::from_vec( dim as usize +1, vertices.vec() ) )
                .hit_merge_by_predicate( OrderOperatorAutoReverse )
                .dedup()
        )
}

/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// an iterator that runs over `dim`-dimensional subsimplices in **ascending** lexicographic 
/// order.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::subsimplices_dim_d_iter_ascend;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let one_simplices       =   subsimplices_dim_d_iter_ascend( &dowker_simplices, 1 ).unwrap();
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![0,1], vec![0,2], vec![0,3], vec![1,2], ] );
/// ```
pub fn subsimplices_dim_d_iter_ascend< Vertex >(
    complex_facets: & Vec< SortedVec< Vertex >>, 
    dim: isize 
) 
    -> 
    Option< 
            Dedup< KMerge<  Combinations<Cloned<Iter<Vertex>>> > > 
        > 
    where Vertex: Ord + Clone
{
    if dim < -1 { return None }    
    Some( 
            complex_facets
                .iter()
                .map( |x| x.vec().iter().cloned().combinations( (dim + 1isize) as usize )  )
                .kmerge()
                .dedup()
        )
}



/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension, then lexicographically, in ascending order.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::subsimplices_dim_0_thru_d_iter_ascend;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let all_simplices       =   subsimplices_dim_0_thru_d_iter_ascend( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( all_simplices, vec![ vec![0], vec![1], vec![2], vec![3], vec![0,1], vec![0,2], vec![0,3], vec![1,2], ] );
/// ```
pub fn subsimplices_dim_0_thru_d_iter_ascend< Vertex >( 
    complex_facets: & Vec< SortedVec< Vertex >>, 
    max_dim: isize 
) 
    -> 
    Flatten< std::vec::IntoIter<Dedup< KMerge<  itertools::Combinations< Cloned< Iter<'_, Vertex>>  > > >>>
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .map(   |dim|
                subsimplices_dim_d_iter_ascend(complex_facets, dim).unwrap()
            )
        .collect_vec()
        .into_iter()
        .flatten()
}

/// Assuming the user-provided facets have vertices sorted in strictly ascending order, returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension (strictly descending order), then lexicographically (strictly descending order).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::subsimplices_dim_0_thru_d_iter_descend;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let all_simplices       =   subsimplices_dim_0_thru_d_iter_descend( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( all_simplices, vec![ vec![1,2], vec![0,3], vec![0,2], vec![0,1], vec![3], vec![2], vec![1], vec![0], ] );
/// ```
pub fn subsimplices_dim_0_thru_d_iter_descend< Vertex >( 
    complex_facets: &Vec< SortedVec< Vertex >>, 
    max_dim: isize 
) 
    -> 
    Flatten< std::vec::IntoIter<
            Dedup< 
                    HitMerge< 
                            CombinationsReverse< Vertex, &Vec< Vertex > >,
                            OrderOperatorAutoReverse,
                        >,  
                >     
        >>
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .rev()
        .map(   |dim|
                subsimplices_dim_d_iter_descend( complex_facets, dim ).unwrap()
            )
        .collect_vec()
        .into_iter()
        .flatten()
}


/// Assuming the user-provided facets have vertices sorted in strictly ascending order, returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension (strictly ascending order), then lexicographically (strictly descending order).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let all_simplices       =   subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( all_simplices, vec![ vec![3], vec![2], vec![1], vec![0], vec![1,2], vec![0,3], vec![0,2], vec![0,1], ] );
/// ```
pub fn subsimplices_dim_0_thru_d_iter_ascend_dim_descend_lex< Vertex >( 
        complex_facets: &Vec< SortedVec< Vertex >>, 
        max_dim: isize 
    ) 
    -> 
    Flatten< std::vec::IntoIter<
            Dedup< 
                    HitMerge< 
                            CombinationsReverse< Vertex, &Vec< Vertex > >,
                            OrderOperatorAutoReverse,
                        >,  
                >     
        >>
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .map(   |dim|
                subsimplices_dim_d_iter_descend( complex_facets, dim ).unwrap()
            )
        .collect_vec()
        .into_iter()
        .flatten()
}



/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// a vector of vectors whose kth entry contains all dimension-k subsimplices,
/// sorted in in lexicographic order.
pub fn subsimplices_dim_0_thru_d_vecvec< Vertex >( 
            complex_facets: & Vec< SortedVec< Vertex >>, 
            max_dim: isize 
        ) 
        -> 
        Vec< Vec< Vec< Vertex >>> 
        where Vertex: Ord + Clone
{
    if max_dim < 0 { return Vec::new() } // return the empty set if max_dim < 0
    let mut seq             =   Vec::with_capacity( max_dim as usize );
    for dim in 0 .. max_dim + 1  {
        let vec: Vec<_>     =   subsimplices_dim_d_iter_ascend(
                                    complex_facets,
                                    dim
                                ).unwrap()
                                .collect();
        seq.push( vec );
    }
    seq
}

/// Similar to `subsimplices_dim_0_thru_d_vecvec`, but stores all simplices in a single vector, sorted first by dimension and second by lexicographic order.
pub fn subsimplices_dim_0_thru_d_concatenated< Vertex >( 
    complex_facets: & Vec< SortedVec< Vertex >>, 
    max_dim: isize 
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
pub fn set_of_subsequences< IterFacet, Vertex >( facets: IterFacet ) -> HashSet< Vec< Vertex > > 
    where   IterFacet:      IntoIterator< Item = Vec< Vertex > >,
            Vertex:    Ord + Hash + Clone
{
    println!("THIS FUNCTION COULD PROBABLY BE MADE MUCH MORE EFFICIENT");    
    let mut faces       =   HashSet::new();
    for facet in facets {
        for seq_length in 1 .. facet.len() {
            for comb in facet.iter().cloned().combinations( seq_length ) {
                // faces.insert( SortedVec::new(comb) );
                faces.insert( comb );
            }
        }
    }
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

/// use oat_rust::topology::simplicial::simplices::vector::FacetIteratorNoReturnAscendingLex;

/// // Create the iterator
/// let mut facet_iterator_noreturn     
///     =   FacetIteratorNoReturnAscendingLex::new( vec![0, 1, 2], None );
///
/// // Test it                                                
/// let mut answers = vec![
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![0, 1], deleted_vertex_index: Some(2) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![0, 2], deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![1, 2], deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![]    , deleted_vertex_index: None    },            
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![0, 1], deleted_vertex_index: Some(2) },                                                                        
/// ];
///
/// for i in 0..5 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }      
//
/// // Re-initialize with a new simplex
///
/// facet_iterator_noreturn.reinitialize_with_simplex( vec![0 ,3] );
///
/// // Test again        
///
/// answers = vec![
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![3], deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![] , deleted_vertex_index: None    },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },                                                                                                      
/// ];    
///
/// for i in 0..4 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }   
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct FacetIteratorNoReturnAscendingLex< Vertex: Ord >
{
    pub simplex: Vec< Vertex > ,
    pub facet: Vec< Vertex >,
    pub deleted_vertex_index: Option<usize>
}

impl < Vertex: Ord > FacetIteratorNoReturnAscendingLex < Vertex > 
    where Vertex: Clone
{

    /// Initialize a no-return facet iterator.  
    /// 
    /// When it is first constructed, the internal state of the iterator does not represent a facet.
    /// The internal state will be updated to represent the first facet when `next()` is called for the first time.
    pub fn new( simplex: Vec< Vertex >, buffer: Option< Vec<Vertex> > ) -> Self {
        let buffer = buffer.unwrap_or( Vec::with_capacity( simplex.len()-1 ) );
        
        FacetIteratorNoReturnAscendingLex {
            simplex,
            facet: buffer,
            deleted_vertex_index: None
        }
    }

    /// Reinitialize a no-return facet iterator with a new simplex.
    /// 
    /// The result is essentially the same as replacing our iterator with `FacetIteratorNoReturnAscendingLex::new( simplex )`.
    pub fn reinitialize_with_simplex( &mut self, simplex: Vec< Vertex > ) {
        // if necessary, expand the capacity of the facet vector
        if simplex.len() > 1 + self.facet.capacity() { 
            self.facet.reserve_exact(
                simplex.len() - 1 - self.facet.capacity()
            ) 
        }
        // replace the old simplex with the new
        self.simplex    =   simplex;
        // update the state to indicate that it does not represent a facet
        self.facet.clear();
        self.deleted_vertex_index = None;
    }
}


impl < Vertex: Ord >
    Iterator for 
    FacetIteratorNoReturnAscendingLex < Vertex >     
    where Vertex : Clone
{
    type Item    =   ();

    fn next( &mut self ) -> Option<()> {

        if let Some( deleted_index ) = self.deleted_vertex_index {

            if deleted_index == 0 {
                // if we start from the facet obtained by deleting vertex 0, then the 
                // next state should **not** represent a facet
                self.deleted_vertex_index   =   None;
                self.facet.clear();
                None
                
            } else {
                // if we start from the facet obtained by deleting vertex k > 0, then 
                // the next state should represent the facet obtained by deleting vertex k-1
                let next_deleted_index  =   deleted_index - 1;
                self.facet[ next_deleted_index ] = self.simplex[ deleted_index ].clone(); // replace the deleted index and remove the next one
                self.deleted_vertex_index = Some( next_deleted_index );
                Some( () )
            }
        
        } else {

            self.facet.clear();
            for i in 0..self.simplex.len()-1 {   // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
                self.facet.push( self.simplex[ i ].clone() ) 
            }      
            // set deleted vertex equal to last
            self.deleted_vertex_index = Some( self.simplex.len()-1 );  // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
            Some(())            
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

        let complex_facets          =   vec![ 
                                                                SortedVec::new(vec![0, 1, 2]).unwrap()
                                                            ];

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
        let mut facet_iterator_noreturn     =   FacetIteratorNoReturnAscendingLex::new(
                                                    vec![0, 1, 2],
                                                    None
                                                );

        // Test it                                                
        let mut answers = vec![
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![0, 1], deleted_vertex_index: Some(2) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![0, 2], deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![1, 2], deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![]    , deleted_vertex_index: None    },            
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![0, 1], deleted_vertex_index: Some(2) },                                                                        
        ];

        for i in 0..5 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }      

        // Re-initialize with a new simplex

        facet_iterator_noreturn.reinitialize_with_simplex( vec![0 ,3] );

        // Test again        

        answers = vec![
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![3], deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![] , deleted_vertex_index: None    },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },                                                                                                      
        ];    

        for i in 0..4 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }     
                
    }       

}