//! 
//! 
//! TO DO: RIGHT NOW WE DUPLICATE DISMAT (ONCE IN THE PARAM STRUCT AND ONCE IN THE MATRIX); TRY TO REMOVE ONE OF THESE COPIES
//! 

use itertools::Itertools;
use ordered_float::OrderedFloat;
use oat_rust::matrices::operations::umatch::row_major::{new_umatchrowmajor_with_clearing, UmatchRowMajor};
use oat_rust::rings::operator_traits::{Semiring, DivisionRing, Ring};
use oat_rust::utilities::partial_order::{OrderComparatorAutoLt, OrderComparatorLtByKey};
use oat_rust::utilities::distances::{minmax};
use std::fmt::Debug;
use std::marker::PhantomData;

use crate::simplicial::simplices::SimplexFiltered;

use super::data_structures::{CliqueBoundaryMatrix, SimplexIter};

/// A struct to hold the parameters of a filtered clique complex.  Constructed via `new(..)`.
/// 
/// - `dismat`: a disimilarity matrix
/// - `maxdim`: maximum homology dimension (will construct simplices of dimension 0 through maxdim+1, inclusive)
/// - `maxdis`: maximum dissimilarity threshold
/// - `ring_operator`: ring operator for the coefficient ring
#[derive(Clone, Debug)]
pub struct CliqueParams
                <RingOperator, RingElement>
    where 
        RingOperator:   Clone + Semiring< RingElement > + Ring< RingElement > + DivisionRing< RingElement >,
        RingElement:    Clone + Debug,
{
    dismat:                 Vec<Vec< OrderedFloat<f64> >>,
    maxdim:                 usize,
    maxdis:                 OrderedFloat<f64>,
    ring_operator:          RingOperator,
    phantom_ringelement:    PhantomData<RingElement>,
}

impl <RingOperator, RingElement>
    
    CliqueParams
        <RingOperator, RingElement>

    where 
        RingOperator:   Clone + Semiring< RingElement > + Ring< RingElement > + DivisionRing< RingElement >,
        RingElement:    Clone + Debug,        
    {
    /// Construct nicely formatted parameters for the clique condstruction.
    pub fn new(
        dismat:                 & Vec<Vec< f64 >>,
        maxdim:                 usize,
        maxdis:                 Option<f64>,
        ring_operator:          RingOperator,
        ) ->
        Self {
        
        let m = dismat.len();


        let mut dismat_ordf64 = vec![ vec![OrderedFloat(0.0); m]; m ]; // preallocate matrix filled with zeros

        for i in 0 .. m {
            for j in i .. m {
                if dismat[i][j] != dismat[j][i] {
                    panic!("Input matrix is not symmetric: entry {:?} doesn't equal entry {:?}", (i,j), (j,i));
                } else {
                    dismat_ordf64[i][j] = OrderedFloat(dismat[i][j]);
                    dismat_ordf64[j][i] = dismat_ordf64[i][j];
                }
            }
        }        

        // convert dismat to ordered floats
        // let dismat =  
        //     (0..m)
        //         .map(   |x| 
        //                 (0..m)
        //                     .map(|y| OrderedFloat(dismat[x][y]) ) 
        //                     .collect_vec()
        //             )
        //         .collect_vec();
    
        // if no maxdis is provided then use the enclosing radius
        let maxdis = match maxdis {
                Some( t ) => OrderedFloat(t),
                None => {  minmax( & dismat_ordf64 ).clone() }      
            }; 
        CliqueParams { dismat: dismat_ordf64, maxdim, maxdis, ring_operator, phantom_ringelement: PhantomData }       
    }

    pub fn dismat_ref(&self) -> & Vec<Vec< OrderedFloat<f64> >> { &self.dismat }
    pub fn maxdim( &self ) -> usize { self.maxdim.clone() }
    pub fn maxdis( &self ) -> OrderedFloat<f64> { self.maxdis.clone() }
    pub fn ring_operator( &self ) -> RingOperator { self.ring_operator.clone() }
}    


/// A vector of simplices sorted first by dimension (ascending `0..params.maxdim`) then by diameter (descending) then by lexicographic order
pub fn cliques_in_order<RingOperator, RingElement>( 
            params: & CliqueParams<RingOperator, RingElement> 
        ) 
        -> 
        Vec< SimplexFiltered<OrderedFloat<f64>> >
    where 
        RingOperator:   Clone + Semiring< RingElement > + Ring< RingElement > + DivisionRing< RingElement >,
        RingElement:    Clone + Debug,        
    {
        let vec: Vec<_> = (0..params.maxdim+1).map(
            |dim|
            {
                let mut vec = SimplexIter::new( 
                            dim,
                            & params.dismat,
                            params.maxdis,                
                        )
                        .collect_vec();
                vec.sort_by(|x,y| y.cmp(x) );
                vec
            }
        )
        .flatten()
        .collect();
        return vec
}


/// Construct the (lazy) boundary matrix of a filtered clique complex.
/// 
/// - `dismat`: a disimilarity matrix
/// - `maxdis`: maximum dissimilarity threshold
/// - `ring_operator`: ring operator for the coefficient ring
pub fn get_clique_boundary_matrix< RingOperator, RingElement >(
        dismat:     & Vec<Vec<f64>>,  // dissimilarity matrix
        maxdis:     Option<f64>,
        ring_operator:  RingOperator
    )
    ->
    CliqueBoundaryMatrix<OrderedFloat<f64>, RingElement, RingOperator> 

    where 
        RingOperator:   Clone + Semiring< RingElement > + Ring< RingElement > + DivisionRing< RingElement >,
        RingElement:    Clone + Debug,
{
    let params = CliqueParams::new( dismat, 0, maxdis, ring_operator ); // max dimension doesn't matter
    
    // let m = dismat.len();

    // // convert dismat to ordered floats
    // let dismat =  
    //     (0..m)
    //         .map(   |x| 
    //                 (0..m)
    //                     .map(|y| OrderedFloat(dismat[x][y]) ) 
    //                     .collect_vec()
    //             )
    //         .collect_vec();
    // // let dismat2 = dismat.clone();

    // // if no maxdis is provided then use the enclosing radius
    // let maxdis = match maxdis {
    //         Some( t ) => OrderedFloat(t),
    //         None => {  minmax( & dismat ).clone() }      
    //     };
    // let cutoff_matrix = get_cutoff_matrix( & dismat, maxdis );          

    // let boundary_matrix = 
    //     CliqueBoundaryMatrix
    //         {
    //             ring_operator: ring_operator.clone(),
    //             /// The "maxdis" value represents the maximum of filtration value
    //             maxdis,
    //             /// A vector representing the dissimilarity matrix by listing all its rows
    //             dismat,
    //             /// A vector representing the neighborhood within "maxdis" of each vertex
    //             cutoff_matrix,
    //             phantomsnzval: PhantomData,
    //         };
    return CliqueBoundaryMatrix::new( & params );
}



/// U-match factorization of a filtered clique complex.
/// 
/// - `boundary_matrix`: a filtered clique boundary matrix
/// - `maxdim`: maximum homology dimension (will construct simplices of dimension 0 through maxdim+1, inclusive)
pub fn umatch_from_clique< 'a, RingOperator, RingElement >(        
        boundary_matrix:    &'a CliqueBoundaryMatrix<OrderedFloat<f64>, RingElement, RingOperator>,
        maxdim:             usize,
    )
    ->
    UmatchRowMajor<
            &'a  CliqueBoundaryMatrix<OrderedFloat<f64>, RingElement, RingOperator>, 
            RingOperator, 
            OrderComparatorLtByKey<SimplexFiltered<OrderedFloat<f64>>, RingElement, (SimplexFiltered<OrderedFloat<f64>>, RingElement), OrderComparatorAutoLt<SimplexFiltered<OrderedFloat<f64>>>>, 
            OrderComparatorLtByKey<SimplexFiltered<OrderedFloat<f64>>, RingElement, (SimplexFiltered<OrderedFloat<f64>>, RingElement), OrderComparatorAutoLt<SimplexFiltered<OrderedFloat<f64>>>>
        >
    where 
        RingOperator:   Clone + Semiring< RingElement > + Ring< RingElement > + DivisionRing< RingElement >,
        RingElement:    Clone + Debug,
{
    // let m = dismat.len();

    // // convert dismat to ordered floats
    // let dismat =  
    //     (0..m)
    //         .map(   |x| 
    //                 (0..m)
    //                     .map(|y| OrderedFloat(dismat[x][y]) ) 
    //                     .collect_vec()
    //             )
    //         .collect_vec();
    // let dismat2 = dismat.clone();

    // // if no maxdis is provided then use the enclosing radius
    // let maxdis = match maxdis {
    //         Some( t ) => OrderedFloat(t),
    //         None => {  minmax( & dismat ).clone() }      
    //     };
    // let cutoff_matrix = get_cutoff_matrix( & dismat, maxdis );          
    
    // let boundary_matrix = 
    //     CliqueBoundaryMatrix
    //         {
    //             ring_operator: ring_operator.clone(),
    //             /// The "maxdis" value represents the maximum of filtration value
    //             maxdis,
    //             /// A vector representing the dissimilarity matrix by listing all its rows
    //             dismat,
    //             /// A vector representing the neighborhood within "maxdis" of each vertex
    //             cutoff_matrix,
    //             phantomsnzval: PhantomData,
    //         };

    let iter_keymaj = 
        (0..maxdim+1).map(
            |dim|
            {
                let mut vec = SimplexIter::new( 
                            dim,
                            & boundary_matrix.dismat,
                            boundary_matrix.maxdis,                
                        )
                        .collect_vec();
                vec.sort_by(|x,y| y.filvalue.cmp(&x.filvalue));
                vec
            }
        )
        .flatten();

    let umatch = new_umatchrowmajor_with_clearing(
                boundary_matrix, 
                iter_keymaj.clone(), 
                boundary_matrix.ring_operator.clone(), 
                OrderComparatorAutoLt::new(), 
                OrderComparatorAutoLt::new(), 
            );          
    
    return umatch
    
}