//! Bars and barcodes for persistence modules

use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::algebra::chains::jordan::JordanBasisMatrix;
use crate::algebra::vectors::entries::{KeyValSet, KeyValNew};
use crate::algebra::matrices::operations::multiply::vector_matrix_multiply_minor_descend_simplified;
use crate::algebra::matrices::operations::umatch::row_major::Umatch;
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};

use crate::utilities::order::JudgePartialOrder;

// use polars::prelude::*;

use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::hash::Hash;

use derive_getters::{Getters, Dissolve};




/// A single bar in a persistent homology barcode (optionally including generators)
/// 
/// This object stores information about
/// 
/// - `id_number`: an integer used to distinguish this bar from all other bars in the barcode
/// - `birth`: the left endpoint of the bar
/// - `death`: the right endpoint of the bar (if it exists)
/// - `dimension`: the dimension of the homology group
/// - `birth_column`: in the special case of simplicial complexes, this is called the "birth simplex"
/// - `death_column`: in the special case of simplicial complexes, this is called the "death simplex" (if it exists)
/// - `cycle_representative`: a cycle representative in homology, which is born when the bar is born
/// - `bounding_chain`: a bounding chain, which is born when the bar dies
/// 
/// The primary motivation for this object is convenience -- reducing the number and complexity of
/// method calls needed to access commonly used properties of the barcode.
#[derive(Getters,Dissolve)]
#[derive(Debug,Clone)]
pub struct Bar< Index, Entry > {
    id_number:              usize,
    #[getter(skip)]
    dimension:              isize,
    birth:                  OrderedFloat< f64 >,    
    birth_column:           Index,
    death:                  Option< OrderedFloat< f64 > >,    
    death_column:           Option< Index >,
    cycle_representative:   Option< Vec< Entry > >,
    bounding_chain:         Option< Vec< Entry > >,
}

impl < Index, Entry >

    Bar
        < Index, Entry > {            
    pub fn dimension( &self ) -> isize { self.dimension }
    /// The birth fitration value, as an `OrderedFloat<f64>`
    pub fn birth_ordf64(&self)  -> OrderedFloat<f64> { self.birth }
    /// The death fitration value, as an `OrderedFloat<f64>`
    pub fn death_ordf64(&self)  -> OrderedFloat<f64> { self.death.unwrap_or( OrderedFloat( f64::INFINITY ) ) }    
    /// The birth fitration value, as an `f64`
    pub fn birth_f64( &self ) -> f64 { self.birth_ordf64().into_inner() }
    /// The death fitration value, as an `f64`
    pub fn death_f64( &self ) -> f64 { self.death_ordf64().into_inner() }
    /// The birth-death pair as `(f64,f64)`
    pub fn interval_f64( &self ) -> (f64,f64) {  (self.birth_f64(), self.death_f64()) }
    /// The birth-death pair as `(OrderedFloat<f64>,OrderedFloat<f64>)`
    pub fn interval_ordf64( &self ) -> (OrderedFloat<f64>,OrderedFloat<f64>) {  (self.birth_ordf64(), self.death_ordf64()) }    
}

/// Wrapper for a vector of [`Bar`]s, with handy utilities.
#[derive(Debug,Clone)]
pub struct Barcode
            < Index, Entry >
    { bars: Vec< Bar< Index, Entry > > }

impl < Index, Entry > 
    
    Barcode 
        < Index, Entry > {

    /// Construct a barcode from a list of bars
    pub fn new< I >( bars: I ) -> Self 
        where 
            I: IntoIterator< Item=Bar<Index,Entry> >
    { Barcode{ bars: bars.into_iter().collect() } }

    /// An iterator that runs over all the bars
    pub fn iter(&self) -> std::slice::Iter<'_, Bar<Index, Entry>> {
        self.bars.iter()
    }

    /// Immutable reference to the internally stored vector of [`Bar`]s
    pub fn bars( &self ) -> & Vec< Bar< Index, Entry > > { &self.bars }

    /// Vector containing references to all bars of dimension `dim`
    pub fn bars_in_dim( &self, dim:isize ) -> Vec< & Bar< Index, Entry > > { 
        self.bars.iter().filter(|x| x.dimension()==dim ).collect() 
    }

    /// Immutable reference to the `bafr_id_number`th internally stored bar.
    pub fn bar( &self, bar_id_number: usize ) -> & Bar< Index, Entry > { &self.bars[bar_id_number] }

    /// Return a sorted list of all endpoints of intervals in the barcode.
    pub fn endpoints_ordf64_within_dim( &self, dim: isize ) -> Vec< OrderedFloat<f64> > {
        let iter_a = self.bars.iter().filter(|x|x.dimension()==dim).map(|x| x.birth_ordf64() );
        let iter_b = self.bars.iter().filter(|x|x.dimension()==dim).map(|x| x.death_ordf64() );        
        let mut finite_endpoints = HashSet::new();
        for val in iter_a { finite_endpoints.insert( val ); }
        for val in iter_b { finite_endpoints.insert( val ); }
        let mut finite_endpoints = finite_endpoints.drain().collect_vec();
        finite_endpoints.sort();
        finite_endpoints
    }

    /// Return a sorted list of all endpoints of intervals in the barcode
    /// 
    /// The search runs over all bars in the barcode
    pub fn endpoints_ordf64( &self ) -> Vec< OrderedFloat<f64> > {
        let iter_a = self.bars.iter().map(|x| x.birth_ordf64() );
        let iter_b = self.bars.iter().map(|x| x.death_ordf64() );        
        let mut finite_endpoints = HashSet::new();
        for val in iter_a { finite_endpoints.insert( val ); }
        for val in iter_b { finite_endpoints.insert( val ); }
        let mut finite_endpoints = finite_endpoints.drain().collect_vec();
        finite_endpoints.sort();
        finite_endpoints
    }    

    /// Either the `Some( maximum of the set { finite endpoints in intervals in the barcode })`,
    /// or `None`, if there exists no finite endpoint.
    /// 
    /// The search runs over bars of every dimension.
    pub fn max_finite_endpoint( &self ) -> Option< OrderedFloat<f64> > {
        let endpoints = self.endpoints_ordf64();
        let limit = &OrderedFloat(f64::INFINITY);
        return endpoints.iter().cloned().filter(|x| x < limit ).max()
    }

    /// Returns a vector of triples `(id, birth, death)`, where `id` is the uniue id of the bar.
    pub fn intervals_f64( &self, dim: isize ) -> Vec< ( f64, f64, usize ) > { 
        self.bars
            .iter()
            .filter_map(|x|
                match x.dimension == dim {
                    true      => Some( ( x.birth_f64(),  x.death_f64(), *x.id_number(), ) ) ,
                    false     => None
                }                
            ) 
            .collect_vec()
    }

    /// A vector of tuples `(t, betti_t)` where `t` is an endpoint of an
    /// interval in the barcode, and `betti_t` is the dimension `dim` betti
    /// number at filtration parameter `t`.
    pub fn betti_curve( &self, dim: isize ) -> Vec< ( OrderedFloat<f64>, usize ) > {
        println!(" SEE ALSO THE OLD IMPLEMENTATION OF BETTI CURVE -BOTTOM OF THIS FILE");

        let mut hash     =   HashMap::new();
        let endpoints = self.endpoints_ordf64_within_dim(dim);
        for endpoint in endpoints {
            for bar in self.bars() {
                if dim != bar.dimension() { continue }
                if endpoint < *bar.birth() { continue }
                if let & Some(death_time) = bar.death() {
                    if *endpoint < *death_time { let val = hash.get_mut(&endpoint).unwrap(); *val += 1 }
                } else {
                    let val = hash.get_mut(&endpoint).unwrap(); *val += 1
                }
            }
        }
        let mut curve = hash.drain().collect_vec();
        curve.sort();
        curve
    }

    // /// Export to a Polars dataframe
    // pub fn into_polars( self ) -> DataFrame {
    //     unzip_n!(pub 8);

    //     // let dissolve = |x: &Bar<Index,Entry>| (
    //     //     x.id_number, x.dimension, x.birth.into_inner(), x.death.map(|x| x.into_inner()),
    //     //     x.birth_column, x.death_column, x.cycle_representative,
    //     //     x.bounding_chain
    //     // ); 
    //     let (id, dim, birth, death, 
    //         birth_column, death_column, 
    //         cycle, bounding ) 
    //         = self.bars.into_iter().map( |x| x.dissolve() ).unzip_n_vec();
    //     // let sid: polars::series::Series = id.iter().map(|x| x.clone() as u64).collect();
    //     df![
    //         "id"=> &id, "dim"=>dim, "birth"=>birth, "death"=>death, 
    //         "birth_column"=>birth_column, "death_column"=>death_column,
    //         "cycle"=>cycle, "bounding"=>bounding
    //     ].ok().unwrap()
    // }
}

/// Extract the barcode of a filtered chain complex from a U-match factorization of its boundary matrix.
/// 
/// The output is formatted as a vector of ['Bar']s.
/// 
/// - `iter_keymaj` is an iterator that runs over every row-index `p` of the boundary matrix
/// such that dim(`p`) < max { dim(`q`) : `q` is a row-index }.
/// - `dim_fn` returns the dimension of each index
/// - `fil_fn` returns the birth time of each index
/// - `return_cycle_representatives` determines wheter cycle repersentatives are calculated
/// - `return_bounding_chains` determines whether bounding chains are returned
/// 
/// **Note** if both bounding chains and cycle representatives are returned, then this function
/// internally checks that the boundary of the bounding chain is indeed the corresponding
/// cycle representative.
pub fn barcode< 
            I, 
            DimFn, 
            FilFn, 
            KeyBoth,
            EntryBoth,
            Mapping, 
            RingOperator, 
            OrderOperatorRowEntries, 
        >
    ( 
        umatch:                         & Umatch<  Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorRowEntries  >, 
        iter_keymaj:                    I, 
        mut dim_fn:                     DimFn, 
        mut fil_fn:                     FilFn, 
        return_cycle_representatives:   bool,
        return_bounding_chains:         bool,
    ) 
    -> Barcode<
                Mapping::RowIndex, 
                Mapping::EntryMinor 
            > 
    where
        I:                                          Iterator<Item=KeyBoth>, 
        DimFn:                                      FnMut( & KeyBoth)-> isize,
        FilFn:                                      FnMut( & KeyBoth)-> OrderedFloat< f64 >,      
        Mapping:                               ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorRowEntries:        Clone + JudgePartialOrder<  Mapping::EntryMajor  > + JudgePartialOrder< Mapping::EntryMinor>,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq    

        
{
    let mut barcode             =   Vec::new();
    let jordan                  =   JordanBasisMatrix::new( umatch );

    let matching = umatch.matching_ref(); 

    for keymaj in iter_keymaj {
        if matching.contains_keymin( & keymaj ) { continue } // in this case we ignore the key, since it doesn't correspond to a cycle

        let death_column = matching.keymaj_to_keymin( & keymaj );
        let death = death_column.as_ref().map( &mut fil_fn );
        let birth = fil_fn( & keymaj );        
        
        if death == Some( birth ) { continue } // if birth = death, don't include the bar

        let cycle_representative   =   match return_cycle_representatives {
            false   =>  { None },
            true    =>  { Some( jordan.view_minor_descend( keymaj.clone() ).collect_vec() )
                            // if death_column.is_none() {
                            //     // the chain never dies, so look it up the chain as a column of the domain comb
                            //     Some(
                            //         umatch.comb_domain()
                            //             .view_minor_descend( keymaj.clone() )
                            //             .into_iter()
                            //             .collect_vec()
                            //         )                                
                            // } else {
                            //     // the chain in question is a boundary, so look it up as a column of the codomain comb
                            //     let boundary_vec        =   umatch.comb_codomain()
                            //                                         .view_minor_descend( keymaj.clone() )
                            //                                         .into_iter()
                            //                                         .collect_vec();
                            //     let scalar              =   umatch.matching_ref().keymaj_to_snzval(&keymaj);
                            //     let boundary_vec_scaled =  boundary_vec.iter().cloned().scale( scalar, umatch.ring_operator() ).collect_vec();
                            //     Some( boundary_vec_scaled )
                            // }
                                                      
                        }
        };

        // compute bounding chain and verify that it does bound the 
        let bounding_chain  =   match return_bounding_chains {
            false   =>  { None },
            true    =>  {
                            if let Some( keymin ) = matching.keymaj_to_keymin( & keymaj ) {
                                let bounding_chain  =   JordanBasisMatrix::new( umatch ).view_minor_descend( keymin ).collect_vec();                           
                                                        // umatch.comb_domain()
                                                        //     .view_minor_descend( keymin )
                                                        //     .into_iter()
                                                        //     .collect_vec();
                                let boundary_vec    =   vector_matrix_multiply_minor_descend_simplified(
                                                                &bounding_chain, 
                                                                umatch.mapping_ref(), 
                                                                umatch.ring_operator(), 
                                                                umatch.order_operator_major(),
                                                            ).collect_vec();
                                if return_bounding_chains && return_cycle_representatives {
                                    assert_eq!( cycle_representative, Some(boundary_vec) ); // check that the bounding chain bounds the right cycle
                                }
                                Some( bounding_chain )
                            } else {
                                None
                            }
                        }
        };

        let dimension = dim_fn( & keymaj );
        let id_number = barcode.len();
        let birth_column    =   keymaj;
        let bar     =   Bar{ id_number, dimension, birth, birth_column, death, death_column, cycle_representative, bounding_chain };

        barcode.push(bar);
    }
    Barcode{ bars: barcode }
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

    use crate::algebra::vectors::entries::KeyValGet;
    use crate::algebra::matrices::display::print_indexed_major_views;
    use crate::algebra::matrices::types::third_party::IntoCSR;
    use crate::algebra::rings::operator_structs::ring_native::FieldRationalSize;
    use crate::algebra::{matrices::query::MatrixEntry, chains::factored::factor_boundary_matrix, };

    use crate::topology::point_cloud::unit_circle;
    use crate::topology::simplicial::simplices::filtered::SimplexFiltered;
    use crate::topology::simplicial::{from::graph_weighted::{ChainComplexVrFiltered, }, };

    use crate::utilities::distances::rowwise_distances;
    use crate::utilities::iterators::general::minmax;
    use crate::utilities::order::{OrderOperatorAuto, };

    use std::sync::Arc;
    use sprs::CsMatBase;




    #[test]
    fn test_barcode_random_symmetric_matrix() {

    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        let npoints = 20;
        let maxdim = 1;

        let dissimilarity_matrix_data = crate::utilities::random::random_symmetric_matrix_zero_diag(npoints);
        let dissimilarity_matrix_sparse = dissimilarity_matrix_data.clone().into_csr( npoints, npoints );        
        let dissimilarity_matrix = & dissimilarity_matrix_sparse;

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


        for i in 0 .. npoints {
            for j in i .. npoints {
                assert_eq!( dissimilarity_matrix.entry_major_at_minor(i,j), dissimilarity_matrix.entry_major_at_minor(j,i) );
            }
        }


        let ring_operator = crate::algebra::rings::operator_structs::ring_native::FieldRationalSize::new();
        let boundary_matrix_data = ChainComplexVrFiltered::new( dissimilarity_matrix, npoints, dissimilarity_value_max, dissimilarity_value_min, ring_operator );
        let boundary_matrix = Arc::new(boundary_matrix_data);
        let keymaj_vec = boundary_matrix.cliques_in_order(maxdim);
    
        let iter_keymaj = keymaj_vec.iter().cloned();    
            
        println!("starting umatch");
        let factored = factor_boundary_matrix(
                    boundary_matrix, 
                    ring_operator, 
                    OrderOperatorAuto::new(), 
                    iter_keymaj.clone(),
                );
    
        println!("setting up to unpack");  
        let dim_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.dimension() as isize;
        let fil_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.filtration();    
        let barcode = barcode( factored.umatch(), iter_keymaj, dim_fn, fil_fn, true , true);
        for bar in barcode.bars() {  println!("{:#?}", bar)}
        for row in 0 .. npoints {
            println!("view: {:?},", dissimilarity_matrix_data[row].iter().cloned().map(|x| x.into_inner()).collect_vec() );
        }
    }    


    #[test]
    fn test_barcode_random_symmetric_matrix_with_theshold() {

    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        let npoints = 20;
        let maxdim = 1;
        let dissimilarity_value_max = OrderedFloat(1.1);
        let dissimilarity_value_min = OrderedFloat(0.0);


        let dissimilarity_matrix_data = crate::utilities::random::random_symmetric_matrix_zero_diag(npoints).into_csr( npoints, npoints );
        let mut tri = sprs::TriMat::new( (npoints, npoints) );
        for (v, (i,j) ) in dissimilarity_matrix_data.iter() { 
            if v.clone() <= dissimilarity_value_max{
                tri.add_triplet(i, j, v.clone()); 
            }            
        }
        let dissimilarity_matrix_data: CsMatBase<_,_,_,_,_> = tri.to_csr();
        let dissimilarity_matrix = & dissimilarity_matrix_data;

        for i in 0 .. npoints {
            for j in i .. npoints {
                assert_eq!( dissimilarity_matrix.entry_major_at_minor(i,j), dissimilarity_matrix.entry_major_at_minor(j,i) );
            }
        }
     
        let ring_operator = crate::algebra::rings::operator_structs::ring_native::FieldRationalSize::new();
        let boundary_matrix_data = ChainComplexVrFiltered::new( dissimilarity_matrix, npoints, dissimilarity_value_max, dissimilarity_value_min, ring_operator );
        let boundary_matrix = Arc::new(boundary_matrix_data);
        let keymaj_vec = boundary_matrix.cliques_in_order(maxdim);
    
        let iter_keymaj = keymaj_vec.iter().cloned();    
            
        println!("starting umatch");
        let factored = factor_boundary_matrix(
                    boundary_matrix, 
                    ring_operator, 
                    OrderOperatorAuto::new(), 
                    iter_keymaj.clone(),
                );
    
        println!("setting up to unpack");  
        let dim_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.dimension() as isize;
        let fil_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.filtration();    
        let barcode = barcode( factored.umatch(), iter_keymaj, dim_fn, fil_fn, true , true);
        for bar in barcode.bars() {  println!("{:#?}", bar)}
    }    


    #[test]
    fn test_barcode_circle() {

    
        use crate::algebra::matrices::types::third_party::IntoCSR;
        
        let npoints = 50;
        let maxdim = 2;

        let cloud = unit_circle(npoints, Some(0.0 .. 0.1) );
        let dissimilarity_matrix_data = rowwise_distances(cloud).into_csr( npoints, npoints );
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


        for i in 0 .. npoints {
            for j in i .. npoints {
                assert_eq!( dissimilarity_matrix.entry_major_at_minor(i,j), dissimilarity_matrix.entry_major_at_minor(j,i) );
            }
        }
     
        let ring_operator = crate::algebra::rings::operator_structs::ring_native::FieldRationalSize::new();
        let boundary_matrix_data = ChainComplexVrFiltered::new( dissimilarity_matrix, npoints, dissimilarity_value_max, dissimilarity_value_min, ring_operator );
        let boundary_matrix = Arc::new(boundary_matrix_data);
        let keymaj_vec = boundary_matrix.cliques_in_order(maxdim);
    
        let iter_keymaj = keymaj_vec.iter().cloned();    
            
        // println!("starting umatch");
        let factored = factor_boundary_matrix(
                    boundary_matrix, 
                    ring_operator, 
                    OrderOperatorAuto::new(), 
                    iter_keymaj.clone(),
                );
    
        // println!("setting up to unpack");  
        let dim_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.dimension() as isize;
        let fil_fn = |x: &SimplexFiltered<OrderedFloat<f64>> | x.filtration();    
        let barcode = barcode( factored.umatch(), iter_keymaj, dim_fn, fil_fn, true , true);
        // for bar in barcode.bars() {  println!("{:#?}", bar)}
        // println!("number of pairs: {:?}", factored.umatch().matching_ref().num_pairs() );
    }        




}









    // /// The betti numbers of a chain complex.
    // /// 
    // /// The betti numbers of a chain complex with boundary matrix `D` can be computed as 
    // /// follows, c.f. [U-match Factorization](https://arxiv.org/abs/2108.08831) (recall
    // /// that `D` has rows and columns for chains in every dimension): (i) obtain
    // /// a U-match factorization of `D` with matching matrix `M`, (ii) the `k`th betti number
    // /// is the number of indices `i` such that `M[i,:]=0` and `M[:,i]=0`, and index `i`
    // /// corrsponds to a chain of dimension `k`.
    // /// 
    // /// This function computes betti numbers according to this formula, assuming that `Self`
    // /// is the matching matrix of some boundary matrix `D`.  Argument `I` is an iterator
    // /// that runs over all the row (equivalently column) indices of `D`.  If you only need
    // /// homology up through dimension `d`, you only need to include indices for chains
    // /// of dimension `d` and below.  Argument `dim_fn`
    // /// is a function that returns the dimension of the chain associated with each index.
    // /// 
    // /// *Remark* under the hood this method is identical to [unmatched_histo].  It is 
    // /// included as a separate function primarily for the purpose of readability.
    // pub fn betti_numbers< I, F >( &self, iter_keyboth: I, dim_fn: F ) 
    //         -> 
    //         Vec< usize > 
    //     where
    //         I: Iterator<Item=KeyBoth>, 
    //         F: FnMut(KeyBoth)->usize  
    //     {
    //     return self.matching_ref().unmatched_histo(iter_keyboth, dim_fn)
    // }   