//! The Jordan canonical form (up to permutation) of a boundary matrix with rows/columns indexed by chains of all dimensions

use crate::algebra::matrices::operations::umatch::row_major::comb::{CombDomainViewMinorDescend, CombCodomainViewMinorDescend};
use crate::utilities::iterators::general::IterTwoType;
use crate::algebra::vectors::entries::{KeyValSet, KeyValNew};
use crate::algebra::matrices::{operations::umatch::row_major::Umatch, query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients}};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};
use crate::utilities::order::JudgePartialOrder;
use crate::algebra::vectors::operations::{VectorOperations, Scale};

use std::fmt::Debug;
use std::hash::Hash;




//  ---------------------------------------------------------------------------------------
//  JORDAN BASIS VECTOR
//  ---------------------------------------------------------------------------------------


/// An element of a Jordan basis
/// 
/// This is a wrapper struct; it's used as a convenient way to reduce the
/// complexity of the type signature of a composite of two iterator types.
pub struct JordanBasisMatrixVector< 
                    'a, 
                    Mapping, 
                    RingOperator, 
                    OrderOperatorRowEntries,                     
                    KeyBoth,
                    EntryBoth,
                >
    where
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
    column:      IterTwoType<     
                        CombDomainViewMinorDescend<'a, Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorRowEntries>,
                        Scale< 
                                CombCodomainViewMinorDescend< Mapping, RingOperator, OrderOperatorRowEntries >,
                                KeyBoth, 
                                RingOperator, 
                                Mapping::Coefficient,
                            > 
                        
                    >,
}

impl < 
            'a, 
            Mapping, 
            RingOperator, 
            OrderOperatorRowEntries,                     
            KeyBoth,
            EntryBoth,
        >

    Iterator for

    JordanBasisMatrixVector< 
                    'a, 
                    Mapping, 
                    RingOperator, 
                    OrderOperatorRowEntries,                     
                    KeyBoth,
                    EntryBoth,
                >
    where
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
    type Item = EntryBoth;

    fn next(&mut self) -> Option<Self::Item> { self.column.next() }
}




//  ================================================================================
//  JORDAN BASIS
//  ================================================================================



//  ---------------------------------------------------------------------------------------
//  MATRIX
//  ---------------------------------------------------------------------------------------



/// The Jordan basis of a boundary matrix decomposed with U-match factorization.
/// 
/// This struct represents a square invertible upper-triangular matrix whose columns form
/// a Jordan basis for a boundary matrix `D`.  Concretely, the construction proceeds as follows
/// 
/// - Obtain a [U-match factorization](https://arxiv.org/abs/2108.08831) of `D` this means a matrix equation `RM = DC`
///   where `R` and `C` are square upper-unitriangular matrices and `M` is a generalized
///   matching matrix.  This factorization is related to the `R=DV` and `RU=D` factorizations.
/// - Let `J = C`.  Then, for each nonzero entry `M[i,j]`, replace column
///   `J[:,i]` with column `D * J[:,j]`.
/// 
/// The resulting matrix, `J`, is upper triangular, and its columns form a Jordan bases for
/// `D`.
pub struct JordanBasisMatrix< 
                    'a, 
                    Mapping, 
                    RingOperator, 
                    OrderOperatorRowEntries, 
                    KeyBoth,
                    EntryBoth,                    
                >
    where
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
    umatch: &'a Umatch<  Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorRowEntries  >,                 
}


impl <  'a,
        KeyBoth,
        EntryBoth,
        Mapping, 
        RingOperator, 
        OrderOperatorRowEntries, 
    >

        JordanBasisMatrix< 
                'a,
                Mapping, 
                RingOperator, 
                OrderOperatorRowEntries, 
                KeyBoth,
                EntryBoth,                
            >
    where
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

    /// Create a new Jordan basis.
    pub fn new( 
                    umatch: &'a Umatch<  Mapping, RingOperator, OrderOperatorRowEntries, OrderOperatorRowEntries  >,                 
            ) -> Self {
        JordanBasisMatrix { umatch }
    }

    // THIS METHOD IS DEPRECATED IN FAVOR OF THE MATRIX ORACLE METHOD view_minor_descend.  THEREFORE OK TO DELETE.
    // /// Extract a column of the Jordan basis associated with the U-match factorization.
    // pub fn keymaj_to_column_vec( &self, keymaj:  &KeyBoth ) -> Vec< EntryBoth > {
    //     match self.umatch.matching_ref().keymaj_to_keymin( keymaj ) {
    //             None => {
    //                 // the chain never becomes a boundary, so look it up as a column of the domain comb
    //                 self.umatch.comb_domain()
    //                     .view_minor_descend( keymaj.clone() )
    //                     .into_iter()
    //                     .collect_vec()
    //             }
    //             _ => {
    //                 // the chain eventually becomes a boundary, so look it up as a column of the codomain comb
    //                 let boundary_vec        =   self.umatch.comb_codomain()
    //                                                     .view_minor_descend( keymaj.clone() )
    //                                                     .into_iter()
    //                                                     .collect_vec();
    //                 // we scale the column to ensure that it's actually in the image of another column
    //                 let scalar              =   self.umatch.matching_ref().keymaj_to_snzval(&keymaj);
    //                 let boundary_vec_scaled =  boundary_vec.iter().cloned().scale( scalar, self.umatch.ring_operator() ).collect_vec();
    //                 boundary_vec_scaled
    //             }
    //     }
    // }
}    


//  ---------------------------------------------------------------------------------------
//  ORACLE IMPLEMNTATION: INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------------------

impl <  'a,
        KeyBoth,
        EntryBoth,
        Mapping, 
        RingOperator, 
        OrderOperatorRowEntries, 
    >

    IndicesAndCoefficients for

    JordanBasisMatrix< 
            'a,
            Mapping, 
            RingOperator, 
            OrderOperatorRowEntries, 
            KeyBoth,
            EntryBoth,                
        >
    where
        Mapping:                               IndicesAndCoefficients< 
                                                            ColIndex=KeyBoth, 
                                                            RowIndex=KeyBoth, 
                                                            EntryMajor=EntryBoth, 
                                                            EntryMinor=EntryBoth 
                                                        > + 
                                                    ViewRowAscend + 
                                                    ViewColDescend,     
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
    type EntryMajor = EntryBoth;
    type EntryMinor = EntryBoth;    
    type RowIndex = KeyBoth;
    type ColIndex = KeyBoth;    
    type Coefficient = Mapping::Coefficient;    
}


//  ---------------------------------------------------------------------------------------
//  ORACLE IMPLEMNTATION: VIEW MINOR DESCEND
//  ---------------------------------------------------------------------------------------


impl <  'a,
        KeyBoth,
        EntryBoth,
        Mapping, 
        RingOperator, 
        OrderOperatorRowEntries, 
    >

    ViewColDescend for

    JordanBasisMatrix< 
            'a,
            Mapping, 
            RingOperator, 
            OrderOperatorRowEntries, 
            KeyBoth,
            EntryBoth,                
        >
    where
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
    type ViewMinorDescend = JordanBasisMatrixVector< 'a, Mapping, RingOperator, OrderOperatorRowEntries, KeyBoth, EntryBoth, >;

    type ViewMinorDescendIntoIter = Self::ViewMinorDescend;

    fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend {
        
        let keymaj = index;
        let column = 
            match self.umatch.matching_ref().keymaj_to_keymin( &keymaj ) {
                None => {
                    // the chain never becomes a boundary, so look it up as a column of the domain comb
                    IterTwoType::Iter1(
                        self.umatch.comb_domain()
                        .view_minor_descend( keymaj.clone() )
                    )

                }
                _ => {             
                    // the chain eventually becomes a boundary, so look it up as a column of the codomain comb
                    let boundary_vec        =   self.umatch.comb_codomain()
                                                        .view_minor_descend( keymaj.clone() );
                    // we scale the column to ensure that it's actually in the image of another column
                    let scalar              =   self.umatch.matching_ref().keymaj_to_snzval(&keymaj);
                    let boundary_vec_scaled =  boundary_vec.into_iter().scale( scalar, self.umatch.ring_operator() );
                    IterTwoType::Iter2( boundary_vec_scaled )
                }
            };
        JordanBasisMatrixVector{ column }
    }
}





mod tests {
    use std::sync::Arc;

    use itertools::Itertools;
    use ordered_float::OrderedFloat;
    use sprs::CsMatBase;

    use crate::algebra::matrices::display::print_indexed_major_views;
    use crate::algebra::matrices::types::third_party::IntoCSR;
    use crate::algebra::rings::operator_structs::ring_native::FieldRationalSize;
    use crate::algebra::{matrices::query::MatrixEntry, chains::factored::factor_boundary_matrix, };
    use crate::topology::point_cloud::unit_circle;
    use crate::topology::simplicial::simplices::filtered::SimplexFiltered;
    use crate::utilities::distances::rowwise_distances;
    use crate::utilities::order::{OrderOperatorAuto, };

    use crate::topology::simplicial::{from::graph_weighted::{ChainComplexVrFiltered,}, };

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_ph_jordan_basis_random_symmetric_matrix() {

        use crate::algebra::vectors::entries::KeyValGet;
        use crate::algebra::matrices::types::third_party::IntoCSR;
        use crate::utilities::iterators::general::minmax;
        
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
        // let boundary_matrix = ChainComplexVrFilteredArc::new( Arc::new(boundary_matrix_data) );
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

        let matching    =   factored.umatch().matching_ref();
        let jordan      =   factored.jordan_basis_matrix_packet();
        let boundary    =   factored.umatch().mapping_ref_packet();

        // check that jordan vectors map to each other
        for keymaj in factored.row_indices() {


            if let Some( keymin ) = matching.keymaj_to_keymin( &keymaj ) {
                let preim   =   factored.jordan_basis_vector( keymin );
                let image = factored.jordan_basis_vector( keymaj );
                let prod = preim.multiply_matrix_packet_minor_descend( boundary.clone() );
                assert!( prod.eq( image.collect_vec() ) );
            } else if let Some( paired_row ) = matching.keymin_to_keymaj( & keymaj ) {
                let preim   =   factored.jordan_basis_vector( keymaj );
                let image = factored.jordan_basis_vector( paired_row );
                let prod = preim.multiply_matrix_packet_minor_descend( boundary.clone() );
                assert!( prod.eq( image.collect_vec() ) );                
            } else {
                let mut prod 
                    =   factored
                            .jordan_basis_vector(keymaj)
                            .multiply_matrix_packet_minor_descend( boundary.clone() );
                assert!( prod.next() == None );
            }
        }


        for keymaj in factored.row_indices() {
            let a = factored.jordan_basis_vector(keymaj.clone());
            let b = jordan.matrix.view_minor_descend( keymaj.clone() );
            assert!( a.collect_vec() == b.collect_vec() );
        }
    }    
}