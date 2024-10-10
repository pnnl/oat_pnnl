//! U-match factorization of boundary matrices, with functionality for (co)homology, cycles, boundaries, etc.



use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
use crate::algebra::vectors::entries::{KeyValSet, KeyValNew};
use crate::algebra::matrices::query::{ViewRowAscend, ViewColDescend, IndicesAndCoefficients};
use crate::algebra::matrices::operations::umatch::row_major::{Umatch, ParetoShortCircuit};
use crate::algebra::rings::operator_traits::{Semiring, Ring, DivisionRing};
use crate::utilities::order::{JudgePartialOrder, OrderOperatorByKeyCutsom};

use std::cmp::Ordering;
use std::collections::{HashMap};
use std::fmt::Debug;
use std::hash::Hash;

use super::jordan::{JordanBasisMatrix, JordanBasisMatrixVector};





/// Factor a boundary matrix
/// 
/// Arguments
/// - `boundary_matrix`: the boundary matrix
/// - `ring_operator`: the ring operator for the matrix coefficients
/// - `order_operator`: represents the linear order on row/column indices
/// - `row_indices`: an iterator that runs over row indices in
///    (first) ascending order of dimension, and (second) descending order,
///    as determined by `order_operator`.
/// 
/// Returns a `FactoredBoundaryMatrix`, which provides access to bases for
/// homology, cycles, boundaries, etc.
pub fn factor_boundary_matrix< 
                Mapping, 
                RingOperator, 
                OrderOperatorKeyBoth, 
                KeyBoth,
                EntryBoth,                    
                I,                     
            > 
        (
            boundary_matrix:    Mapping,
            ring_operator:      RingOperator,
            order_operator:   OrderOperatorKeyBoth,
            row_indices:        I,
        ) 
        ->  FactoredBoundaryMatrix< 
                    Mapping, 
                    RingOperator, 
                    OrderOperatorByKeyCutsom< KeyBoth, Mapping::Coefficient, EntryBoth, OrderOperatorKeyBoth >, 
                    KeyBoth,
                    EntryBoth,                    
                    I,                     
                > 
    where
        I:                                          Clone + IntoIterator<Item=KeyBoth>,
        Mapping:                                    ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:                   IntoIterator + ParetoShortCircuit< EntryBoth >,
        // Mapping::ViewMajorAscendIntoIter:      Clone,      
        Mapping::EntryMajor:                        Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:                  IntoIterator,        
        Mapping::EntryMinor:                            Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorKeyBoth:                       Clone + JudgePartialOrder<  Mapping::RowIndex  > + JudgePartialOrder< Mapping::ColIndex >,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq,

{
    let umatch = Umatch::factor_with_clearing(
            // the matrix we want to factor
            boundary_matrix, 
            // the relevant rows of the matrix, sorted first by dimension and then in reverse lexicographic order
            // computing homology of dimension d, we only need to iterate over simplices of dimension d and below.
            row_indices.clone().into_iter(), 
            // the operator for the coefficient ring
            ring_operator, 
            // a struct that specifies the order on column indices
            order_operator.clone(),
            // a struct that specifies the order on row indices  (should be same as order on column indices)               
            order_operator,
        );        
    FactoredBoundaryMatrix { umatch, row_indices }
}




//  ==============================================================================
//  FACTORED BOUNDARY MATRIX STRUCT
//  ==============================================================================


/// A factored boundary matrix
/// 
/// Holds compressed represenation of (i) the Jordan canonical form of the matrix, 
/// and its corresponding Jordan basis, (ii) a U-match factorization.
/// 
/// These can be used to compute bases for homology, cycle space, boundary spaces,
/// and more.  See the methods on this struct for details.
#[derive(Debug,Clone)]
pub struct FactoredBoundaryMatrix< 
                    Mapping, 
                    RingOperator, 
                    OrderOperatorBoth, 
                    KeyBoth,
                    EntryBoth,                    
                    I,                     
                > 
    where
        I:                                          Clone + IntoIterator<Item=KeyBoth>, 
        Mapping:                               ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorBoth:        Clone + JudgePartialOrder<  Mapping::EntryMajor  > + JudgePartialOrder< Mapping::EntryMinor>,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq                  
{
    umatch:                         Umatch<  Mapping, RingOperator, OrderOperatorBoth, OrderOperatorBoth  >, 
    row_indices:                    I, 
    // safe_homology_dims:             HashSet< isize >, // indicates the homology dimensions where it is safe to compute betti numbers, cycle representatives, etc.
}


impl < 
        Mapping, 
        RingOperator, 
        OrderOperatorBoth, 
        KeyBoth,
        EntryBoth,                    
        I, 
    > 

    FactoredBoundaryMatrix< 
                    Mapping, 
                    RingOperator, 
                    OrderOperatorBoth, 
                    KeyBoth,
                    EntryBoth,                    
                    I, 
                > 
    where
        I:                                          Clone + IntoIterator<Item=KeyBoth>, 
        Mapping:                               ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorBoth:                        Clone + JudgePartialOrder<  Mapping::EntryMajor  > + JudgePartialOrder< Mapping::EntryMinor>,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq                  
{    






    // Jordan matrix
    // ---------------

    /// Invertible upper triangular matrix whose columns form a Jordan basis of the boundary matrix.
    /// 
    /// The matrix is indexed by the same keys as the boundary matrix.
    /// 
    /// For details, see the documenation for [JordanBasisMatrix](super::jordan_basis::JordanBasisMatrix)
    pub fn jordan_basis_matrix( &self ) -> JordanBasisMatrix< '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, >  
    { JordanBasisMatrix::new( & self.umatch ) }

    /// Invertible upper triangular matrix whose columns form a Jordan basis of the boundary matrix, wrapping in a convenient [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket)
    /// 
    /// The matrix is indexed by the same keys as the boundary matrix.
    /// 
    /// For details, see the documenation for [JordanBasisMatrix](super::jordan_basis::JordanBasisMatrix)
    pub fn jordan_basis_matrix_packet( &self ) 
        -> 
        MatrixAlgebraPacket< 
                JordanBasisMatrix< '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, >,
                RingOperator,
                OrderOperatorBoth,
                OrderOperatorBoth,
            >
    { MatrixAlgebraPacket{ matrix: JordanBasisMatrix::new( & self.umatch ), ring: self.umatch.ring_operator(), row_entry_order: self.umatch.order_operator_major(), col_entry_order: self.umatch.order_operator_minor() } }

    /// A Jordan basis vector
    /// 
    /// This function is equivalent to  `self.jordan_basis_matrix().view_minor_descend(column_index)` 
    pub fn jordan_basis_vector( &self, column_index: KeyBoth ) -> JordanBasisMatrixVector< '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, >  
    {  self.jordan_basis_matrix().view_minor_descend(column_index)  }

    /// A reference to the internally stored U-match factorization
    pub fn umatch(&self) -> & Umatch<  Mapping, RingOperator, OrderOperatorBoth, OrderOperatorBoth  >
    { & self.umatch }

    // Index iterators
    // ---------------

    /// Runs over the row indices of the boundary matrix in reverse order
    pub fn row_indices( &self ) -> I { self.row_indices.clone() }    

    /// Harmonic indices of the Jordan matrix
    /// 
    /// Includes only indices which are visited by `self.row_indices()`
    pub fn indices_harmonic( &self ) -> 
        JordanIndices    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanIndices{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::harmonic() } }
    
    /// Cycle indices of the Jordan matrix
    /// 
    /// Includes only indices which are visited by `self.row_indices()`
    pub fn indices_cycle( &self ) -> 
        JordanIndices    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanIndices{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::cycle() } }
    
    /// Cocycle indices of the Jordan matrix
    /// 
    /// Includes only indices which are visited by `self.row_indices()`
    pub fn indices_cocycle( &self ) -> 
        JordanIndices    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanIndices{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::cocycle() } }
    
    /// Boundary indices of the Jordan matrix
    /// 
    /// Includes only indices which are visited by `self.row_indices()`
    pub fn indices_boundary( &self ) -> 
        JordanIndices    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanIndices{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::boundary() } }
    
    /// Coboundary indices of the Jordan matrix
    /// 
    /// Includes only indices which are visited by `self.row_indices()`
    pub fn indices_coboundary( &self ) -> 
        JordanIndices    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanIndices{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::coboundary() } }

    /// Non-harmonic indices of the Jordan matrix
    /// 
    /// Includes only indices which are visited by `self.row_indices()`
    pub fn indices_non_harmonic( &self ) -> 
        JordanIndices    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanIndices{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::non_harmonic() } }   


    // Chain iterators
    // ---------------

    /// Harmonic elements in the Jordan basis
    /// 
    /// Iterates over `self.indices_harmonic()`, converting each index to the corresponding column of the Jordan basis
    /// matrix.
    /// 
    /// The set of dimension-`d` chains returned by this iterator is a basis of cycle representatives for homology in dimension `d`,
    /// provided that `self.row_indices()` runs over all indices of the boundary matrix of dimension-`d` and dimension-`d-1`.
    pub fn basis_harmonic( &self ) -> 
        JordanColumns    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanColumns{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::harmonic() } }
    
    /// Cycles in the Jordan basis
    /// 
    /// Iterates over `self.indices_cycle()`, converting each index to the corresponding column of the Jordan basis
    /// matrix.
    /// 
    /// The set of dimension-`d` chains returned by this iterator is a basis for the space of `d`-dimensional cycles,
    /// provided that `self.row_indices()` runs over all indices of the boundary matrix of dimension-`d`
    pub fn basis_cycle( &self ) -> 
        JordanColumns    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanColumns{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::cycle() } }
    
    /// Jordan basis vectors indexed by cocycle indices (**not** cocycles; see details)
    /// 
    /// Iterates over `self.indices_cocycle()`, converting each index to the corresponding column of the Jordan basis matrix.
    /// 
    /// This is **not** a basis for the cocycle space of the chain complex; for that we need to iterate over elements of
    /// a Jordan *cobasis* (that is, we would need to iterate over the rows of the *inverse* of the Jordan basis matrix, not 
    /// over columns of the Jordan basis matrix).
    pub fn basis_cocycle( &self ) -> 
        JordanColumns    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanColumns{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::cocycle() } }
    
    /// Boundaries in the Jordan basis
    /// 
    /// Iterates over `self.indices_boundary()`, converting each index to the corresponding column of the Jordan basis
    /// matrix.
    /// 
    /// The set of dimension-`d` chains returned by this iterator is a basis for the space of `d`-dimensional boundaries,
    /// provided that `self.row_indices()` runs over all indices of the boundary matrix of dimension-`d+1`
    pub fn basis_boundary( &self ) -> 
        JordanColumns    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanColumns{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::boundary() } }
    
    /// Jordan basis vectors indexed by coboundary indices (**not** coboundaries; see details)
    /// 
    /// Iterates over `self.indices_coboundary()`, converting each index to the corresponding column of the Jordan basis matrix.
    /// 
    /// This is **not** a basis for the coboundary space of the chain complex; for that we need to iterate over elements of
    /// a Jordan *cobasis* (that is, we would need to iterate over the rows of the *inverse* of the Jordan basis matrix, not 
    /// over columns of the Jordan basis matrix).
    pub fn basis_coboundary( &self ) -> 
        JordanColumns    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanColumns{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::coboundary() } }

    /// Non-harmonic indices of the Jordan matrix
    /// 
    /// Iterates over `self.indices_non_harmonic()`, converting each index to the corresponding column of the Jordan basis matrix.
    pub fn basis_non_harmonic( &self ) -> 
        JordanColumns    < '_, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, I, >  
    { JordanColumns{ umatch: &self.umatch, iter_keymaj: self.row_indices().into_iter(), criteria: JordanCriteria::non_harmonic() } }  


    // Betti numbers
    // -------------
    /// The betti numbers of the factored complex (but see details)
    /// 
    /// These are correct for every dimension `d` such that `self.row_indices()` contains every
    /// index of dimension `d-1` and `d`.
    /// 
    /// If `d` satisfies this condition and the returned hashmap
    /// has no key `d`, then the value for `d` is zero.
    pub fn betti_numbers< F: FnMut(KeyBoth)-> isize >( &self, mut dim_fn: F ) -> HashMap< isize, isize > {
        let mut curve = HashMap::new();
        for key in self.indices_harmonic() {
            let dim = dim_fn(key);
            * curve.entry( dim ).or_insert(0) += 1;
        }
        curve
    }

    // Cycle numbers
    // -------------
    /// The dimensions of the cycle spaces
    /// 
    /// These are correct for every dimension `d` such that `self.row_indices()` contains every
    /// index of dimension `d-1` and `d`.
    /// 
    /// If `d` satisfies this condition and the returned hashmap
    /// has no key `d`, then the value for `d` is zero.
    pub fn cycle_numbers< F: FnMut(KeyBoth)-> isize >( &self, mut dim_fn: F ) -> HashMap< isize, isize > {
        let mut curve = HashMap::new();
        for key in self.indices_cycle() {
            let dim = dim_fn(key);
            * curve.entry( dim ).or_insert(0) += 1;
        }
        curve
    }    

    // Boundary numbers
    // ----------------
    /// The dimensions of the boundary spaces
    /// 
    /// These are correct for every dimension `d` such that `self.row_indices()` contains every
    /// index of dimension `d-1` and `d`.  
    /// 
    /// If `d` satisfies this condition and the returned hashmap
    /// has no key `d`, then the value for `d` is zero.
    pub fn boundary_numbers< F: FnMut(KeyBoth)-> isize >( &self, mut dim_fn: F ) -> HashMap< isize, isize > {
        let mut curve = HashMap::new();
        for key in self.indices_boundary() {
            let dim = dim_fn(key);
            * curve.entry( dim ).or_insert(0) += 1;
        }
        curve
    }      







    /// Columns of cycle minimization constraint matrix in Escolar Hiraoka 2015 (with extended functionality).
    /// 
    /// Let `birth_column` be the index of a column in the boundary matrix.  This function returns
    /// every column index `i` such that 
    /// - `i` has the same dimension in homology as `birth_column`
    /// - `i` strictly precedes `birth_column` in the linear order on column indices
    /// and element `i` of the Jordan basis is
    /// - a cycle
    /// - bounded no later than `birth_column` (if `birth_column` indexes an essential class then this )
    /// 
    /// # Arguments
    /// 
    /// - `birth_column`: the column we wish to compare
    /// - `dim_fn`: returns the homology dimension of each column index of the boundary matrix
    /// 
    /// # Implications for computation
    /// 
    /// Adding any linear combination of these columns to a the jordan basis vector of `birth_column` will yield
    /// a cycle representative `v` for a persistent homology class with the same birth and death time as the cycle `u`
    /// of `birth_column` (namely, `u` is the column of the Jordan basis indexed by `birth_simplex`).  Moreover,
    /// swapping `u` for `v` will not create any new linear dependencies in the basis; that is, it will produce a
    /// basis for cycle space which includes a basis for persistent homology. We can perform this swap on every column
    /// in the jordan basis, simultaneously, and the result will be a valid basis for persistent homology.

    pub fn indices_escolar_hiraoka< DimFn >( &self, birth_column: & KeyBoth, mut dim_fn: DimFn, ) 
        -> Vec< KeyBoth > 
        where
            DimFn:  FnMut( & KeyBoth)-> isize,       
    {        
        // define some parameters
        let order_operator = self.umatch().order_operator_major();
        let matching_ref      =   self.umatch.matching_ref();
        
        // a convenience function; we want to use our order operator, but it technically only operates
        // on matrix entries, not on indices; so we embed indices in entries
        let zero = RingOperator::zero();
        let to_entry = |i: KeyBoth| EntryBoth::new( i, zero.clone() );        
        
        // calculate some baselines for comparison
        let objective_dimension     =   dim_fn( birth_column );
        let death_column_opt         =   matching_ref                              //  death filtration value
                                            .keymaj_to_keymin( birth_column );
        
        // now we can enumerate the desired column indices
        let mut columns = Vec::new();
        for keymaj in self.row_indices() {

            // must have equal dimension
            if dim_fn( & keymaj ) != objective_dimension { continue }
            // must be born strictly before birth_column
            if order_operator.judge_ge( &to_entry(keymaj.clone()), &to_entry(birth_column.clone()) ) { continue }
            // must be a cycle
            if matching_ref.contains_keymin( & keymaj ) { continue }
            // cannot die strictly later than birth_column
            let comp_death_opt          =   matching_ref                                      
                                            .keymaj_to_keymin( birth_column );
            if let Some( death_column ) = death_column_opt.as_ref() {
                if let Some( comp_death ) = comp_death_opt {
                    if order_operator.judge_gt( &to_entry(comp_death.clone()), &to_entry(death_column.clone()) ) {
                        continue
                    }
                } else {
                    continue
                }
            }
            
            columns.push(keymaj);
        }        
        columns
    }


    /// Columns of cycle minimization constraint matrix in Escolar Hiraoka 2015 (**relaxed** and with extended functionality).
    /// 
    /// Similar to [FactoredBoundaryMatrix::indices_escolar_hiraoka], but we relax the constraint on which indices to include,
    /// by allowing some indices with equal filtration value.
    /// 
    /// Let `birth_column` be the index of a column in the boundary matrix.  This function returns
    /// every column index `i` such that 
    /// - column `i` of the jordan basis is a cycle with the same dimension in homology as `birth_column`
    /// - birth time of the PH class represented by `i` ≤ birth time of the PH class represented by `birth_column`
    /// - death time of the PH class represented by `i` ≤ death time of the PH class represented by `birth_column`
    /// 
    ///
    /// # Implications for computation
    /// 
    /// Adding any linear combination of these columns to a the jordan basis vector of `birth_column` will yield
    /// a cycle representative `v` for a persistent homology class with the same birth and death time as the cycle `u`
    /// of `birth_column` (namely, `u` is the column of the Jordan basis indexed by `birth_simplex`).  Moreover,
    /// swapping `u` for `v` will not create any new linear dependencies in the basis; that is, it will produce a
    /// basis for cycle space which includes a basis for persistent homology. **HOWEVER**, this gaurantee only holds
    /// for a single substitution; were we to optimize two different columns and swap them both into the Jordan basis,
    /// we might create linear dependency.
    /// 
    /// # Arguments
    /// 
    /// - `birth_column`: the column we wish to compare
    /// - `dim_fn`: returns the homology dimension of each column index of the boundary matrix
    /// - `filtration_order`: a closure operator that compares two indices on the basis of their **filtration alone**
    pub fn indices_escolar_hiraoka_relaxed< DimFn, FiltrationOrder >( &self, birth_column: & KeyBoth, mut dim_fn: DimFn, mut filtration_order: FiltrationOrder ) 
        -> Vec< KeyBoth > 
        where
            DimFn:  FnMut( & KeyBoth)-> isize,
            FiltrationOrder:  FnMut( & KeyBoth, & KeyBoth )-> Ordering,
    {        
        // define some parameters
        let matching_ref      =   self.umatch.matching_ref();      
        
        // calculate some baselines for comparison
        let objective_dimension     =   dim_fn( birth_column );
        let death_column_opt         =   matching_ref                              //  death filtration value
                                            .keymaj_to_keymin( birth_column );
        
        // now we can enumerate the desired column indices
        let mut columns = Vec::new();
        for keymaj in self.row_indices() {

            if keymaj == *birth_column { continue }

            // must have equal dimension
            if dim_fn( & keymaj ) != objective_dimension { continue }
            // must be born strictly before birth_column
            if filtration_order( &keymaj, birth_column) == Ordering::Greater { continue }
            // must be a cycle
            if matching_ref.contains_keymin( & keymaj ) { continue }
            // cannot die strictly later than birth_column
            let comp_death_opt          =   matching_ref                                      
                                            .keymaj_to_keymin( birth_column );
            if let Some( death_column ) = death_column_opt.as_ref() {
                if let Some( comp_death ) = comp_death_opt {
                    if filtration_order( &comp_death, death_column) == Ordering::Greater { continue }
                } else {
                    continue
                }
            }
            
            columns.push(keymaj);
        }        
        columns
    }    

}


//  ----------------------------------------------------------------------------------------
//  ITERATOR OF INDICES
//  ----------------------------------------------------------------------------------------


/// A flag used to indicate which indices we should keep.
#[derive(Debug,Clone)]
pub struct JordanCriteria{ boundary: bool, bounding: bool, harmonic: bool}

impl JordanCriteria{
    /// Harmonic indices
    fn harmonic() -> Self { JordanCriteria { boundary: false, bounding: false, harmonic: true } }    
    /// Cycle indices
    fn cycle() -> Self { JordanCriteria { boundary: true, bounding: false, harmonic: true } }
    /// Cocycle indices
    fn cocycle() -> Self { JordanCriteria { boundary: false, bounding: true, harmonic: true } }    
    /// Boundary indices
    fn boundary() -> Self { JordanCriteria { boundary: true, bounding: false, harmonic: false } }
    /// Coboundary indices
    fn coboundary() -> Self { JordanCriteria { boundary: false, bounding: true, harmonic: false } }
    /// Non-harmonic indices
    fn non_harmonic() -> Self { JordanCriteria { boundary: true, bounding: true, harmonic: false } }
}

/// Indices of (non)zero rows and columns
/// 
/// This struct is an iterator that runs over a subset of the rows/column indices of the matching
/// matrix `M` of a U-match factorization `RM=DC`, where `D` is a boundary matrix.  This matrix 
/// is equal (up to nonzero scaling and permutation) to the Jordan canonical form of `D`, hence
/// the name.
/// 
/// Depending on initialization, the iterator will run over the (non)zero rows and columns of `M`.
/// It is constructed by methods on the [`FactoredBoundaryMatrix`] struct.
#[derive(Debug,Clone)]
pub struct JordanIndices< 
                    'a,
                    Mapping, 
                    RingOperator, 
                    OrderOperatorBoth, 
                    KeyBoth,
                    EntryBoth,                    
                    I, 
                > 
    where
        I:                                          Clone + IntoIterator<Item=KeyBoth>, 
        Mapping:                               ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorBoth:        Clone + JudgePartialOrder<  Mapping::EntryMajor  > + JudgePartialOrder< Mapping::EntryMinor>,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq  
{
    umatch:         &'a Umatch<  Mapping, RingOperator, OrderOperatorBoth, OrderOperatorBoth  >,
    iter_keymaj:    I::IntoIter,
    criteria:       JordanCriteria,
}        

impl < 
        'a,
        Mapping, 
        RingOperator, 
        OrderOperatorBoth, 
        KeyBoth,
        EntryBoth,                    
        I, 
    > 

    Iterator for 

    JordanIndices< 
                    'a,
                    Mapping, 
                    RingOperator, 
                    OrderOperatorBoth, 
                    KeyBoth,
                    EntryBoth,                    
                    I, 
                > 
    where
        I:                                          Clone + IntoIterator<Item=KeyBoth>, 
        Mapping:                               ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorBoth:        Clone + JudgePartialOrder<  Mapping::EntryMajor  > + JudgePartialOrder< Mapping::EntryMinor>,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq    
{
    type Item = KeyBoth;

    fn next(&mut self) -> Option<Self::Item> {
        let matching            =   self.umatch.matching_ref();
        let criterion           =   &self.criteria;
        self.iter_keymaj
            .find(  |x|  
                        (   criterion.bounding  &&  matching.contains_keymin(x) )
                    ||
                        (   criterion.boundary  &&  matching.contains_keymaj(x) )
                    ||
                        (   criterion.harmonic  &&  matching.lacks_key(x)       )            
                )                
    }
}


//  ----------------------------------------------------------------------------------------
//  ITERATOR OF VECTORS
//  ----------------------------------------------------------------------------------------




/// Jordan basis vectors in the (co)kernel of the boundary matrix
/// 
/// This struct is an iterator that runs over a subset of the rows/column indices of the matching
/// matrix `M` of a U-match factorization `RM=DC`, where `D` is a boundary matrix.  This matrix 
/// is equal (up to nonzero scaling and permutation) to the Jordan canonical form of `D`, hence
/// the name.
/// 
/// Depending on initialization, the iterator will run over the (non)zero rows and columns of `M`.
/// It is constructed by methods on the [`FactoredBoundaryMatrix`] struct.
#[derive(Debug,Clone)]
pub struct JordanColumns< 
                    'a,
                    Mapping, 
                    RingOperator, 
                    OrderOperatorBoth, 
                    KeyBoth,
                    EntryBoth,                    
                    I, 
                > 
    where
        I:                                          Clone + IntoIterator<Item=KeyBoth>, 
        Mapping:                               ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorBoth:        Clone + JudgePartialOrder<  Mapping::EntryMajor  > + JudgePartialOrder< Mapping::EntryMinor>,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq  
{
    umatch:         &'a Umatch<  Mapping, RingOperator, OrderOperatorBoth, OrderOperatorBoth  >,
    iter_keymaj:    I::IntoIter,
    criteria:       JordanCriteria,
}        

impl < 
        'a,
        Mapping, 
        RingOperator, 
        OrderOperatorBoth, 
        KeyBoth,
        EntryBoth,                    
        I, 
    > 

    Iterator for 

    JordanColumns< 
                    'a,
                    Mapping, 
                    RingOperator, 
                    OrderOperatorBoth, 
                    KeyBoth,
                    EntryBoth,                    
                    I, 
                > 
    where
        I:                                          Clone + IntoIterator<Item=KeyBoth>, 
        Mapping:                               ViewRowAscend<EntryMajor = EntryBoth> + 
                                                    ViewColDescend<EntryMinor = EntryBoth> + 
                                                    IndicesAndCoefficients< ColIndex=KeyBoth, RowIndex=KeyBoth >,     
        KeyBoth:                                    Clone + Debug + Hash + std::cmp::Eq, // required for the hashing performed by the generalized matching array   // !!!! try deleting debug
        Mapping::Coefficient:                       Clone + Debug,
        Mapping::ViewMajorAscend:              IntoIterator,        
        Mapping::EntryMajor:         Clone + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >, 
        Mapping::ViewMinorDescend:             IntoIterator,        
        Mapping::EntryMinor:        Clone + Debug + KeyValSet< Mapping::RowIndex, Mapping::Coefficient > + KeyValNew< Mapping::RowIndex, Mapping::Coefficient >,  // !!!! try to delete debug
        OrderOperatorBoth:        Clone + JudgePartialOrder<  Mapping::EntryMajor  > + JudgePartialOrder< Mapping::EntryMinor>,
        RingOperator:                               Clone + Semiring< Mapping::Coefficient > + Ring< Mapping::Coefficient > + DivisionRing< Mapping::Coefficient >,        
        EntryBoth:                                  std::cmp::PartialEq    
{
    type Item = JordanBasisMatrixVector< 'a, Mapping, RingOperator, OrderOperatorBoth, KeyBoth, EntryBoth, >;

    fn next(&mut self) -> Option<Self::Item> {
        let matching            =   self.umatch.matching_ref();
        let criterion           =   &self.criteria;
        self.iter_keymaj
            .find(  |x|  
                        (   criterion.bounding  &&  matching.contains_keymin(x) )
                    ||
                        (   criterion.boundary  &&  matching.contains_keymaj(x) )
                    ||
                        (   criterion.harmonic  &&  matching.lacks_key(x)       )            
                )
            .map(
                |index|
                JordanBasisMatrix::new( self.umatch ).view_minor_descend(index)
            )            
    }
}