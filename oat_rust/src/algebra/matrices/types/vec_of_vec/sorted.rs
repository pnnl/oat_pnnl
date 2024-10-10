//! Internal vectors are stored in sorted order
//! 
//! See the [parent module](super) for details on what vector-of-vectors format means.

use crate::algebra::vectors::entries::{KeyValGet, KeyValTypes};
use crate::algebra::matrices::query::{   ViewRow,
                                        ViewRowAscend,
                                        ViewRowDescend, ViewColDescend, IndicesAndCoefficients, ViewColAscend, ViewCol, MatrixEntry, MatrixOracle,
                                        // ViewCol, 
                                        // ViewColAscend,
                                        // ViewColDescend,
                                        // WhichMajor,
                                        // MajorDimension
                                    };

                                    use crate::utilities::binary_search::{find_sorted_binary_oracle};
                                    use crate::utilities::order::{JudgePartialOrder, is_sorted_strictly, OrderOperatorByKey, };
                                    use crate::utilities::statistics::histogram;


use rand::Rng;                                          // we use this module to generate random elements
use rand::distributions::{Bernoulli, Distribution};     // we use this module to generate random elements

use std::iter::{Rev, Cloned};
use std::marker::PhantomData;
use std::slice::Iter;
use itertools::Itertools;

use super::super::bimajor::{MatrixBimajor, MatrixBimajorData};



/// Vector-of-vectors sparse matrix format, with internal vectors stored in sorted order
/// 
/// See the [parent module](super) for details on what vector-of-vectors format means.
/// 
/// - Entries in each internal vector are stored in *strictly* ascending order of index (no repeat indices).
/// - Order of incides is determined by Rust's `PartialOrd` trait.
/// 
/// # See also
/// 
/// If you need a customized order for entries, see the [variant for custom orders](super::sorted_custom).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// 
/// // Create a new matrix.
/// let matrix  =   VecOfVec::new(
///                     vec![ vec![(1,1.)], vec![], vec![(2,2.)]  ],
///                 );
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VecOfVec
            < ColIndex, Coefficient >

{
    vec_of_vec:         Vec< Vec< ( ColIndex, Coefficient ) > >,
    order_operator:     OrderOperatorByKey < ColIndex, Coefficient, ( ColIndex, Coefficient ) >,
    // phantom_lifetime:   PhantomData< &'a ColIndex >,
}

impl < ColIndex, Coefficient >
        
        VecOfVec
            < ColIndex, Coefficient > 

{
    /// Make a new `VecOfVec`.  Throws an error if one of the internal vectors is not sorted
    /// in strictly asecending order, by index.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// // This line should not panic.
    /// let matrix = VecOfVec::new( vec![ vec![ (0,5), (1,6) ] ] ); // no panic
    /// 
    /// // This line should panic.
    /// // let matrix = VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ] );
    /// // And here is the test that confirms it.
    /// let result = std::panic::catch_unwind(|| {VecOfVec::new( vec![ vec![ (1,6), (0,5) ] ] )} );
    /// assert!(result.is_err()); 
    /// ```
    pub fn new( vecvec: Vec < Vec < ( ColIndex, Coefficient ) > > ) -> Self  
        where   ColIndex: PartialOrd,
                ( ColIndex, Coefficient ):     KeyValGet< ColIndex, Coefficient >, // if we comment this out then  we get an error sayint that `OrderOperatorByKey` doesn't implement the `JudgePartialOrder` trait; probably this has something to do with the fact that `OrderOperatorByKey` autoimplements `OrderOperatorByKey` for structs that implement `KeyValGet`
    {
        for vec in vecvec.iter() {
            if ! is_sorted_strictly( vec, &mut OrderOperatorByKey::new() ) {
                panic!("Attempt to construct `VecOfVec` from a vector of unsorted vectors.  To proceed, each internal vector must be sorted in *strictly* ascending order, by index.") 
            }
        }        

        VecOfVec{   
                    vec_of_vec:         vecvec,   
                    order_operator:   OrderOperatorByKey::new(),                 
                    // phantom_lifetime:   PhantomData,                  
                }
    }

    /// Returns a clone of the internally stored order comparator.
    pub fn clone_order_operator( &self ) -> OrderOperatorByKey < ColIndex, Coefficient, ( ColIndex, Coefficient ) > 
        where OrderOperatorByKey < ColIndex, Coefficient, ( ColIndex, Coefficient ) > :   Clone
    { self.order_operator.clone() }

    /// Returns an immutable view of the `Vec< Vec< (usize, Coefficient) > >` that stores the entries of of the matrix, internally.
    pub fn vec_of_vec_borrowed( &self ) -> & Vec< Vec< (ColIndex, Coefficient) > > { & self.vec_of_vec }


    /// Append a vector to the internally stored `Vec<Vec<EntryType>>` struct.  The new vector must be sorted in strictly ascending order; the function panics, otherwise.
    pub fn push_vec( &mut self, new_sorted_vec: Vec < (ColIndex, Coefficient) > ) 
        where 
            OrderOperatorByKey<ColIndex, Coefficient, (ColIndex, Coefficient)>: JudgePartialOrder< (ColIndex, Coefficient)>
    {
        if ! is_sorted_strictly( & new_sorted_vec, &mut self.order_operator ) {
            panic!("Attempt to append a non-strictly-sorted vector to `VecOfVec`.  To proceed, the new vector must be sorted in *strictly* ascending order, by index.") 
        }        
        self.vec_of_vec.push( new_sorted_vec );
    }    

    /// Ensures that `self` has at least `m` major views.
    /// 
    /// If `self` has fewer than `m` major slices, then
    /// - (if necessary) extend the capacity of `self.vec_of_vec` to `m` using `reserve_exact`
    /// - push as many empty vectors to `self.vec_of_vec` as necessary, to ensure there are `m` major views
    pub fn ensure_number_of_rows( &mut self, m: usize ) {
        if m > self.vec_of_vec.len() {
            if m > self.vec_of_vec.capacity() {
                self.vec_of_vec.reserve_exact( m-self.vec_of_vec.capacity() ) // reserve the exact right amount of extra capacity, if needed
            }
        }
        for p in 0 .. (m-self.vec_of_vec.len()) {
            self.vec_of_vec.push( Vec::with_capacity(0) );
        }

    }

    /// Creates a [`VecOfVec`] from an iterable that runs over iterables that run over tuples.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::query::ViewRowAscend;
    /// 
    /// let iter = (0..2).map( |x| vec![(x,x)] );
    /// let vec_of_vec = VecOfVec::from_iterable( iter );
    /// itertools::assert_equal( (& vec_of_vec).view_major_ascend(0), vec![(0,0)]);
    /// itertools::assert_equal( (& vec_of_vec).view_major_ascend(1), vec![(1,1)]) 
    /// ```
    pub fn from_iterable< I >( iter: I ) -> VecOfVec< ColIndex, Coefficient > 
        where
            I:          IntoIterator,
            I::Item:    IntoIterator< Item = (ColIndex, Coefficient) >,
            ColIndex:     Clone + PartialOrd,
            Coefficient:     Clone,
    {
        let vecvec =    iter.into_iter().map( |x| x.into_iter().collect_vec() ).collect_vec();
        VecOfVec::new( vecvec )
    }

    pub fn vec_of_vec( &self ) -> &Vec< Vec< ( ColIndex, Coefficient ) > > { &self.vec_of_vec }


}



impl < Coefficient >

    VecOfVec< usize, Coefficient > 
{

    /// Maximum column index of an explicitly stored entry (or None, if no entries are stored)
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// 
    /// let a: Vec<Vec<(usize,f64)>>    =   vec![vec![ (0,1.) ]];   // matrix with 1 entry
    /// let b: Vec<Vec<(usize,f64)>>    =   Vec::new();             // empty matrix
    /// 
    /// assert_eq!( VecOfVec::new(a).max_column_index(), Some(0) );
    /// assert_eq!( VecOfVec::new(b).max_column_index(), None    );
    /// ```
    pub fn max_column_index( &self ) -> Option< usize > {
        let mut max_col = None;
        for row in self.vec_of_vec.iter() {
            for col_index in row.iter().map(|x| x.0) {
                max_col = max_col.max( Some(col_index) )
            }
        }     
        max_col   
    }


    /// Prints a dense form of the matrix
    /// 
    /// - Entries that are not explicitly stored will appear as `filler_coefficient`.
    /// - The number of columns will equal the maximum column index of any explicitly stored entry.
    /// - If no entry is explicitly stored, then this function prints `Matrix contains no explicitly stored entries`
    pub fn print_dense( &self, filler_coefficient: Coefficient )
            where
                Coefficient:    std::fmt::Debug + Clone,
    {
        match self.max_column_index() {
            None => { println!("Matrix contains no explicitly stored entries") }
            Some( max_column_index ) => {
                let mut to_print = vec![ filler_coefficient.clone(); max_column_index + 1 ]; // we need +1 becuase of zero indexing
                for row_data in self.vec_of_vec.iter() {
                    for entry in to_print.iter_mut() { *entry = filler_coefficient.clone() }; // clear the entries
                    for (j,v) in row_data.iter().cloned() { to_print[j]=v } // fill with the nonzero entries
                    println!("{:?}", &to_print);
                }
            }
        }
    }

    /// Generates a new `VecOrVecSimple` representing the transpose, with copied (not borrowed) data.
    /// 
    /// Compare this method with [transpose](crate::algebra::matrices::operations::MatrixOperations::transpose), which
    /// generates a *lazy* transpose.
    /// 
    /// # Arguments
    /// 
    /// If `num_rows_of_transpose` is less than 1 + the maximum column index of a
    /// any explicitly stored entry in `self`, then returns `None`.  You can check the maximum column index with the method
    /// [VecOfVec::max_column_index]. Otherwise returns the transpose of a
    /// matrix of `self`, where `self` is regarded as a matrix of size `(self.vec_of_vec.len()) x num_rows_of_transpose`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::MatrixOperations;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::query::ViewRowAscend;
    /// 
    /// use itertools::Itertools;
    ///  
    /// let matrix     
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"0,0"), (1,"0,1"), ],
    ///                 vec![ (0,"1,0"), (1,"1,1"), ],
    ///             ]
    ///         );
    /// 
    /// let transpose   
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"0,0"), (1,"1,0"), ],
    ///                 vec![ (0,"0,1"), (1,"1,1"), ],
    ///             ]
    ///         );
    /// 
    /// let transpose_lazy  =   matrix.transpose();
    /// let transpose_deep  =   & matrix.transpose_deep(2).unwrap();
    ///
    ///  
    /// for p in 0..2 {
    ///     let a   =   transpose_deep.view_major_ascend(p).collect_vec();
    ///     let b   =   transpose_lazy.view_major_ascend(p).collect_vec();
    ///     let c   =   transpose.view_major_ascend(p).collect_vec();
    ///     assert_eq!(&a, &b);
    ///     assert_eq!(&a, &c);
    /// }
    /// 
    /// // check that matrix returns None when too few rows are specified
    /// assert_eq!( matrix.antitranspose_deep(1), None );
    /// ```
    pub fn transpose_deep( &self, num_rows_of_transpose: usize ) -> Option< Self >
    
        where
            Coefficient:  Clone
    {
        // count the explicit entries in each column
        let mut hist = vec![0; num_rows_of_transpose]; // a histogram that counts the number of nonzero entries in each column
        for row in self.vec_of_vec.iter() {
            for col_index in row.iter().map(|x| x.0) {
                if col_index + 1 > num_rows_of_transpose {
                    return None // in this case we haven't reserved enough space
                } else {
                    hist[ col_index ] +=1; // otherwise count the column in the historgram
                }
            }
        }

        // allocate exactly the right amount of space in the new matrix
        let mut transpose_data = Vec::with_capacity(num_rows_of_transpose);
        for num_indices in hist.into_iter() {
            transpose_data.push( Vec::with_capacity(num_indices) ); 
        }

        // fill the new matrix
        for row in 0 .. self.vec_of_vec.len() {
            for (col,val) in self.vec_of_vec[row].iter() {
                transpose_data[*col].push( (row,val.clone()) )
            }
        }

        return Some( VecOfVec { vec_of_vec: transpose_data, order_operator: OrderOperatorByKey::new() } )
    }


    /// Generates a new `VecOrVecSimple` representing the anti-transpose, with copied (not borrowed) data.
    /// 
    /// # Arguments
    /// 
    /// If `num_rows_of_antitranspose` is less than 1 + the maximum column index of a
    /// any explicitly stored entry in `self`, then `None` is returned.
    /// You can check the maximum column index of `self` with the method [VecOfVec::max_column_index].
    /// Otherwise returns the antitranspose of a
    /// matrix of `self`, where `self` is regarded as a matrix of size `(self.vec_of_vec.len()) x num_rows_of_antitranspose`.
    /// 
    /// # Caution
    /// 
    /// There are two important differences between this function and [antitranspose](crate::algebra::matrices::operations::MatrixOperations::transpose):
    /// 
    /// - [matrix.antitranspose()](crate::algebra::matrices::operations::MatrixOperations::transpose) is a lazy method that does not generate new data
    /// - the set of (key,val) pairs that appear in a major (respectively, minor) view of [matrix.antitranspose()](crate::algebra::matrices::operations::MatrixOperations::transpose)
    ///   are the *same* as the entries in a minor (respectively, major) view of `matrix`; only the sequence in which those entries appear is different.
    ///   By contrast, the keys in the (key,val) pairs of [matrix.antitranspose_deep()](VecOfVec::antitranspose_deep) are different;
    ///   they are obtained by subtracting the original keys from (# rows in the antitransposed matrix - 1).
    /// - For this reason, [matrix.antitranspose_deep()](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) is only available for
    ///   very specific types of matrices; [AntiTranspose] is available for a much broader class.
    /// 
    /// # Detail
    /// 
    /// This is equlvalent to calling [VecOfVec::transpose_deep], producing a struct of form `[ inner0, .., innerm ]`,
    /// then replacing this struct with `[ innerm, .., inner0 ]`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::operations::MatrixOperations;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::query::ViewRowAscend;
    /// 
    /// use itertools::Itertools;
    /// 
    /// let matrix     
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"0,0"), (1,"0,1"), (2,"0,2"), ],
    ///                 vec![ (0,"1,0"), (1,"1,1"), (2,"1,2"), ],
    ///             ]
    ///         );
    ///                 
    /// let antitranspose   
    ///     =   & VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"1,2"), (1,"0,2"), ],
    ///                 vec![ (0,"1,1"), (1,"0,1"), ],
    ///                 vec![ (0,"1,0"), (1,"0,0"), ],                        
    ///             ]
    ///         );
    ///                     
    /// let antitranspose_deep  =  matrix.antitranspose_deep(3).unwrap();
    /// let antitranspose_deep  =  & antitranspose_deep;
    ///                     
    /// // check that calculation is correct
    /// for p in 0..3 {
    ///     let a   =   antitranspose.view_major_ascend(p).collect_vec();
    ///     let b   =   antitranspose_deep.view_major_ascend(p).collect_vec();
    ///     assert_eq!(a, b)
    /// }
    /// 
    /// // check that matrix returns None when too few rows are specified
    /// assert_eq!( matrix.antitranspose_deep(2), None );
    /// ```
    pub fn antitranspose_deep( &self, num_rows_of_antitranspose: usize ) -> Option< Self >
    
        where
            Coefficient:  Clone + std::fmt::Debug,
    {
        let num_rows_original = self.vec_of_vec.len();
        
        let mut transpose   =   self.transpose_deep(num_rows_of_antitranspose)?;
        transpose.vec_of_vec.reverse();
        // println!("vec_vec len: {:?}", transpose.vec_of_vec.len() );
        // for p in 0..transpose.vec_of_vec.len() {
        //     println!("raw row before: {:?}", &transpose.vec_of_vec[p] );
        //     transpose.vec_of_vec[p].reverse();
        //     println!("raw row before: {:?}", &transpose.vec_of_vec[p] );
        // }        
        for row in &mut transpose.vec_of_vec {
            row.reverse(); // reverse order of entries
            for (j,v) in row.iter_mut() { 
                let jj  = num_rows_original - *j -1; // have to subtract 1, due to zero-indexing 
                *j = jj;  // replace each column index in the new matrix with (# columns in the new matrix) - column_index
            }
        }

        return Some(transpose)
    }    

    /// Returns a [MatrixBiajorData](crate::algebra::matrices::types::bimajor::MatrixBimajorData) that contains both row-major and colum-major copies of `self`
    /// 
    /// This doubles memory usage (it is not lazy); however it allows efficient access to **both rows and columns**.
    /// 
    /// # Arguments
    /// 
    /// If `num_rows_of_transpose` is less than 1 + the maximum column index of a
    /// any explicitly stored entry in `self`, then `None` is returned.  Otherwise returns the antitranspose of a
    /// matrix of `self`, where `self` is regarded as a matrix of size `(self.vec_of_vec.len()) x num_rows_of_transpose`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    /// use oat_rust::algebra::matrices::query::{ViewRow, ViewRowAscend, ViewRowDescend};
    /// use oat_rust::algebra::matrices::query::{ViewCol, ViewColAscend, ViewColDescend};        
    /// 
    /// use itertools::Itertools;
    /// 
    /// let matrix   
    ///     =   VecOfVec::new(
    ///             vec![
    ///                 vec![ (0,"a"), (1,"b"), (2,"c"), ],
    ///                 vec![ (0,"d"), (1,"e"), (2,"f"), ],
    ///             ]
    ///         );
    ///     
    /// let bimajor =   matrix.clone().bimajor(3).unwrap();
    ///     
    /// for p in 0..2 {
    ///     assert_eq!(     ( & matrix  ).view_major(p).collect_vec(),
    ///                     ( & bimajor ).view_major(p).collect_vec(),            );
    ///     assert_eq!(     ( & matrix  ).view_major_ascend(p).collect_vec(),
    ///                     ( & bimajor ).view_major_ascend(p).collect_vec(),     );
    ///     assert_eq!(     ( & matrix  ).view_major_descend(p).collect_vec(),
    ///                     ( & bimajor ).view_major_descend(p).collect_vec(),    );
    /// }
    /// 
    /// for p in 0..3 {
    ///     assert_eq!(     ( & matrix  ).view_minor(p).collect_vec(),
    ///                     ( & bimajor ).view_minor(p).collect_vec(),            );
    ///     assert_eq!(     ( & matrix  ).view_minor_ascend(p).collect_vec(),
    ///                     ( & bimajor ).view_minor_ascend(p).collect_vec(),     );
    ///     assert_eq!(     ( & matrix  ).view_minor_descend(p).collect_vec(),
    ///                     ( & bimajor ).view_minor_descend(p).collect_vec(),    );
    /// }     
    /// ```
    pub fn bimajor( self, num_rows_of_transpose: usize ) -> Option< MatrixBimajorData< Self, Self > >

        where
            Coefficient:  Clone + std::fmt::Debug,
    {        
        self.transpose_deep(num_rows_of_transpose).map(
            |transpose|
            MatrixBimajorData{ matrix_major_data: self, matrix_minor_data: transpose }
        )
    }        




}


impl VecOfVec< usize, usize > {

    //  RANDOM MATRICES
    //  ===========================================================================    

    //! Some functions that return matrices, e.g. `f(m) = random_m_by_m_sparse_matrix`.



    /// Random invertible upper-triangular row-major matrix with coefficients mod p
    /// 
    /// - The returned matrix has size `matrix_size` x `matrix_size`
    /// - Diagonal entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
    /// - Entries strictly above the diagonal are equally likely to be (i) structurally zero (and therefore algebraically) zero, 
    /// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
    /// In case (iii), entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
    pub fn random_mod_p_upper_triangular_invertible( 
            matrix_size: usize, 
            modulus: usize 
        ) ->
            VecOfVec< usize, usize >
    {
        let mut rng = rand::thread_rng(); // this generates random integers
        let mut vec_of_vec = vec![];
        for majkey in 0 .. matrix_size {
            let coefficient_leading         =   rng.gen_range( 1 .. modulus );
            let mut new_vec     =   vec![ (majkey, coefficient_leading) ]; // start with a nonzero diagonal element            
            for q in majkey+1 .. matrix_size { // fill out the rest of this row of the matrix
                let coefficient   = rng.gen_range( 1 .. modulus );
                let flag = rng.gen_range(0usize .. 3usize);
                if      flag == 0 { new_vec.push( ( q, 0 )           ) }
                else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
                else              { continue }
            }
            vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        }

        VecOfVec::new(vec_of_vec)// formally wrap the matrix in a VecOfVec struct
    }

    /// Generate a random invertible upper triangular matrix with integer coefficients.
    /// 
    /// - Entries on the diagonal are equal to 1.
    /// - Each entry strictly above the diagonal is either (i) strucutrually zero, 
    /// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
    /// Algebraically nonzero entries are drawn iid from the uniform distribution on [1, .., `modulus`).
    /// 
    /// Nonzero entries are drawn in an iid fasion from [1, .., `modulus`).
    pub fn random_mod_p_upper_unitriangular( 
            matrix_size: usize, 
            modulus: usize
        ) -> 
            VecOfVec< usize, usize > {

        let mut rng = rand::thread_rng(); // this generates random integers
        let mut vec_of_vec = vec![];
        for majkey in 0 .. matrix_size {
            let coefficient_leading         =   1;
            let mut new_vec     =   vec![ (majkey, coefficient_leading) ]; // start with a nonzero diagonal element            
            for q in majkey+1 .. matrix_size { // fill out the rest of this row of the matrix
                let coefficient   = rng.gen_range( 0 .. modulus );
                let flag = rng.gen_range(0usize .. 3usize);
                if      flag == 0 { new_vec.push( ( q, 0 )           ) }
                else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
                else              { continue }
            }
            vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        }

        VecOfVec::new( vec_of_vec ) // formally wrap the matrix in a VecOfVec struct
    }


    /// Generate a random row-major matrix of size `num_rows` x `num_cols` with coefficients
    /// in `{ 0, .., p-1 }`.
    /// 
    /// Each entry is equally likely to be (i) structurally zero, 
    /// (ii) structurally nonzero but algebraically zero, or (iii) structurally and algebraically nonzero.
    /// In case (iii), entries are drawn independently from the uniform distribution on `{ 1, .., p-1 }`.
    pub fn random_mod_p( 
            num_rows: usize, 
            num_cols: usize, 
            modulus: usize 
        ) ->
            VecOfVec< usize, usize >
    {
        let mut rng = rand::thread_rng(); // this generates random integers
        let mut vec_of_vec = vec![];
        for majkey in 0 .. num_rows {
            let mut new_vec     =   vec![]; // start with an empty vector
            for q in majkey+1 .. num_cols { // fill it in
                let coefficient   = rng.gen_range( 1 .. modulus );
                let flag = rng.gen_range(0usize .. 3usize);
                if      flag == 0 { new_vec.push( ( q, 0 )           ) }
                else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
                else              { continue }
            }
            vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        }

        VecOfVec::new(vec_of_vec)// formally wrap the matrix in a VecOfVec struct
    }

    /// Generate a random `VecOfVec< usize, usize >`
    /// 
    /// Nonzero entries are drawn in an iid fasion from [0, .., `modulus`) if `allow_nonstructural_zero == true` and from [1, .., `modulus`)
    /// if `allow_nonstructural_zero == false`.
    pub fn random_mod_p_with_density( 
                num_indices_major:          usize, 
                num_indices_minor:          usize, 
                approximate_density:        f64, 
                modulus:                    usize,
                allow_nonstructural_zero:   bool,
            ) 
            -> 
            VecOfVec< usize, usize > 
    {

        let mut rng = rand::thread_rng(); // this generates random integers

        let d = Bernoulli::new( approximate_density ).unwrap(); // bernoulli random variable returning `true` with probability `approximate_density`
        let v = d.sample(&mut rand::thread_rng());
        println!("{} is from a Bernoulli distribution", v);

        let mut vecvec = Vec::new(); // initialize empty vector of vectors
        for keymaj in 0 .. num_indices_major {
            vecvec.push( Vec::new() ); // push a new vector representing a major view
            for keymin in 0 .. num_indices_minor { // for each minor index
                if d.sample( &mut rand::thread_rng() ) { // add a structural nonzero entry with probability `approximate_density`
                    let coefficient   = match allow_nonstructural_zero{ 
                        true => { rng.gen_range( 0 .. modulus ) },
                        false => { rng.gen_range( 1 .. modulus ) }
                    };
                    vecvec[ keymaj ].push( (keymin, coefficient) );
                }
            }
            vecvec[ keymaj ].shrink_to_fit();
        }

        VecOfVec::new( vecvec )

    }    
}


//  ORACLE IMPLEMENTATIONS FOR &'a VecOfVec
//  ---------------------------------------------

///  The Matrix Oracle trait
/// 
///  # Details
/// 
/// For details see the documentation for the MatrixOracle trait.
impl < 'a, ColIndex, Coefficient > 

    MatrixOracle for 
    &'a VecOfVec < ColIndex, Coefficient >

    where   ColIndex:       Clone + Ord,    
            Coefficient:    Clone,    

{ 
    type Coefficient        =   Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   usize;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   ColIndex;       // The type of column indices
    
    type RowEntry           =   ( ColIndex, Coefficient );          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   ( usize, Coefficient );       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Cloned< Iter< 'a, (ColIndex, Coefficient) > >;  // What you get when you ask for a row.
    type RowIter            =   Cloned< Iter< 'a, (ColIndex, Coefficient) > >;  // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse         =   Cloned<Rev<std::slice::Iter<'a, (ColIndex, Coefficient)>>>;  // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter     =   Cloned<Rev<std::slice::Iter<'a, (ColIndex, Coefficient)>>>;  // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    
    type Column             =   VecOfVecMatrixColumn< 'a, ColIndex, Coefficient >; // What you get when you ask for a column
    type ColumnIter         =   VecOfVecMatrixColumn< 'a, ColIndex, Coefficient >; // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse      =   VecOfVecMatrixColumnReverse< 'a, ColIndex, Coefficient >; // What you get when you ask for a column with the order of entries reversed 
    type ColumnReverseIter  =   VecOfVecMatrixColumnReverse< 'a, ColIndex, Coefficient >; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    // entry lookup
    fn entry( & self, row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient > {
        let row = & self.vec_of_vec[ row ];
        find_sorted_binary_oracle( 
                    0, 
                    row.len() as isize - 1, 
                    |p| column.cmp( & row[p as usize].0 ) 
                ).map(|x| row[ x as usize ].1.clone() )        
    }  

    // row lookup
    fn row(                     & self,  index: Self::RowIndex    )   -> Self::Row   { 
        self.vec_of_vec[index].iter().cloned() 
    }
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>  {
        if index >= self.vec_of_vec.len() { return None }
        else { return Some( self.row( index ) ) }
    }
    fn row_reverse(             & self,  index: Self::RowIndex    )       -> Self::RowReverse  { 
        self.vec_of_vec[index].iter().rev().cloned()
    }
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>  {
        if index >= self.vec_of_vec.len() { return None }
        else { return Some( self.row_reverse( index ) ) }        
    }  
    
    // column lookup
    fn column(                  & self,  index: Self::ColumnIndex )       -> Self::Column {
        VecOfVecMatrixColumn{
            vec_of_vec:             self,
            keymaj:                 0,
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column> {
        Some( self.column(index) )
    }   
    fn column_reverse(          & self,  index: Self::ColumnIndex )       -> Self::ColumnReverse  { 
        VecOfVecMatrixColumnReverse{
            vec_of_vec:             self,
            keymaj:                 self.vec_of_vec.len(),
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }            
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse>{
        Some( self.column_reverse(index) )
    }
}   



//  IndicesAndCoefficients
impl < 'a, ColIndex, Coefficient > 

    IndicesAndCoefficients for 
    &'a VecOfVec < ColIndex, Coefficient >

{ 
    type EntryMajor =   (ColIndex, Coefficient);    
    type EntryMinor = (usize, Coefficient);    
    type ColIndex = ColIndex; 
    type RowIndex = usize; 
    type Coefficient = Coefficient; }   



// MatrixEntry
impl < 'a, ColIndex, Coefficient > 
    
    MatrixEntry for 
    
    &'a VecOfVec 
        < ColIndex, Coefficient > 

    where   ColIndex:     Clone + Ord,    
            Coefficient:     Clone,
{      
    /// Retrieve the entry at the location specified by `keymaj` and `keymin`.
    fn entry_major_at_minor( &self, keymaj: Self::RowIndex, keymin: Self::ColIndex, ) -> Option< Self::Coefficient > {
        let view = & self.vec_of_vec[ keymaj ];
        find_sorted_binary_oracle( 
                    0, 
                    view.len() as isize - 1, 
                    |p| keymin.cmp( & view[p as usize].0 ) 
                ).map(|x| view[ x as usize ].1.clone() )
    }
}


// ViewRow
impl < 'a, ColIndex, Coefficient > 
    
    ViewRow for 
    
    &'a VecOfVec 
        < ColIndex, Coefficient > 

    where   ColIndex:     Clone,    
            Coefficient:     Clone,
{      
    type ViewMajor          =   Cloned< Iter< 'a, (ColIndex, Coefficient) > >;
    type ViewMajorIntoIter  =   Cloned< Iter< 'a, (ColIndex, Coefficient) > >;

    fn view_major( & self, index: usize ) -> Cloned< Iter< 'a, (ColIndex, Coefficient) > > {
        return self.vec_of_vec[index].iter().cloned()
    } 
}

// impl < 'a, ColIndex, Coefficient > OracleRefInherit for &'a VecOfVec < ColIndex, Coefficient > {}

impl < 'a, ColIndex, Coefficient > 
    
    ViewRowAscend  for 
    
    &'a VecOfVec < ColIndex, Coefficient > 

    where   ColIndex:     Clone,    
            Coefficient:     Clone,    
{
    type ViewMajorAscend            =   Cloned< Iter< 'a, (ColIndex, Coefficient) > >;
    type ViewMajorAscendIntoIter    =   Cloned< Iter< 'a, (ColIndex, Coefficient) > >;

    /// Assumes that entries in each vector are sorted in ascending order.
    fn view_major_ascend( & self, index: usize ) -> Cloned< Iter< 'a, (ColIndex, Coefficient) > > {
        self.view_major( index )
    } 
}

impl < 'a, ColIndex, Coefficient > 
    
    ViewRowDescend  for 
    
    &'a VecOfVec < ColIndex, Coefficient > 

    where   ColIndex:     Clone,    
            Coefficient:     Clone,  
{
    type ViewMajorDescend           =   Cloned<Rev<std::slice::Iter<'a, (ColIndex, Coefficient)>>>;
    type ViewMajorDescendIntoIter   =   Cloned<Rev<std::slice::Iter<'a, (ColIndex, Coefficient)>>>;

    /// Assumes that entries in each vector are sorted in ascending order.    
    fn view_major_descend( & self, index: usize ) -> Cloned<Rev<std::slice::Iter<'a, (ColIndex, Coefficient)>>> {
        return self.vec_of_vec[index].iter().rev().cloned()
    } 
}

/// Represents a minor view of a `VecOfVec`, with entries appearing in descending order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::{VecOfVec, VecOfVecMatrixColumnReverse};
/// use oat_rust::algebra::matrices::query::ViewColDescend;
///         
/// // define a row-major sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_descend( 0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_descend( 1 ) );
/// itertools::assert_equal( vec![ (1,7), (0,5) ], (& matrix).view_minor_descend( 2 ) );  
/// ```
pub struct VecOfVecMatrixColumnReverse< 'a, ColIndex, Coefficient > 
    where   ColIndex:     Clone,    
            Coefficient:     Clone,  
{
    vec_of_vec:             &'a VecOfVec< ColIndex, Coefficient >,
    keymaj:                 usize,
    keymin:                 ColIndex,
    phantom_keymin:         PhantomData< ColIndex >,
}

// Iterator
impl < 'a, ColIndex, Coefficient > 

        Iterator for 

        VecOfVecMatrixColumnReverse< 'a, ColIndex, Coefficient >

    where   ColIndex:     Clone + PartialEq,    
            Coefficient:     Clone,          
{
    type Item = (usize, Coefficient);

    fn next( &mut self ) -> Option< Self::Item > {
        // extract the underlying vector of vectors
        let vecvec_data =  self.vec_of_vec.vec_of_vec_borrowed();                
        // (assuming the matrix is row-major) scan the rows of the matrix from bottom to top
        while self.keymaj > 0 {
            // drop the row number by 1
            self.keymaj -= 1;            
            // get the row
            let view_major = & vecvec_data[ self.keymaj ];
            // scan the row to see if it contains an entry of form ( my_keymin, snzval )
            for ( keymin, snzval ) in view_major {
                // if it does, then return ( row_number, snzval )
                if keymin != & self.keymin { continue }
                else { return Some( ( self.keymaj, snzval.clone() ) ) }
            }
        }
        None
    }
}        

/// Represents a minor view of a `VecOfVec`, with entries appearing in ascending order of index.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::{VecOfVec, VecOfVecMatrixColumn};
/// use oat_rust::algebra::matrices::query::ViewColAscend;
///         
/// // define a row-major sparse matrix representing the array
/// //  0   5   5
/// //  0   0   7
/// let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
/// itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_ascend( 0 ) );
/// itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_ascend( 1 ) );
/// itertools::assert_equal( vec![ (0,5), (1,7) ], (& matrix).view_minor_ascend( 2 ) );  
/// ```
pub struct VecOfVecMatrixColumn< 'a, ColIndex, Coefficient > 
    where   ColIndex:     Clone,    
            Coefficient:     Clone,  
{
    vec_of_vec:             &'a VecOfVec< ColIndex, Coefficient >,
    keymaj:                 usize,
    keymin:                 ColIndex,
    phantom_keymin:         PhantomData< ColIndex >,
}

// Iterator
impl < 'a, ColIndex, Coefficient > 

        Iterator for 

        VecOfVecMatrixColumn< 'a, ColIndex, Coefficient >

    where   ColIndex:     Clone + PartialEq,    
            Coefficient:     Clone,          
{
    type Item = (usize, Coefficient);

    fn next( &mut self ) -> Option< Self::Item > {
        // extract the underlying vector of vectors
        let vecvec_data =  self.vec_of_vec.vec_of_vec_borrowed();                
        // (assuming the matrix is row-major) scan the rows of the matrix from bottom to top
        while self.keymaj < vecvec_data.len() {
            // get the row
            let view_major = & vecvec_data[ self.keymaj ];
            // scan the row to see if it contains an entry of form ( my_keymin, snzval )
            for ( keymin, snzval ) in view_major {
                // if it does, then return ( row_number, snzval )
                if keymin != & self.keymin { continue }
                else { 
                    // grow the row number by 1 (this is one of two branches where we do this)
                    self.keymaj += 1;
                    
                    // return the entry
                    return Some( ( self.keymaj - 1, snzval.clone() ) ) 
                }
            }
            // grow the row number by 1 (this is one of two branches where we do this)
            self.keymaj += 1;              
        }      

        // in this case the iterator is exhausted
        None
    }
}            


//  ViewColDescend
impl < 'a, ColIndex, Coefficient > 

    ViewColDescend for 
    
    &'a VecOfVec
        < ColIndex, Coefficient >

    where   ColIndex:     Clone + PartialEq,    
            Coefficient:     Clone,          
    
{
    type ViewMinorDescend = VecOfVecMatrixColumnReverse< 'a, ColIndex, Coefficient >;
    type ViewMinorDescendIntoIter = Self::ViewMinorDescend;

    fn view_minor_descend( &self, index: ColIndex ) -> VecOfVecMatrixColumnReverse< 'a, ColIndex, Coefficient > {
        VecOfVecMatrixColumnReverse{
            vec_of_vec:             self,
            keymaj:                 self.vec_of_vec.len(),
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }
}

//  ViewColAscend
impl < 'a, ColIndex, Coefficient > 

    ViewColAscend for 
    
    &'a VecOfVec
        < ColIndex, Coefficient >

    where   ColIndex:     Clone + PartialEq,    
            Coefficient:     Clone,          
    
{
    type ViewMinorAscend = VecOfVecMatrixColumn< 'a, ColIndex, Coefficient >;
    type ViewMinorAscendIntoIter = Self::ViewMinorAscend;

    fn view_minor_ascend( &self, index: ColIndex ) -> VecOfVecMatrixColumn< 'a, ColIndex, Coefficient > {
        VecOfVecMatrixColumn{
            vec_of_vec:             self,
            keymaj:                 0,
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }
}

//  ViewCol
impl < 'a, ColIndex, Coefficient > 

    ViewCol for 
    
    &'a VecOfVec
        < ColIndex, Coefficient >

    where   ColIndex:     Clone + PartialEq,    
            Coefficient:     Clone,          
    
{
    type ViewMinor = VecOfVecMatrixColumn< 'a, ColIndex, Coefficient >;
    type ViewMinorIntoIter = Self::ViewMinor;

    fn view_minor( &self, index: ColIndex ) -> VecOfVecMatrixColumn< 'a, ColIndex, Coefficient > {
        VecOfVecMatrixColumn{
            vec_of_vec:             self,
            keymaj:                 0,
            keymin:                 index,
            phantom_keymin:         PhantomData,            
        }
    }
}







#[cfg(test)]
mod tests {    

    use crate::algebra::matrices::display::print_indexed_major_views;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_vec_of_vec_from_iterable() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewRowAscend;
        
        let iter = (0..2).map( |x| vec![(x,x)] );
        let vec_of_vec = VecOfVec::from_iterable( iter );
        itertools::assert_equal( (& vec_of_vec).view_major_ascend(0), vec![(0,0)]);
        itertools::assert_equal( (& vec_of_vec).view_major_ascend(1), vec![(1,1)])   
    }

    #[test]    
    fn test_vec_of_vec_simple_descending_minor_view() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewColDescend;
        
        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
        itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_descend( 0 ) );
        itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_descend( 1 ) );
        itertools::assert_equal( vec![ (1,7), (0,5) ], (& matrix).view_minor_descend( 2 ) );        
    }

    #[test]     
    fn test_vec_of_vec_simple_ascending_minor_view() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewColAscend;

        // define a row-major sparse matrix representing the array
        //  0   5   5
        //  0   0   7
        let matrix = VecOfVec::new( vec![ vec![ (1,5), (2,5) ], vec![ (2,7) ] ] );
        println!("waytpoint 1");
        itertools::assert_equal(Vec::<(usize,i32)>::new() , (& matrix).view_minor_ascend( 0 ) );
        println!("waytpoint 2");        
        itertools::assert_equal( vec![ (0,5) ], (& matrix).view_minor_ascend( 1 ) );
        println!("waytpoint 3");        
        itertools::assert_equal( vec![ (0,5), (1,7) ], (& matrix).view_minor_ascend( 2 ) ); 
    }

    #[test]    
    fn test_vec_of_vec_simple_antitranspose_deep() {
        use crate::algebra::matrices::operations::MatrixOperations;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::ViewRowAscend;
        
        use itertools::Itertools;
        
        let matrix     
            =   & VecOfVec::new(
                    vec![
                        vec![ (0,"0,0"), (1,"0,1"), (2,"0,2"), ],
                        vec![ (0,"1,0"), (1,"1,1"), (2,"1,2"), ],
                    ]
                );

        let antitranspose   
            =   & VecOfVec::new(
                    vec![
                        vec![ (0,"1,2"), (1,"0,2"), ],
                        vec![ (0,"1,1"), (1,"0,1"), ],
                        vec![ (0,"1,0"), (1,"0,0"), ],                        
                    ]
                );

        let antitranspose_deep  =  matrix.antitranspose_deep(3).unwrap();
        let antitranspose_deep  =  & antitranspose_deep;

        // check that calculation is correct
        for p in 0..3 {
            let a   =   antitranspose.view_major_ascend(p).collect_vec();
            let b   =   antitranspose_deep.view_major_ascend(p).collect_vec();
            assert_eq!(a, b)
        }

        // check that matrix returns None when too few rows are specified
        assert_eq!( matrix.antitranspose_deep(2), None );
    }



    #[test]
    fn test_vec_of_vec_simple_bimajor_data() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::{ViewRow, ViewRowAscend, ViewRowDescend};
        use crate::algebra::matrices::query::{ViewCol, ViewColAscend, ViewColDescend};        
        
        use itertools::Itertools;
        
        let matrix   
            =   VecOfVec::new(
                    vec![
                        vec![ (0,"a"), (1,"b"), (2,"c"), ],
                        vec![ (0,"d"), (1,"e"), (2,"f"), ],
                    ]
                );
            
        let bimajor =   matrix.clone().bimajor(3).unwrap();
            
        for p in 0..2 {
            assert_eq!(     ( & matrix  ).view_major(p).collect_vec(),
                            ( & bimajor ).view_major(p).collect_vec(),            );
            assert_eq!(     ( & matrix  ).view_major_ascend(p).collect_vec(),
                            ( & bimajor ).view_major_ascend(p).collect_vec(),     );
            assert_eq!(     ( & matrix  ).view_major_descend(p).collect_vec(),
                            ( & bimajor ).view_major_descend(p).collect_vec(),    );
        }
        for p in 0..3 {
            assert_eq!(     ( & matrix  ).view_minor(p).collect_vec(),
                            ( & bimajor ).view_minor(p).collect_vec(),            );
            assert_eq!(     ( & matrix  ).view_minor_ascend(p).collect_vec(),
                            ( & bimajor ).view_minor_ascend(p).collect_vec(),     );
            assert_eq!(     ( & matrix  ).view_minor_descend(p).collect_vec(),
                            ( & bimajor ).view_minor_descend(p).collect_vec(),    );
        }        
    }

}
