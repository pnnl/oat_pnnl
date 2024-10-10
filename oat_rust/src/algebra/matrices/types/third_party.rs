//! Wrappers for third-party matrix types
//! 
//! # Tips on construction
//! 
//! The best reference for third-party matrix structs is their documentation.  However,
//! here are a few common construction methods
//! 
//! Build a dense 2 dimensional ndarray from a vector of vectors
//! 
//! ```
//! let mat = vec![ vec![1;2]; 3 ]; // a vector of vectors
//! let mat = ndarray::Array2::from_shape_vec((3, 2), mat);
//! ```
//! 
//! Build a sparse 2x2 identity matrix two different ways
//! 
//! ```
//! use oat_rust::algebra::matrices::types::third_party::IntoCSR;
//! 
//! let eye =   vec![ vec![(0,1)], vec![(1,1)] ].into_csr(2,2);
//! 
//! let eye =   vec![ (0,0,1), (1,1,1) ].into_csr(2,2);
//! ```

use crate::algebra::matrices::query::{ViewRowAscend, IndicesAndCoefficients, MatrixEntry, MatrixOracle, column_helper::SparseColumn};

use std::iter::Rev;
use std::ops::{Deref, Range};
use std::sync::Arc;

use sprs::{ TriMat, CsMatBase, SpIndex, TriMatBase };
use sprs::vec::{VectorIterator, IntoSparseVecIter};


//  =====================================================================================
//  & sprs::CsMatBase
//  =====================================================================================


//  VIEW TYPE
//  ------------------------

/// A wrapper for the `sprs` struct `VectorIterator<'a, N, I >`, which returns
/// `(usize, N)` instead of `(usize, &N)`
pub struct SprsRowIteratorCloned< 'a, N, I > {
    vector_iterator: VectorIterator<'a, N, I >
}

impl < 'a, N, I >

    Iterator for 

    SprsRowIteratorCloned
        < 'a, N, I >

    where
        I: SpIndex,
        N: Clone,          
{
    type Item = ( usize, N );
    fn next(&mut self) -> Option<Self::Item> {
        self.vector_iterator.next().map(|(k,v)| (k,v.clone()))
    }
}  


/// Runs over the entries of an outer view of a `sprs::CsMatBase`
/// 
/// **Note** "Outer views" of a CSR matrix are rows.  Outer views of a
/// CSC matrix are columns.
/// 
/// **In OAT we currently regard every instance of `sprs::CsMatBase` as
/// a CSR matrix, even if the `sprs` module regards it as a CSC matrix
/// (that is, even if `matrix.storage()` returns `sprs::CompressedStorage::CSR`)**.
/// This is equivalent to saying that we always regard "outer views" of a `CsMatBase`
/// as rows.
pub struct SprsRowIterator< 'a, N, I, > {
    ordinals:           Range< usize >,
    coefficients:       &'a [N],
    column_indices:     &'a [I],
}

impl < 'a, N, I, >

    Iterator for 

    SprsRowIterator< 'a, N, I, >

    where 
        I:      SpIndex,
        N:      Clone,
{
    type Item   =   (usize, N);
    fn next( &mut self ) -> Option< Self::Item > {
        self.ordinals.next().map(
            |i| 
            ( self.column_indices[i].index(), self.coefficients[i].clone() )
        )
    }
}    

/// Runs over the entries of an outer view of a `sprs::CsMatBase` **in reverse order**
/// 
/// **Note** "Outer views" of a CSR matrix are rows.  Outer views of a
/// CSC matrix are columns.
/// 
/// **In OAT we currently regard every instance of `sprs::CsMatBase` as
/// a CSR matrix, even if the `sprs` module regards it as a CSC matrix
/// (that is, even if `matrix.storage()` returns `sprs::CompressedStorage::CSR`)**.
/// This is equivalent to saying that we always regard "outer views" of a `CsMatBase`
/// as rows.
pub struct SprsRowReverseIterator< 'a, N, I > {
    ordinals:           Rev<Range<usize>>,
    coefficients:       &'a [N],
    column_indices:     &'a [I],
}

impl < 'a, N, I >

    Iterator for 

    SprsRowReverseIterator< 'a, N, I >

    where 
        I:      SpIndex,
        N:      Clone,    
{
    type Item   =   (usize, N);
    fn next( &mut self ) -> Option< Self::Item > {
        self.ordinals.next().map(
            |i| {
                let i = i.index(); // convert the index to usize
                ( self.column_indices[i].index(), self.coefficients[i].clone() ) // return the column index and coefficient
            }
        )
    }
}  


//  -------------------------------------------------------------------------------------
//  & sprs::CsMatBase
//  -------------------------------------------------------------------------------------


//  MATRIX ORACLE
//  -------------

impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >

    MatrixOracle for

    &'a CsMatBase
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    where
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,
        N: Clone,
{           

    type Coefficient            =   N; // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   usize; // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   usize; // The type of column indices    
    type RowEntry               =   (usize,N); // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   (usize,N); // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =   SprsRowIterator< 'a, N, I, >;  // What you get when you ask for a row.
    type RowIter                =   SprsRowIterator< 'a, N, I, >;  // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse             =   SprsRowReverseIterator< 'a, N, I, >;  // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter         =   SprsRowReverseIterator< 'a, N, I, >;  // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    type Column                 =   SparseColumn< // What you get when you ask for a row.
                                        &'a CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr >,
                                        Range< usize >,
                                    >;
    type ColumnIter             =   SparseColumn< // What you get when you call `column.into_iter()`, where `column` is a column
                                        &'a CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr >,
                                        Range< usize >,
                                    >;  
    type ColumnReverse          =   SparseColumn<  // What you get when you ask for a column with the order of entries reversed                             
                                        &'a CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr >,
                                        Rev< Range< usize > >,
                                    >;  
    type ColumnReverseIter      =   SparseColumn<  // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)
                                        &'a CsMatBase < N, I, IptrStorage, IndStorage, DataStorage, Iptr >,
                                        Rev< Range< usize > >,
                                    >; 

    fn entry(                   &   self, row: Self::RowIndex, col: Self::ColumnIndex ) ->  Option< Self::Coefficient > { 
        if self.is_csc() {
            panic!("A program called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats.");  
        }
        self.get( row, col ).cloned()
    }
    fn row(                     &   self, index: Self::RowIndex   )   -> Self::Row {   
        if index >= self.rows() {
            panic!("A program called `self.row( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows.", index, self.rows() );
        }    
        self.row_opt( index ).unwrap()  
    } 
    fn row_opt(                 &   self, index: Self::RowIndex   )   -> Option<Self::Row> {         
        if self.is_csc() { 
            panic!("A program called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if index >= self.rows() {
            None
        } else {
            let ordinals_raw: Range<Iptr>   =   self.indptr().outer_inds(index);
            let ordinals: Range< usize >    =   ordinals_raw.start.index() .. ordinals_raw.end.index();
            let column_indices: &[I]        =   self.indices();
            let coefficients: &[N]          =   self.data();
            Some( SprsRowIterator{ ordinals, coefficients, column_indices } )
        }
    }
    fn row_reverse(             &   self, index: Self::RowIndex   )   -> Self::RowReverse {  
        if index >= self.rows() {
            panic!("A program called `self.row( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows.", index, self.rows() );
        }         
        self.row_reverse_opt( index ).unwrap()  
    }     
    fn row_reverse_opt(         &   self, index: Self::RowIndex   )   -> Option<Self::RowReverse>
    { 
        if self.is_csc() { 
            panic!("A program called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if index >= self.rows() {
            None // the index is too big
        } else {
            let ordinals_raw: Range<Iptr>       =   self.indptr().outer_inds(index);
            let ordinals: Rev< Range< usize > > =   ( ordinals_raw.start.index() .. ordinals_raw.end.index() ).rev();
            let column_indices: &[I]            =   self.indices();
            let coefficients: &[N]              =   self.data();
            Some( SprsRowReverseIterator{ ordinals, coefficients, column_indices } )
        }
    }    
    /// Get a column
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::query::{MatrixOracle};
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; 
    /// 
    /// 
    /// // data for a matrix:
    /// // |  1   0 |
    /// // | -1   0 |
    /// // |  0   0 |
    /// // |  0   0 |           
    /// let matrix = vec![ 
    ///         vec![ (0, 1.)  ], 
    ///         vec![ (0,-1.)  ], 
    ///         vec![          ], 
    ///         vec![          ],                 
    ///     ];
    /// 
    /// // wrap the data in a VecOfVec sparse matrix struct
    /// let matrix = VecOfVec::new(matrix);
    /// 
    /// // get some sparse columns
    /// let mut column_0        =   (& matrix ).column(0);
    /// let mut column_0_rev    =   (& matrix ).column_reverse(0);
    /// let mut column_1        =   (& matrix ).column(1);
    /// let mut column_1_rev    =   (& matrix ).column_reverse(1);
    /// 
    /// // check the columns are correct
    /// assert_eq!( column_0.next(),        Some( (0,  1.) )    );
    /// assert_eq!( column_0.next(),        Some( (1, -1.) )    );
    /// assert_eq!( column_0.next(),        None                );    
    /// 
    /// assert_eq!( column_0_rev.next(),    Some( (1, -1.) )    );
    /// assert_eq!( column_0_rev.next(),    Some( (0,  1.) )    );
    /// assert_eq!( column_0_rev.next(),    None                );   
    /// 
    /// assert_eq!( column_1.next(),        None        );
    /// assert_eq!( column_1_rev.next(),    None        );   
    /// ```
    fn column(                  &   self, index: Self::ColumnIndex)   -> Self::Column    { 
        if index > self.cols() {
            panic!("A program called `self.col( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows", index, self.rows() );
        }
        self.column_opt(index).unwrap() 
    }
    fn column_opt(              &   self, index: Self::ColumnIndex)   -> Option<Self::Column> {
        if self.is_csc() { 
            panic!("A program called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if index >= self.cols() {
            return None // if the index is too great, then return None
        } else {
            let row_index_iterator  =   0 .. self.rows();
            let matrix              =   self;
            let column_index        =   index;
            Some( SparseColumn{ matrix, row_index_iterator, column_index } )
        }
    }   
    /// Get a column with entries in reverse order
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::query::{MatrixOracle};
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; 
    /// 
    /// 
    /// // data for a matrix:
    /// // |  1   0 |
    /// // | -1   0 |
    /// // |  0   0 |
    /// // |  0   0 |           
    /// let matrix = vec![ 
    ///         vec![ (0, 1.)  ], 
    ///         vec![ (0,-1.)  ], 
    ///         vec![          ], 
    ///         vec![          ],                 
    ///     ];
    /// 
    /// // wrap the data in a VecOfVec sparse matrix struct
    /// let matrix = VecOfVec::new(matrix);
    /// 
    /// // get some sparse columns
    /// let mut column_0        =   (& matrix ).column(0);
    /// let mut column_0_rev    =   (& matrix ).column_reverse(0);
    /// let mut column_1        =   (& matrix ).column(1);
    /// let mut column_1_rev    =   (& matrix ).column_reverse(1);
    /// 
    /// // check the columns are correct
    /// assert_eq!( column_0.next(),        Some( (0,  1.) )    );
    /// assert_eq!( column_0.next(),        Some( (1, -1.) )    );
    /// assert_eq!( column_0.next(),        None                );    
    /// 
    /// assert_eq!( column_0_rev.next(),    Some( (1, -1.) )    );
    /// assert_eq!( column_0_rev.next(),    Some( (0,  1.) )    );
    /// assert_eq!( column_0_rev.next(),    None                );   
    /// 
    /// assert_eq!( column_1.next(),        None        );
    /// assert_eq!( column_1_rev.next(),    None        );   
    /// ```    
    fn column_reverse(          &   self, index: Self::ColumnIndex)   -> Self::ColumnReverse       {
        if index > self.cols() {
            panic!("A program called `self.col( {:?} )` on a `sprs:CsMatBase` object.  However, the object only has {:?} rows", index, self.rows() );
        }
        self.column_reverse_opt(index).unwrap() 
    }        
    fn column_reverse_opt(      &   self, index: Self::ColumnIndex)   -> Option<Self::ColumnReverse> {
        if self.is_csc() { 
            panic!("A program called a method from the `MatrixOracle` trait on a `sprs:CsMatBase` object.  However, the object is in CSC format.  `MatrixOracle` is only available for CSR format.  The `sprs` library offers convenient methods for converting between these two formats."); 
        } else if index >= self.cols() {
            return None // if the index is too great, then return None
        } else {
            let row_index_iterator  =   ( 0 .. self.rows() ).rev();
            let matrix              =   self;
            let column_index        =   index;
            Some( SparseColumn{ matrix, row_index_iterator, column_index } )        
        }
    }

} 


//  INDICES AND COEFFICIENTS
//  ------------------------

impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >

    IndicesAndCoefficients for

    &'a CsMatBase
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    where
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,     
{
    type RowIndex=usize; type ColIndex=usize; type Coefficient=N; type EntryMajor = ( usize, N ); type EntryMinor = ( usize, N );
}   

//  MAJOR ASCEND
//  ------------------------

impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >
    
    ViewRowAscend for 

    &'a CsMatBase
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    where
        N:                              Clone,
        I:                              SpIndex,
        Iptr:                           SpIndex,
        IptrStorage: Deref<Target   =   [Iptr]>,
        IndStorage: Deref<Target    =   [I]>,
        DataStorage: Deref<Target   =   [N]>,                 
{
    type ViewMajorAscend            =   SprsRowIteratorCloned< 'a, N, I >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymaj: usize ) -> Self::ViewMajorAscend {
        SprsRowIteratorCloned{ vector_iterator: self.outer_view( keymaj ).unwrap().into_sparse_vec_iter() }
    }
}


//  ENTRY
//  ------------------------

impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >
    
    MatrixEntry for 

    &'a CsMatBase
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    where
        N: Clone,
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,  
        N: Clone,               
{
    fn entry_major_at_minor( &self, keymaj: usize, keymin: usize ) -> Option< Self::Coefficient > {
        self.get( keymaj, keymin ).cloned()
    }
}


//  -------------------------------------------------------------------------------------
//  Arc< sprs::CsMatBase >
//  -------------------------------------------------------------------------------------


//  VIEW TYPE
//  ------------------------

/// A wrapper for the `sprs` struct `VectorIterator<'a, N, I >`, which returns
/// `(usize, N)` instead of `(usize, &N)`
pub struct VectorIteratorArc
            < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 
    where
        N: Clone,
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,  
{
    range:  Range< usize >,
    matrix: Arc< CsMatBase< N, I, IptrStorage, IndStorage, DataStorage, Iptr > >
}

impl < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    Iterator for 

    VectorIteratorArc
        < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

    where
        N: Clone,
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,      
{
    type Item = ( usize, N );
    fn next(&mut self) -> Option<Self::Item> {
        self.range.next().map( 
                |x| 
                ( self.matrix.indices()[ x ].index(), self.matrix.data()[ x ].clone() ) 
            )
    }
}  

//  INDICES AND COEFFICIENTS
//  ------------------------

impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >

    IndicesAndCoefficients for

    Arc< 
                CsMatBase 
                    < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 
            >

    where
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,     
{
    type EntryMajor = (usize,N);
    type EntryMinor = (usize,N);    
    type RowIndex=usize; 
    type ColIndex=usize; 
    type Coefficient=N;
}   

//  MAJOR ASCEND
//  ------------------------

impl < N, I, IptrStorage, IndStorage, DataStorage, Iptr >
    
    ViewRowAscend for 

    Arc< 
                CsMatBase 
                    < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 
            >

    where
        N:                              Clone,
        I:                              SpIndex,
        Iptr:                           SpIndex,
        IptrStorage: Deref<Target   =   [Iptr]>,
        IndStorage: Deref<Target    =   [I]>,
        DataStorage: Deref<Target   =   [N]>,                 
{
    type ViewMajorAscend            =   VectorIteratorArc
                                            < N, I, IptrStorage, IndStorage, DataStorage, Iptr >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( &self, keymaj: usize ) -> Self::ViewMajorAscend {
        let range = self.indptr().outer_inds_sz( keymaj );
        VectorIteratorArc{ matrix: self.clone(), range }
    }
}


//  ENTRY
//  ------------------------

impl < N, I, IptrStorage, IndStorage, DataStorage, Iptr >
    
    MatrixEntry for 

    Arc< 
                CsMatBase 
                    < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 
            >

    where
        N: Clone,
        I: SpIndex,
        Iptr: SpIndex,
        IptrStorage: Deref<Target = [Iptr]>,
        IndStorage: Deref<Target = [I]>,
        DataStorage: Deref<Target = [N]>,  
        N: Clone,               
{
    fn entry_major_at_minor( &self, keymaj: usize, keymin: usize ) -> Option< Self::Coefficient > {
        self.get( keymaj, keymin ).cloned()
    }
}

      


//  -------------------------------------------------------------------------------------
//  Conversion to sparsemat CSR matrices
//  -------------------------------------------------------------------------------------

/// Provides a `.into_csr()` method to convert objects into CSR matrices
pub trait IntoCSR< I, N > 
    where
        N: Clone + num::Num,
        I: SpIndex,
{
    fn into_csr( self, nrows: usize, ncols: usize ) -> CsMatBase< N, I, Vec<I>, Vec<I>, Vec<N> >;
}


//  Implement for Vec< Vec< N > > 
//  ---------------------------------------- 
impl < N >
    
    IntoCSR
        < usize, N > for
    
    Vec< Vec< N > > 

    where
        N: Clone + num::Num,

{
    /// Converts Vec< Vec< N > > into a CSR matrix 
    /// 
    /// If `v` has type `Vec< Vec< N > >` and `v[i][j]=z` then 
    /// `v.into_csr(p,q)[i,j]=z`.
    /// 
    /// Throws an error if `p` and `q` aren't large enough.
    fn into_csr( self, nrows: usize, ncols: usize ) 
        -> 
        CsMatBase<N, usize, Vec<usize>, Vec<usize>, Vec<N>>
    {
        let mut mat = TriMat::new( ( nrows, ncols ) );
        for ( rownum, row ) in self.into_iter().enumerate() {
            for ( colnum,val ) in row.into_iter().enumerate() {
                mat.add_triplet( rownum, colnum, val );
            }
        }
        mat.to_csr()
    }
}



//  Implement for Vec< Vec< (usize, N) > > 
//  ---------------------------------------- 
impl < N >
    
    IntoCSR
        < usize, N > for
    
    Vec< Vec< (usize, N) > > 

    where
        N: Clone + num::Num,

{
    /// Converts Vec< Vec< (usize, N) > > into a CSR matrix 
    /// 
    /// If `v` has type `Vec< Vec< (usize, N) > >` and `v[i]` contains
    /// `(j,z)` then `v.into_csr(p,q)[i,j]=z`.
    /// 
    /// 
    /// Throws an error if `p` and `q` aren't large enough.
    fn into_csr( self, nrows: usize, ncols: usize ) 
        -> 
        // CsMatBase<N, usize, Vec<_>, Vec<usize>, Vec<N>, _>
        // CsMatI< N, usize, usize >
        CsMatBase<N, usize, Vec<usize>, Vec<usize>, Vec<N>>
    {
        let mut mat = TriMat::new( ( nrows, ncols ) );
        for ( rownum, row ) in self.into_iter().enumerate() {
            for ( colnum,val ) in row {
                mat.add_triplet( rownum, colnum, val );
            }
        }
        mat.to_csr()
    }
}


//  Implement for Vec< (I, I, N) > 
//  ---------------------------------------- 
impl < I, N >
    
    IntoCSR
        < I,N > for
    
    Vec< (I, I, N) >

    where
        N: Clone + num::Num,
        I: SpIndex,

{
    /// Converts Vec< Vec< (I, I, N) > > into a CSR matrix
    /// 
    /// If `v` has type `Vec< (I, I, N) >` then for each
    /// `(i,j,z)` in `v` we have `v.into_csr(p,q)[i,j]=z`.  If `v` contains
    /// multiple triplets with the equal values for `i` and `j` then 
    /// one will be overwritten.
    /// 
    /// 
    /// Throws an error if `p` and `q` aren't large enough.
    /// 
    /// This is just a wrapper for the method `from_triplets` in the 
    /// `sprs` module.
    fn into_csr( self, nrows: usize, ncols: usize ) 
        -> 
        CsMatBase< N, I, Vec<I>, Vec<I>, Vec<N> >
    {
        let capacity = self.len();
        let mut iv = Vec::with_capacity( capacity );
        let mut jv = Vec::with_capacity( capacity );
        let mut zv = Vec::with_capacity( capacity );                
        for (i,j,z) in self {
            iv.push(i); jv.push(j); zv.push(z);
        }        
        TriMatBase::from_triplets( (nrows, ncols), iv, jv, zv ).to_csr()
    }
}



// let mut a = TriMat::new((4, 4));
// a.add_triplet(0, 0, 3.0_f64);
// a.add_triplet(1, 2, 2.0);
// a.add_triplet(3, 0, -2.0);

// // This matrix type does not allow computations, and must to
// // converted to a compatible sparse type, using for example
// let b: CsMat<_> = a.to_csr();




//  -------------------------------------------------------------------------------------
//  ndarray
//  -------------------------------------------------------------------------------------


// impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >

//     IndicesAndCoefficients for

//     &'a CsMatBase
//         < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

//     where
//         I: SpIndex,
//         Iptr: SpIndex,
//         IptrStorage: Deref<Target = [Iptr]>,
//         IndStorage: Deref<Target = [I]>,
//         DataStorage: Deref<Target = [N]>,     
// {
//     type RowIndex=usize; type ColIndex=usize; type Coefficient=&'a N;
// }   


// impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >
    
//     ViewRowAscend for 

//     &'a CsMatBase
//         < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

//     where
//         I: SpIndex,
//         Iptr: SpIndex,
//         IptrStorage: Deref<Target = [Iptr]>,
//         IndStorage: Deref<Target = [I]>,
//         DataStorage: Deref<Target = [N]>,                 
// {
//     type ViewMajorAscend            =   VectorIterator<'a, N, I >;
//     type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;
//     type EntryMajor =   ( usize, &'a N );

//     fn view_major_ascend( &self, keymaj: usize ) -> Self::ViewMajorAscend {
//         self.outer_view( keymaj ).unwrap().into_sparse_vec_iter()
//     }
// }


// impl < 'a, N, I, IptrStorage, IndStorage, DataStorage, Iptr >
    
//     MatrixEntry for 

//     &'a CsMatBase
//         < N, I, IptrStorage, IndStorage, DataStorage, Iptr > 

//     where
//         I: SpIndex,
//         Iptr: SpIndex,
//         IptrStorage: Deref<Target = [Iptr]>,
//         IndStorage: Deref<Target = [I]>,
//         DataStorage: Deref<Target = [N]>,  
//         N: Clone,               
// {
//     fn entry( &self, keymaj: usize, keymin: usize ) -> Option< Self::Coefficient > {
//         self.get( keymin, keymaj )
//     }
// }





#[cfg(test)]
mod tests {

    #[test]
    fn test_matrix_bimajor_data() {   

        use crate::algebra::matrices::query::{MatrixOracle};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; 


        // data for a matrix:
        // |  1   0 |
        // | -1   0 |
        // |  0   0 |
        // |  0   0 |           
        let matrix = vec![ 
                vec![ (0, 1.)  ], 
                vec![ (0,-1.)  ], 
                vec![          ], 
                vec![          ],                 
            ];
        
        // wrap the data in a VecOfVec sparse matrix struct
        let matrix = VecOfVec::new(matrix);

        // get some sparse columns
        let mut column_0        =   (& matrix ).column(0);
        let mut column_0_rev    =   (& matrix ).column_reverse(0);
        let mut column_1        =   (& matrix ).column(1);
        let mut column_1_rev    =   (& matrix ).column_reverse(1);

        // check the columns are correct
        assert_eq!( column_0.next(),        Some( (0,  1.) )    );
        assert_eq!( column_0.next(),        Some( (1, -1.) )    );
        assert_eq!( column_0.next(),        None                );    

        assert_eq!( column_0_rev.next(),    Some( (1, -1.) )    );
        assert_eq!( column_0_rev.next(),    Some( (0,  1.) )    );
        assert_eq!( column_0_rev.next(),    None                );   

        assert_eq!( column_1.next(),        None        );
        assert_eq!( column_1_rev.next(),    None        );                        
    }  





}

