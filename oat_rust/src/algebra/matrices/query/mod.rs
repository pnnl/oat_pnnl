//! Look up the (structural) nonzero entries in a row or column.
//! 
//! - [views](#views)
// //! - [non-integer index types](#non-integer-index-types)
//! - [entries](#entries)
//! - [rows and columns](#views-of-rows-and-columns)
// //! - [indexing](#indexing)
//! 
//! # Views
//! 
//! To talk about look-up operations, we use the language of views:
//! 
//! - a *major view* is either a row of a row-major matrix or a column of a column-major matrix
//! - a *minor view* is either a column of a row-major matrix or a row of a column-major matrix
//! 
//! Indices are named after views
//! 
//! - a *major key* is either a row index of a row-major matrix or a column index of a column-major matrix
//! - a *minor key* is either a column index of a row-major matrix or a row index of a column-major matrix
//! 
//! Most matrix algebra libraries use integers as indices.  But indices in OAT can be strings, floats, lists, or *anything else*!  You can also mix and match indices on rows and columns.  For example, you could define a matrix that indexes [rows with integers and columns with strings](crate::algebra::matrices#build-your-own-matrix).
//! 
//! 
//! # Entries
//! 
//! Here's an example of how to look up a single entry.  The terms "major" and "minor" are explained in [views](#views).  For a different example, see the [matrix module](crate::algebra::matrices#build-your-own-matrix).
//! 
//! ```
//! // Import some tools
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; // a matrix data structure
//! use oat_rust::algebra::matrices::query::MatrixEntry; // trait to retreive entries
//!  
//! // Create a sparse 2x2 upper triangular matrix.
//! let matrix_data =    VecOfVec::new(
//!                         vec![   vec![   (0, 0.),   (1, 1.)  ], 
//!                                 vec![              (1, 2.)  ] 
//!                             ], 
//!                      ); 
//! 
//! // For technical reasons, the matrix oracle traits aren't defined on `VecOfVec` directly.  
//! // Instead, the oracle traits are implemented on references to `VecOfVec`.  References
//! // in Rust are denoted by an `&` symbol.
//! let matrix  =   & matrix_data; 
//! 
//! // Get entries with the `.entry` method of the `MatrixEntry` trait
//! assert_eq!( matrix.entry_major_at_minor( 0, 0 ), Some( 0. ) );
//! assert_eq!( matrix.entry_major_at_minor( 0, 1 ), Some( 1. ) );
//! assert_eq!( matrix.entry_major_at_minor( 1, 0 ), None       );
//! assert_eq!( matrix.entry_major_at_minor( 1, 1 ), Some( 2. ) );
//! ```
//! 
//! 
//! 
//! # Rows and columns
//! 
//! 
//! Let's construct a matrix, then look up the entries in its rows and columns.  The terms "major" and "minor" are explained in [views](#views).  For a different example, see the [matrix module](crate::algebra::matrices#build-your-own-matrix).
//! 
//!    
//! ```
//! // Import some tools
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; // a matrix data structure
//! use oat_rust::algebra::matrices::query::{ViewRow, ViewCol}; // viewing traits
//!  
//! // Create a sparse 2x2 upper triangular matrix with 1's above the diagonal. OAT has many
//! // different pre-defined types for sparse matrices.  One example is `VecOfVec`.  
//! let matrix_data =    VecOfVec::new(
//!                         vec![   vec![   (0, 1.),   (1, 1.)  ], 
//!                                 vec![              (1, 1.)  ] 
//!                             ], 
//!                      );
//! 
//! // For technical reasons, the matrix oracle traits aren't defined on `VecOfVec` directly.  
//! // Instead, the oracle traits are implemented on references to `VecOfVec`.  References 
//! // in Rust are denoted by an `&` symbol.
//! let matrix  =   & matrix_data; 
//!  
//! // Look up the nonzero entries in row 0.
//! // To do so, we use the `view_major` method associated with the `ViewRow` trait.
//! let view    =   matrix.view_major( 0 ); 
//! itertools::assert_equal( view, vec![(0, 1.), (1, 1.)] );
//! 
//! // Look up the nonzero entries in column 1.
//! // To do so, we use the `view_minor` method associated with the `ViewCol` trait.
//! let view    =   matrix.view_minor( 1 ); 
//! itertools::assert_equal( view, vec![(0, 1.), (1, 1.)] );
//! ```
//! 
//! 
//! 
// //! # Views and major dimensions
// //! 
// //! To understand matrix oracles in OAT, you need just two basic concepts: major dimension
// //! and view.  These are easiest to understand in terms of an example:
// //! 
// //! **Example: Vector-of-Vectors**
// //! 
// //! Consider the following 2x2 matrix:
// //! 
// //! ```text
// //! 5 6
// //! 7 0
// //! ```
// //! 
// //! We can represent this matrix as a vector-of-vectors, where each internal 
// //! vector represents an ordered list of the nonzero entries in a given row
// //! (we call this type of representation *row-major*):
// //! 
// //! ```
// //! vec![
// //!     vec![ (0, 5), (1, 6) ],
// //!     vec![ (0, 7) ],
// //! ];
// //! ```
// //! 
// //! Alternatively, we can represent the matrix as a vector-of-vectors where each internal 
// //! vector represents an ordered list of the nonzero entries in a given column (we call
// //! this type of repreentation *column-major*):
// //! 
// //! ```
// //! vec![
// //!     vec![ (0, 5), (1, 7) ],
// //!     vec![ (0, 6) ],
// //! ];
// //! ```
// //! 
// //! If you want to reconstruct a full row of the matrix from one of these vector-of-vector representations,
// //! then you will probably use a different procedure than you would use to
// //! reconstruct a full column of the matrix.
// //! For example, given a *row-major* representation of a matrix,  
// //! - you can reconstruct row `i` by iterating over the entries in the `i`th internal vector
// //! - you can reconstruct column `i`, by searching through each and every internal vector,
// //! to check if that vector contains an entry of form `(i, *)`.
// //! 
// //! Many sparse matrix data structures share this characteristic: they use different procedures to reconstruct
// //! rows versus columns.  Typically, one of the two procedures will be significantly faster / more 
// //! efficient than the other; in the case of vector-of-vectors, for example, it's much faster to iterate
// //! over a single internal vector than it is to search across all of them, so reconstructing rows is
// //! easier in row-major vector-of-vector representations, while reconstructing columns is easier in 
// //! column-major vector-of-vector representations.
// //! 
// //! For the purpose of writing code, it's often irrelevant whether the internal vectors of a vector-of-vectors
// //! represent rows versus columns.  All that matters is that 
// //! - you have a data structure that can produce sparse vectors in one of two ways, and 
// //! - one way is much faster / more efficient than the other.
// //! 
// //! A sparse vector constructed via the efficient method is called a **major view**.  A sparse vector 
// //! constructed via the less efficient method is called a **minor view**.  A data structure that
// //! represents a matrix is **row-major** if major views represent rows, and **column-major** if 
// //! major-views represent columns.
// //! 
// //! That's it!  Now you understand about major views and row-major versus column-major reprsentations!



pub mod column_helper;





//  =======================================================


use std::fmt::Debug;
use std::iter::IntoIterator;

use derive_new::new;

use crate::algebra::rings::operator_traits::Semiring;
use crate::algebra::vectors::entries::KeyValGet;
use crate::utilities::order::JudgeOrder;
 // auto-implement a trait on references to objects that implement the trait



// /// An enum with two values: `Row` and `Col`.
// #[derive(Clone, Debug)]
// pub enum MajorDimension{
//     Row,
//     Col
// }

/// Access entries, rows, and columns of a matrix
/// 
/// 
/// The matrix [Oracle](crate::algebra::matrices::query::Oracle) trait provides a simple, consistent set of commands to retreive information about a matrix.
/// We designed this trait because OAT uses many different types of matrices, and we want to use the same set of commands for all of them.
/// 
/// - [Use this trait](how-to-use-this-trait)
/// - [Make this trait available for new types of objects](make-this-trait-available-for-new-types-of-objects)
/// 
// /// The trait itself has two main ingredients:
// /// 
// /// - a list of [lookup  methods](#look-up-methods) to access entries, rows, and columns
// /// - a list of [input/ouput types](#types) that define the inputs and outputs of the lookup functions
/// 
/// # Examples
/// 
/// #### Entries
/// 
/// To get the coefficient in a given row and column, use [`oracle.entry( row, column )`](crate::algebra::matrices::query::Oracle::entry)
/// 
/// #### Rows
/// 
/// To access a row, use one of the following 
/// 
/// - [`oracle.row(index)`](crate::algebra::matrices::query::Oracle::row): Returns a Rust iterable that runs over the entries in the row. Throws an error if the `index` is invalid.  See the [documentation](crate::algebra::matrices::query::Oracle::row) for full details.
/// - [`oracle.row_reverse(index)`](crate::algebra::matrices::query::Oracle::row_reverse): Reverses order of entries
/// - [`oracle.row_opt(index)`](crate::algebra::matrices::query::Oracle::row): Returns `None` if the index is invalid (no error)
/// - [`oracle.row_reverse_opt(index)`](crate::algebra::matrices::query::Oracle::row_reverse): Combines the `_opt` and `_reverse` methods
/// 
/// #### Columns
/// 
/// To access a row, use one of the following 
/// 
/// - [`oracle.column(index)`](crate::algebra::matrices::query::Oracle::row): throws an error if the index is invalid
/// - [`oracle.column_reverse(index)`](crate::algebra::matrices::query::Oracle::row_reverse): reverses order of entries
/// - [`oracle.column_opt(index)`](crate::algebra::matrices::query::Oracle::row): returns `None` if the index is invalid (no error)
/// - [`oracle.column_reverse_opt(index)`](crate::algebra::matrices::query::Oracle::row_reverse): combines the `_opt` and `_reverse` methods
/// 
/// 
/// # Make this trait available for new types of objects
/// 
/// There are a lot of different look-up functions, so there are lots of different input and output types.
/// The Oracle trait keeps a list of these types, which turns out to useful in some situations.
/// To help explain what each type means, we've added some "fake code" below.
/// Each line has format 
/// 
/// `oracle.method_name( InputType) -> OutputType`
/// 
/// Here are some examples:
/// 
/// ```ignore
/// oracle.entry( ColIndex, ColumnIndex     )   ->  Option< Coefficient >
/// 
/// oracle.row(                 ColIndex    )   ->  Row
/// oracle.row_opt(             ColIndex    )   ->  Option<Row>
/// oracle.row_reverse(         ColIndex    )   ->  RowReverse
/// oracle.row_reverse_opt(     ColIndex    )   ->  Option<RowReverse>
/// 
/// oracle.column(              ColumnIndex )   ->  Column
/// oracle.column_opt(          ColumnIndex )   ->  Option<Column>
/// oracle.column_reverse(      ColumnIndex )   ->  ColumnReverse
/// oracle.column_reverse_opt(  ColumnIndex )   ->  Option<ColumnReverse>
/// 
/// Row.into_iterator()                         ->  RowIter;
/// RowReverse.into_iterator()                  ->  RowReverseIter;
/// Column.into_iterator()                      ->  ColumnIter;
/// ColumnReverse.into_iterator()               ->  ColumnReverseIter;            
/// 
/// RowEntry.key()                              ->   ColumnIndex;
/// RowEntry.val()                              ->   Coefficient;
/// ColumnEntry.key()                           ->   ColIndex;
/// ColumnEntry.val()                           ->   Coefficient;
/// ```
/// 
///  - `Coefficient`: the entries of the matrix (integers/floats/strings/etc.)
///  - `ColIndex`: the X in the statment "I want row X"
///  - `ColumnIndex`: the X in the statement "I want column X"
///  - `RowEntry`: the X in the statement "the first entry in my row is X"
///  - `ColumnEntry`: the X in the statement "the first entry in my column is X"
///  - `Row`: the oracle gives you an object of this type when you call `oracle.row(index)`
///  -  RowIter: 
///  -  RowReverse;
///  -  RowReverseIter;
///  -  Column;
///  -  ColumnIter;
///  -  ColumnReverse;
///  -  ColumnReverseIter;
/// 
pub trait MatrixOracle{

    // type Row:               // What you get when you ask for a row.
    //                         IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowIter           >;
    // type RowIter:           // What you get when you call `row.into_iter()`, where `row` is a row
    //                         Iterator<       Item    =   Self::RowEntry                                              >;
    // type RowReverse:        // What you get when you ask for a row with the order of entries reversed
    //                         IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowReverseIter    >;
    // type RowReverseIter:    // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    //                         Iterator<       Item    =  Self::RowEntry                                               >;
    // type Column:            // What you get when you ask for a column
    //                         IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnIter        >;
    // type ColumnIter:        // What you get when you call `column.into_iter()`, where `column` is a column
    //                         Iterator<       Item    =   Self::ColumnEntry                                           >;
    // type ColumnReverse:     // What you get when you ask for a column with the order of entries reversed                             
    //                         IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnReverseIter >;  
    // type ColumnReverseIter: // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)
    //                         Iterator<       Item    =   Self::ColumnEntry                                           >;

    type Coefficient            ;    // The type of coefficient stored in each entry of the matrix    
    type RowIndex               ;    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            ;    // The type of column indices    
    type RowEntry               :   KeyValGet<  Self::ColumnIndex,  Self::Coefficient   >;  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            :   KeyValGet<  Self::RowIndex,     Self::Coefficient   >;  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    :   IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowIter           >;  // What you get when you ask for a row.
    type RowIter                :   Iterator<       Item    =   Self::RowEntry                                              >;  // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse             :   IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowReverseIter    >;  // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter         :   Iterator<       Item    =   Self::RowEntry                                              >;  // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    type Column                 :   IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnIter        >;  // What you get when you ask for a column
    type ColumnIter             :   Iterator<       Item    =   Self::ColumnEntry                                           >;  // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse          :   IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnReverseIter >;  // What you get when you ask for a column with the order of entries reversed                             
    type ColumnReverseIter      :   Iterator<       Item    =   Self::ColumnEntry                                           >;  // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    fn entry(                   &   self, row: Self::RowIndex, column: Self::ColumnIndex ) ->  Option< Self::Coefficient >;
    fn row(                     &   self, index: Self::RowIndex   )   -> Self::Row                 { self.row_opt(index).unwrap() }
    fn row_opt(                 &   self, index: Self::RowIndex   )   -> Option<Self::Row>;    
    fn row_reverse(             &   self, index: Self::RowIndex   )   -> Self::RowReverse          { self.row_reverse_opt(index).unwrap() }    
    fn row_reverse_opt(         &   self, index: Self::RowIndex   )   -> Option<Self::RowReverse>;    
    fn column(                  &   self, index: Self::ColumnIndex)   -> Self::Column              { self.column_opt(index).unwrap() }
    fn column_opt(              &   self, index: Self::ColumnIndex)   -> Option<Self::Column>;    
    fn column_reverse(          &   self, index: Self::ColumnIndex)   -> Self::ColumnReverse       { self.column_reverse_opt(index).unwrap() }            
    fn column_reverse_opt(      &   self, index: Self::ColumnIndex)   -> Option<Self::ColumnReverse>;    

} 



impl < 'a, T: MatrixOracle >
    
    MatrixOracle for 
    &'a T

{   
    type Coefficient        =   T::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   T::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   T::ColumnIndex;       // The type of column indices
    
    type RowEntry           =   T::RowEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   T::ColumnEntry;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   T::Row;               // What you get when you ask for a row.
    type RowIter                =   T::RowIter;           // What you get when you call `row.into_iter()`, where `row` is a row
    type RowReverse             =   T::RowReverse;        // What you get when you ask for a row with the order of entries reversed
    type RowReverseIter         =   T::RowReverseIter;    // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
    type Column                 =   T::Column;            // What you get when you ask for a column   
    type ColumnIter             =   T::ColumnIter;        // What you get when you call `column.into_iter()`, where `column` is a column
    type ColumnReverse          =   T::ColumnReverse;     // What you get when you ask for a column with the order of entries reversed                                
    type ColumnReverseIter      =   T::ColumnReverseIter; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

    fn entry(                   & self,  row: Self::RowIndex, column: Self::ColumnIndex )   ->  Option< Self::Coefficient >
        { (*self).entry( row, column ) }
    
    fn row(                     & self,  index: Self::RowIndex    )       -> Self::Row
        { (*self).row( index ) }    
    fn row_opt(                 & self,  index: Self::RowIndex    )   -> Option<Self::Row>
        { (*self).row_opt( index ) }        
    fn row_reverse(             & self,  index: Self::RowIndex    )       -> Self::RowReverse
        { (*self).row_reverse( index ) }
    fn row_reverse_opt(         & self,  index: Self::RowIndex    )   -> Option<Self::RowReverse>
        { (*self).row_reverse_opt( index ) }    
    
    fn column(                  & self,  index: Self::ColumnIndex )       -> Self::Column
        { (*self).column( index ) }
    fn column_opt(              & self,  index: Self::ColumnIndex )   -> Option<Self::Column>
        { (*self).column_opt( index ) }
    fn column_reverse(          & self,  index: Self::ColumnIndex )       -> Self::ColumnReverse
        { (*self).column_reverse( index ) }    
    fn column_reverse_opt(      & self,  index: Self::ColumnIndex )   -> Option<Self::ColumnReverse>
        { (*self).column_reverse_opt( index ) }        
} 


// pub trait Oracle{

//     /// The type of coefficient stored in each entry of the matrix
//     type Coefficient;   
//     /// The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
//     type KeyRow;
//     /// The type of column indices
//     type KeyColumn;    
//     /// The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
//     type EntryRow;
//     /// The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
//     type EntryColumn;
//     /// The type of object that represents a row, e.g. a list of key-value pairs
//     type Row;
//     /// The type of object that represents a row, e.g. a list of key-value pairs
//     type RowIncreasing;
//     type RowDecreasing;
//     /// The type of object that represents a column, e.g. a list of key-value pairs
//     type Column;
//     type ColumnIncreasing;
//     type ColumnDecreasing;

//     fn entry(               & self,  row: Self::KeyRow, column: Self::KeyColumn )   ->  Option< Self::Coefficient >;

//     fn row(                     & self,  index: Self::KeyRow    )       -> Self::Row                    { self.row_opt(index).unwrap() }
//     fn row_increasing(          & self,  index: Self::KeyRow    )       -> Self::RowIncreasing          { self.row_increasing_opt(index).unwrap() }    
//     fn row_decreasing(          & self,  index: Self::KeyRow    )       -> Self::RowDecreasing          { self.row_decreasing_opt(index).unwrap() }    
//     fn column(                  & self,  index: Self::KeyColumn )       -> Self::Column                 { self.column_opt(index).unwrap() }
//     fn column_increasing(       & self,  index: Self::KeyColumn )       -> Self::ColumnIncreasing       { self.column_increasing_opt(index).unwrap() }        
//     fn column_decreasing(       & self,  index: Self::KeyColumn )       -> Self::ColumnDecreasing       { self.column_decreasing_opt(index).unwrap() }        
//     fn row_opt(                 & self,  index: Self::KeyRow    )   -> Option<Self::Row>;
//     fn row_increasing_opt(      & self,  index: Self::KeyRow    )   -> Option<Self::RowIncreasing>;
//     fn row_decreasing_opt(      & self,  index: Self::KeyRow    )   -> Option<Self::RowDecreasing>;
//     fn column_opt(              & self,  index: Self::KeyColumn )   -> Option<Self::Column>;
//     fn column_increasing_opt(   & self,  index: Self::KeyColumn )   -> Option<Self::ColumnIncreasing>;
//     fn column_decreasing_opt(   & self,  index: Self::KeyColumn )   -> Option<Self::ColumnDecreasing>;    
// } 



//  =======================================================


/// Provides look-up commands for information relevant to matrix algebra
/// 
/// The `MatrixOracle` trait only allows the user to look up entries of a matrix.
/// This isn't enough to perform matrix algebra on it's own.  For example, if we
/// have two matrices with entries of type `usize`, and we want to multiply over
/// the integers modulo 4, then we have to tell the computer that we want to
/// work modulo 4.
/// 
/// The `MatrixAlgebra` trait provides a language to specify this information, via several look-up
/// commands: 
/// 
/// - `matrix.ring_operator()` returns an object called a called an [ring operator](crate::algebra::rings#terminology)
///   that performs algebraic operations on
///   the matrix coefficients
/// - `matrix.row_index_order()`  returns an
///   object called an [order operator](crate::utilities::order#terminology) that compares the order of row indices.
// ///   For example, if indices are weighted simplices then sometimes you might want to order
// ///   simplices lexicographically, and other times you might want to order them according to weight.
// ///   Whatever your preference, you can define a struct that returns the right value for
// ///   for `matrix.row_index_order()`. 
///   Similarly, `matrix.column_index_order()` returns [order operator](crate::utilities::order#terminology)
///   that compares the order of column indices.
/// - `matrix.row_entry_order()`  returns an [order operator](crate::utilities::order#terminology)
///   for row entries (which are different from row indices). The order on row entries should always agree with column 
///   index order, in the sense that `(i, x) < (j, y)` whenever `i < j`. **NB** because `i` and `j` are
///   column indices, we compare them with the `matrix.column_index_order()` order operator.
///   You might ask why we need and order operator for entries, if it is supposed to agree with the
///   order operator for column indices.  The answer is that sometimes
///   it can be more efficient to compare two entries directly than to first extract their indices,
///   and then compare those. The `matrix.column_entry_order()` function is similar, but for
///   column indices.
/// 
// /// compare first
// ///   dimension, then by weight, then by lexicographic order, or some other custom ordering.
// ///   At other times you might want to choose a different order. You can define different types
// ///   of matrix structs that return 
// /// - `matrix.row_index_order()` is similar to `matrix.row_index_order`, but for column indices
// /// - `row_entry_order`, and `column_entry_order`.  Each
// /// of these commands will return an object that says something about how to work
// /// with the matrix, algebraically.  The `ring_operator` object performs basic
// /// ring operations (addition and multiplication). The `row_entry_order` object
// /// can be used to compare the order in which different row entries should appear.
// /// 
// ///   The following pseudocode will not run, but it
// /// gives the general idea:
// /// 
// /// ```ignore
// /// // define a matrix
// /// let matrix              =   vec![
// ///                                 vec![ (0,0.), (1,1.) ],
// ///                                 vec![ (2,2.), (3,3.) ]
// ///                             ];
// /// 
// /// // perform some ring operations
// /// let ring_operator       =   matrix.ring_operator();
// /// let x                   =   ring_operator.add( 2, 2 );          // add 2 and 2
// /// let y                   =   ring_operator.multiply( 2, 2 );     // multiply 2 and 2
// /// 
// /// // perform some order operations
// /// let row_entry_order     =   matrix.row_entry.order();
// /// // this will assign z value of true, indicating that entry (1,1.) precedes entry (2,2.)
// /// let z                   =   row_entry_order.lt( (1,1.), (2,2.) );   
// /// // this will assign w value of false, indicating that entry (1,1.) doesn't precede entry (2,2.)
// /// let w                   =   row_entry_order.lt( (2,2.), (1,1.) );   
// /// ```
/// 
/// # Where is this used?
/// 
/// OAT uses the `MatrixAlgebra` look-up commands for many of its basic operations: 
/// multiplying matrices, computing dot-products, etc.
/// 
/// # I havea  matrix oracle but it doesn't implement this trait; what should I do?
/// 
/// You don't have to implement this trait.  You can just wrap your matrix oracle in
/// a `MatrixAlgebraPacket`!  
pub trait MatrixAlgebra: MatrixOracle {
    type RingOperator:      Semiring<   Self::Coefficient >;
    type RowEntryOrder:     JudgeOrder< Self::RowEntry >;
    type RowIndexOrder:     JudgeOrder< Self::RowIndex >;
    type ColumnEntryOrder:  JudgeOrder< Self::ColumnEntry >;
    type ColumnIndexOrder:  JudgeOrder< Self::ColumnIndex >;

    fn ring_operator( &self ) -> Self::RingOperator;
    fn row_entry_order( &self ) -> Self::RowEntryOrder;
    fn row_index_order( &self ) -> Self::RowIndexOrder;
    fn column_entry_order( &self ) -> Self::ColumnEntryOrder;
    fn column_index_order( &self ) -> Self::ColumnIndexOrder;
}


impl < 'a, T: MatrixAlgebra >
    
    MatrixAlgebra for 
    &'a T
{
    type RingOperator       =   T::RingOperator;
    type RowEntryOrder      =   T::RowEntryOrder;
    type RowIndexOrder      =   T::RowIndexOrder;
    type ColumnEntryOrder   =   T::ColumnEntryOrder;
    type ColumnIndexOrder   =   T::ColumnIndexOrder;

    fn ring_operator( &self ) -> Self::RingOperator { (*self).ring_operator() }
    fn row_entry_order( &self ) -> Self::RowEntryOrder { (*self).row_entry_order() }
    fn row_index_order( &self ) -> Self::RowIndexOrder { (*self).row_index_order() }    
    fn column_entry_order( &self ) -> Self::ColumnEntryOrder { (*self).column_entry_order() }
    fn column_index_order( &self ) -> Self::ColumnIndexOrder { (*self).column_index_order() }    
}








//  ---------------------------------------------------------------------------
//  MAJOR DIMENSION 
//  ---------------------------------------------------------------------------

// /// Specifies a major dimension (row or column).
// pub trait WhichMajor{ fn major_dimension( &self ) -> MajorDimension; }


//  ---------------------------------------------------------------------------
//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------

/// Specifies the major indices, minor indices, and coefficients of a sparse matrix.
/// 
/// - `Coefficient`: the coefficients of the matrix
/// - `RowIndex`: keys used to look up major views
/// - `ColIndex`: keys used to look u pminor views
/// - `EntryMajor`: an entry of a major view
/// - `EntryMinor`: an entry of a minor view
pub trait IndicesAndCoefficients {
    type Coefficient;   
    type RowIndex;
    type ColIndex;
    type EntryMajor;
    type EntryMinor;
    // type Kmaj;
    // type Kmin;
    // type Emaj;
    // type Emin;
    // type KMaj;
    // type KMin;
    // type EMaj;
    // type EMin;    
}

impl < 'a, T: IndicesAndCoefficients >
    
    IndicesAndCoefficients for 
    &'a T
{   
    type Coefficient=T::Coefficient;     
    type ColIndex=T::ColIndex; 
    type RowIndex=T::RowIndex; 
    type EntryMajor=T::EntryMajor; 
    type EntryMinor=T::EntryMinor;    
}    


// pub trait MatrixTypeParameters {
//     type ColIndex;
//     type RowIndex;
//     type Coefficient;    
//     type KeyValPairMin;
//     type KeyValPairMaj;
// }

// impl < 'a, T: MatrixTypeParameters >
    
//     MatrixTypeParameters for 
//     &'a T
// {   type ColIndex = T::ColIndex; type RowIndex = T::RowIndex; type Coefficient = T::Coefficient; type KeyValPairMaj = T::KeyValPairMaj; type KeyValPairMin = T::KeyValPairMin;   }    


//  ---------------------------------------------------------------------------
//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------

/// Look up an a single entry by row and column
pub trait MatrixEntry: IndicesAndCoefficients {

    /// Look up an entry by major key (first) and minor key (second)
    /// 
    /// Return the coefficient stored in position `keymin` of major view `keymaj`, or `None`, if no coefficient is stored there.
    /// 
    /// # Comments
    /// 
    /// The verbose method name `entry_major_at_minor` was chosen over the shorter `entry`, because the former
    /// makes the order of arguments more explicit.    
    fn entry_major_at_minor( &self, keymaj: Self::RowIndex, keymin: Self::ColIndex, ) -> Option< Self::Coefficient >;
}

impl < 'a, T: IndicesAndCoefficients + MatrixEntry >
    
    MatrixEntry for 
    &'a T
{   
    /// Look up an entry by major key (first) and minor key (second)
    /// 
    /// Return the coefficient stored in position `keymin` of major view `keymaj`, or `None`, if no coefficient is stored there.
    /// 
    /// # Comments
    /// 
    /// The verbose method name `entry_major_at_minor` was chosen over the shorter `entry`, because the former
    /// makes the order of arguments more explicit.
    fn entry_major_at_minor( &self, keymaj: Self::RowIndex, keymin: Self::ColIndex, ) -> Option< Self::Coefficient > { (*self).entry_major_at_minor(keymaj, keymin,) }    
}    


//  ---------------------------------------------------------------------------
//  ORACLE REFERENCE INHERIT
//  ---------------------------------------------------------------------------

// pub trait OracleRefInherit {}

// /// If `T` implements `ViewColDesecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: IndicesAndCoefficients > IndicesAndCoefficients for &'a T {
//     type ColIndex     =   T::ColIndex;
//     type RowIndex     =   T::RowIndex;        
//     type Coefficient     =   T::Coefficient;                
// }

// /// If `T` implements `OracleInherit`, then so does `&'a T`.
// impl < 'a, T: OracleRefInherit > OracleRefInherit for & 'a T {}

// /// If `T` implements `ViewRow` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + ViewRow< RowIndex, ViewMajor >, RowIndex, ViewMajor:IntoIterator > ViewRow< RowIndex, ViewMajor > for &'a T {
//     fn view_major(&self, index: Self::RowIndex) -> Self::ViewMajor {
//         (*self).view_major(index)
//     }
// }

// /// If `T` implements `ViewRowAsecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + ViewRowAscend< RowIndex, ViewMajorAscend >, RowIndex, ViewMajorAscend:IntoIterator > ViewRowAscend< RowIndex, ViewMajorAscend > for &'a T {
//     fn view_major_ascend(&self, index: Self::RowIndex) -> Self::ViewMajorAscend {
//         (*self).view_major_ascend(index)
//     }
// }

// /// If `T` implements `ViewRowDesecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + ViewRowDescend< RowIndex, ViewMajorDescend >, RowIndex, ViewMajorDescend:IntoIterator > ViewRowDescend< RowIndex, ViewMajorDescend > for &'a T {
//     fn view_major_descend(&self, index: Self::RowIndex) -> Self::ViewMajorDescend {
//         (*self).view_major_descend(index)
//     }
// }

// /// If `T` implements `ViewCol` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + ViewCol< RowIndex, ViewMinor >, RowIndex, ViewMinor > ViewCol< RowIndex, ViewMinor > for &'a T {
//     fn view_minor(&self, index: Self::RowIndex) -> Self::ViewMinor {
//         (*self).view_minor(index)
//     }
// }

// /// If `T` implements `ViewColAsecnd` and `OracleRefInherit` then so does `&'a T`.
// impl < 'a, T: OracleRefInherit + ViewColAscend< RowIndex, ViewMinorAscend >, RowIndex, ViewMinorAscend:IntoIterator > ViewColAscend< RowIndex, ViewMinorAscend > for &'a T {
//     fn view_minor_ascend(&self, index: Self::RowIndex) -> Self::ViewMinorAscend {
//         (*self).view_minor_ascend(index)
//     }
// }

// impl < 'a, T: IndicesAndCoefficients + ViewColDescend > ViewColDescend for &'a T {
//     type ViewMinorDescend       =   T::ViewMinorDescend;
//     type ViewMinorDescendIntoIter   =   T::ViewMinorDescendIntoIter;        
//     type EntryMinor  =   T::EntryMinor;                
//     fn view_minor_descend(&self, index: T::ColIndex) -> T::ViewMinorDescend {
//         (*self).view_minor_descend(index)
//     }
// }


//  ===========================================================================
//  MAJOR VIEWS
//  ===========================================================================


//  ---------------------------------------------------------------------------
//  UNORDERED
//  ---------------------------------------------------------------------------

//  ViewRowInherit IS DEPRECATED IN FAVOR OF OracleRefInherit
//  
// /// A trait with no methods, used as a criterion for auto-implementation of the ViewRow trait on references.
// pub trait ViewRowInherit {}

// /// If `T` implements `ViewRowInherit`, the so does `&'a T`.
// impl < 'a, T: ViewRowInherit > 

//     ViewRowInherit for 
    
//     & 'a T 
//     {}

/// Entries may appear in unsorted order (but the same order, every time).
// #[auto_impl(&)] 
pub trait ViewRow: IndicesAndCoefficients
{
    type ViewMajor:             IntoIterator< IntoIter = Self::ViewMajorIntoIter, Item = Self::EntryMajor >;
    type ViewMajorIntoIter:     Iterator< Item = Self::EntryMajor >;

    /// Get a major view.
    ///
    /// The order in which terms appear should be the same every time the
    /// function is called; however, the order need not be sorted.
    fn   view_major( &self, index: Self::RowIndex ) -> Self::ViewMajor;

    // ----------------
    // provided methods
    // ----------------

    /// Get an iterator of views indexed by a collection of `indices`
    fn   views_major< J >( self, indices: J ) 
        -> 
        ViewsMajor< Self, J::IntoIter > 
        where 
            J:              IntoIterator,
            J::IntoIter:    Iterator< Item=Self::RowIndex >,
            Self:   Sized     
        {
            ViewsMajor { matrix: self, indices: indices.into_iter() }
        }    
}

impl < 'a, T > 

    ViewRow for 
    
    &'a T

    where
        T:                  ViewRow + IndicesAndCoefficients,
        &'a T:              IndicesAndCoefficients< 
                                    EntryMajor=T::EntryMajor,        
                                    EntryMinor=T::EntryMinor,
                                    ColIndex=T::ColIndex, 
                                    RowIndex=T::RowIndex, 
                                    Coefficient=T::Coefficient 
                                >,        
    {
        type ViewMajor            =   T::ViewMajor;
        type ViewMajorIntoIter    =   T::ViewMajorIntoIter;              

        fn   view_major( &self, index: T::RowIndex ) -> T::ViewMajor { (*self).view_major( index ) }
    }    

/// Contains a matrix and an iterable; returns one view of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new)]
pub struct ViewsMajor
                < T, I > 
    where
        T:      IndicesAndCoefficients + ViewRow,
        I:      Iterator< Item = T::RowIndex >        
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    ViewsMajor
        < T, I > 

    where
        T:      IndicesAndCoefficients + ViewRow,
        I:      Iterator< Item = T::RowIndex >  
{
    type Item = T::ViewMajor;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.view_major(x) )
    }
}   

//  ---------------------------------------------------------------------------
//  ASCENDING ORDER
//  ---------------------------------------------------------------------------


/// Entries appear in **strictly** ascending order, according to index.
/// 
/// Consecutive entries must have distinct indices.
// #[auto_impl(&)] 
pub trait ViewRowAscend:    IndicesAndCoefficients
{
    type ViewMajorAscend: IntoIterator< IntoIter = Self::ViewMajorAscendIntoIter, Item = Self::EntryMajor >;
    type ViewMajorAscendIntoIter: Iterator< Item = Self::EntryMajor >;

    /// Get a major view with entries sorted in ascending order of index.
    fn   view_major_ascend( &self, index: Self::RowIndex ) -> Self::ViewMajorAscend;

    // ----------------
    // provided methods
    // ----------------

    /// Get an iterator of major views
    fn   views_major_ascend< J >( self, indices: J ) 
        -> 
        ViewsMajorAscend< Self, J::IntoIter > 
        where
            J:              IntoIterator,
            J::IntoIter:    Iterator< Item=Self::RowIndex >,
            Self:   Sized     
        {
            ViewsMajorAscend { matrix: self, indices: indices.into_iter() }
        }    
}

impl < 'a, T > 

    ViewRowAscend for 
    
    &'a T

    where
        T:                  ViewRowAscend + IndicesAndCoefficients,
        &'a T:              IndicesAndCoefficients< 
                                    EntryMajor=T::EntryMajor,        
                                    EntryMinor=T::EntryMinor,
                                    ColIndex=T::ColIndex, 
                                    RowIndex=T::RowIndex, 
                                    Coefficient=T::Coefficient 
                                >,

    {
        type ViewMajorAscend            =   T::ViewMajorAscend;
        type ViewMajorAscendIntoIter    =   T::ViewMajorAscendIntoIter;                

        fn   view_major_ascend( &self, index: T::RowIndex ) -> T::ViewMajorAscend { (*self).view_major_ascend( index ) }
    }

/// Contains a matrix and an iterable; returns one view of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new)]
pub struct ViewsMajorAscend
                < T, I > 
    where
        T:      IndicesAndCoefficients + ViewRowAscend,
        I:      Iterator< Item = T::RowIndex >        
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    ViewsMajorAscend
        < T, I > 

    where
        T:      IndicesAndCoefficients + ViewRowAscend,
        I:      Iterator< Item = T::RowIndex >  
{
    type Item = T::ViewMajorAscend;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.view_major_ascend(x) )
    }
}            


//  ---------------------------------------------------------------------------
//  DESCENDING ORDER
//  ---------------------------------------------------------------------------


/// Entries appear in **strictly** descending order, according to index.  
/// 
/// Consecutive entries must have distinct indices.
pub trait ViewRowDescend: IndicesAndCoefficients
{
    type ViewMajorDescend: IntoIterator< IntoIter = Self::ViewMajorDescendIntoIter, Item = Self::EntryMajor >;
    type ViewMajorDescendIntoIter: Iterator< Item = Self::EntryMajor >;
        
    /// Get a major view with entries sorted in descending order of index.
    fn   view_major_descend( &self, index: Self::RowIndex ) -> Self::ViewMajorDescend;

    // ----------------
    // provided methods
    // ----------------

    /// Get an iterator of views indexed by a collection of `indices`
    fn   views_major_descend< J >( self, indices: J ) 
        -> 
        ViewsMajorDescend< Self, J::IntoIter > 
        where 
            J:              IntoIterator,
            J::IntoIter:    Iterator< Item=Self::RowIndex >,
            Self:   Sized     
        {
            ViewsMajorDescend { matrix: self, indices: indices.into_iter() }
        }    
}

impl < 'a, T > 

    ViewRowDescend for 
    
    &'a T

    where
        T:                  ViewRowDescend + IndicesAndCoefficients,
        &'a T:              IndicesAndCoefficients< 
                                    EntryMajor=T::EntryMajor,        
                                    EntryMinor=T::EntryMinor,
                                    ColIndex=T::ColIndex, 
                                    RowIndex=T::RowIndex, 
                                    Coefficient=T::Coefficient 
                                >,

    {
        type ViewMajorDescend            =   T::ViewMajorDescend;
        type ViewMajorDescendIntoIter    =   T::ViewMajorDescendIntoIter;            

        fn   view_major_descend( &self, index: T::RowIndex ) -> T::ViewMajorDescend { (*self).view_major_descend( index ) }
    }

/// Contains a matrix and an iterable; returns one view of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new)]
pub struct ViewsMajorDescend
                < T, I > 
    where
        T:      IndicesAndCoefficients + ViewRowDescend,
        I:      Iterator< Item = T::RowIndex >        
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    ViewsMajorDescend
        < T, I > 

    where
        T:      IndicesAndCoefficients + ViewRowDescend,
        I:      Iterator< Item = T::RowIndex >  
{
    type Item = T::ViewMajorDescend;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.view_major_descend(x) )
    }
}   

// /// Entries appear in **non-strictly** ascending order, according to index.  
// /// 
// /// Consecutive entries may have identical indices.
// // #[auto_impl(&)] 
// pub trait ViewRowAscendLax: IndicesAndCoefficients
// {
//     type ViewMajorAscendLax: IntoIterator< IntoIter = Self::ViewMajorAscendLaxIntoIter, Item = Self::EntryMajor >;
//     type ViewMajorAscendLaxIntoIter: Iterator< Item = Self::EntryMajor >;

//     /// Get a major view with entries sorted in ascending order of index.
//     fn   view_major_ascend( &self, index: Self::RowIndex ) -> Self::ViewMajorAscendLax;
// }

// /// Entries appear in **non-strictly** descending order, according to index.   Consecutive entries may have identical indices.
// // #[auto_impl(&)] 
// pub trait ViewRowDescendLax: IndicesAndCoefficients
// {
//     type ViewMajorDescendLax: IntoIterator< IntoIter = Self::ViewMajorDescendLaxIntoIter, Item = Self::EntryMajor >;
//     type ViewMajorDescendLaxIntoIter: Iterator< Item = Self::EntryMajor >;

//     /// Get a major view with entries sorted in descending order of index.
//     fn   view_major_descend( &self, index: Self::RowIndex ) -> Self::ViewMajorDescendLax;
// }

// FOR FUTURE CONSIDERATION
// pub trait ViewRowAscendScoped < ColIndex, RowIndex, Coefficient>
// {
//     type PairMajorAscendScoped: KeyValGet < ColIndex, Coefficient, >;
//     type ViewMajorAscendScoped: IntoIterator < Item = PairMajorAscendScoped >;
//     /// Get a major view with entries sorted in ascending order of index, clipped to range [min,
//     /// max).
//     fn   view_major_ascend_scoped( &self, index: Self::RowIndex, min: ColIndex, max: ColIndex ) -> Self::ViewMajorAscendScoped;
// }


//  ===========================================================================
//  MINOR VIEWS
//  ===========================================================================


//  ---------------------------------------------------------------------------
//  UNORDERED
//  ---------------------------------------------------------------------------


/// Entries may not appear in sorted order.
// #[auto_impl(&)] 
pub trait ViewCol: IndicesAndCoefficients
{
    type ViewMinor: IntoIterator< IntoIter = Self::ViewMinorIntoIter, Item = Self::EntryMinor >;
    type ViewMinorIntoIter: Iterator< Item = Self::EntryMinor >;

    /// Get a minor view.
    ///
    /// The order in which terms appear should be the same every time the
    /// function is called; however, the order need not be sorted.
    fn   view_minor( &self, index: Self::ColIndex ) -> Self::ViewMinor;

    // ----------------
    // provided methods
    // ----------------

    /// Get an iterator of views indexed by a collection of `indices`
    fn   views_minor< J >( self, indices: J ) 
        -> 
        ViewsMinor< Self, J::IntoIter > 
        where 
            J:              IntoIterator,
            J::IntoIter:    Iterator< Item=Self::ColIndex >,
            Self:   Sized     
        {
            ViewsMinor { matrix: self, indices: indices.into_iter() }
        }    
}

impl < 'a, T > 

    ViewCol for 
    
    &'a T

    where
        T:                  ViewCol + IndicesAndCoefficients,
        &'a T:              IndicesAndCoefficients< 
                                    EntryMajor=T::EntryMajor,        
                                    EntryMinor=T::EntryMinor,
                                    ColIndex=T::ColIndex, 
                                    RowIndex=T::RowIndex, 
                                    Coefficient=T::Coefficient 
                                >,

    {
        type ViewMinor            =   T::ViewMinor;
        type ViewMinorIntoIter    =   T::ViewMinorIntoIter;              

        fn   view_minor( &self, index: T::ColIndex ) -> T::ViewMinor { (*self).view_minor( index ) }
    }

/// Contains a matrix and an iterable; returns one view of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new)]
pub struct ViewsMinor
                < T, I > 
    where
        T:      IndicesAndCoefficients + ViewCol,
        I:      Iterator< Item = T::ColIndex >
{
    /// A matrix oracle
    matrix:     T,
    /// A collection of indices
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    ViewsMinor
        < T, I > 

    where
        T:      IndicesAndCoefficients + ViewCol,
        I:      Iterator< Item = T::ColIndex >  
{
    type Item = T::ViewMinor;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.view_minor(x) )
    }
}   


//  ---------------------------------------------------------------------------
//  ASCENDING ORDER
//  ---------------------------------------------------------------------------

/// Entries appear in **strictly** ascending order, according to index.  
/// 
/// Consecutive entries have distinct indices.
// #[auto_impl(&)] 
pub trait ViewColAscend: IndicesAndCoefficients
{
    type ViewMinorAscend: IntoIterator< IntoIter = Self::ViewMinorAscendIntoIter, Item = Self::EntryMinor >;
    type ViewMinorAscendIntoIter: Iterator< Item = Self::EntryMinor >;

    /// Get a minor view with entries sorted in ascending order of index.
    fn   view_minor_ascend( &self, index: Self::ColIndex ) -> Self::ViewMinorAscend;

    // ----------------
    // provided methods
    // ----------------

    /// Get an iterator of views indexed by a collection of `indices`
    fn   views_minor_ascend< J >( self, indices: J ) 
        -> 
        ViewsMinorAscend< Self, J::IntoIter > 
        where 
            J:              IntoIterator,
            J::IntoIter:    Iterator< Item=Self::ColIndex >,
            Self:   Sized     
        {
            ViewsMinorAscend { matrix: self, indices: indices.into_iter() }
        }    
}

impl < 'a, T > 

    ViewColAscend for 
    
    &'a T

    where
        T:                  ViewColAscend + IndicesAndCoefficients,
        &'a T:              IndicesAndCoefficients< 
                                    EntryMajor=T::EntryMajor,        
                                    EntryMinor=T::EntryMinor,
                                    ColIndex=T::ColIndex, 
                                    RowIndex=T::RowIndex, 
                                    Coefficient=T::Coefficient 
                                >,

    {
        type ViewMinorAscend            =   T::ViewMinorAscend;
        type ViewMinorAscendIntoIter    =   T::ViewMinorAscendIntoIter;              

        fn   view_minor_ascend( &self, index: T::ColIndex ) -> T::ViewMinorAscend { (*self).view_minor_ascend( index ) }
    }

/// Contains a matrix and an iterable; returns one view of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new)]
pub struct ViewsMinorAscend
                < T, I > 
    where
        T:      IndicesAndCoefficients + ViewColAscend,
        I:      Iterator< Item = T::ColIndex >        
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    ViewsMinorAscend
        < T, I > 

    where
        T:      IndicesAndCoefficients + ViewColAscend,
        I:      Iterator< Item = T::ColIndex >  
{
    type Item = T::ViewMinorAscend;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.view_minor_ascend(x) )
    }
}   

//  ---------------------------------------------------------------------------
//  DESCENDING ORDER
//  ---------------------------------------------------------------------------

/// Entries appear in **strictly** descending order, according to index.  
/// 
/// Consecutive entries have distinct indices.
// #[auto_impl(&)] 
pub trait ViewColDescend: IndicesAndCoefficients
{
    type ViewMinorDescend: IntoIterator< IntoIter = Self::ViewMinorDescendIntoIter, Item = Self::EntryMinor >;
    type ViewMinorDescendIntoIter: Iterator< Item = Self::EntryMinor >;

    /// Get a minor view with entries sorted in descending order of index.
    fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend;

    // ----------------
    // provided methods
    // ----------------

    /// Get an iterator of views indexed by a collection of `indices`
    fn   views_minor_descend< J >( self, indices: J ) 
        -> 
        ViewsMinorDescend< Self, J::IntoIter > 
        where 
            J:              IntoIterator,
            J::IntoIter:    Iterator< Item=Self::ColIndex >,
            Self:   Sized     
        {
            ViewsMinorDescend { matrix: self, indices: indices.into_iter() }
        }    
}

impl < 'a, T > 

    ViewColDescend for 
    
    &'a T

    where
        T:                  ViewColDescend + IndicesAndCoefficients,
        &'a T:              IndicesAndCoefficients< 
                                    EntryMajor=T::EntryMajor,        
                                    EntryMinor=T::EntryMinor,
                                    ColIndex=T::ColIndex, 
                                    RowIndex=T::RowIndex, 
                                    Coefficient=T::Coefficient 
                                >,

    {
        type ViewMinorDescend            =   T::ViewMinorDescend;
        type ViewMinorDescendIntoIter    =   T::ViewMinorDescendIntoIter;                

        fn   view_minor_descend( &self, index: T::ColIndex ) -> T::ViewMinorDescend { (*self).view_minor_descend( index ) }
    }

/// Contains a matrix and an iterable; returns one view of the matrix for each index returned by the iterable.
#[derive(Clone,Copy,Debug,new)]
pub struct ViewsMinorDescend
                < T, I > 
    where
        T:      IndicesAndCoefficients + ViewColDescend,
        I:      Iterator< Item = T::ColIndex >        
{
    matrix:     T,
    indices:    I,
}        

impl < T, I >

    Iterator for 
    
    ViewsMinorDescend
        < T, I > 

    where
        T:      IndicesAndCoefficients + ViewColDescend,
        I:      Iterator< Item = T::ColIndex >  
{
    type Item = T::ViewMinorDescend;
    fn next(&mut self) -> Option<Self::Item> {
        self.indices.next().map(|x| self.matrix.view_minor_descend(x) )
    }
}   


//  ---------------------------------------------------------------------------
//  ASCENDING LAX ORDER
//  ---------------------------------------------------------------------------


// /// Entries appear in **non-strictly** ascending order, according to index.  
// /// 
// /// Consecutive entries may have identical indices.
// // #[auto_impl(&)] 
// pub trait ViewColAscendLax: IndicesAndCoefficients
// {
//     type ViewMinorAscendLax: IntoIterator< IntoIter = Self::ViewMinorAscendLaxIntoIter, Item = Self::EntryMinor >;
//     type ViewMinorAscendLaxIntoIter: Iterator< Item = Self::EntryMinor >;

//     /// Get a minor view with entries sorted in ascending order of index.
//     fn   view_minor_ascend( &self, index: Self::ColIndex ) -> Self::ViewMinorAscendLax;
// }

// //  ---------------------------------------------------------------------------
// //  DESCENDING LAX ORDER
// //  ---------------------------------------------------------------------------

// /// Entries appear in **non-strictly** descending order, according to index.   
// /// 
// /// Consecutive entries may have identical indices.
// // #[auto_impl(&)] 
// pub trait ViewColDescendLax: IndicesAndCoefficients
// {
//     type ViewMinorDescendLax: IntoIterator< IntoIter = Self::ViewMinorDescendLaxIntoIter, Item = Self::EntryMinor >;
//     type ViewMinorDescendLaxIntoIter: Iterator< Item = Self::EntryMinor >;

//     /// Get a minor view with entries sorted in descending order of index.
//     fn   view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescendLax;
// }




















//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    use crate::algebra::matrices::query::IndicesAndCoefficients;


    
    #[test] 
    fn test_trait_implementation_demos() {
        // import crates
        use crate::algebra::matrices::query::ViewRow;          

        
        //  -----------------------------------------------------------------------------------------
        //  implement ViewRow on a WRAPPER FOR &'a Vec< Vec< usize > > 
        //  -----------------------------------------------------------------------------------------

        // NB:  (1) the salient features of this segment are that 
        //          (a) it uses an immutable reference to a vector-of-vectors.  This is because the
        //              `lifetime` of the reference can help us to overcome some issues concerning
        //              lifetimes and major views.
        //          (b) it introduces and uses a *wrapper* struct that contains a `&'a Vec< Vec< usize > >` 
        //              This is beecause Rust's "orphan rule" prevents the implementation of a
        //              foreign trait on a foreign type (see the following segment for details)

        // This struct is just a wrapper around `&'a Vec< Vec< (usize,usize)`
        struct VecOfVecReferenceWrapper< 'a > { vec_of_vec: &'a Vec< Vec< (usize,usize) > > }

        impl < 'a > IndicesAndCoefficients for VecOfVecReferenceWrapper< 'a > { 
            type EntryMajor = &'a (usize,usize);
            type EntryMinor = &'a (usize,usize);                    
            type ColIndex = usize;  
            type RowIndex = usize;  
            type Coefficient = usize;
        }

        impl < 'a >

                ViewRow for 

                VecOfVecReferenceWrapper< 'a >

        {
            type ViewMajor          = &'a [(usize,usize)];
            type ViewMajorIntoIter  = std::slice::Iter< 'a, (usize,usize) >;

            fn view_major( & self, index: Self::RowIndex ) -> Self::ViewMajor {
                return self.vec_of_vec[index].as_slice()
            } 
        }

        // Get a major view
        let matrix = VecOfVecReferenceWrapper{ vec_of_vec: &vec![ vec![ (1,1), (2,2) ]  ] };
        let row = matrix.view_major( 0 );
        itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );

        // Get a major view from a *nested* immutable reference.
        let _row = (&(& matrix)).view_major( 0 );      


        //  -----------------------------------------------------------------------------------------
        //  implement ViewRow on &'a Vec< Vec< usize > > 
        //  -----------------------------------------------------------------------------------------


        // NB:  SO FAR AS WE ARE CURRENTLY AWARE it is not possible to implement ViewRow on 
        //      &'a Vec< Vec< usize > > **except when the code resides in the same file where this
        //      trait is defined.**  The reason for this is Rust's orphan rule 
        //      (c.f. https://doc.rust-lang.org/error-index.html#E0117), which prevents the
        //      implementation of a foreign trait on a foreign type (where "foreign" means a trait or
        //      type that is defined outside the file where you want to implement the trait).
        //      Therfore the following code compiles in the unit tests for the file
        //      src/oat/matrices/entry_lookup.rs, but it does not compile in doc string tests or
        //      unit tests defined in other files.

        impl < 'a > IndicesAndCoefficients for &'a Vec < Vec < (usize,usize) > > { 
            type EntryMajor = &'a ( usize, usize ); 
            type EntryMinor = &'a ( usize, usize ); 
            type ColIndex = usize;  
            type RowIndex = usize;  
            type Coefficient = usize;
        }

        impl < 'a >

                ViewRow for 

                &'a Vec < Vec < (usize,usize) > > 

        {
            type ViewMajor          =   &'a [(usize,usize)];
            type ViewMajorIntoIter  =   std::slice::Iter< 'a, (usize,usize) >;

            fn view_major( & self, index: Self::RowIndex ) -> Self::ViewMajor {
                return self[index].as_slice()
            } 
        }

        // Get a major view from an immutable reference.     
        let mut matrix = vec![ vec![ (1,1), (2,2) ]  ];
        let row = (& matrix).view_major( 0 );
        itertools::assert_equal( row, vec![ (1,1), (2,2) ].iter() );

        // Get a major view from a *nested* immutable reference.
        let _row = (&(& matrix)).view_major( 0 );      
        
        // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        let _row_a = (& matrix).view_major( 0 );
        matrix.push( vec![ (4,4), (5,5)] );
        let row_b = (& matrix).view_major( 1 );   

        itertools::assert_equal(row_b, vec![ (4,4), (5,5)].iter() );
        // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );        
        
        
        //  -----------------------------------------------------------------------------------------
        //  implement ViewRow on &'a mut Vec< Vec< usize > > 
        //  -----------------------------------------------------------------------------------------

        // NB:  SO FAR AS WE ARE CURRENTLY AWARE it is not possible to implement ViewRow on &'a mut Vec< Vec< usize > > 
        //      in a manner similar to the example above.  Indeed, the following code is a facsimile of the 
        //      section above, with `&'a mut` replaceing `&'a`.  It does not seem to compile, due to 
        //      a lifetime conflict.

        // impl < 'a, 'b, (usize,usize) >

        //         ViewRow
        //         <  usize,  &'a [(usize,usize)]  > for 

        //         &'a mut Vec < Vec < (usize,usize) > > 

        //         where   (usize,usize): Clone,
        //                 'b: 'a,
        // {
        //     fn view_major( & self, index: usize ) -> &'a [(usize,usize)] {
        //         return self[index].as_slice()
        //     } 
        // }

        // // Get a major view from an immutable reference.     
        // let mut matrix = vec![ vec![ 1, 2, 3 ]  ];
        // let row = (& matrix).view_major( 0 );
        // itertools::assert_equal( row, vec![ vec![ (1,1), (2,2) ]  ].iter() );

        // // Get a major view from a *nested* immutable reference.
        // let _row = (&(& matrix)).view_major( 0 );      
        
        // // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        // let _row_a = (& matrix).view_major( 0 );
        // matrix.push( vec![ (4,4), (5,5) ] );
        // let row_b = (& matrix).view_major( 0 );   

        // itertools::assert_equal(row_b.iter().cloned() , vec![ 1,2,3 ] );
        // // itertools::assert_equal( _row_a.iter().cloned() , vec![ 1,2,3] );           


        //  -----------------------------------------------------------------------------------------
        //  implement ViewRow on a struct with a lifetime generic parameter
        //  -----------------------------------------------------------------------------------------

        //  NB: SO FAR AS WE ARE AWARE it is not straightforward to implement `ViewRow` in a way
        //      that leverages a generic lifetime parameter associated with a struct.  This is in 
        //      contrast to some of the example shown below, which *do* leverage the lifetime
        //      parameter of a struct (it is key to note that different traits are implemented, and
        //      in the later examples the trait explicitly involves a lifetime parameter).  
        //      
        //  The following is a record of our "best attempt" at an implementation of `ViewRow` that
        //  does make use of a lifetime parameter.  The attempt failed, but others may have more
        //  success in the future.

        // use std::marker::PhantomData;   

        // // Define the struct
        // pub struct VecOfVecWithLifetime
        // < 'a, IndexCoeffPair >
        // {
        //     pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
        //     phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
        // }

        // // Implement the trait
        // impl < 'a >

        //         ViewRow
        //         <  usize, Cloned<Rev<std::slice::Iter<'a, (usize,usize)>>>  > for 

        //         VecOfVecWithLifetime< 'a, (usize,usize) >

        //         where   (usize,usize): Clone,
        // {
        //     fn view_major( & self, index: usize ) -> Cloned<Rev<std::slice::Iter<'a, (usize,usize)>>> {
        //         return self.vec_of_vec[index].iter().rev().cloned()
        //     } 
        // }


        //  -----------------------------------------------------------------------------------------
        //  -----------------------------------------------------------------------------------------
        //  Documentation of alternate oracle trait: with lifetimes
        //  -----------------------------------------------------------------------------------------
        //  -----------------------------------------------------------------------------------------        
        use std::marker::PhantomData;

        /// Entries may not appear in sorted order.
        pub trait ViewRowWithLifetime< 'a, MajKey, ViewMajor >
        {
            /// Get a major view.
            ///
            /// The order in which terms appear should be the same every time the
            /// function is called; however, the order need not be sorted.
            fn   view_major_with_life<'b: 'a>( &'b self, index: MajKey ) -> ViewMajor;
        }

        //  -----------------------------------------------------------------------------------------
        //  implement ViewRowWithLifetime on &'a mut Vec< Vec< (usize,usize) > >
        //  -----------------------------------------------------------------------------------------

        impl < 'a >

                ViewRowWithLifetime
                <  'a, usize,  &'a [(usize,usize)]  > for 

                &'a mut Vec < Vec < (usize,usize) > > 

                where   (usize,usize): Clone,
        {
            fn view_major_with_life<'b: 'a>( &'b self, index: usize ) -> &'a [(usize,usize)] {
                return self[index].as_slice()
            } 
        }

        // Define a mutable reference to a matrix
        let mut matrix_core = vec![ vec![ (1,1), (2,2) ]  ];
        let matrix_ref = &mut matrix_core;

        // Get a major view from the mutable reference.     
        let row = matrix_ref.view_major_with_life( 0 );
        itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );

        // Get a major view from a *nested* immutable reference.
        let _row = (& matrix_ref).view_major_with_life( 0 );      
        
        // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` 
        // below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        let _row_a = matrix_ref.view_major_with_life( 0 );
        matrix_ref.push( vec![ (4,4), (5,5) ] );
        let row_b = matrix_ref.view_major_with_life( 1 );   

        itertools::assert_equal(row_b.iter().cloned() , vec![ (4,4), (5,5) ] );
        // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] ); 
        
        
        //  -----------------------------------------------------------------------------------------
        //  implement ViewRowWithLifetime on VecOfVecWithLifetime< 'a, (usize,usize) >
        //  -----------------------------------------------------------------------------------------

        // DEFINE THE STRUCT
        pub struct VecOfVecWithLifetime
                    < 'a, IndexCoeffPair >
        
        {
            // pub major_dimension: MajorDimension, 
            pub vec_of_vec: Vec< Vec< IndexCoeffPair > >,
            pub phantom_kvpair: PhantomData< &'a IndexCoeffPair >,
        }

        // DEFINE ITS METHODS
        impl    < 'a, IndexCoeffPair >
                
                VecOfVecWithLifetime 
                < 'a, IndexCoeffPair > 

        {
            // Make new (empty) VecOfVec. 
            pub fn new( vecvec: Vec < Vec < IndexCoeffPair > > ) -> Self  
            {
                VecOfVecWithLifetime{   
                            // major_dimension: major_dimension,
                            vec_of_vec:     vecvec,                    
                            phantom_kvpair: PhantomData,                   
                        }
            }
        }

        // IMPLEMENT ViewRowWithLifetime
        impl < 'a >

                ViewRowWithLifetime
                <  'a, usize,  &'a Vec< (usize,usize) >  > for 

                VecOfVecWithLifetime< 'a, (usize,usize) >

                where   (usize,usize): Clone,
        {
            fn view_major_with_life<'b: 'a>( &'b self, index: usize ) -> &'a Vec< (usize,usize) > {
                &self.vec_of_vec[index]
            } 
        }        
        
        // Define a mutable reference to a matrix
        let mut matrix = VecOfVecWithLifetime::new( vec![ vec![ (1,1), (2,2) ]  ] );

        // Get a major view from the mutable reference.     
        let row = matrix.view_major_with_life( 0 );
        itertools::assert_equal( row.iter().cloned(), vec![ (1,1), (2,2) ] );

        // Get a major view from a *nested* immutable reference.
        let _row = matrix.view_major_with_life( 0 );      
        
        // NB: If you try to modify the vec-of-vec you may (depending on circumstances) run into lifetime errors
        // Example: uncommenting the line `itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );` below will result in an error at line `matrix.push( vec![ (4,4), (5,5) ] );`
        let _row_a = matrix.view_major_with_life( 0 );
        matrix.vec_of_vec.push( vec![ (4,4), (5,5) ] );
        let row_b = matrix.view_major_with_life( 1 );   

        itertools::assert_equal(row_b.iter().cloned() , vec![ (4,4), (5,5) ] );
        // itertools::assert_equal( _row_a.iter().cloned() , vec![ (1,1), (2,2) ] );         


    }

}        
