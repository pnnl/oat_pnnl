//! Wrapper for several pieces of data commonly used together with a matrix

use derive_getters::Dissolve;

use crate::{algebra::{matrices::{query::{ViewRowAscend, IndicesAndCoefficients, ViewColDescend, MatrixOracle, MatrixAlgebra}, operations::multiply::{vector_matrix_multiply_major_ascend_simplified, vector_matrix_multiply_minor_descend_simplified}}, vectors::{entries::{KeyValGet, KeyValSet}, operations::LinearCombinationSimplified}, rings::operator_traits::Semiring}, utilities::order::{ReverseOrder, JudgePartialOrder}};

use super::transpose::{Transpose, AntiTranspose};





/// Wrapper for several pieces of data commonly used together with a matrix
#[derive(Clone,Debug,Dissolve)]
pub struct MatrixAlgebraPacket
                < Matrix, RingOperator, OrderOperatorRowEntriesRight, OrderOperatorColumnEntriesLeft > 
{
    pub matrix:               Matrix,
    pub ring:                 RingOperator,
    pub row_entry_order:      OrderOperatorRowEntriesRight,
    pub col_entry_order:      OrderOperatorColumnEntriesLeft,
}

impl < Matrix, RingOperator, OrderOperatorRowEntriesRight, OrderOperatorColumnEntriesLeft >

    MatrixAlgebraPacket
        < Matrix, RingOperator, OrderOperatorRowEntriesRight, OrderOperatorColumnEntriesLeft >

    // where
    //     Matrix:                 Clone,
    //     RingOperator:           Clone,
    //     OrderOperatorRowEntriesRight:     Clone,
    //     OrderOperatorColumnEntriesLeft:     Clone,        
{
    /// A reference to the matrix 
    pub fn matrix_ref(&self) -> & Matrix { & self.matrix }
    /// A clone of the ring operator
    pub fn ring(&self) -> RingOperator where RingOperator: Clone { self.ring.clone() }    
    /// A clone of the order operator on row entries
    pub fn row_entry_order(&self) -> OrderOperatorRowEntriesRight where OrderOperatorRowEntriesRight: Clone { self.row_entry_order.clone() }        
    /// A clone of the order operator on column entries
    pub fn col_entry_order(&self) -> OrderOperatorColumnEntriesLeft where OrderOperatorColumnEntriesLeft: Clone { self.col_entry_order.clone() }    

    /// AntiTranspose the matrix
    pub fn antitranspose( self )
         ->     
        MatrixAlgebraPacket
                < AntiTranspose<Matrix>, RingOperator, ReverseOrder<OrderOperatorColumnEntriesLeft>, ReverseOrder<OrderOperatorRowEntriesRight>, >
    {
        MatrixAlgebraPacket { 
            matrix: AntiTranspose::new(self.matrix), 
            ring: self.ring, 
            row_entry_order: ReverseOrder::new( self.col_entry_order ), 
            col_entry_order: ReverseOrder::new( self.row_entry_order ),
        }
    }  

    /// Transpose the matrix
    pub fn transpose( self )
         ->     
        MatrixAlgebraPacket
                < Transpose<Matrix>, RingOperator, OrderOperatorColumnEntriesLeft, OrderOperatorRowEntriesRight, >
    {
        MatrixAlgebraPacket { 
            matrix: Transpose::new(self.matrix), 
            ring: self.ring, 
            row_entry_order: self.col_entry_order, 
            col_entry_order: self.row_entry_order,
        }
    }


    /// Left-multiplication with another matrix
    /// 
    /// Returns `self * other`
    pub fn multiply_left< MatrixRight, OrderOperatorRowEntries2, OrderOperatorMinor2 >( 
                self, 
                other:  MatrixAlgebraPacket<
                                MatrixRight,
                                RingOperator,
                                OrderOperatorRowEntries2,
                                OrderOperatorMinor2,                                
                            > 
            )
            ->
            ProductPacketDEPRECATEFORNEWORACLEVERSION< Matrix, MatrixRight, RingOperator, OrderOperatorRowEntries2, OrderOperatorColumnEntriesLeft, >
    {
        ProductPacketDEPRECATEFORNEWORACLEVERSION{ 
            matrix_left:        self.matrix,
            matrix_right:        other.matrix,
            ring:  self.ring,
            row_entry_order:    other.row_entry_order,
            col_entry_order:    self.col_entry_order,
        }
    }

    /// Right-multiplication with another matrix
    /// 
    /// Returns `other * self`
    pub fn multiply_right< MatrixRight, OrderOperatorRowEntries2, OrderOperatorMinor2 >( 
                self, 
                other:  MatrixAlgebraPacket<
                                MatrixRight,
                                RingOperator,
                                OrderOperatorRowEntries2,
                                OrderOperatorMinor2,                                
                            > 
            )
            ->
            ProductPacketDEPRECATEFORNEWORACLEVERSION< MatrixRight, Matrix, RingOperator, OrderOperatorRowEntriesRight, OrderOperatorMinor2, >
    {
        ProductPacketDEPRECATEFORNEWORACLEVERSION{ 
            matrix_left:        other.matrix,
            matrix_right:        self.matrix,
            ring:  self.ring,
            row_entry_order:    self.row_entry_order,
            col_entry_order:    other.col_entry_order,
        }
    }    




}





/// Product of two matrices
/// 
/// Unlike the [ProductMatrix](crate::algebra::matrices::operations::multiply::ProductMatrix) struct,
/// this one contains enough information to return minor views.
#[derive(Clone,Debug,Dissolve)]
pub struct ProductPacketDEPRECATEFORNEWORACLEVERSION< 
        MatrixLeft, 
        MatrixRight, 
        RingOperator,
        OrderOperatorRowEntriesRight,
        OrderOperatorColumnEntriesLeft,
    > 
{
    pub matrix_left:        MatrixLeft,
    pub matrix_right:       MatrixRight,
    pub ring:               RingOperator,
    pub row_entry_order:    OrderOperatorRowEntriesRight,
    pub col_entry_order:    OrderOperatorColumnEntriesLeft,    
    // packet1:    MatrixAlgebraPacket< MatrixLeft, RingOperator, OrderOperatorRowEntriesEntryMatrixL, OrderOperatorRowEntriesEntryMatrixR >,
    // packet2:    MatrixAlgebraPacket< MatrixRight, RingOperator, OrderOperatorRowEntriesEntryMatrixR, OrderOperatorMinorEntryMatrixR >,
}


//  MATRIX ORACLE
//  ---------------------------------------------------------------------------


// impl  < 
//         MatrixLeft, 
//         MatrixRight, 
//         RingOperator,
//         OrderOperatorRowEntriesRight,
//         OrderOperatorColumnEntriesLeft,
//     > 

//     MatrixOracle for
    
//     ProductPacketDEPRECATEFORNEWORACLEVERSION< 
//         MatrixLeft, 
//         MatrixRight, 
//         RingOperator,
//         OrderOperatorRowEntriesRight,
//         OrderOperatorColumnEntriesLeft,
//     >    
//     where
//         MatrixLeft:                    MatrixOracle,
//         MatrixRight:                    MatrixOracle< Coefficient = MatrixLeft::Coefficient, RowIndex = MatrixLeft::ColumnIndex >, 
//         MatrixRight::RowEntry:          KeyValSet < MatrixRight::ColumnIndex, MatrixRight::Coefficient >,  
//         MatrixRight::ColumnIndex:       Clone + PartialEq, // PartialEq is required by the struct that simplifies sparse vector iterators; it has to be able to compare the indices of different entries
//         MatrixRight::Coefficient:       Clone,                          
//         RingOperator:               Clone + Semiring< MatrixLeft::Coefficient >,
//         OrderOperatorRowEntriesRight:    Clone + JudgePartialOrder<  MatrixRight::RowEntry >,                 
// {   

//     type Coefficient            =   MatrixLeft::Coefficient;    // The type of coefficient stored in each entry of the matrix    
//     type RowIndex               =   MatrixLeft::RowIndex; // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
//     type ColumnIndex            =   MatrixRight::ColumnIndex; // The type of column indices    
//     type RowEntry               =   MatrixRight::RowEntry;  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
//     type ColumnEntry            =   MatrixLeft::ColumnEntry;  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
//     type Row                    =   LinearCombinationSimplified
//                                         < MatrixRight::RowIter, MatrixRight::ColumnIndex, MatrixRight::Coefficient, RingOperator, OrderOperatorRowEntriesRight >;
//     type RowIter                =   LinearCombinationSimplified
//                                         < MatrixRight::RowIter, MatrixRight::ColumnIndex, MatrixRight::Coefficient, RingOperator, OrderOperatorRowEntriesRight >;
//     type RowReverse             =   LinearCombinationSimplified
//                                         < MatrixRight::RowReverseIter, MatrixRight::ColumnIndex, MatrixRight::Coefficient, RingOperator, OrderOperatorRowEntriesRight >;
//     type RowReverseIter         =   LinearCombinationSimplified
//                                         < MatrixRight::RowReverseIter, MatrixRight::ColumnIndex, MatrixRight::Coefficient, RingOperator, OrderOperatorRowEntriesRight >;
//     type Column                 =   LinearCombinationSimplified
//                                         < MatrixLeft::ColumnIter, MatrixLeft::ColumnIndex, MatrixLeft::Coefficient, RingOperator, OrderOperatorColumnEntries >;
//     type ColumnIter             =   LinearCombinationSimplified
//                                         < MatrixLeft::ColumnIter, MatrixLeft::ColumnIndex, MatrixLeft::Coefficient, RingOperator, OrderOperatorColumnEntries >;
//     type ColumnReverse          =   LinearCombinationSimplified
//                                         < MatrixLeft::ColumnReverseIter, MatrixLeft::ColumnIndex, MatrixLeft::Coefficient, RingOperator, OrderOperatorColumnEntries >;
//     type ColumnReverseIter      =   LinearCombinationSimplified
//                                         < MatrixLeft::ColumnReverseIter, MatrixLeft::ColumnIndex, MatrixLeft::Coefficient, RingOperator, OrderOperatorColumnEntries >;


//     fn entry(                   &   self, row: Self::RowIndex, column: Self::ColumnIndex ) ->  Option< Self::Coefficient > {
//         println!("Add a test for this.");
//         let mut return_value    =   None;
//         for entry in self.matrix_left.row( row ) {
//             let scale       =   entry.val();
//             let row2        =   entry.key(); // pull out a row index
//             let _       =   self.matrix_right
//                                     .row( row2.clone() ) // look up the corresponding row
//                                     .into_iter()
//                                     .find(|x| x.key()==column ) // find an entry with the correct column
//                                     .map(   
//                                         |x| 
//                                         {
//                                             // multiply the entry in the correct column with the scalar by which we multiply the row
//                                             let term    =   self.ring.multiply( scale, x.val() );  
//                                             if let Some( sum ) = return_value {
//                                                 // if our running some already has a nonzero entry, add the new term
//                                                 return_value    =   Some(  self.ring.add( sum, term ) );
//                                             } else {
//                                                 // otherwise create the new term
//                                                 return_value    =   Some( term );
//                                             }
//                                         }
//                                     );
//             return_value
//         }
//     }
//     fn row(                     &   self, index: Self::RowIndex   )   -> Self::Row { 
//         vector_matrix_multiply_major_ascend_simplified( 
//             self.matrix_left.view_major_ascend( index ),
//             & self.matrix_right,
//             self.ring.clone(),
//             self.row_entry_order.clone(),
//         )
//     }
//     fn row_opt(                 &   self, index: Self::RowIndex   )   -> Option<Self::Row> {

//     }   
//     fn row_reverse(             &   self, index: Self::RowIndex   )   -> Self::RowReverse {
//         vector_matrix_multiply_major_ascend_simplified( 
//             self.matrix_left.view_major_ascend( index ),
//             self.matrix_right.reverse_ref(),
//             self.ring.clone(),
//             self.row_entry_order.clone().reverse(),
//         )        
//     }
//     fn row_reverse_opt(         &   self, index: Self::RowIndex   )   -> Option<Self::RowReverse>;    
//     fn column(                  &   self, index: Self::ColumnIndex)   -> Self::Column              { self.column_opt(index).unwrap() }
//     fn column_opt(              &   self, index: Self::ColumnIndex)   -> Option<Self::Column>;    
//     fn column_reverse(          &   self, index: Self::ColumnIndex)   -> Self::ColumnReverse       { self.column_reverse_opt(index).unwrap() }            
//     fn column_reverse_opt(      &   self, index: Self::ColumnIndex)   -> Option<Self::ColumnReverse>;    

// } 


//  INDICES AND COEFFICIENTS
//  ---------------------------------------------------------------------------
impl  < 
        MatrixLeft, 
        MatrixRight, 
        RingOperator,
        OrderOperatorRowEntriesRight,
        OrderOperatorColumnEntriesLeft,
    > 

    IndicesAndCoefficients for
    
    ProductPacketDEPRECATEFORNEWORACLEVERSION< 
        MatrixLeft, 
        MatrixRight, 
        RingOperator,
        OrderOperatorRowEntriesRight,
        OrderOperatorColumnEntriesLeft,
    >    
    where
        MatrixLeft:                            IndicesAndCoefficients,
        MatrixRight:                            IndicesAndCoefficients< Coefficient = MatrixLeft::Coefficient >,
{   
    type Coefficient = MatrixLeft::Coefficient;
    type EntryMajor = MatrixRight::EntryMajor;
    type EntryMinor = MatrixLeft::EntryMinor;
    type RowIndex = MatrixLeft::RowIndex;
    type ColIndex = MatrixRight::ColIndex;
}




// ViewRowAscend
impl  < 
        MatrixLeft, 
        MatrixRight, 
        RingOperator,
        OrderOperatorRowEntriesRight,
        OrderOperatorColumnEntriesLeft,
    > 

    ViewRowAscend for
    
    ProductPacketDEPRECATEFORNEWORACLEVERSION< 
        MatrixLeft, 
        MatrixRight, 
        RingOperator,
        OrderOperatorRowEntriesRight,
        OrderOperatorColumnEntriesLeft,      
    >    
    where
        MatrixLeft:                            ViewRowAscend + IndicesAndCoefficients,
        MatrixRight:                            ViewRowAscend + IndicesAndCoefficients< Coefficient = MatrixLeft::Coefficient, RowIndex = MatrixLeft::ColIndex >, 
        MatrixLeft::ViewMajorAscend:           IntoIterator,
        MatrixRight::ViewMajorAscend:           IntoIterator,
        MatrixLeft::EntryMajor:                KeyValGet < MatrixLeft::ColIndex, MatrixLeft::Coefficient >,
        MatrixRight::EntryMajor:                KeyValGet < MatrixRight::ColIndex, MatrixRight::Coefficient > + KeyValSet < MatrixRight::ColIndex, MatrixRight::Coefficient >,  
        MatrixRight::ColIndex:                  Clone + PartialEq, // PartialEq is required by the struct that simplifies sparse vector iterators; it has to be able to compare the indices of different entries
        MatrixRight::Coefficient:              Clone,                          
        RingOperator:                       Clone + Semiring< MatrixLeft::Coefficient >,
        OrderOperatorRowEntriesRight:             Clone + JudgePartialOrder<  MatrixRight::EntryMajor >,                 

{   
    type ViewMajorAscend            =   LinearCombinationSimplified
                                            < MatrixRight::ViewMajorAscendIntoIter, MatrixRight::ColIndex, MatrixRight::Coefficient, RingOperator, OrderOperatorRowEntriesRight >;
    type ViewMajorAscendIntoIter    =   Self::ViewMajorAscend;

    fn view_major_ascend( & self, index: Self::RowIndex ) 
        -> 
        LinearCombinationSimplified
            < MatrixRight::ViewMajorAscendIntoIter, MatrixRight::ColIndex, MatrixRight::Coefficient, RingOperator, OrderOperatorRowEntriesRight >
    {
        vector_matrix_multiply_major_ascend_simplified( 
                self.matrix_left.view_major_ascend( index ),
                & self.matrix_right,
                self.ring.clone(),
                self.row_entry_order.clone(),
            )   
    }
}


// ViewColDescend
impl  < 
        MatrixLeft, 
        MatrixRight, 
        RingOperator,
        OrderOperatorRowEntriesRight,
        OrderOperatorColumnEntriesLeft,
    > 

    ViewColDescend for
    
    ProductPacketDEPRECATEFORNEWORACLEVERSION< 
        MatrixLeft, 
        MatrixRight, 
        RingOperator,
        OrderOperatorRowEntriesRight,
        OrderOperatorColumnEntriesLeft,      
    >    
    where
        MatrixLeft:                            ViewColDescend + IndicesAndCoefficients,
        MatrixRight:                            ViewColDescend + IndicesAndCoefficients< Coefficient = MatrixLeft::Coefficient, RowIndex = MatrixLeft::ColIndex, >, 
        MatrixLeft::ViewMinorDescend:          IntoIterator,
        MatrixRight::ViewMinorDescend:          IntoIterator,
        MatrixLeft::EntryMinor:                KeyValGet < MatrixLeft::RowIndex, MatrixLeft::Coefficient > + KeyValSet < MatrixLeft::RowIndex, MatrixLeft::Coefficient >,
        MatrixRight::EntryMinor:                KeyValGet < MatrixRight::RowIndex, MatrixRight::Coefficient >,  
        MatrixLeft::RowIndex:                  Clone + PartialEq, // PartialEq is required by the struct that simplifies sparse vector iterators; it has to be able to compare the indices of different entries
        MatrixRight::Coefficient:              Clone,                          
        RingOperator:                       Clone + Semiring< MatrixLeft::Coefficient >,
        OrderOperatorColumnEntriesLeft:             Clone + JudgePartialOrder<  MatrixLeft::EntryMinor >,                 

{   
    type ViewMinorDescend           =   LinearCombinationSimplified< 
                                                MatrixLeft::ViewMinorDescendIntoIter, 
                                                MatrixLeft::RowIndex, 
                                                MatrixLeft::Coefficient, 
                                                RingOperator, 
                                                ReverseOrder< OrderOperatorColumnEntriesLeft >,
                                            >;
    type ViewMinorDescendIntoIter   =   Self::ViewMinorDescend;

    fn view_minor_descend( & self, index: MatrixRight::ColIndex ) 
        -> 
        LinearCombinationSimplified< 
                MatrixLeft::ViewMinorDescendIntoIter, 
                MatrixLeft::RowIndex, 
                MatrixLeft::Coefficient, 
                RingOperator, 
                ReverseOrder< OrderOperatorColumnEntriesLeft >,
            >
    {
        vector_matrix_multiply_minor_descend_simplified( 
                self.matrix_right.view_minor_descend( index ),
                & self.matrix_left,
                self.ring.clone(),
                self.col_entry_order.clone(),
            )   
    }
}





// /// Product of two matrices
// /// 
// /// Unlike the [ProductMatrix](crate::algebra::matrices::operations::multiply::ProductMatrix) struct,
// /// this one contains enough information to return minor views.
// #[derive(Clone,Debug,Dissolve)]
// pub struct ProductPacketDEPRECATEFORNEWORACLEVERSION< 
//         MatrixLeft,
//         MatrixRight,
//         RingOperator,
//         OrderOperatorRowEntriesLeft,
//         OrderOperatorRowEntriesRight,
//         OrderOperatorColumnEntriesLeft,
//         OrderOperatorColumnEntriesLeft,                  
//     > 
// {
//     pub matrix_left:        MatrixLeft,
//     pub matrix_right:       MatrixRight,
//     pub ring:               RingOperator,
//     pub row_entry_order:    OrderOperatorRowEntriesRight,
//     pub col_entry_order:    OrderOperatorColumnEntriesLeft,    
//     // packet1:    MatrixAlgebraPacket< MatrixLeft, RingOperator, OrderOperatorRowEntriesEntryMatrixL, OrderOperatorRowEntriesEntryMatrixR >,
//     // packet2:    MatrixAlgebraPacket< MatrixRight, RingOperator, OrderOperatorRowEntriesEntryMatrixR, OrderOperatorMinorEntryMatrixR >,
// }


// // Matrix Oracle
// impl  < 
//         MatrixLeft, 
//         MatrixRight, 
//         RingOperator,
//         OrderOperatorColumnEntriesLeft,
//         OrderOperatorRowEntriesRight,        
//     > 

//     MatrixOracle for
    
//     ProductPacket< 
//         MatrixLeft, 
//         MatrixRight, 
//         RingOperator,
//         OrderOperatorRowEntriesLeft,         
//         OrderOperatorColumnEntriesLeft,      
//         OrderOperatorRowEntriesRight,                
//     >    
//     where
//         MatrixLeft:                         MatrixOracle,
//         MatrixRight:                        MatrixOracle< 
//                                                 Coefficient     =   MatrixLeft::Coefficient, 
//                                                 RowIndex        =   MatrixLeft::ColumnIndex, 
//                                             >, 
//         // MatrixLeft::EntryMinor:                KeyValGet < MatrixLeft::RowIndex, MatrixLeft::Coefficient > + KeyValSet < MatrixLeft::RowIndex, MatrixLeft::Coefficient >,
//         // MatrixRight::EntryMinor:                KeyValGet < MatrixRight::RowIndex, MatrixRight::Coefficient >,  
//         MatrixLeft::RowIndex:               Clone + PartialEq, // PartialEq is required by the struct that simplifies sparse vector iterators; it has to be able to compare the indices of different entries
//         MatrixLeft::Coefficient:            Clone,                          
//         RingOperator:                       Clone + Semiring< MatrixLeft::Coefficient >,
//         OrderOperatorColumnEntriesLeft:     Clone + JudgePartialOrder<  MatrixLeft::ColumnEntry >,

// {   

//     // type Row:               // What you get when you ask for a row.
//     //                         IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowIter           >;
//     // type RowIter:           // What you get when you call `row.into_iter()`, where `row` is a row
//     //                         Iterator<       Item    =   Self::RowEntry                                              >;
//     // type RowReverse:        // What you get when you ask for a row with the order of entries reversed
//     //                         IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowReverseIter    >;
//     // type RowReverseIter:    // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
//     //                         Iterator<       Item    =  Self::RowEntry                                               >;
//     // type Column:            // What you get when you ask for a column
//     //                         IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnIter        >;
//     // type ColumnIter:        // What you get when you call `column.into_iter()`, where `column` is a column
//     //                         Iterator<       Item    =   Self::ColumnEntry                                           >;
//     // type ColumnReverse:     // What you get when you ask for a column with the order of entries reversed                             
//     //                         IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnReverseIter >;  
//     // type ColumnReverseIter: // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)
//     //                         Iterator<       Item    =   Self::ColumnEntry                                           >;

//     type Coefficient            =   MatrixLeft::Coefficient;// The type of coefficient stored in each entry of the matrix    
//     type RowIndex               =   MatrixLeft::RowIndex;// The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
//     type RowEntry               =   MatrixRight::RowEntry;// The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
//     type ColumnIndex            =   MatrixRight::ColumnIndex;// The type of column indices    
//     type ColumnEntry            =   MatrixLeft::ColEntry;// The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
//     type Row                    =   LinearCombinationSimplified< 
//                                         MatrixRight::Row, 
//                                         MatrixRight::ColumnIndex, 
//                                         MatrixRight::Coefficient, 
//                                         RingOperator, 
//                                         OrderOperatorRowEntriesRight,
//                                     >; // What you get when you ask for a row.
//     type RowIter                =   Self::Row; // What you get when you call `row.into_iter()`, where `row` is a row
//     type RowReverse             =   LinearCombinationSimplified< 
//                                         MatrixRight::RowReverse, 
//                                         MatrixRight::ColumnIndex, 
//                                         MatrixRight::Coefficient, 
//                                         RingOperator, 
//                                         ReverseOrder< OrderOperatorRowEntriesRight >,
//                                     >; // What you get when you ask for a row with the order of entries reversed
//     type RowReverseIter         =   Self::RowReverse; // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
//     type Column                 =   LinearCombinationSimplified< 
//                                         MatrixLeft::Column, 
//                                         MatrixLeft::RowIndex, 
//                                         MatrixLeft::Coefficient, 
//                                         RingOperator, 
//                                         OrderOperatorColumnEntriesLeft,
//                                     >; // What you get when you ask for a column
//     type ColumnIter             =   Self::Column; // What you get when you call `column.into_iter()`, where `column` is a column
//     type ColumnReverse          =   LinearCombinationSimplified< 
//                                         MatrixLeft::ColumnReverse, 
//                                         MatrixLeft::RowIndex, 
//                                         MatrixLeft::Coefficient, 
//                                         RingOperator, 
//                                         ReverseOrder< OrderOperatorColumnEntriesLeft >,
//                                     >;// What you get when you ask for a column with the order of entries reversed                             
//     type ColumnReverseIter      =   Self::ColumnReverse; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

//     fn entry(                   &   self, row: Self::RowIndex, column: Self::ColumnIndex ) ->  Option< Self::Coefficient >;
//     fn row(                     &   self, index: Self::RowIndex   )   -> Self::Row                 { self.row_opt(index).unwrap() }
//     fn row_opt(                 &   self, index: Self::RowIndex   )   -> Option<Self::Row>;    
//     fn row_reverse(             &   self, index: Self::RowIndex   )   -> Self::RowReverse          { self.row_reverse_opt(index).unwrap() }    
//     fn row_reverse_opt(         &   self, index: Self::RowIndex   )   -> Option<Self::RowReverse>;    
//     fn column(                  &   self, index: Self::ColumnIndex)   -> Self::Column              { self.column_opt(index).unwrap() }
//     fn column_opt(              &   self, index: Self::ColumnIndex)   -> Option<Self::Column>;    
//     fn column_reverse(          &   self, index: Self::ColumnIndex)   -> Self::ColumnReverse       { self.column_reverse_opt(index).unwrap() }            
//     fn column_reverse_opt(      &   self, index: Self::ColumnIndex)   -> Option<Self::ColumnReverse>;    

// } 



// /// Product of two matrices
// /// 
// /// Unlike the [ProductMatrix](crate::algebra::matrices::operations::multiply::ProductMatrix) struct,
// /// this one contains enough information to return minor views.
// #[derive(Clone,Debug,Dissolve)]
// pub struct ProductPacketDEPRECATEFORNEWORACLEVERSION< 
//         MatrixLeft,
//         MatrixRight,
//         RingOperator,
//         OrderOperatorRowEntriesLeft,
//         OrderOperatorRowEntriesRight,
//         OrderOperatorColumnEntriesLeft,
//         OrderOperatorColumnEntriesLeft,                  
//     > 
// {
//     pub matrix_left:        MatrixLeft,
//     pub matrix_right:       MatrixRight,
//     pub ring:               RingOperator,
//     pub row_entry_order:    OrderOperatorRowEntriesRight,
//     pub col_entry_order:    OrderOperatorColumnEntriesLeft,    
//     // packet1:    MatrixAlgebraPacket< MatrixLeft, RingOperator, OrderOperatorRowEntriesEntryMatrixL, OrderOperatorRowEntriesEntryMatrixR >,
//     // packet2:    MatrixAlgebraPacket< MatrixRight, RingOperator, OrderOperatorRowEntriesEntryMatrixR, OrderOperatorMinorEntryMatrixR >,
// }




/// Product of two matrices
/// 
/// Unlike the [ProductMatrix](crate::algebra::matrices::operations::multiply::ProductMatrix) struct,
/// this one contains enough information to return minor views.
#[derive(Clone,Debug,Dissolve)]
pub struct ProductMatrix< 
        MatrixLeft,
        MatrixRight,                 
    > 
{
    pub matrix_left:        MatrixLeft,
    pub matrix_right:       MatrixRight,
}


// // Matrix Oracle
// impl  < 
//         MatrixLeft, 
//         MatrixRight,       
//     > 

//     MatrixOracle for
    
//     ProductPacket< 
//         MatrixLeft, 
//         MatrixRight,               
//     >    
//     where
//         MatrixLeft:                         MatrixOracle + MatrixAlgebra,
//         MatrixRight:                        MatrixOracle< 
//                                                 Coefficient         =   MatrixLeft::Coefficient, 
//                                                 RowIndex            =   MatrixLeft::ColumnIndex,                                              
//                                             > +
//                                             MatrixAlgebra<
//                                                 RingOperator        =   MatrixLeft::RingOperator,
//                                                 ColumnEntryOrder    =   MatrixLeft::RowEntryOrder,
//                                             >, 
//         // MatrixLeft::EntryMinor:                KeyValGet < MatrixLeft::RowIndex, MatrixLeft::Coefficient > + KeyValSet < MatrixLeft::RowIndex, MatrixLeft::Coefficient >,
//         // MatrixRight::EntryMinor:                KeyValGet < MatrixRight::RowIndex, MatrixRight::Coefficient >,  
//         MatrixLeft::RowIndex:               Clone + PartialEq, // PartialEq is required by the struct that simplifies sparse vector iterators; it has to be able to compare the indices of different entries
//         MatrixLeft::Coefficient:            Clone,                          
//         RingOperator:                       Clone + Semiring< MatrixLeft::Coefficient >,
//         OrderOperatorColumnEntriesLeft:     Clone + JudgePartialOrder<  MatrixLeft::ColumnEntry >,

// {   

//     // type Row:               // What you get when you ask for a row.
//     //                         IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowIter           >;
//     // type RowIter:           // What you get when you call `row.into_iter()`, where `row` is a row
//     //                         Iterator<       Item    =   Self::RowEntry                                              >;
//     // type RowReverse:        // What you get when you ask for a row with the order of entries reversed
//     //                         IntoIterator<   Item    =   Self::RowEntry,     IntoIter    =   Self::RowReverseIter    >;
//     // type RowReverseIter:    // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
//     //                         Iterator<       Item    =  Self::RowEntry                                               >;
//     // type Column:            // What you get when you ask for a column
//     //                         IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnIter        >;
//     // type ColumnIter:        // What you get when you call `column.into_iter()`, where `column` is a column
//     //                         Iterator<       Item    =   Self::ColumnEntry                                           >;
//     // type ColumnReverse:     // What you get when you ask for a column with the order of entries reversed                             
//     //                         IntoIterator<   Item    =   Self::ColumnEntry,  IntoIter    =   Self::ColumnReverseIter >;  
//     // type ColumnReverseIter: // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)
//     //                         Iterator<       Item    =   Self::ColumnEntry                                           >;

//     type Coefficient            =   MatrixLeft::Coefficient;// The type of coefficient stored in each entry of the matrix    
//     type RowIndex               =   MatrixLeft::RowIndex;// The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
//     type RowEntry               =   MatrixRight::RowEntry;// The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
//     type ColumnIndex            =   MatrixRight::ColumnIndex;// The type of column indices    
//     type ColumnEntry            =   MatrixLeft::ColEntry;// The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
//     type Row                    =   LinearCombinationSimplified< 
//                                         MatrixRight::Row, 
//                                         MatrixRight::ColumnIndex, 
//                                         MatrixRight::Coefficient, 
//                                         RingOperator, 
//                                         OrderOperatorRowEntriesRight,
//                                     >; // What you get when you ask for a row.
//     type RowIter                =   Self::Row; // What you get when you call `row.into_iter()`, where `row` is a row
//     type RowReverse             =   LinearCombinationSimplified< 
//                                         MatrixRight::RowReverse, 
//                                         MatrixRight::ColumnIndex, 
//                                         MatrixRight::Coefficient, 
//                                         RingOperator, 
//                                         ReverseOrder< OrderOperatorRowEntriesRight >,
//                                     >; // What you get when you ask for a row with the order of entries reversed
//     type RowReverseIter         =   Self::RowReverse; // What you get when you call `row_reverse.into_iter()`, where `row_reverse` is a reversed row (which is a row with order of entries reversed)
//     type Column                 =   LinearCombinationSimplified< 
//                                         MatrixLeft::Column, 
//                                         MatrixLeft::RowIndex, 
//                                         MatrixLeft::Coefficient, 
//                                         RingOperator, 
//                                         OrderOperatorColumnEntriesLeft,
//                                     >; // What you get when you ask for a column
//     type ColumnIter             =   Self::Column; // What you get when you call `column.into_iter()`, where `column` is a column
//     type ColumnReverse          =   LinearCombinationSimplified< 
//                                         MatrixLeft::ColumnReverse, 
//                                         MatrixLeft::RowIndex, 
//                                         MatrixLeft::Coefficient, 
//                                         RingOperator, 
//                                         ReverseOrder< OrderOperatorColumnEntriesLeft >,
//                                     >;// What you get when you ask for a column with the order of entries reversed                             
//     type ColumnReverseIter      =   Self::ColumnReverse; // What you get when you call `column_reverse.into_iter()`, where `column_reverse` is a reversed column (which is a column with order of entries reversed)

//     fn entry(                   &   self, row: Self::RowIndex, column: Self::ColumnIndex ) ->  Option< Self::Coefficient >;
//     fn row(                     &   self, index: Self::RowIndex   )   -> Self::Row                 { self.row_opt(index).unwrap() }
//     fn row_opt(                 &   self, index: Self::RowIndex   )   -> Option<Self::Row>;    
//     fn row_reverse(             &   self, index: Self::RowIndex   )   -> Self::RowReverse          { self.row_reverse_opt(index).unwrap() }    
//     fn row_reverse_opt(         &   self, index: Self::RowIndex   )   -> Option<Self::RowReverse>;    
//     fn column(                  &   self, index: Self::ColumnIndex)   -> Self::Column              { self.column_opt(index).unwrap() }
//     fn column_opt(              &   self, index: Self::ColumnIndex)   -> Option<Self::Column>;    
//     fn column_reverse(          &   self, index: Self::ColumnIndex)   -> Self::ColumnReverse       { self.column_reverse_opt(index).unwrap() }            
//     fn column_reverse_opt(      &   self, index: Self::ColumnIndex)   -> Option<Self::ColumnReverse>;    

// } 