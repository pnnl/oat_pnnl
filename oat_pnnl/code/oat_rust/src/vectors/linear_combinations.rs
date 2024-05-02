//! Structs that represent  lienear combinations of other vectors, in a **memory efficient manner**.
//!
//! Suppose you have two iterators that each representa  sparse vector: `vec_a` and `vec_b`.
//! How do you represent `vec_a + vec_b`?
//! 
//! One option is to place all the entries from `vec_a` and `vec_b` in to a Rust vector, then combine entries
//! that have equal indices.  However, this may require a great deal of memory if `vec_a` and `vec_b` are long.
//! 
//! Instead, we could create a new iterator that only extracts entries from `vec_a` and `vec_b` when needed.
//! There is an object in the `Itertools` library that does something very similar: ['merge'](https://docs.rs/itertools/0.7.2/itertools/trait.Itertools.html#method.merge).
//! The `LinearCombinationUnsimplified` and `LinearCombinationSimplified` structs do essentially this.

use crate::{utilities::{iterators::merge::heap_of_iterators::HitMerge, partial_order::StrictlyLess}, rings::operator_traits::Semiring, entries::{KeyValSet, KeyValGet}};

use super::operations::{Scale, Simplify};





//  UNSIMPLIFIED
//  -------------------------------------------------------------------------------------

/// A wrapper struct representing a (simplified) linear combination of vectors.  If each
/// vector in the combination returns elements in (strictly or nonstrictly) ascending 
/// order of index, then this struct returns entries in **non-strictly** ascending order
/// of index (repeat indices may occurr).
/// 
/// This struct is intended simply to simplify the expression of a more complicated type
pub struct LinearCombinationUnsimplified
                < EntryIter, Index, SnzVal, RingOperator, OrderComparator >
    where
        EntryIter:          Iterator,
        EntryIter::Item:    KeyValSet< Index, SnzVal >,
        SnzVal:             Clone,
        RingOperator:       Semiring< SnzVal >,
        OrderComparator:    StrictlyLess<  EntryIter::Item >,
{
    pub linear_combination_unsimplified:    HitMerge<   
                                                    // type parameter #1 = a scalar multiple of a vector
                                                    Scale<
                                                        EntryIter,
                                                        Index,
                                                        RingOperator,
                                                        SnzVal
                                                    >,
                                                    // type parameter #2 = the object that determines order on entries
                                                    OrderComparator,
                                                >,
}

// implement unwrap
impl < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    LinearCombinationUnsimplified
        < EntryIter, Index, SnzVal, RingOperator, OrderComparator >    

    where
        EntryIter:          Iterator,
        EntryIter::Item:    KeyValSet< Index, SnzVal >,
        SnzVal:             Clone,
        RingOperator:       Semiring< SnzVal >,
        OrderComparator:    StrictlyLess<  EntryIter::Item >,     
    {
        /// Returns the merged entry iterator contained in the wrapper.
        pub fn unwrap_unsimplifed_lin_comb( self ) 
                -> 
                HitMerge<   
                        // type parameter #1 = a scalar multiple of a vector
                        Scale<
                            EntryIter,
                            Index,
                            RingOperator,
                            SnzVal
                        >,
                        // type parameter #2 = the object that determines order on entries
                        OrderComparator,
                    > {
                self.linear_combination_unsimplified
            }
    }     

// implement Clone
impl < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    Clone for 
    
    LinearCombinationUnsimplified  
        < EntryIter, Index, SnzVal, RingOperator, OrderComparator >    

    where
        EntryIter:          Clone + Iterator,
        EntryIter::Item:    Clone + KeyValSet< Index, SnzVal >,
        SnzVal:             Clone,
        RingOperator:       Clone + Semiring< SnzVal >,
        OrderComparator:    Clone + StrictlyLess<  EntryIter::Item >,        
        Scale<EntryIter, Index, RingOperator, SnzVal>: Clone,
{
    fn clone( &self ) -> Self { LinearCombinationUnsimplified{ linear_combination_unsimplified: self.linear_combination_unsimplified.clone() } }
}

// implement the iterator trait on the wrapper
impl < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    Iterator
        for

    LinearCombinationUnsimplified
        < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    where
        EntryIter:          Iterator,
        EntryIter::Item:    KeyValGet< Index, SnzVal > + KeyValSet< Index, SnzVal >,
        SnzVal:             Clone,
        RingOperator:       Semiring< SnzVal >,
        OrderComparator:    StrictlyLess<  EntryIter::Item >,  
    {
        type Item = EntryIter::Item;

        fn next( &mut self ) -> Option< Self::Item > {
            self.linear_combination_unsimplified.next()
        }
    }



//  SIMPLIFIED
//  -------------------------------------------------------------------------------------

/// A wrapper struct representing a (simplified) linear combination of vectors.  If each
/// vector in the combination returns elements in (strictly or nonstrictly) ascending 
/// order of index, then this struct returns entries in **strictly** ascending order
/// of index (no repeat indices).
/// 
/// This struct is intended simply to simplify the expression of a more complicated type
pub struct LinearCombinationSimplified
                < EntryIter, Index, SnzVal, RingOperator, OrderComparator >
    where
        EntryIter:          Iterator,
        EntryIter::Item:    KeyValGet< Index, SnzVal > + KeyValSet< Index, SnzVal >,
        SnzVal:             Clone,
        RingOperator:       Semiring< SnzVal >,
        OrderComparator:    StrictlyLess<  EntryIter::Item >,        
{
    pub linear_combination_simplified:  Simplify<
                                                HitMerge<   
                                                        // type parameter #1 = a scalar multiple of a vector
                                                        Scale<
                                                            EntryIter,
                                                            Index,
                                                            RingOperator,
                                                            SnzVal
                                                        >,
                                                        // type parameter #2 = the object that determines order on entries
                                                        OrderComparator,
                                                    >,
                                                Index,
                                                RingOperator,
                                                SnzVal,
                                            > 
}

// implement unwrap
impl < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    LinearCombinationSimplified
        < EntryIter, Index, SnzVal, RingOperator, OrderComparator >    

    where
        EntryIter:          Iterator,
        EntryIter::Item:    KeyValGet< Index, SnzVal > + KeyValSet< Index, SnzVal >,
        SnzVal:             Clone,
        RingOperator:       Semiring< SnzVal >,
        OrderComparator:    StrictlyLess<  EntryIter::Item >,         
    {

        /// Returns the simplified, merged entry iterator contained in the wrapper.        
        pub fn unwrap_simplified_lin_comb( self ) 
                -> 
                Simplify<
                    HitMerge<   
                            // type parameter #1 = a scalar multiple of a vector
                            Scale<
                                EntryIter,
                                Index,
                                RingOperator,
                                SnzVal
                            >,
                            // type parameter #2 = the object that determines order on entries
                            OrderComparator,
                        >,
                    Index,
                    RingOperator,
                    SnzVal,
                > {
                self.linear_combination_simplified
            }
    }        

// implement Clone
impl < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    Clone for 
    
    LinearCombinationSimplified  
        < EntryIter, Index, SnzVal, RingOperator, OrderComparator >    

    where
        EntryIter:          Clone + Iterator,
        EntryIter::Item:    Clone + KeyValSet< Index, SnzVal >,
        Index:             Clone,
        SnzVal:             Clone,
        RingOperator:       Clone + Semiring< SnzVal >,
        OrderComparator:    Clone + StrictlyLess<  EntryIter::Item >,        
        Scale<EntryIter, Index, RingOperator, SnzVal>: Clone,
{
    fn clone( &self ) -> Self { LinearCombinationSimplified{ linear_combination_simplified: self.linear_combination_simplified.clone() } }
}

// implement the iterator trait on the wrapper
impl < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    Iterator
        for

    LinearCombinationSimplified
        < EntryIter, Index, SnzVal, RingOperator, OrderComparator >

    where
        EntryIter:          Iterator,
        EntryIter::Item:    KeyValGet< Index, SnzVal > + KeyValSet< Index, SnzVal >,
        SnzVal:             Clone,
        RingOperator:       Semiring< SnzVal >,
        OrderComparator:    StrictlyLess<  EntryIter::Item >,  
        Index:             PartialEq,
    {
        type Item = EntryIter::Item;

        fn next( &mut self ) -> Option< Self::Item > {
            self.linear_combination_simplified.next()
        }
    }
