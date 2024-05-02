//! Operations on sparse vectors: [`Gather`] , [`Scale`], [`DropZeros`], [`Simplify`], etc.
//!
// //! There are a variety of operations one might wish to apply to a sparse vector or collection of sparse vectors.
//! This module provides functionality for a variety of operations on sparse vectors.  The following examples are not exhaustive; scroll to the bottom for a complete list.
//! 
//! 
//! # Examples
//! 
//! 
//! Let `iter_a`, `iter_b`, and `iter_c` be sparse vectors, i.e. iterators that run over 
//! sparse matrix entries.  For example, we could define `iter_a`, `iter_b`, `iter_c` as follows
//! 
//! ```
//! // First define the entries
//! // Note that `vec!` creates a standard rust vector, which is different 
//! // from the sort of vector we care about)
//! let entries_a   =   vec![ (1, 1.), (4, 4.) ];
//! let entries_b   =   vec![ (2, 2.), (3, 3.) ];
//! let entries_c   =   vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
//! 
//! // Now define the sparse vector iterators
//! // Note that `iter()` creates an interator, and `cloned()` reformats 
//! // the entries of each iterator.
//! let iter_a      =   entries_a.iter().cloned(); 
//! let iter_b      =   entries_b.iter().cloned();
//! let iter_c      =   entries_c.iter().cloned();
//! ```
//! 
//! Let's also define the operator of the coefficient ring we want to work with.
//!     
//! ```
//! // Load the module that allows us to define our coefficient ring.
//! use oat_rust::rings::operator_structs::ring_native::*;

//! // Define the operator of a coefficient ring (in this case, floating point real numbers)
//! let ring_operator = DivisionRingNative::<f64>::new();   
//! ```
//!     
//! * **Scale, drop zeros, gather, simplify**
//! 
//!     We can scale, drop zero entries, and gather terms as follows
//!     
//!     ```
//!     use oat_rust::vectors::operations::*;
//!     use oat_rust::rings::operator_structs::ring_native::*;
//! 
//!     # // Define the vector
//!     # let entries_a   =   vec![ (1, 1.), (4, 4.) ];
//!     # let entries_c   =   vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
//!     # let iter_a      =   entries_a.iter().cloned();
//!     # let iter_c      =   entries_c.iter().cloned();
//!     #
//!     # // Define the operator of the coefficient ring (in this case, floating point real numbers)
//!     # let ring_operator = DivisionRingNative::<f64>::new();        
//!       
//!     // SCALE A VECTOR BY 2.
//!     // Example: convert [ (1, 1.), (4, 4.) ] into [ (1, 2.), (4, 8.) ]
//!     let scaled : Vec<_> = iter_a
//!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
//!                             .scale( 2., ring_operator.clone() )
//!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
//!     assert_eq!( scaled, vec![ (1, 2.), (4, 8.) ]);
//!       
//!     // DROP ZERO ENTRIES
//!     // Example: convert [ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ] into [ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]
//!     let dropped : Vec<_> = iter_c
//!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
//!                             .drop_zeros( ring_operator.clone() )
//!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
//!     
//!     assert_eq!( dropped, vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]);
//!       
//!     // GATHER CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX
//!     // The resulting vector has no repeating consecutive indices; each index gets 
//!     // the sum of the corresponding coefficients.
//!     // Example: convert [(1,1.), (1,0.5), (2,0.), (1,0.)] into [(1,1.5), (2,0.), (1,0.)]
//!     let gathered : Vec<_> = iter_c
//!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
//!                             .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
//!                             .gather( ring_operator.clone() )
//!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
//!     assert_eq!( gathered, vec![ (1, 1.), (2, 2.), (3, 6.), (4, 0.) ]);   
//! 
//!     // GATHER CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX, AND DROP RESULTING ZEROS
//!     let simplified : Vec<_> = iter_c
//!                                 .clone() // this makes a copy of the iterator, so the original stays unchanged
//!                                 .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
//!                                 .simplify( ring_operator.clone() )
//!                                 .collect(); // this collects the entries of the iterator into a standard Rust vector
//!     assert_eq!( simplified, vec![ (1, 1.), (2, 2.), (3, 6.), ]);  
//!     ```
//! * **Combine iterators in sorted order** (basic)
//! 
//!   We can combine two iterators, `A` and `B`, into a single iterator `C` using the 
//!     [merge](https://docs.rs/itertools/0.7.2/itertools/fn.merge.html) function from [itertools](https://docs.rs/itertools/latest/itertools/). 
//!     The resulting iterator, `C`, is a 
//!     [Merge struct](https://docs.rs/itertools/0.7.8/itertools/structs/struct.Merge.html).
//!     Iterator `C` will iterate over all the entries in `A` and `B`.
//!     If the items of `A` and `B` appear in sorted order, then the items of `C` will also 
//!     appear in sorted order.
//!     ```
//!     # let entries_a   =   vec![ (1, 1.), (4, 4.) ];
//!     # let entries_b   =   vec![ (2, 2.), (3, 3.) ];
//!     # let iter_a      =   entries_a.iter().cloned(); 
//!     # let iter_b      =   entries_b.iter().cloned();
//! 
//!     use itertools::merge;
//!     use std::iter::FromIterator;
//!     
//!     // Merge [ (1, 1.), (4, 4.) ] and [ (2, 2.), (3, 3.) ].
//!     // The entries in these vectors are in sorted order, so the resulting iterator will 
//!     // also produce items in sorted order.
//!     let iter_merged   =   merge(iter_a, iter_b);
//!     let entries_mrgd  =   Vec::from_iter(iter_merged);
//!     assert_eq!( entries_mrgd, vec![ (1, 1.), (2, 2.), (3, 3.), (4, 4.) ])
//!     ```
//! 
//!     We can merge any `k` iterators of the same kind using the 
//!     [merge](https://docs.rs/itertools/0.7.2/itertools/fn.merge.html) 
//!     function from 
//!     [itertools](https://docs.rs/itertools/latest/itertools/).  
//! 
//! * **Combine iterators in sorted order** (advanced)
//!  
//!     For advanced usage (eg matrix reduction), we also provide a
//!     customized merge process in the [hit_merge](crate::utilities::iterators::merge::heap_of_iterators) module.
//! 
//! * **Add**
//! 
//!     We can add the vectors represented by `iter_a` and `iter_b` by 
//!     first combining (e.g., with the `merge` function discussed above), 
//!     then applying the `gather` method.
//! 
//! * **Subtract**
//! 
//!     We can subtract  `iter_a` from `iter_b` by first scaling `iter_a` by `-1`, then adding.

use serde_json::map::Entry;

use crate::utilities::functions::evaluate::EvaluateFunction;
use crate::utilities::iterators::general::{PeekUnqualified};
use crate::entries::{KeyValGet, KeyValSet, KeyValNew};
use crate::rings::operator_traits::{Semiring, Ring};
use crate::utilities::sets::MapHasKeyOrSequenceHasElement;
use std::collections::HashMap;
use std::fmt::{Debug};
use std::marker::PhantomData;



//  ---------------------------------------------------------------------------
//  ITERATOR WRAPPER / TRANSFORMATIONS 
//  ---------------------------------------------------------------------------




//  ---------------------------------------------------------------------------
//  DROP ZEROS


/// Iterates over the same items as `self.undropped`, skipping any items with coefficient 0.
/// 
/// Formally, we skip any item `x` such that `self.ring.is_0( x.val() )==true`.
/// 
#[derive(Debug, Clone)]
pub struct DropZeros 
    < EntryIter, Index, RingOperator, RingElement > 
    where   EntryIter:          Iterator,
            EntryIter::Item:    KeyValGet < Index, RingElement>,
            RingOperator:       Semiring < RingElement >,
{
    undropped:              EntryIter,
    ring_operator:          RingOperator,
    phantom_index:          PhantomData< Index >,
    phantom_ringelement:     PhantomData< RingElement >,    
}

impl    < EntryIter, Index, RingOperator, RingElement > 
        
        Iterator for DropZeros 
        
        < EntryIter, Index, RingOperator, RingElement > 

        where   EntryIter:           Iterator,
                EntryIter::Item:     KeyValGet < Index, RingElement>,
                RingOperator:   Semiring < RingElement >,

{
    type Item = EntryIter::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        let mut next = self.undropped.next();

        while let Some( ref x ) = next {
            if self.ring_operator.is_0( x.val() ) { next = self.undropped.next(); }
            else {break} 
        }
        return next 
    }
}

//  ---------------------------------------------------------------------------
//  SCALE


/// Iterates over the same items as `self.unscaled`, with all coefficients scaled by `self.scale`.
#[derive(Debug, Clone)]
pub struct Scale      
    
    < EntryIter, Index, RingOperator, RingElement > 
    
    where   EntryIter:           Iterator,
            EntryIter::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
            RingOperator:        Semiring < RingElement >,
            // RingElement: Debug + Clone,          
{
    unscaled:       EntryIter,
    scale:          RingElement,    
    ring_operator:  RingOperator,
    phantom_index:  PhantomData< Index >,
}

impl    < EntryIter, Index, RingOperator, RingElement > 
        
        Iterator for Scale
        
        < EntryIter, Index, RingOperator, RingElement > 
   
        where   EntryIter:          Iterator,
                EntryIter::Item:    KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
                RingOperator:       Semiring < RingElement >,
                // Index:          Debug + Clone,
                RingElement:        Clone,

{
    type Item = EntryIter::Item;

    fn next( &mut self) -> Option< Self::Item > 
        {
            if let Some( mut x ) = self.unscaled.next() { 
                x.set_val( 
                    self.ring_operator.multiply( 
                        x.val(), 
                        self.scale.clone(), 
                    )
                );
                Some(x)
            }
            else { None }
        }
}


//  ---------------------------------------------------------------------------
//  GATHER COEFFICIENTS 


/// Iterates over the same items as `self.ungathered`, except that 
/// consecutive entries with equal indices are merged into a single entry whose
/// coefficient is the sum of the coefficients.
#[derive(Debug, Clone)]
pub struct Gather
    
    < EntryIter, Index, RingOperator, RingElement > 

    where   EntryIter:           Iterator + PeekUnqualified,
            EntryIter::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
            RingOperator:   Semiring < RingElement >,
            // <EntryIter::Item as KeyValGet>::Key: Debug + Clone,
            // RingElement: Debug + Clone,    

{
    ungathered: EntryIter,
    ring_operator: RingOperator,
    phantom_index:          PhantomData< Index >,
    phantom_ringelement:     PhantomData< RingElement >,     
}



impl    < EntryIter, Index, RingOperator, RingElement > 

        Iterator for Gather
    
        < EntryIter, Index, RingOperator, RingElement > 
   
        where   EntryIter:           Iterator + PeekUnqualified,
                EntryIter::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
                RingOperator:   Semiring < RingElement >,
                Index:          PartialEq,
                // <EntryIter::Item as KeyValGet>::Key: Debug + Clone + PartialEq,
                // RingElement: Debug + Clone,  
{
    type Item = EntryIter::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        if let Some( mut x ) = self.ungathered.next() {
            while let Some( peek ) = self.ungathered.peek_unqualified() {
                if peek.key() == x.key() { 
                    x.set_val(
                        self.ring_operator.add( 
                            x.val(), 
                            peek.val() 
                        )
                    );
                    let _ = self.ungathered.next(); // we have already gotten what we need
                }
                else { break }
            }
            return Some( x )
        }
        else 
        { None }
    }
}


//  ---------------------------------------------------------------------------
//  SIMPLIFY


/// Combines consecutive elements with equal index into a single entry; drops this entry, if it is zero.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::rings::operator_structs::ring_native::DivisionRingNative;
/// use oat_rust::vectors::operations::Simplify;
/// use itertools::Itertools;
/// 
/// 
/// // Define the operator of the coefficient ring_operator.
/// let ring_operator = DivisionRingNative::<f64>::new();  
/// 
/// // Simplify some iterators:
/// let vec = vec![ (0, 0.), (0, 1.), (1, 1.), (1, -1.), (2, 1.)];
/// let unsimplified = vec.iter().cloned().peekable();
/// let simplified = Simplify::new( unsimplified, ring_operator.clone() );       
/// assert_eq!( simplified.collect_vec(), vec![ (0,1.), (2, 1.) ] );
/// 
/// let vec: Vec<(usize, f64)> = vec![];
/// let unsimplified = vec.iter().cloned().peekable();            
/// let simplified = Simplify::new( unsimplified, ring_operator.clone() );       
/// assert_eq!( simplified.collect_vec(), vec![] );
/// 
/// let vec = vec![ (0, 0.), (0, 0.) ];
/// let unsimplified = vec.iter().cloned().peekable();            
/// let simplified = Simplify::new( unsimplified, ring_operator.clone() );              
/// assert_eq!( simplified.collect_vec(), vec![] );
/// ```
#[derive(Debug, Clone)]
pub struct Simplify
    
    < EntryIter, Index, RingOperator, RingElement > 

    where   EntryIter:          Iterator + PeekUnqualified,
            EntryIter::Item:    KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
            RingOperator:       Semiring< RingElement >,

{
    pub unsimplified:           EntryIter,
    pub ring_operator:          RingOperator,
    pub phantom_index:          PhantomData< Index >,
    pub phantom_ringelement:    PhantomData< RingElement >,
}


impl    < EntryIter, Index, RingOperator, RingElement > 
        
        Simplify
        < EntryIter, Index, RingOperator, RingElement > 

        where   EntryIter:           Iterator + PeekUnqualified,
        EntryIter::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
        RingOperator:   Semiring< RingElement >,        
{
    pub fn new( unsimplified: EntryIter, ring_operator: RingOperator ) 
        -> 
        Simplify < EntryIter, Index, RingOperator, RingElement >   {
            
        Simplify{
                unsimplified:           unsimplified,
                ring_operator:          ring_operator,
                phantom_index:          PhantomData,
                phantom_ringelement:    PhantomData,
            }
        
    }
}


impl    < EntryIter, Index, RingOperator, RingElement > 

        Iterator for Simplify
    
        < EntryIter, Index, RingOperator, RingElement > 
   
        where   EntryIter:           Iterator + PeekUnqualified,
                EntryIter::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
                RingOperator:   Semiring< RingElement >,
                Index:          PartialEq,
{
    type Item = EntryIter::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        while let Some( mut x ) = self.unsimplified.next() {
            // ----------!!!
            // let mut counter = 0;
            // ---------- !!!

            while let Some( peek ) = self.unsimplified.peek_unqualified() {
                // ----------!!!
                // println!("SIMPLIFY is trying to add a {:?}th term", counter);
                // counter +=1;
                // ----------!!!

                if peek.key() == x.key() { 
                    x.set_val(
                        self.ring_operator.add( 
                            x.val(), 
                            peek.val() 
                        )
                    );
                    let _ = self.unsimplified.next(); // we have already gotten what we need
                }
                else { break }
            }
            match self.ring_operator.is_0( x.val() ){
                false   => { return Some( x ) },    //  return the gathered entry, if its coefficient is nonzero                
                true    => { continue },         //  otherwise skip to the next entry
            }
        }
        // else 
        // { println!("returning a simplified term"); None }        
        // println!("returning a simplified term"); 
        return None
    }
}










//  ---------------------------------------------------------------------------
//  NEGATE


/// Negates the coefficient of each entry, i.e., swaps `( i, a )` with `(i, -a)`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::vectors::operations::Negate;
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
/// 
/// let vec = vec![ (0, 0usize), (1, 1usize) ];
/// let ring_operator = PrimeOrderFieldOperator::new(7);
/// let negated = Negate::new( vec.iter().cloned(), ring_operator );
/// itertools::assert_equal( negated, vec![ (0, 0usize), (1, 6usize) ] );
/// ```
#[derive(Debug, Clone)]
pub struct Negate 
    < EntryIter, Index, RingOperator, RingElement > 
    where   EntryIter:          Iterator,
            EntryIter::Item:    KeyValGet < Index, RingElement> + KeyValSet < Index, RingElement>,
            RingOperator:       Semiring < RingElement >,
{
    unnegated:              EntryIter,
    ring_operator:          RingOperator,
    phantom_index:          PhantomData< Index >,
    phantom_ringelement:    PhantomData< RingElement >,    
}

impl    < EntryIter, Index, RingOperator, RingElement > 

        Negate 
            < EntryIter, Index, RingOperator, RingElement > 

    where   EntryIter:          Iterator,
            EntryIter::Item:    KeyValGet < Index, RingElement> + KeyValSet < Index, RingElement>,
            RingOperator:       Semiring < RingElement >,            
{
    pub fn new( unnegated: EntryIter, ring_operator: RingOperator ) 
            ->         
            Negate < EntryIter, Index, RingOperator, RingElement > 
        { Negate{ unnegated, ring_operator, phantom_index: PhantomData, phantom_ringelement: PhantomData } }
}            

impl    < EntryIter, Index, RingOperator, RingElement > 
        
        Iterator for 
        
        Negate 
            < EntryIter, Index, RingOperator, RingElement > 

        where   EntryIter:           Iterator,
                EntryIter::Item:     KeyValGet < Index, RingElement> + KeyValSet < Index, RingElement>,
                RingOperator:        Semiring < RingElement > + Ring < RingElement >,

{
    type Item = EntryIter::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        match self.unnegated.next() {
            Some( mut next_item ) => {
                next_item.set_val( 
                        self.ring_operator.negate( next_item.val() )
                    );
                return Some( next_item )
            },
            None => { return None }
        }
    }
}




//  ---------------------------------------------------------------------------
//  MASK



/// A "mask" over the iterator; it returns an index-coefficient pair iff the index is a key in the
/// gien hashmap.
/// 
/// # Design notes
/// 
/// It would be possible to construct an object with the same behavior via Rust's [filter](https://doc.rust-lang.org/stable/std/iter/trait.Iterator.html#method.filter)
/// method.  However, to achieve the desired result one needs access to the hashmap; the most natural way
/// to provide such acceess to a [`Filter`](https://doc.rust-lang.org/stable/std/iter/struct.Filter.html) object is
/// to use a closure.  However, closures give rise to a host of troubles vis-a-vis type inferrence.  The object
/// that we construct here is intended to serve as a basic primitive; a building block that will be used 
/// throughout the library and combined in complicated ways with other structs.  Therefore we prefer to define
/// a new object with better properties vis-a-vis type inference.
pub struct OnlyIndicesInsideCollection<
        EntryIter, CollectionToInclude, Index, RingElement,
    >
    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < Index, RingElement >,
            CollectionToInclude:    MapHasKeyOrSequenceHasElement< Index >,
{
    collection_to_include:      CollectionToInclude, 
    entry_iter:                 EntryIter,
    phantom_index:              PhantomData< Index >,
    phantom_ringelement:        PhantomData< RingElement >,    
}

impl < EntryIter, CollectionToInclude, Index, RingElement, > 

    OnlyIndicesInsideCollection
        < EntryIter, CollectionToInclude, Index, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < Index, RingElement >,
            CollectionToInclude:    MapHasKeyOrSequenceHasElement< Index >,        
{
    /// Create a new `OnlyIndicesInsideCollection`.  See documentation for [`OnlyOutsideCollection`].
    pub fn new
        ( entry_iter: EntryIter, collection_to_include: CollectionToInclude ) 
        -> 
        OnlyIndicesInsideCollection
            < EntryIter, CollectionToInclude, Index, RingElement > 
        {
            OnlyIndicesInsideCollection{ collection_to_include, entry_iter, phantom_index: PhantomData, phantom_ringelement: PhantomData }
        }        
}

impl < EntryIter, Index, RingElement, CollectionToInclude, > 
    
    Iterator 
        for

    OnlyIndicesInsideCollection
        < EntryIter, CollectionToInclude, Index, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < Index, RingElement >,
            CollectionToInclude:    MapHasKeyOrSequenceHasElement< Index >,
{
    type Item = EntryIter::Item;

    fn next( &mut self ) -> Option< Self::Item > {
        // println!("OnlyIndicesInsideCollection: next(): BEGINNING");
        while let Some( entry ) = self.entry_iter.next() {
            // println!("OnlyIndicesInsideCollection: next(): BEGINNING WHILE LOOP");
            // println!("IT APPEARS THAT THE ERROR OCCURS WHEN WE TRY TO CALL self.collection_to_include.map_has_key_or_sequence_has_element ");
            match self.collection_to_include.map_has_key_or_sequence_has_element( & entry.key() ) {
                true    =>  { return Some( entry ) },
                false   =>  { continue },
            }
        }
        // println!("OnlyIndicesInsideCollection: next(): ENDING");        
        return None
    }
}


//  ---------------------------------------------------------------------------
//  OTHER MASK



/// A "mask" over an iterator; it returns an index-coefficient pair iff the index is **not** a key in the
/// given hashmap or an element of the given vector.
/// 
/// # Design notes
/// 
/// It would be possible to construct an object with the same behavior via Rust's [filter](https://doc.rust-lang.org/stable/std/iter/trait.Iterator.html#method.filter)
/// method.  However, to achieve the desired result one needs access to the hashmap; the most natural way
/// to provide such acceess to a [`Filter`](https://doc.rust-lang.org/stable/std/iter/struct.Filter.html) object is
/// to use a closure.  However, closures give rise to a host of troubles vis-a-vis type inferrence.  The object
/// that we construct here is intended to serve as a basic primitive; a building block that will be used 
/// throughout the library and combined in complicated ways with other structs.  Therefore we prefer to define
/// a new object with better properties vis-a-vis type inference.
pub struct OnlyIndicesOutsideCollection<
        EntryIter, CollectionToExclude, Index, RingElement,
    >
{
    collection_to_exclude:      CollectionToExclude, 
    entry_iter:                 EntryIter,
    phantom_index:              PhantomData< Index >,
    phantom_ringelement:        PhantomData< RingElement >,
}

impl < EntryIter, CollectionToExclude, Index, RingElement, > 

    OnlyIndicesOutsideCollection
        < EntryIter, CollectionToExclude, Index, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < Index, RingElement >,
            CollectionToExclude:    MapHasKeyOrSequenceHasElement< Index >,        
{
    /// Create a new `OnlyOutsideCollection`.  See documentation for [`OnlyOutsideCollection`].
    pub fn new
        ( entry_iter: EntryIter, collection_to_exclude: CollectionToExclude ) 
        -> 
        OnlyIndicesOutsideCollection
            < EntryIter, CollectionToExclude, Index, RingElement, > 
        {
            OnlyIndicesOutsideCollection{ entry_iter: entry_iter, collection_to_exclude: collection_to_exclude, phantom_index: PhantomData, phantom_ringelement: PhantomData }
        }        
}


impl < EntryIter, CollectionToExclude, Index, RingElement, > 
    
    Iterator 
        for

    OnlyIndicesOutsideCollection
        < EntryIter, CollectionToExclude, Index, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < Index, RingElement >,
            CollectionToExclude:    MapHasKeyOrSequenceHasElement< Index >,
{
    type Item = EntryIter::Item;

    fn next( &mut self ) -> Option< Self::Item > {
        while let Some( entry ) = self.entry_iter.next() {
            match self.collection_to_exclude.map_has_key_or_sequence_has_element( & entry.key() ) {
                false   =>  return Some( entry ),
                true    =>  continue,
            }
        }
        return None
    }
}







//  ---------------------------------------------------------------------------
//  MAP BY INDEX


/// Similar to map for iterators, however the only thing we map is the index.
/// 
/// This object always iteratos over tuples.  If you want to interate over some other kind of objec,
/// it may be worth considering the [`MapByTransform`](oat_rust::utilities::iterators::general) struct.
///
/// # Examples
/// 
/// ```
/// use std::collections::HashMap;
/// use oat_rust::vectors::operations::ChangeIndexSimple;
/// 
/// // create iterator
/// let entry_iter_data: Vec< (usize, usize) > = vec![(1,1), (2,2)];
/// let entry_iter = entry_iter_data.iter();
/// 
/// // hashmaps implement `EvaluateFunction` automatically; see documentaion for `EvaluateFunction` for details
/// let mut hash = HashMap::new();
/// hash.insert( 1usize, -1i32 );
/// hash.insert( 2usize, -2i32 );
/// 
/// // create a filter/mapped iterator        
/// let iter_filtermapped_by_index: ChangeIndexSimple<_,_,_,i32,usize> = ChangeIndexSimple::new( entry_iter, &hash );
/// 
/// // check that it is correct        
/// itertools::assert_equal( iter_filtermapped_by_index, vec![(-1,1usize), (-2,2usize)] );
/// ```
pub struct ChangeIndexSimple<
        EntryIter, IndexChanger, IndexOld, IndexNew, RingElement,
    >
{
    index_changer:              IndexChanger, 
    entry_iter:                 EntryIter,
    phantom_indexold:           PhantomData< IndexOld >,
    phantom_indexnew:           PhantomData< IndexNew >,    
    phantom_ringelement:        PhantomData< RingElement >,
}

impl < EntryIter, IndexChanger, IndexOld, IndexNew, RingElement, > 

    ChangeIndexSimple
        < EntryIter, IndexChanger, IndexOld, IndexNew, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < IndexOld, RingElement >,
{
    /// Create a new `ChangeIndexSimple`.  See documentation for [`ChangeIndexSimple`].
    pub fn new
        ( entry_iter: EntryIter, index_changer: IndexChanger ) 
        -> 
        ChangeIndexSimple
            < EntryIter, IndexChanger, IndexOld, IndexNew, RingElement, > 
        {
            ChangeIndexSimple{ entry_iter, index_changer, phantom_indexold: PhantomData, phantom_indexnew: PhantomData, phantom_ringelement: PhantomData }
        }        
}


impl < EntryIter, IndexChanger, IndexOld, IndexNew, RingElement, > 
    
    Iterator 
        for

    ChangeIndexSimple
        < EntryIter, IndexChanger, IndexOld, IndexNew, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < IndexOld, RingElement >,
            IndexChanger:            EvaluateFunction< IndexOld, IndexNew >,
{
    type Item = ( IndexNew, RingElement );

    fn next( &mut self ) -> Option< Self::Item > {
        match self.entry_iter.next() {
            Some( entry ) => {
                return Some(  (  self.index_changer.evaluate_function( entry.key() ), entry.val()  )  )
            },
            None => { return None }
        }
    }
}



//  ---------------------------------------------------------------------------
//  CHANGE ENTRY TYPE


/// Wrapper around an entry iterator; the wrapper iterates over the same sequence of index-coefficient pairs, but changes the type of each entry to `EntryNew`
pub struct ChangeEntryType
                <   EntryIter, EntryNew, Index, RingElement, >
{
    entry_iter:                 EntryIter,
    phantom_entrynew:           PhantomData< EntryNew >,
    phantom_index:              PhantomData< Index >,    
    phantom_ringelement:        PhantomData< RingElement >,
}

impl < EntryIter, EntryNew, Index, RingElement, >

    ChangeEntryType
        < EntryIter, EntryNew, Index, RingElement, >

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < Index, RingElement >,
            EntryNew:               KeyValGet < Index, RingElement > + KeyValNew < Index, RingElement >,            
{
    /// Create a new `ChangeEntryType` struct.  See documentation for [`ChangeEntryType`].
    pub fn new
        ( entry_iter: EntryIter ) 
        -> 
        ChangeEntryType
            < EntryIter, EntryNew, Index, RingElement, >
        {
            ChangeEntryType{ entry_iter, phantom_entrynew: PhantomData, phantom_index: PhantomData, phantom_ringelement: PhantomData }
        }        
}


impl < EntryIter, EntryNew, Index, RingElement, >
    
    Iterator 
        for

    ChangeEntryType
        < EntryIter, EntryNew, Index, RingElement, >

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < Index, RingElement >,
            EntryNew:               KeyValGet < Index, RingElement > + KeyValNew < Index, RingElement >,
{
    type Item = EntryNew;

    fn next( &mut self ) -> Option< Self::Item > {
        match self.entry_iter.next() {
            Some( entry ) => {
                return Some(  EntryNew::new(  entry.key(), entry.val()  )  )
            },
            None => { return None }
        }
    }
}


//  ---------------------------------------------------------------------------
//  FILTER-MAP BY INDEX


/// Similar to filter-map for iterators, however the only thing we filter and map is the index.
///
/// # Examples
/// 
/// ```
/// use std::collections::HashMap;
/// use oat_rust::vectors::operations::FilterChangeIndex;
/// 
/// // create iterator
/// let entry_iter_data: Vec< (usize, usize) > = vec![(1,1), (2,2), (3,3)];
/// let entry_iter = entry_iter_data.iter();
/// 
/// // hashmaps implement `EvaluateFunction` automatically; see documentaion for `EvaluateFunction` for details
/// let mut hash = HashMap::new();
/// hash.insert( 1usize, -1i32 );
/// hash.insert( 2usize, -2i32 );
/// 
/// // create a filter/mapped iterator        
/// let iter_filtermapped_by_index: FilterChangeIndex<_,_,_,i32,usize> = FilterChangeIndex::new( entry_iter, &hash );
/// 
/// // check that it is correct        
/// itertools::assert_equal( iter_filtermapped_by_index, vec![(-1,1usize), (-2,2usize)] );
/// ```
pub struct FilterChangeIndex<
        EntryIter, IndexFilterChanger, IndexOld, IndexNew, RingElement,
    >
{
    index_filter_changer:       IndexFilterChanger, 
    entry_iter:                 EntryIter,
    phantom_indexold:           PhantomData< IndexOld >,
    phantom_indexnew:           PhantomData< IndexNew >,    
    phantom_ringelement:        PhantomData< RingElement >,
}

impl < EntryIter, IndexFilterChanger, IndexOld, IndexNew, RingElement, > 

    FilterChangeIndex
        < EntryIter, IndexFilterChanger, IndexOld, IndexNew, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < IndexOld, RingElement >,
{
    /// Create a new `OnlyOutsideCollection`.  See documentation for [`OnlyOutsideCollection`].
    pub fn new
        ( entry_iter: EntryIter, index_filter_changer: IndexFilterChanger ) 
        -> 
        FilterChangeIndex
            < EntryIter, IndexFilterChanger, IndexOld, IndexNew, RingElement, > 
        {
            FilterChangeIndex{ entry_iter, index_filter_changer, phantom_indexold: PhantomData, phantom_indexnew: PhantomData, phantom_ringelement: PhantomData }
        }        
}


impl < EntryIter, IndexFilterChanger, IndexOld, IndexNew, RingElement, > 
    
    Iterator 
        for

    FilterChangeIndex
        < EntryIter, IndexFilterChanger, IndexOld, IndexNew, RingElement, > 

    where   EntryIter:              Iterator,
            EntryIter::Item:        KeyValGet < IndexOld, RingElement >,
            IndexFilterChanger:      EvaluateFunction< IndexOld, Option<IndexNew> >,
{
    type Item = ( IndexNew, RingElement );

    fn next( &mut self ) -> Option< Self::Item > {
        while let Some( entry ) = self.entry_iter.next() {
            match self.index_filter_changer.evaluate_function( entry.key() ) {
                None    =>  return None,
                Some( index_new)    => return Some( (index_new, entry.val() ) ),
            }
        }
        return None
    }
}



















//  ---------------------------------------------------------------------------
//  TRANSFORMATION TRAIT
//  ---------------------------------------------------------------------------


/// Provides some convenient methods for transforming vectors that are auto-implemented on most vector iterators.
///
/// The methods in this trait are implemented automatically for structs that
/// implement `Iterator< Item : KeyValGet >`.  They are not complicated
/// in general; rather, they're provided as a convenience for 
/// chained transformations.
/// 
/// Most methods require the entry iterator to be sized; *however* recall that
/// this trait is just a convenience.  You can usually perform the desired
/// transformation on you iterator directly, even if you cannot perform it
/// with this trait because the iterator violates the `Sized` criterion.
pub trait Transforms< Index, RingElement>

    where   Self:           Iterator,
            Self::Item:     KeyValGet < Index, RingElement>,

{

    /// Returns an interator that iterates over the same items as `self`, 
    /// skipping any items with coefficient 0.
    fn  drop_zeros< RingOperator >
        ( self, ring_operator: RingOperator ) 
        -> 
        DropZeros< Self, Index, RingOperator, RingElement > 
        
        where   Self:           Iterator + Sized,
                Self::Item:     KeyValGet < Index, RingElement>,
                RingOperator:   Semiring < RingElement >,
                // <Self::Item as KeyValGet>::Key: Debug + Clone,
                // <Self::Item as KeyValGet>::Val: Debug + Clone,

    {
        DropZeros{ undropped: self, ring_operator: ring_operator, phantom_index: PhantomData, phantom_ringelement: PhantomData }
    }

    /// Returns an interator that iterates over the same items as `self`, 
    /// with all coefficients scaled by `scalar`.
    fn  scale 
        < RingOperator > 
        ( self, scalar: RingElement, ring_operator: RingOperator )
        -> 
        Scale < Self, Index, RingOperator, RingElement >
        
        where   Self:           Iterator + Sized,
                Self::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
                RingOperator:   Semiring < RingElement >,
                // <Self::Item as KeyValGet>::Key: Debug + Clone,
                // RingElement: Clone,
    {
        Scale{ unscaled: self, scale: scalar, ring_operator: ring_operator, phantom_index: PhantomData }
    }

    /// Returns an interator that iterates over the same items as `self`, except that 
    /// consecutive entries with equal indices are merged into a single entry whose
    /// coefficient is the sum of the coefficients.  
    fn gather < RingOperator > ( self, ring_operator: RingOperator )
        -> Gather < Self, Index, RingOperator, RingElement >

        where   Self:           Iterator + Sized + PeekUnqualified,
                Self::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
                RingOperator:   Semiring < RingElement >,
                Index:          PartialEq,
                // <Self::Item as KeyValGet>::Key: Debug + Clone,
                // <Self::Item as KeyValGet>::Val: Debug + Clone,               
    {
        Gather{ ungathered: self, ring_operator: ring_operator, phantom_index: PhantomData, phantom_ringelement: PhantomData  } 
    }

    /// Returns an interator that iterates over the same items as `self`, except that 
    /// (i) consecutive entries with equal indices are merged into a single entry whose
    /// coefficient is the sum of the coefficients, (ii) if this entry has coefficient 0,
    /// then it is omitted.  
    fn simplify < RingOperator > ( self, ring_operator: RingOperator )
        -> Simplify < Self, Index, RingOperator, RingElement >

        where   Self:           Iterator + Sized + PeekUnqualified,
                Self::Item:     KeyValGet < Index, RingElement > + KeyValSet < Index, RingElement >,
                RingOperator:   Semiring < RingElement >,
                Index:          PartialEq,
                // <Self::Item as KeyValGet>::Key: Debug + Clone,
                // <Self::Item as KeyValGet>::Val: Debug + Clone,               
    {
        Simplify{ unsimplified: self, ring_operator: ring_operator, phantom_index: PhantomData, phantom_ringelement: PhantomData  } 
    }        

    fn negate < RingOperator > (self, ring_operator: RingOperator )
        -> Negate< Self, Index, RingOperator, RingElement >
        where 
            Self:           Sized,
            RingOperator:   Ring< RingElement >,
            Self::Item:     KeyValSet< Index, RingElement >
    {
        Negate::new( self, ring_operator )
    }


    /// Include a collection.
    fn include_collection
        < CollectionToInclude: MapHasKeyOrSequenceHasElement<Index> >
        ( self, collection_to_include: CollectionToInclude ) 
        -> 
        OnlyIndicesInsideCollection< Self, CollectionToInclude, Index, RingElement > 
            where
                CollectionToInclude:    MapHasKeyOrSequenceHasElement<Index>,
                Self:                   Sized,
                Self::Item:             KeyValGet< Index, RingElement >,
    {
        OnlyIndicesInsideCollection::new( self, collection_to_include ) 
    }

    /// Exclude a collection.
    fn exclude_collection
        < CollectionToExclude: MapHasKeyOrSequenceHasElement<Index> >
        ( self, collection_to_exclude: CollectionToExclude ) 
        -> 
        OnlyIndicesOutsideCollection< Self, CollectionToExclude, Index, RingElement > 
            where
                Self:   Sized
    {
        OnlyIndicesOutsideCollection::new( self, collection_to_exclude ) 
    }
        
}

// We implement this trait automatically on all iterators.
impl    < EntryIter, Index, RingElement >

        Transforms < Index, RingElement >
        
        for EntryIter  

        where   EntryIter:           Iterator,
                EntryIter::Item:     KeyValGet < Index, RingElement>,
                // <EntryIter::Item as KeyValGet>::Key: Debug + Clone,
                // RingElement: Debug + Clone,          
{} // everything implemented automatically


//     DEFINE SPARSE VECTOR ITERATOR AS
//
//     struct Svi< Iter, Index, RingOperator, Coeff > 
//         where   Iter: Iterator< Item = KeyValItem<Index,Coeff> >,
//                 RingOperator: Semiring < Coeff >,
// 
//         {    
// }
// 
// 
//     svi.transform( DropZeros::y,
//                    Scale::Coeff(3)
//                    Merge::y,
//                    ring
//                 )
// 
// 
//     M0 + M1 + M4 + sum(M[2..4,:]) + sum(2 * M[4..6,:]) + 
// 
//     R = gen_ring();
//     C = gen_comparator();
// 
//     let svi0  = M.maj(majkeys[0]).simplify(&R);
//     let svi1  = M.maj(majkeys[1]).simplify(&R);
//     let svi2  = (2..4)
//                     .map(|x| 
//                          M.maj(majkeys[x].simplify(&R)
//                     )
//     let svi3  = (4..6)
//                     .map(|x| 
//                          M.maj(majkeys[x].simplify(&R)
//                         )
//                     .kmerge_by(&C)
//                     .simplify(&R)
//                     .map(|y|
//                             (y.0, R.multiply( &y.1, &2)
//                     )
// 
//     let svi4  = M.maj(majkeys[6]).simplify(&R);
// 
//     let agg = (svi0, svi1, svi4)
//                 .hit_merge(             &C )
//                 .hit_bulk_insert( svi2, &C )
//                 .merge_by( svi3,        &C )
//                 .simplify()
// 
// 
//     !!!!!!
//     drain_monomials
//     drain_monomials_ordered
// 
//     let agg = LciSimplified::combine( ( (svi0, t0), (svi1, t1), (svi3, t3) ), R, C )
//                 .add_combination(( (svi4, t4), (svi5, t5) ))// universal input format
//                 .add_svi( svi4, None    )
//                 .add_svi( svi5, Some(2) )
//                 .add_lci( lci2, None )                            // option for other lci
//                 .add_lci( lci2, Some(3) )                            // option for other lci
//                 .add_sum( (svi6, svi7,  svi8 ), None    )            // avoid 1; avoid nested parentheses 
//                 .add_sum( (svi9, svi10, svi11), Some(2) )  // option to scale
// 
// 
//                 .add_svi( svi4 )                            // avoid 1; avoid nested parentheses 
//                 .add_svi_scaled( svi5, 2 )                  // option to scale
//               
//                 .add_sum( (svi6, svi7, svi8) )            // avoid 1; avoid nested parentheses 
//                 .add_sum_scaled( (svi9, svi10, svi11), 2 )  // option to scale
//                 
//                 .add_lci( lci2 )                            // option for other lci
//                 .add_lci_scaled( lci3, 4 )
// 
//                 .add_svi( svi4 )                            // avoid 1; avoid nested parentheses 
//                 .add_lci( lci2 )                            // option for other lci
//                 .add_sum( (svi6, svi7, svi8) )              // avoid 1 
//                 .add_svi_scaled( svi5, 2 )                  // option to scale
//                 .add_lci_scaled( lci3, 4 )                  // option to scale
//                 .add_sum_scaled( (svi9, svi10, svi11), 2 )  // option to scale
//                 
// 
// 
//                 .add_svi( svi4, R.one() )
//                 .add_svi( svi7, 2       )
//                 .add_lci( lci0, R.one() )
//                 .add_lci( lci8, 3
//                 .add_sum( (svi6, svi7,  svi8 ), R.one() )
//                 .add_sum( (svi9, svi10, svi11), 4       )
// 
//                 .add( Term::Svi(svi4) )
//                 .add( Term::Lci(lci5) )
// 
// 
//                 .add( Term::Svi(svi4), STerm::Svi(svi5, 2), STerm::Lci(lci, 4) )
//                 .add_svi(svi1).add_lci(lci2).add_sum_scaled( (svi2, svi3), 2)
//                 .add_svi(svi1, None).add_svi(svi
// 
// 
//                 .add_combination( ((svi0, 1), (svi2, 2), (svi3, 3)) ).add_sum( (svi4,) ).add
//                 .add_lci( lci, 4 )
//                 .add_sum(
// 
//                 X = LciSimplified::sum( (svi0, svi1), R, C ).add_scaled( (svi
// 
//                 X = LciSimplified::new().add( (
//                         Term::Svi( svi0 ),
//                         Term::Svi( svi1 ).x(2),
//                         Term::Lci( lci0 ),
//                         Term::Lci( lci1 ).x(3),
//                         Term::Sum( (svi2, svi3, svi4) ),
//                         Term::Sum( (svi5, svi6, svi7) ).x(4),
//                         Term::Com( ( (svi10, 3), (svi11, 4) ),
//                         Term::Com( ( (svi12, 3), (svi13, 4) ).x(7),
//                         Term::Com( iter1 ).x(7),
//                         Term::Sum( iter0 ).x(4),
//                         )
//                     )
// 
//                 X.add( Term::Sum(( svi0, svi1, svi2 )), Term::Svi( svi3 ).x(2) );
// 
//                 x.change_add( y, C);
//                 x.change_add_simplify( y, C);
// 
// 
//     enum LTerms< I, Index, Coeff >
//         where I: IntoIterator< Item = KeyValItem<Index,Coeff> >
// 
//     {
//         Scale(   I       ),
//         Unscaled( I, Coeff),
//         Combination( 
//     }
// 
//     lc.add( LTerm::Unscaled( svi    ) )
//     lc.add( LTerm::Scale(   svi, x ) )
//    
//     lc.add_n( LTerms::Unscaled( ( svi0, svi1, svi2 )    ) )
//     lc.add_n( Lterms::Scale(   ( svi0, svi1, svi2 ), y ) )
//   
//     lc.add_n( Lterms::Combination( ( (svi0, t0), (svi1, t1), (svi2, t2) ) )
// 
// 
//     
//     lc.add(        svi0    )
//     lc.add_scaled( svi1, 2 )
// 
//     lc.add_n(        (svi2, svi3)    )
//     lc.add_n_scaled( (svi4, svi5), 2 )
//     
//     lc.add_n_multiscaled(( (svi6, 2), (svi7, 3) ))
//     lc.add_lci( lci )
// //--------------------------------------------------------------
// 
//     lc.add_1( svi );
//     lc.add_n( (svi0, svi1) );
//     lc.add_1_times( svi, 2 );
//     lc.add_n_times( ( svi0, svi1 ), 2 )
//     
//     lc.add_combination( ( (svi0, 2), (svi1, 3) ) );
// 
// 
//     lc.add( LTerm::times_1( svi0    ) )
//     lc.add( LTerm::times_x( svi0, 2 ) )
// 
//     lc.add( LTerms::times_1( (svi0, svi1)    ) )
//     lc.add( LTerms::times_x( (svi0, svi1), 3 ) )
// 
//     lc.add( LTerms::combine( ( (svi0, 2), (svi1, 3) ) )
//     lc.add( lci.into_forgetful_terms() )
// 
// //--------------------------------------------------------------
//     lc.add( svi )
//     lc.add_scaled(svi, 2)
// 
// 
//     lc.add( LTerm::new(svi).x(2) )
//     lc.add_n( (svi0, svi1).map(|x| (x, 2) ) )
// 
// 
//     lc.add( LTerm::scaled( svi, 2 ) )
//     lc.add( LTerm::unscaled( x )
// 
// 
//     lc.add( LTerm::scaled(   svi    ) )
//     lc.add( LTerm::unscaled( svi, x ) )
//    
//     lc.add( LTerms::scaled(   ( svi0, svi1, svi2 )    ) )
//     lc.add( LTerms::unscaled( ( svi0, svi1, svi2 ), y ) )
//   
//     lc.add( LTerms::combination( ( (svi0, t0), (svi1, t1), (svi2, t2) ) )
//     lc.add( LTerms::lci( lci0 ) )
// 
//     lc.add( LTerm::scaled(svi) )
//     lc.add( svi )
// 
//     lc.add_mn( ..) // FOR THESE COMPLEX LARGE SCALE OPERATIONS -- JUST EXPLAIN HOW TO DO A BULK
//     INSERT AND HEAPIFY
// 
// 
//     X.add( LTerm::Scale(   (svi0,), y )
//     X.add( Lterm::Scale(   (svi0, svi1, svi2), y )
//     X.add( LTerm::Unscaled( (svi0,) )
//     X.add( LTerm::Unscaled( (svi0, svi1, svi2) ) )
//     X.add( Lterm::Combination( (svi0, t0), (svi1, t1), (svi2, t2) )
// 
// 
//     !!!!! NOPE -- JUST .merge(y, C).simplify();
// 
//     heter_add( svi0, svi1, C )
//     heter_add_simplify( svi0, svi1, C, R)
// 
//     
// 
// 
//     // Generally speaking, better to build/add in bulk.
//     let agg = simple_lc::sum( (svi0, svi1, svi2), R, C)
//                 .add_vector( svi3    );
//                 .add_scaled( svi3, 2 );
//                 .add_combination( svi4, Coeff::Scalars((3, 2, 1)) )
//                 .add_combination( svi4, Coeff::Uniform(3) )
//                 .add_combination( svi4, Coeff::Unit )
//                 .add_combination( svi4, Coeff::Raw( Some(1), None, Some(2) )
//                 .add_combination( ( (svi5, Some(3)), (svi6, Some(2)), (svi7, None) )
// 
//     let agg = 
// 
// 
//                 .add_vector(      svi3, Coeff::unit )
//                 .add_vector(      svi3, Coeff::scalar(3))
//     let agg = svi2.bulk_insert( (svi0, svi1) ). 
// 
// (svi,)kk
//                 .hit_merge( (svi2,), C )
//                 .simplify(&R)
//                 .merge(  hit_merge(svii, C).simplify(&R)  )
// 
//     svi
//         .gather_terms(R.clone())
//         .drop(R.zero())
//         .hit_merge( 
//     
// 
//     svi + svi2 + sum(svi_iter)
//     svi 
//         .add( svi2, &R)
//         .gather_terms(R.clone())
//         .drop_zeros(R.clone())
//         .add( 
//             hit_merge(svi_iter)
//             .gather_terms(R.clone())
//             );
//         .gather()
//         .drop_zeros(&R)
//    
// 
//     spi
//         .into_svi(R)
//         .gather_terms()
//         .
// 
//     kk
// 
//     svir.add( svi );
//     svir.add_simple( svi );
// 
// 
//     vec.add( vec2, ring ).drop_zero(ring).
// 
//     add_1_simple
//     add_k_simple
//     
//     vec.add_1_simple( vec2, ring )
//     vec.add_k_simple(
// 
//     simple_add
//     simple_sum_several
// 
//     simple_add_k_other
//     simple_add_k_same
// 
//     
//     // UNDER CONSTRUCTION
//    //  fn add < RingOperator, Svi2> (mut self, ring_operator: RingOperator, svi2: Svi2) -> 
// 
// 
//    //  svi.add( svi2, ring ).gather(ring)
//    //  svi.add_several( svi_set, ring ).gather(ring)
// 



//  ---------------------------------------------------------------------------
//  TEST
//  ---------------------------------------------------------------------------


#[cfg(test)]
mod tests {
    use itertools::Itertools;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::rings::operator_structs::ring_native::DivisionRingNative;

    #[test]
    fn test_basic_transforms() {

        // Define the coefficient ring
        let ring_operator = DivisionRingNative::<f64>::new();        

        // Define a sequence of vector entries.
        let entry_data = vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
    
        // Define a sparse vector iterator (i.e., an iterator that runs over entries)       
        let sparse_vec = entry_data.iter().cloned();
    
        // SCALE THE VECTOR BY 2.
        let scaled : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .scale( 2., ring_operator.clone() )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        assert_eq!( scaled, vec![ (1, 2.), (2, 4.), (3, 6.), (3, 6.), (4, 0.) ]);
    
        // DROP ZERO ENTRIES
        let dropped : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .drop_zeros( ring_operator.clone() )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        
        assert_eq!( dropped, vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]);
    
        // MERGE CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX
        let gathered : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
                                .gather( ring_operator.clone() )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        assert_eq!( gathered, vec![ (1, 1.), (2, 2.), (3, 6.), (4, 0.) ]);  
        
        // MERGE CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX, AND DROP RESULTING ZEROS
        let gathered : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
                                .simplify( ring_operator.clone() )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        assert_eq!( gathered, vec![ (1, 1.), (2, 2.), (3, 6.), ]);          
        
    }     
    

    #[test]
    fn doctest_draft_test_filter_map() {    
        use std::collections::HashMap;
        use crate::vectors::operations::{FilterChangeIndex, };

        // create iterator
        let entry_iter_data: Vec< (usize, usize) > = vec![(1,1), (2,2), (3,3)];
        let entry_iter = entry_iter_data.iter();
        
        // hashmaps implement `EvaluateFunction` automatically; see documentaion for `EvaluateFunction` for details
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1i32 );
        hash.insert( 2usize, -2i32 );

        // create a filter/mapped iterator        
        let iter_filtermapped_by_index: FilterChangeIndex<_,_,_,i32,usize> = FilterChangeIndex::new( entry_iter, &hash );

        // check that it is correct        
        itertools::assert_equal( iter_filtermapped_by_index, vec![(-1,1usize), (-2,2usize)] );
    }

    #[test]
    fn doctest_draft_test_negate() {
        use crate::vectors::operations::Negate;
        use crate::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;

        let vec = vec![ (0, 0usize), (1, 1usize) ];
        let ring_operator = PrimeOrderFieldOperator::new(7);
        let negated = Negate::new( vec.iter().cloned(), ring_operator );
        itertools::assert_equal( negated, vec![ (0, 0usize), (1, 6usize) ] );
    }



    // Simplify entry-iterators
    // =====================================================================================================================

    // * SEE DOC TESTS    

}


