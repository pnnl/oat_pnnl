//! (Generalized) matching matrices

use serde_json::map::IntoIter;

use crate::entries::{KeyValGet, KeyValNew, ExtractKey, ExtractVal};
use crate::matrices::matrix_oracle_traits::{OracleMajor, OracleMajorAscend, OracleMajorDescend, IndicesAndCoefficients};
use crate::utilities::functions::compose::ComposeFunctions;
use crate::utilities::functions::evaluate::EvaluateFunctionDisambiguator;
use crate::utilities::iterators::general::IterTwoType;
use crate::utilities::sequences_and_ordinals::BiMapSequential;
use std::collections::HashMap;
use std::hash::Hash;
use std::iter::{Empty, FromIterator, Cloned};
use std::marker::PhantomData;


//  ================================================================================
//  MATCHING ARRAY TRAIT
//  ================================================================================


//  TO BEC COMPLETED:
//
//  A TRAIT FOR MATCHING ARRAYS; METHODS ARE TO INCLUDE
//       - KEYMAJ TO KEYMIN
//       - KEYMIN TO KEYMAJ
//       - KEYMAJ TO SNZVAL
//       - KEYMIN TO SNZVAL
//       - PUSH TRIPLE
//       - RETURN OBJECTS THAT EVALUATE THE 4 FUNCTIONS
//       - ? ITERATE OVER MATCHED ELEMENTS IN ORDER


//  ================================================================================
//  GENERALIZED MATCHING ARRAY STRUCT WITH ORDINALS
//  ================================================================================

/// Encodes two pieces of data: (i) a generalized matching array, and (ii) a bijection between 
/// the set of major keys and {0, .., N}.
/// 
/// An unordered matching array can be encoded as a bijection {major keys} <-> {minor keys}
/// together with a function {major keys} -> structural nonzero values.  By contrast, a 
/// `GeneralizedMatchingArrayWithMajorOrdinals` encodes bijections
/// {major keys} <-> {0, .., N} <-> {minor keys} PLUS a function {0, ..., N} -> {structural nonzero valus}.
///
/// Concretely, the struct stores these data in three different private fields
/// - a [BiMapSequential](`crate::utilities::sequences_and_ordinals::BiMapSequential`) 
///  object representing a bijection from major keys to {0, .., N}.
/// - a [BiMapSequential](`crate::utilities::sequences_and_ordinals::BiMapSequential`) 
///  object representing the corresponding bijection from minor keys to {0, .., N}.
/// - a vector representing the sequence of structural nonzero values
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::matrices::matrix_types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
/// 
/// // initialize new matching array
/// let mut matching_array  =   GeneralizedMatchingArrayWithMajorOrdinals::new();
/// 
/// // push some nonzero entries
/// matching_array.push( 0, 1, 3 );
/// matching_array.push( 1, 2, 4 );
/// matching_array.push( 2, 0, 5 );    
/// 
/// assert_eq!( matching_array.vec_snzval_ref(), &vec![ 3, 4, 5 ] );
/// 
/// assert_eq!( matching_array.keymin_to_keymaj( &0 ), Some(1) );
/// assert_eq!( matching_array.keymaj_to_keymin( &0 ), Some(2) );     
/// 
/// assert_eq!( matching_array.keymin_to_ordmaj( &1 ), Some(1) );
/// assert_eq!( matching_array.keymaj_to_ordmaj( &1 ), Some(0) ); 
/// 
/// assert_eq!( matching_array.ordmaj_to_keymin( 0 ), 0 );
/// assert_eq!( matching_array.ordmaj_to_keymaj( 0 ), 1 );     
/// 
/// assert_eq!( matching_array.keymin_to_snzval( &0 ), 3 );
/// assert_eq!( matching_array.keymaj_to_snzval( &0 ), 5 );       
/// 
/// assert_eq!( matching_array.ordmaj_to_snzval( 0 ), 3 );
/// ```
#[derive(Debug, Clone)]
pub struct GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal >
    where   KeyMin:     Clone + Hash + std::cmp::Eq,
            KeyMaj:     Clone + Hash + std::cmp::Eq,
{  
    bimap_min:  BiMapSequential< KeyMin >,
    bimap_maj:  BiMapSequential< KeyMaj >,
    vec_snzval: Vec< SnzVal >,
    // keymaj_2_ordmaj: HashMap< KeyMaj, usize >,
    // keymin_2_ordmaj: HashMap< KeyMin, usize >,
    // ordmaj_2_keymaj: Vec< KeyMaj >,
    // ordmaj_2_keymin: Vec< KeyMin >,    
    // ordmaj_2_snzval: Vec< SnzVal >,
    // phantom: PhantomData< &'a usize >,
}

impl < 'a, KeyMin, KeyMaj, SnzVal > 

    GeneralizedMatchingArrayWithMajorOrdinals
    
    < KeyMin, KeyMaj, SnzVal >
    where   KeyMin:     Clone + Hash + std::cmp::Eq, // without cloning, it is difficult to write a method that returns the minkey at a given index (one gets "cannot move" errors)
            KeyMaj:     Clone + Hash + std::cmp::Eq, // without cloning, it is difficult to write a method that returns the majkey at a given index (one gets "cannot move" errors)            
    
{
    /// Create a new, empty matching array.
    /// 
    /// See [`GeneralizedMatchingArrayWithMajorOrdinals`] for usage.
    pub fn new() -> GeneralizedMatchingArrayWithMajorOrdinals< KeyMin, KeyMaj, SnzVal > 
    {
        GeneralizedMatchingArrayWithMajorOrdinals{    
            bimap_maj:  BiMapSequential::new(),
            bimap_min:  BiMapSequential::new(),
            vec_snzval: Vec::new(),
        }        
    }

    /// Push a new structural nonzero entry to the array.    
    pub fn push( 
                    &mut self,
                    keymin:     KeyMin, 
                    keymaj:     KeyMaj, 
                    snzval:     SnzVal 
                )
            {
                self.bimap_min.push( keymin ); 
                self.bimap_maj.push( keymaj );               
                self.vec_snzval.push( snzval );
            }

    /// Reverse the order of the ordinals
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::matrices::matrix_types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
    /// 
    /// let mut matching = GeneralizedMatchingArrayWithMajorOrdinals::new();
    /// matching.push(1, 1, 1);
    /// matching.push(2, 2, 2);
    /// assert_eq!( matching.bimap_min_ref().vec_ord_to_val(), &vec![1, 2] ); // the vector sending major ordinals to minor keys
    /// assert_eq!( matching.bimap_maj_ref().vec_ord_to_val(), &vec![1, 2] ); // the vector sending major ordinals to major keys
    /// assert_eq!( matching.vec_snzval_ref(), &vec![1, 2] ); // the vector sending major ordinals to linear coefficients
    /// matching.reverse();
    /// assert_eq!( matching.bimap_min_ref().vec_ord_to_val(), &vec![2, 1] ); // the vector sending major ordinals to minor keys
    /// assert_eq!( matching.bimap_maj_ref().vec_ord_to_val(), &vec![2, 1] ); // the vector sending major ordinals to major keys
    /// assert_eq!( matching.vec_snzval_ref(), &vec![2, 1] ); // the vector sending major ordinals to linear coefficients
    /// ```
    pub fn reverse( &mut self ) {
        self.bimap_min.reverse();    
        self.bimap_maj.reverse();
        self.vec_snzval.reverse();             
    }     
    
    /// The number of pairs stored in the array
    pub fn num_pairs(&self) -> usize { self.bimap_maj.len() }

    /// Returns an iterator that runs over all `(&minky, &majkey)` pairs; pairs occur in the usual `ord_to_val` order.
    pub fn minmajpairs( &self ) -> std::iter::Zip<std::slice::Iter<'_, KeyMin>, std::slice::Iter<'_, KeyMaj>> {
        let a = self.bimap_min.vec_ord_to_val().iter();                
        let b = self.bimap_maj.vec_ord_to_val().iter();
        return a.zip(b)
    }

    /// Returns a vector of all `(minky, majkey)` pairs; pairs occur in the usual `ord_to_val` order.
    pub fn minmajpairs_vec( &self ) -> Vec<(KeyMin, KeyMaj)> {
        let a = self.bimap_min.vec_ord_to_val().iter().cloned();                
        let b = self.bimap_maj.vec_ord_to_val().iter().cloned();
        return a.zip(b).collect()
    }    
        
    /// Returns an immutable reference to a `BiMapSequential< KeyMin >` object that stores the 
    /// bijection between minor keys and major ordinals.
    pub fn bimap_min_ref( &self ) -> & BiMapSequential< KeyMin > { & self.bimap_min }

    // DEPRECATED AS UNSAFE        
    // /// Returns a mutable reference to a `BiMapSequential< KeyMin >` object that stores the 
    // /// bijection between minor keys and major ordinals.
    // pub fn bimap_min_refmut( &mut self ) -> &mut BiMapSequential< KeyMin > { &mut self.bimap_min }    

    /// Returns an immutable reference to a `BiMapSequential< KeyMaj >` object that stores the 
    /// bijection between major keys and major ordinals.
    pub fn bimap_maj_ref( &self ) -> & BiMapSequential< KeyMaj > { & self.bimap_maj }    

    // DEPRECATED AS UNSAFE
    // /// Returns a mutable reference to a `BiMapSequential< KeyMaj >` object that stores the 
    // /// bijection between major keys and major ordinals.
    // pub fn bimap_maj_refmut( &mut self ) -> &mut BiMapSequential< KeyMaj > { &mut self.bimap_maj }        

    /// Returns an immutable reference to the sequence of structural nonzero values associated with the array.
    pub fn vec_snzval_ref( &self ) -> & Vec< SnzVal > { & self.vec_snzval }  

    /// A struct that implements `EvalutateFunction` mapping matched minor keys bijectively onto matched major keys.
    pub fn bijection_keymin_to_keymaj( &self ) -> ComposeFunctions< KeyMin, usize, KeyMaj, &HashMap<KeyMin, usize>, &Vec<KeyMaj> > {
        ComposeFunctions::new( 
                self.bimap_min_ref().hashmap_val_to_ord(),
                self.bimap_maj_ref().vec_ord_to_val()
            )
    }

    /// A struct that implements `EvalutateFunction` mapping matched major keys bijectively onto matched minor keys.
    pub fn bijection_keymaj_to_keymin( &self ) -> ComposeFunctions< KeyMaj, usize, KeyMin, &HashMap<KeyMaj, usize>, &Vec<KeyMin> > {
        ComposeFunctions::new( 
                self.bimap_maj_ref().hashmap_val_to_ord(),
                self.bimap_min_ref().vec_ord_to_val()
            )
    }    

    // DEPRECATED AS UNSAFE    
    // /// Returns an immutable reference to the sequence of structural nonzero values associated with the array.
    // pub fn vec_snzval_refmut( &mut self ) -> &mut Vec< SnzVal > { &mut self.vec_snzval }      

    /// Returns the major ordinal of the major key, or returns `None` if the major key is unmatched.
    pub fn keymaj_to_ordmaj( &self, keymaj: &KeyMaj ) -> Option< usize > { 
        self.bimap_maj.ord( keymaj )
    }

    /// Returns the major ordinal of the minor key, or returns `None` if the minor key is unmatched.
    pub fn keymin_to_ordmaj( &self, keymin: &KeyMin ) -> Option< usize > { 
        self.bimap_min.ord( keymin )
    }       

    /// Returns the major key corresponding to the major ordinal.
    pub fn ordmaj_to_keymaj( &self, ordmaj: usize ) -> KeyMaj { 
        self.bimap_maj.val( ordmaj )
    }

    /// Returns the minor key corresponding to the major ordinal.
    pub fn ordmaj_to_keymin( &self, ordmaj: usize ) -> KeyMin { 
        self.bimap_min.val( ordmaj )
    }           

    /// Returns the minor key of the major key, or returns `None` if the major key is unmatched.
    pub fn keymaj_to_keymin( &self, keymaj: &KeyMaj ) -> Option< KeyMin > { 
        match self.keymaj_to_ordmaj( keymaj ) {
            Some( ordmaj ) => { Some( self.ordmaj_to_keymin( ordmaj ) ) },
            None => None
        }
    }

    /// Returns the major key of the minor key, or returns `None` if the major key is unmatched.
    pub fn keymin_to_keymaj( &self, keymin: &KeyMin ) -> Option< KeyMaj > { 
        match self.keymin_to_ordmaj( keymin ) {
            Some( ordmaj ) => { Some( self.ordmaj_to_keymaj( ordmaj ) ) },
            None => None
        }
    }    
 
    /// The nonzero snzval associated with this `keymaj` index.
    pub fn keymaj_to_snzval( &self, keymaj: & KeyMaj ) -> SnzVal 
        where SnzVal: Clone    
    {
        self.vec_snzval[ 
                self.keymaj_to_ordmaj( keymaj ).unwrap()                
            ]            
            .clone()
    }

    /// The nonzero snzval associated with this `keymin` index.
    pub fn keymin_to_snzval( &self, keymin: & KeyMin ) -> SnzVal 
        where SnzVal: Clone    
    {
        self.vec_snzval[ 
                self.keymin_to_ordmaj( keymin ).unwrap() 
            ]
            .clone()
    }    

    /// The nonzero snzval associated with the `ordmaj`th major key.
    pub fn ordmaj_to_snzval( &self, ordmaj: usize ) -> SnzVal 
        where SnzVal: Clone
    {
        self.vec_snzval[ ordmaj ]
            .clone()
    }    

    /// Returns true iff the given minor key is matched.
    pub fn contains_keymin(  &self, keymin: &KeyMin ) -> bool {
        self.bimap_min.contains_key( keymin )
    }

    /// Returns true iff the given major key is matched.
    pub fn contains_keymaj(  &self, keymaj: &KeyMaj ) -> bool {
        self.bimap_maj.contains_key( keymaj )
    }   
    
    /// Returns true iff the given minor key is unmatched.
    pub fn lacks_keymin(  &self, keymin: &KeyMin ) -> bool {
        ! self.bimap_min.contains_key( keymin )
    }
 
    /// Returns true iff the given major key is unmatched.
    pub fn lacks_keymaj(  &self, keymaj: &KeyMaj ) -> bool {
        ! self.bimap_maj.contains_key( keymaj )
    }       

    // THE FOLLOWING METHODS ARE DEPRECATED; HOWEVER WE MAY REVIVE THEM IF WE EVER NEED AN ORDINAL MATRCHING ARRAY WITH ADDITIONAL STRUCTURE        

    // /// Returns the ordinal of the minor key associated with the given major key.
    // /// Panicks if the minor key is not associated with a major key.    
    // pub fn keymaj_to_ordmin( &self, keymaj: &KeyMaj ) -> usize { 
    //     self.ordmaj_to_ordmin[
    //             * self.keymaj_to_ordmaj.get( keymaj ).unwrap()
    //         ]
    // }

    // /// Returns the ordinal of the major key associated with the given minor key.
    // /// Panicks if the major key is not associated with a minor key.
    // pub fn keymin_to_ordmaj( &self, keymin: &KeyMin ) -> usize { 
    //     self.ordmin_to_ordmaj[
    //             * self.keymin_to_ordmin.get( keymin ).unwrap()
    //         ]
    // }    

    // /// Returns either `Some(k)`, where `k` is the ordinal of the major key associated with the given minor key,
    // /// or `None`, if there exists no associated major key.
    // pub fn keymin_to_ordmaj_opt( &self, keymin: &KeyMin ) -> Option< usize > { 
    //     match self.keymin_to_ordmin.get( keymin ) {
    //         Some( & ordmin ) => { Some( self.ordmin_to_ordmaj[ ordmin ].clone() ) },
    //         None => None
    //     }
    // } 
    
    // /// Returns either `Some(k)`, where `k` is the ordinal of the major key associated with the given minor key,
    // /// or `None`, if there exists no associated major key.
    // pub fn keymaj_to_ordmin_opt( &self, keymaj: &KeyMaj ) -> Option< usize > { 
    //     match self.keymaj_to_ordmaj.get( keymaj ) {
    //         Some( & ordmaj ) => { Some( self.ordmaj_to_ordmin[ ordmaj ].clone() ) },
    //         None => None
    //     }
    // }    

    // pub fn keymaj_to_keymin( &self, keymaj: & KeyMaj ) -> & KeyMin { 
    //     & self.ordmin_to_keymin[ 
    //             self.keymaj_to_ordmin( keymaj )
    //         ]        
    // }

    // pub fn keymin_to_keymaj( &self, keymin: & KeyMin ) -> & KeyMaj { 
    //     & self.ordmaj_to_keymaj[ 
    //             self.keymin_to_ordmaj( keymin )
    //         ]
    // }    

    // /// The nonzero snzval associated with this `keymaj` index.
    // pub fn keymaj_to_snzval( &self, keymaj: & KeyMaj ) -> & SnzVal {
    //     & self.ordmaj_to_snzval[ * self.keymaj_to_ordmaj.get( keymaj ).unwrap() ]
    // }

    // /// The nonzero snzval associated with this `keymin` index.
    // pub fn keymin_to_snzval( &self, keymin: & KeyMin ) -> & SnzVal {
    //     & self.ordmaj_to_snzval[ self.keymin_to_ordmaj( keymin ) ]
    // }    
}


impl < 'a, KeyBoth, SnzVal > 

    GeneralizedMatchingArrayWithMajorOrdinals
    
    < KeyBoth, KeyBoth, SnzVal >
    where   KeyBoth:     Clone + Hash + std::cmp::Eq, // without cloning, it is difficult to write a method that returns the minkey at a given index (one gets "cannot move" errors)
{
    /// Returns `true` iff `key` is a matched major or minor key.
    pub fn contains_key( &self, key: &KeyBoth ) -> bool { self.contains_keymin(key) || self.contains_keymaj(key) }

    /// Returns `true` iff `key` is a matched major or minor key.
    pub fn lacks_key( &self, key: &KeyBoth ) -> bool { self.lacks_keymin(key) && self.lacks_keymaj(key) } 
    
    /// A binned count of unmatched items.
    /// 
    /// Concretely, this means vector `v` such that `v[i]` is the number of umatched items 
    /// in the iterator with rank `i` (where rank is determined by `rank_fn`).
    pub fn unmatched_histo< I, F >( &self, iter_keyboth: I, mut rank_fn: F ) 
            -> 
            Vec< usize > 
        where
            I: Iterator<Item=KeyBoth>, 
            F: FnMut(KeyBoth)->usize        
            
        {
        let mut histo = Vec::new();
        let mut rank;
        for key in iter_keyboth { 
            if self.lacks_key( &key ) { 
                rank = rank_fn(key);
                while histo.len() < rank + 1 { histo.push(0) }
                histo[rank] +=1
            } 
        }
        return histo
    }

}


//  EXTEND TRAIT IMPLEMENTATION
//  ------------------------------------------------------------------------------

impl    < 'a, KeyMin, KeyMaj, SnzVal > 
        Extend< ( KeyMin, KeyMaj, SnzVal ) > 
        for 
        GeneralizedMatchingArrayWithMajorOrdinals < KeyMin, KeyMaj, SnzVal > 
    where   KeyMin:     Clone + Hash + std::cmp::Eq,
            KeyMaj:     Clone + Hash + std::cmp::Eq,
            SnzVal:     Clone,
{
    /// Add a nonzero entry to a matching array.
    fn extend< I: IntoIterator<Item= ( KeyMin, KeyMaj, SnzVal ) >>( &mut self, iter: I) {
        for x in iter.into_iter() { self.push( x.0, x.1, x.2 ) }
    }
}


//  FROM_ITER TRAIT IMPLEMENTATION
//  ------------------------------------------------------------------------------

impl    < 'a, KeyMin, KeyMaj, SnzVal > 
        
        FromIterator
        < ( KeyMin, KeyMaj, SnzVal ) > for 
        
        GeneralizedMatchingArrayWithMajorOrdinals 
        < KeyMin, KeyMaj, SnzVal > 

    where   KeyMin:     Clone + Hash + std::cmp::Eq,
            KeyMaj:     Clone + Hash + std::cmp::Eq,
            SnzVal:     Clone,
{
    fn from_iter< I: IntoIterator<Item= ( KeyMin, KeyMaj, SnzVal ) >>(iter: I) -> Self {

        let mut matching_array  =   GeneralizedMatchingArrayWithMajorOrdinals::new();

        for x in iter.into_iter() { matching_array.push( x.0, x.1, x.2 ) }
        
        matching_array
    }
}


//  ORACLE TRAIT IMPLEMENTATIONS
//  ------------------------------------------------------------------------------


// impl < 'a, KeyMin, KeyMaj, SnzVal >

//     OracleRefInherit for
    
//     &'a GeneralizedMatchingArrayWithMajorOrdinals
//         < KeyMin, KeyMaj, SnzVal >

//     where   KeyMin:     Clone + Hash + std::cmp::Eq,
//             KeyMaj:     Clone + Hash + std::cmp::Eq,
//             SnzVal:     Clone,        
// {}


impl < 'a, KeyMin, KeyMaj, SnzVal >

    IndicesAndCoefficients for 

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < KeyMin, KeyMaj, SnzVal >

    where   KeyMin:     Clone + Hash + std::cmp::Eq,
            KeyMaj:     Clone + Hash + std::cmp::Eq,
            SnzVal:     Clone,
{ type KeyMin = KeyMin; type KeyMaj = KeyMaj; type SnzVal = SnzVal; }


impl < 'a, KeyMin, KeyMaj, SnzVal >
     
    OracleMajor for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < KeyMin, KeyMaj, SnzVal >
    
    where   KeyMin:     Clone + Hash + std::cmp::Eq,
            KeyMaj:     Clone + Hash + std::cmp::Eq,
            SnzVal:     Clone,
{    
    type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;
    type ViewMajorEntry     =   ( KeyMin, SnzVal );    
    type ViewMajor          =   IterTwoType<
                                        Empty< ( KeyMin, SnzVal ) >,
                                        std::iter::Once< ( KeyMin, SnzVal ) >,
                                    > ;

    fn view_major( &self, index: Self::KeyMaj ) -> Self::ViewMajor
    {
        match self.keymaj_to_ordmaj( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ordmaj ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ordmaj_to_keymin( ordmaj ),
                                                    self.ordmaj_to_snzval( ordmaj ),  
                                                )          
                            )
                    )
            }
        }
    }
}

impl < 'a, KeyMin, KeyMaj, SnzVal >
     
    OracleMajorAscend for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < KeyMin, KeyMaj, SnzVal >
    
    where   KeyMin:     Clone + Hash + std::cmp::Eq,
            KeyMaj:     Clone + Hash + std::cmp::Eq,
            SnzVal:     Clone,
{    
    type ViewMajorAscendIntoIter = < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscendEntry = ( KeyMin, SnzVal );    
    type ViewMajorAscend =  IterTwoType<
                                    Empty< ( KeyMin, SnzVal ) >,
                                    std::iter::Once< ( KeyMin, SnzVal ) >,
                                > ;

    fn view_major_ascend( &self, index: Self::KeyMaj ) -> Self::ViewMajorAscend
    {
        match self.keymaj_to_ordmaj( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ordmaj ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ordmaj_to_keymin( ordmaj ),
                                                    self.ordmaj_to_snzval( ordmaj ),  
                                                )          
                            )
                    )
            }
        }
    }
}

impl < 'a, KeyMin, KeyMaj, SnzVal >
     
    OracleMajorDescend for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < KeyMin, KeyMaj, SnzVal >
    
    where   KeyMin:     Clone + Hash + std::cmp::Eq,
            KeyMaj:     Clone + Hash + std::cmp::Eq,
            SnzVal:     Clone,
{    
    type ViewMajorDescendIntoIter   =     < Self::ViewMajorDescend as IntoIterator >::IntoIter;
    type ViewMajorDescendEntry      = ( KeyMin, SnzVal );    
    type ViewMajorDescend           =  IterTwoType<
                                                Empty< ( KeyMin, SnzVal ) >,
                                                std::iter::Once< ( KeyMin, SnzVal ) >,
                                            > ;

    fn view_major_descend( &self, index: Self::KeyMaj ) -> Self::ViewMajorDescend
    {
        match self.keymaj_to_ordmaj( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ordmaj ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ordmaj_to_keymin( ordmaj ),
                                                    self.ordmaj_to_snzval( ordmaj ),  
                                                )          
                            )
                    )
            }
        }
    }
}







//  ================================================================================
//  GENERALIZED MATCHING ARRAY
//  ================================================================================


#[derive(Debug, Clone)]
pub struct GeneralizedMatchingArray
                < KeyMin, KeyMaj, EntryMaj, EntryMin, SnzVal >
    where   KeyMin:     Clone + Hash + std::cmp::Eq,  // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            KeyMaj:     Clone + Hash + std::cmp::Eq,  // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            EntryMin:   Clone + KeyValGet< KeyMin, SnzVal > + KeyValNew< KeyMin, SnzVal >, // clone is ncessary in order for HashMap< KeyMaj, EntryMin > to implement `EvaluateFunction< KeyMaj, EntryMin >`
            EntryMaj:   Clone + KeyValGet< KeyMaj, SnzVal > + KeyValNew< KeyMaj, SnzVal >, // clone is ncessary in order for HashMap< KeyMin, EntryMaj > to implement `EvaluateFunction< KeyMin, EntryMaj >`
            SnzVal:     Clone,                       
{  
    keymin_to_entry:    HashMap< KeyMin, EntryMaj >,
    keymaj_to_entry:    HashMap< KeyMaj, EntryMin >,
    phantom_snzval:     PhantomData< SnzVal >
}


impl < 'a, KeyMin, KeyMaj, EntryMin, EntryMaj, SnzVal > 

    GeneralizedMatchingArray
        < KeyMin, KeyMaj, EntryMaj, EntryMin, SnzVal >

    where   KeyMin:     Clone + Hash + std::cmp::Eq, // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            KeyMaj:     Clone + Hash + std::cmp::Eq, // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            EntryMin:   Clone + KeyValGet< KeyMin, SnzVal > + KeyValNew< KeyMin, SnzVal >,
            EntryMaj:   Clone + KeyValGet< KeyMaj, SnzVal > + KeyValNew< KeyMaj, SnzVal >,   
            SnzVal:     Clone, 
{
    /// Create a new, empty matching array.
    /// 
    /// See [`GeneralizedMatchingArrayWithMajorOrdinals`] for usage.
    pub fn new() -> GeneralizedMatchingArray< KeyMin, KeyMaj, EntryMaj, EntryMin, SnzVal > 
    {
        GeneralizedMatchingArray{    
            keymin_to_entry:    HashMap::new(),
            keymaj_to_entry:    HashMap::new(),
            phantom_snzval:     PhantomData,
        }        
    }

    /// Push a new structural nonzero entry to the array.    
    pub fn push( 
                    &mut self,
                    keymin:     KeyMin, 
                    keymaj:     KeyMaj, 
                    snzval:     SnzVal 
                )
    {
        let a = self.keymin_to_entry.insert( keymin.clone(), EntryMaj::new( keymaj.clone(), snzval.clone() ) ); 
        let b = self.keymaj_to_entry.insert( keymaj, EntryMin::new( keymin, snzval ) );               
        if ! a.is_none()  { panic!("attempted to over-write the entry matched to a minor key in a generalized matching matrix")  }
        if ! b.is_none()  { panic!("attempted to over-write the entry matched to a major key in a generalized matching matrix")  }
    }

        
    /// Returns an object `X` such that `X.evaluate_function( keymin ) = keymaj`, where `keymin` is a matched minor key
    /// and `keymaj` is the associated major key; panics when the user calls `X.evaluate_function( keymin )` on
    /// an unmatched index `keymin`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_keymaj_to_keymin( &self )` works by cloning the entry associated with `keymaj` then
    /// extracting the associated index.  This is not always the most efficient way to get the minor key matched
    /// to a major key (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_keymin_to_keymaj( &self ) 
            -> 
            ComposeFunctions< 
                    KeyMin, &EntryMaj, KeyMaj, 
                    & HashMap< KeyMin, EntryMaj >,                     
                    ExtractKey< KeyMaj, SnzVal > 
                > 
    {
        ComposeFunctions::new( & self.keymin_to_entry, ExtractKey::new() )
    }

    /// Returns an object `X` such that `X.evaluate_function( keymaj ) = keymin`, where `keymaj` is a matched major key
    /// and `keymin` is the associated minor key; panics when the user calls `X.evaluate_function( keymaj )` on
    /// an unmatched index `keymaj`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_keymaj_to_keymin( &self )` works by cloning the entry associated with `keymaj` then
    /// extracting the associated index.  This is not always the most efficient way to get the minor key matched
    /// to a major key (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_keymaj_to_keymin( &self ) 
            -> 
            ComposeFunctions< 
                    KeyMaj, & EntryMin, KeyMin, 
                    & HashMap< KeyMaj, EntryMin >,                     
                    ExtractKey< KeyMin, SnzVal > 
                > 
    {
        ComposeFunctions::new( & self.keymaj_to_entry, ExtractKey::new() )
    }    

    /// Returns an object `X` such that `X.evaluate_function( keymin ) = snzval`, where `keymin` is a matched minor key
    /// and `snzval` is the scalar value of the associated nonzero entry; panics when the user calls `X.evaluate_function( keymin )` on
    /// an unmatched index `keymin`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_keymaj_to_keymin( &self )` works by cloning the entry associated with `keymaj` then
    /// extracting the associated scalar.  This is not always the most efficient way to get the minor key matched
    /// to a major key (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_keymin_to_snzval( &self ) 
            -> 
            ComposeFunctions< 
                    KeyMin, & EntryMaj, SnzVal, 
                    & HashMap< KeyMin, EntryMaj >,                     
                    ExtractVal< KeyMaj, SnzVal > 
                > 
    {
        ComposeFunctions::new( & self.keymin_to_entry, ExtractVal::new() )
    }

    /// Returns an object `X` such that `X.evaluate_function( keymaj ) = snzval`, where `keymaj` is a matched major key
    /// and `snzval` is the scalar value of the associated nonzero entry; panics when the user calls `X.evaluate_function( keymaj )` on
    /// an unmatched index `keymaj`.    
    /// 
    /// # Design note
    /// 
    /// The function `f = function_keymaj_to_keymin( &self )` works by cloning the entry associated with `keymaj` then
    /// extracting the associated scalar.  This is not always the most efficient way to get the scalar associated
    /// to a major key (some data structures allow more efficient work-arounds), but (i) it is as efficient as we have
    /// thus far been able to make it, without making limiting assumptions about the data structures involved, and (ii) 
    /// we don't currently foresee that this function will be used heavily.
    pub fn function_keymaj_to_snzval( &self ) 
            -> 
            ComposeFunctions< 
                    &KeyMaj, & EntryMin, SnzVal, 
                    &HashMap< KeyMaj, EntryMin >,                     
                    ExtractVal< KeyMin, SnzVal > 
                > 
    {
        ComposeFunctions::new( & self.keymaj_to_entry, ExtractVal::new() )
    }    




    
    /// Returns an immutable reference `& HashMap< KeyMin, EntryMaj >` representing the partial bijection from 
    /// minor keys to pairs of form `(associated_major_key, associated_nonzero_scalar)`.
    pub fn hashmap_ref_keymin_to_entryymaj( &self ) -> & HashMap< KeyMin, EntryMaj > { & self.keymin_to_entry }

    /// Returns an immutable reference `& HashMap< KeyMaj, EntryMaj >` representing the partial bijection from 
    /// major keys to pairs of form `(associated_minor_key, associated_nonzero_scalar)`.
    pub fn hashmap_ref_keymaj_to_entryymin( &self ) -> & HashMap< KeyMaj, EntryMin > { & self.keymaj_to_entry }




    /// Returns the nonzero entry that corresponds to the major key, or returns `None` if the major key is unmatched.
    pub fn keymaj_to_entrymin_refopt( &self, keymaj: &KeyMaj ) -> Option< & EntryMin > { 
        self.keymaj_to_entry.get( keymaj )
    }

    /// Returns the nonzero entry that corresponds to the minor key, or returns `None` if the minor key is unmatched.
    pub fn keymin_to_entrymaj_refopt( &self, keymin: &KeyMin ) -> Option< & EntryMaj > { 
        self.keymin_to_entry.get( keymin )
    }    
 
    /// The nonzero snzval associated with this `keymaj` index.
    pub fn keymaj_to_snzval_opt( &self, keymaj: & KeyMaj ) -> Option< SnzVal >
        // where SnzVal: Clone    
    {
        self.keymaj_to_entry.get( keymaj ).map( |x| x.val() )
    }

    /// The nonzero snzval associated with this `keymin` index.
    pub fn keymin_to_snzval_opt( &self, keymin: & KeyMin ) -> Option< SnzVal >
        // where SnzVal: Clone    
    {
        self.keymin_to_entry.get( keymin ).map( |x| x.val() )
    }    

    /// Returns true iff the given minor key is matched.
    pub fn contains_keymin(  &self, keymin: &KeyMin ) -> bool {
        self.keymin_to_entry.contains_key( keymin )
    }

    /// Returns true iff the given major key is matched.
    pub fn contains_keymaj(  &self, keymaj: &KeyMaj ) -> bool {
        self.keymaj_to_entry.contains_key( keymaj )
    }    



}











//  ================================================================================
//  WRAPPER STRUCT FOR MAPPING KEYMIN TO KEYMAJ
//  ================================================================================


//  UNDER CONSTRUCTION; THIS HAS NOT YET BEEN NECESSARY















//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Docstring tests (draft)
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod docstring_tests_draft {
    
    #[test]
    fn test_matching_array() {
        use super::GeneralizedMatchingArrayWithMajorOrdinals;

        let mut matching_array  =   GeneralizedMatchingArrayWithMajorOrdinals::new();
        matching_array.push( 0, 1, 3 );
        matching_array.push( 1, 2, 4 );
        matching_array.push( 2, 0, 5 );    
        
        assert_eq!( matching_array.vec_snzval_ref(), &vec![ 3, 4, 5 ] );

        assert_eq!( matching_array.keymin_to_keymaj( &0 ), Some(1) );
        assert_eq!( matching_array.keymaj_to_keymin( &0 ), Some(2) );     

        assert_eq!( matching_array.keymin_to_ordmaj( &1 ), Some(1) );
        assert_eq!( matching_array.keymaj_to_ordmaj( &1 ), Some(0) ); 

        assert_eq!( matching_array.ordmaj_to_keymin( 0 ), 0 );
        assert_eq!( matching_array.ordmaj_to_keymaj( 0 ), 1 );     

        assert_eq!( matching_array.keymin_to_snzval( &0 ), 3 );
        assert_eq!( matching_array.keymaj_to_snzval( &0 ), 5 );       
        
        assert_eq!( matching_array.ordmaj_to_snzval( 0 ), 3 );
    }    
}



