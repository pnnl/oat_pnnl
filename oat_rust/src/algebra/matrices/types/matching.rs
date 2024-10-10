//! (Generalized) matching matrices



use crate::algebra::vectors::entries::{KeyValGet, KeyValNew, ExtractKey, ExtractVal};
use crate::algebra::matrices::query::{ViewRow, ViewRowAscend, ViewRowDescend, IndicesAndCoefficients, ViewCol, ViewColAscend, ViewColDescend};
use crate::utilities::functions::compose::ComposeFunctions;

use crate::utilities::functions::evaluate::{LogicalNot, Map};
use crate::utilities::iterators::general::{IterTwoType, Filter, FilterOutMembers, FilterOutNonmembers,};
use crate::utilities::sequences_and_ordinals::BiMapSequential;
use crate::utilities::sets::MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper;
use std::collections::HashMap;
use std::hash::Hash;
use std::iter::{Empty, FromIterator};
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
/// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
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
/// assert_eq!( matching_array.keymaj_to_keymin( &0 ), Some(1) );  
/// assert_eq!( matching_array.keymin_to_keymaj( &0 ), Some(2) );   
/// 
/// assert_eq!( matching_array.keymaj_to_ord( &1 ), Some(1) ); 
/// assert_eq!( matching_array.keymin_to_ord( &1 ), Some(0) );
/// 
/// assert_eq!( matching_array.ord_to_keymaj( 0 ), 0 );     
/// assert_eq!( matching_array.ord_to_keymin( 0 ), 1 );
/// 
/// assert_eq!( matching_array.keymaj_to_snzval( &0 ), 3 );       
/// assert_eq!( matching_array.keymin_to_snzval( &0 ), 5 );
/// 
/// assert_eq!( matching_array.ord_to_snzval( 0 ), 3 );
/// ```
#[derive(Debug, Clone)]
pub struct GeneralizedMatchingArrayWithMajorOrdinals< ColIndex, RowIndex, Coefficient >
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
{  
    bimap_min:  BiMapSequential< ColIndex >,
    bimap_maj:  BiMapSequential< RowIndex >,
    vec_snzval: Vec< Coefficient >,
    // keymaj_2_ord: HashMap< RowIndex, usize >,
    // keymin_2_ord: HashMap< ColIndex, usize >,
    // ord_2_keymaj: Vec< RowIndex >,
    // ord_2_keymin: Vec< ColIndex >,    
    // ord_2_snzval: Vec< Coefficient >,
    // phantom: PhantomData< &'a usize >,
}

impl < 'a, ColIndex, RowIndex, Coefficient > 

    GeneralizedMatchingArrayWithMajorOrdinals
    
    < ColIndex, RowIndex, Coefficient >
    where   ColIndex:     Clone + Hash + std::cmp::Eq, // without cloning, it is difficult to write a method that returns the minkey at a given index (one gets "cannot move" errors)
            RowIndex:     Clone + Hash + std::cmp::Eq, // without cloning, it is difficult to write a method that returns the majkey at a given index (one gets "cannot move" errors)            
    
{
    /// Create a new, empty matching array.
    /// 
    /// See [`GeneralizedMatchingArrayWithMajorOrdinals`] for usage.
    pub fn new() -> GeneralizedMatchingArrayWithMajorOrdinals< ColIndex, RowIndex, Coefficient > 
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
                    keymaj:     RowIndex,                     
                    keymin:     ColIndex, 
                    snzval:     Coefficient 
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
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
    /// 
    /// let mut matching = GeneralizedMatchingArrayWithMajorOrdinals::new();
    /// matching.push(1, 1, 1.0);
    /// matching.push(2, 2, 2.0);
    /// assert_eq!( matching.bimap_min_ref().ord_to_val_vec(), &vec![1, 2] ); // the vector sending major ordinals to minor keys
    /// assert_eq!( matching.bimap_maj_ref().ord_to_val_vec(), &vec![1, 2] ); // the vector sending major ordinals to major keys
    /// assert_eq!( matching.vec_snzval_ref(), &vec![1.0, 2.0] ); // the vector sending major ordinals to linear coefficients
    /// matching.reverse();
    /// assert_eq!( matching.bimap_min_ref().ord_to_val_vec(), &vec![2, 1] ); // the vector sending major ordinals to minor keys
    /// assert_eq!( matching.bimap_maj_ref().ord_to_val_vec(), &vec![2, 1] ); // the vector sending major ordinals to major keys
    /// assert_eq!( matching.vec_snzval_ref(), &vec![2.0, 1.0] ); // the vector sending major ordinals to linear coefficients
    /// ```
    pub fn reverse( &mut self ) {
        self.bimap_min.reverse();    
        self.bimap_maj.reverse();
        self.vec_snzval.reverse();             
    }     
    
    /// The number of pairs stored in the array
    pub fn num_pairs(&self) -> usize { self.bimap_maj.len() }

    /// Returns an iterator that runs over all `(&RowIndex, &ColIndex)` pairs; pairs occur in the usual `ord_to_val` order.
    pub fn support( &self ) -> std::iter::Zip<std::slice::Iter<'_, RowIndex>, std::slice::Iter<'_, ColIndex>> {
        let b = self.bimap_min.ord_to_val_vec().iter();                
        let a = self.bimap_maj.ord_to_val_vec().iter();
        a.zip(b)
    }

    /// Filters out matched minor indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
    /// 
    /// let mut matching = GeneralizedMatchingArrayWithMajorOrdinals::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_out_matched_minors( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![0,2] )  );
    /// ```
    pub fn filter_out_matched_minors< I: IntoIterator< Item = ColIndex > >( &self, iter: I ) 
        -> 
        FilterOutMembers< I::IntoIter, & HashMap< ColIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
            LogicalNot::new(
                    MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                        self.bimap_min.val_to_ord_hashmap()
                    )
                ) 
            )
    }    

    /// Filters out matched major indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
    /// 
    /// let mut matching = GeneralizedMatchingArrayWithMajorOrdinals::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_out_matched_majors( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![1,3] )  );
    /// ```
    pub fn filter_out_matched_majors< I: IntoIterator< Item = RowIndex > >( &self, iter: I ) 
        -> 
        FilterOutMembers< I::IntoIter, & HashMap< RowIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
            LogicalNot::new(
                        MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                                self.bimap_maj.val_to_ord_hashmap()
                            )
                    ) 
                )
    }   


    /// Filters out unmatched minor indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
    /// 
    /// let mut matching = GeneralizedMatchingArrayWithMajorOrdinals::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_out_unmatched_minors( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![1,3] )  );
    /// ```
    pub fn filter_out_unmatched_minors< I: IntoIterator< Item = ColIndex > >( &self, iter: I ) 
        ->
        FilterOutNonmembers< I::IntoIter, & HashMap< ColIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
            MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                        self.bimap_min.val_to_ord_hashmap()
                    )
                )
    }    

    /// Filters out unmatched major indices
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::matching::GeneralizedMatchingArrayWithMajorOrdinals;
    /// 
    /// let mut matching = GeneralizedMatchingArrayWithMajorOrdinals::new();
    /// matching.push(0, 1, 5.0);
    /// matching.push(2, 3, 7.0);
    /// 
    /// let filtered = matching.filter_only_matched_majors( vec![0,1,2,3] );
    /// assert!(  filtered.eq( vec![0,2] )  );
    /// ```
    pub fn filter_only_matched_majors< I: IntoIterator< Item = RowIndex > >( &self, iter: I ) 
        -> 
        FilterOutNonmembers< I::IntoIter, & HashMap< RowIndex, usize > >
    {
        Filter::new( 
            iter.into_iter(), 
                MapHasKeyOrSequenceHasElementEvaluateFunctionWrapper::new(
                    self.bimap_maj.val_to_ord_hashmap()
                ) 
            )
    }           
        
    /// Returns an immutable reference to a `BiMapSequential< ColIndex >` object that stores the 
    /// bijection between minor keys and major ordinals.
    pub fn bimap_min_ref( &self ) -> & BiMapSequential< ColIndex > { & self.bimap_min }

    // DEPRECATED AS UNSAFE        
    // /// Returns a mutable reference to a `BiMapSequential< ColIndex >` object that stores the 
    // /// bijection between minor keys and major ordinals.
    // pub fn bimap_min_refmut( &mut self ) -> &mut BiMapSequential< ColIndex > { &mut self.bimap_min }    

    /// Returns an immutable reference to a `BiMapSequential< RowIndex >` object that stores the 
    /// bijection between major keys and major ordinals.
    pub fn bimap_maj_ref( &self ) -> & BiMapSequential< RowIndex > { & self.bimap_maj }    

    // DEPRECATED AS UNSAFE
    // /// Returns a mutable reference to a `BiMapSequential< RowIndex >` object that stores the 
    // /// bijection between major keys and major ordinals.
    // pub fn bimap_maj_refmut( &mut self ) -> &mut BiMapSequential< RowIndex > { &mut self.bimap_maj }        

    /// Returns an immutable reference to the sequence of structural nonzero values associated with the array.
    pub fn vec_snzval_ref( &self ) -> & Vec< Coefficient > { & self.vec_snzval }  

    /// A struct that implements `EvalutateFunction` mapping matched minor keys bijectively onto matched major keys.
    pub fn bijection_keymin_to_keymaj( &self ) -> ComposeFunctions< ColIndex, usize, RowIndex, &HashMap<ColIndex, usize>, &Vec<RowIndex> > {
        ComposeFunctions::new( 
                self.bimap_min_ref().val_to_ord_hashmap(),
                self.bimap_maj_ref().ord_to_val_vec()
            )
    }

    /// A struct that implements `EvalutateFunction` mapping matched major keys bijectively onto matched minor keys.
    pub fn bijection_keymaj_to_keymin( &self ) -> ComposeFunctions< RowIndex, usize, ColIndex, &HashMap<RowIndex, usize>, &Vec<ColIndex> > {
        ComposeFunctions::new( 
                self.bimap_maj_ref().val_to_ord_hashmap(),
                self.bimap_min_ref().ord_to_val_vec()
            )
    } 

    /// A struct that implements `EvalutateFunction` mapping matched minor keys bijectively onto matched major keys.
    pub fn bijection_keymin_to_keymaj_opt( &self ) 
        -> ComposeFunctions< 
                ColIndex, 
                Option< usize >, 
                Option< RowIndex >, 
                &HashMap<ColIndex, usize>, 
                Map< &Vec<RowIndex> >,
            > 
    {
        ComposeFunctions::new( 
                self.bimap_min_ref().val_to_ord_hashmap(),
                Map{ mapping_rule: self.bimap_maj_ref().ord_to_val_vec() }
            )
    }

    /// A struct that implements `EvalutateFunction` mapping matched major keys bijectively onto matched minor keys.
    pub fn bijection_keymaj_to_keymin_opt( &self ) 
        -> ComposeFunctions< 
                RowIndex, 
                Option< usize >, 
                Option< ColIndex >, 
                &HashMap<RowIndex, usize>, 
                Map< &Vec<ColIndex> >,
            > {
        ComposeFunctions::new( 
                self.bimap_maj_ref().val_to_ord_hashmap(),
                Map{ mapping_rule: self.bimap_min_ref().ord_to_val_vec() }
            )
    }        

    // DEPRECATED AS UNSAFE    
    // /// Returns an immutable reference to the sequence of structural nonzero values associated with the array.
    // pub fn vec_snzval_refmut( &mut self ) -> &mut Vec< Coefficient > { &mut self.vec_snzval }      

    /// Returns the major ordinal of the major key, or returns `None` if the major key is unmatched.
    pub fn keymaj_to_ord( &self, keymaj: &RowIndex ) -> Option< usize > { 
        self.bimap_maj.ord( keymaj )
    }

    /// Returns the major ordinal of the minor key, or returns `None` if the minor key is unmatched.
    pub fn keymin_to_ord( &self, keymin: &ColIndex ) -> Option< usize > { 
        self.bimap_min.ord( keymin )
    }       

    /// Returns the major key corresponding to the major ordinal.
    pub fn ord_to_keymaj( &self, ord: usize ) -> RowIndex { 
        self.bimap_maj.val( ord )
    }

    /// Returns the minor key corresponding to the major ordinal.
    pub fn ord_to_keymin( &self, ord: usize ) -> ColIndex { 
        self.bimap_min.val( ord )
    }           

    /// Returns the minor key of the major key, or returns `None` if the major key is unmatched.
    pub fn keymaj_to_keymin( &self, keymaj: &RowIndex ) -> Option< ColIndex > { 
        self.keymaj_to_ord( keymaj ).map(|ord| self.ord_to_keymin( ord ))
    }

    /// Returns the major key of the minor key, or returns `None` if the major key is unmatched.
    pub fn keymin_to_keymaj( &self, keymin: &ColIndex ) -> Option< RowIndex > { 
        self.keymin_to_ord( keymin ).map(|ord| self.ord_to_keymaj( ord ))
    }    
 
    /// The nonzero snzval associated with this `keymaj` index.
    pub fn keymaj_to_snzval( &self, keymaj: & RowIndex ) -> Coefficient 
        where Coefficient: Clone    
    {
        self.vec_snzval[ 
                self.keymaj_to_ord( keymaj ).unwrap()                
            ]            
            .clone()
    }

    /// The nonzero snzval associated with this `keymin` index.
    pub fn keymin_to_snzval( &self, keymin: & ColIndex ) -> Coefficient 
        where Coefficient: Clone    
    {
        self.vec_snzval[ 
                self.keymin_to_ord( keymin ).unwrap() 
            ]
            .clone()
    }    

    /// The nonzero snzval associated with the `ord`th major key.
    pub fn ord_to_snzval( &self, ord: usize ) -> Coefficient 
        where Coefficient: Clone
    {
        self.vec_snzval[ ord ]
            .clone()
    }    

    /// Returns true iff the given minor key is matched.
    pub fn contains_keymin(  &self, keymin: &ColIndex ) -> bool {
        self.bimap_min.contains_key( keymin )
    }

    /// Returns true iff the given major key is matched.
    pub fn contains_keymaj(  &self, keymaj: &RowIndex ) -> bool {
        self.bimap_maj.contains_key( keymaj )
    }   
    
    /// Returns true iff the given minor key is unmatched.
    pub fn lacks_keymin(  &self, keymin: &ColIndex ) -> bool {
        ! self.bimap_min.contains_key( keymin )
    }
 
    /// Returns true iff the given major key is unmatched.
    pub fn lacks_keymaj(  &self, keymaj: &RowIndex ) -> bool {
        ! self.bimap_maj.contains_key( keymaj )
    }       

    // THE FOLLOWING METHODS ARE DEPRECATED; HOWEVER WE MAY REVIVE THEM IF WE EVER NEED AN ORDINAL MATRCHING ARRAY WITH ADDITIONAL STRUCTURE        

    // /// Returns the ordinal of the minor key associated with the given major key.
    // /// Panicks if the minor key is not associated with a major key.    
    // pub fn keymaj_to_ordmin( &self, keymaj: &RowIndex ) -> usize { 
    //     self.ord_to_ordmin[
    //             * self.keymaj_to_ord.get( keymaj ).unwrap()
    //         ]
    // }

    // /// Returns the ordinal of the major key associated with the given minor key.
    // /// Panicks if the major key is not associated with a minor key.
    // pub fn keymin_to_ord( &self, keymin: &ColIndex ) -> usize { 
    //     self.ordmin_to_ord[
    //             * self.keymin_to_ordmin.get( keymin ).unwrap()
    //         ]
    // }    

    // /// Returns either `Some(k)`, where `k` is the ordinal of the major key associated with the given minor key,
    // /// or `None`, if there exists no associated major key.
    // pub fn keymin_to_ord_opt( &self, keymin: &ColIndex ) -> Option< usize > { 
    //     match self.keymin_to_ordmin.get( keymin ) {
    //         Some( & ordmin ) => { Some( self.ordmin_to_ord[ ordmin ].clone() ) },
    //         None => None
    //     }
    // } 
    
    // /// Returns either `Some(k)`, where `k` is the ordinal of the major key associated with the given minor key,
    // /// or `None`, if there exists no associated major key.
    // pub fn keymaj_to_ordmin_opt( &self, keymaj: &RowIndex ) -> Option< usize > { 
    //     match self.keymaj_to_ord.get( keymaj ) {
    //         Some( & ord ) => { Some( self.ord_to_ordmin[ ord ].clone() ) },
    //         None => None
    //     }
    // }    

    // pub fn keymaj_to_keymin( &self, keymaj: & RowIndex ) -> & ColIndex { 
    //     & self.ordmin_to_keymin[ 
    //             self.keymaj_to_ordmin( keymaj )
    //         ]        
    // }

    // pub fn keymin_to_keymaj( &self, keymin: & ColIndex ) -> & RowIndex { 
    //     & self.ord_to_keymaj[ 
    //             self.keymin_to_ord( keymin )
    //         ]
    // }    

    // /// The nonzero snzval associated with this `keymaj` index.
    // pub fn keymaj_to_snzval( &self, keymaj: & RowIndex ) -> & Coefficient {
    //     & self.ord_to_snzval[ * self.keymaj_to_ord.get( keymaj ).unwrap() ]
    // }

    // /// The nonzero snzval associated with this `keymin` index.
    // pub fn keymin_to_snzval( &self, keymin: & ColIndex ) -> & Coefficient {
    //     & self.ord_to_snzval[ self.keymin_to_ord( keymin ) ]
    // }    
}


                          


impl < 'a, KeyBoth, Coefficient > 

    GeneralizedMatchingArrayWithMajorOrdinals
        < KeyBoth, KeyBoth, Coefficient > where   
    
    KeyBoth:     Clone + Hash + std::cmp::Eq, // without cloning, it is difficult to write a method that returns the minkey at a given index (one gets "cannot move" errors)
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
        histo
    }

}


//  EXTEND TRAIT IMPLEMENTATION
//  ------------------------------------------------------------------------------

impl    < 'a, ColIndex, RowIndex, Coefficient > 
        Extend< ( RowIndex, ColIndex, Coefficient ) > 
        for 
        GeneralizedMatchingArrayWithMajorOrdinals < ColIndex, RowIndex, Coefficient > 
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{
    /// Add a nonzero entry to a matching array.
    fn extend< I: IntoIterator<Item= ( RowIndex, ColIndex, Coefficient ) >>( &mut self, iter: I) {
        for x in iter.into_iter() { self.push( x.0, x.1, x.2 ) }
    }
}


//  FROM_ITER TRAIT IMPLEMENTATION
//  ------------------------------------------------------------------------------

impl    < 'a, ColIndex, RowIndex, Coefficient > 
        
        FromIterator
        < ( RowIndex, ColIndex, Coefficient ) > for 
        
        GeneralizedMatchingArrayWithMajorOrdinals 
        < ColIndex, RowIndex, Coefficient > 

    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{
    fn from_iter< I: IntoIterator<Item= ( RowIndex, ColIndex, Coefficient ) >>(iter: I) -> Self {

        let mut matching_array  =   GeneralizedMatchingArrayWithMajorOrdinals::new();

        for x in iter.into_iter() { matching_array.push( x.0, x.1, x.2 ) }
        
        matching_array
    }
}


//  ORACLE TRAIT IMPLEMENTATIONS
//  ------------------------------------------------------------------------------


// impl < 'a, ColIndex, RowIndex, Coefficient >

//     OracleRefInherit for
    
//     &'a GeneralizedMatchingArrayWithMajorOrdinals
//         < ColIndex, RowIndex, Coefficient >

//     where   ColIndex:     Clone + Hash + std::cmp::Eq,
//             RowIndex:     Clone + Hash + std::cmp::Eq,
//             Coefficient:     Clone,        
// {}


impl < 'a, ColIndex, RowIndex, Coefficient >

    IndicesAndCoefficients for 

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < ColIndex, RowIndex, Coefficient >

    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{   type ColIndex = ColIndex; 
    type RowIndex = RowIndex; 
    type Coefficient = Coefficient; 
    type EntryMajor = ( ColIndex, Coefficient ); 
    type EntryMinor = ( RowIndex, Coefficient ); }


impl < 'a, ColIndex, RowIndex, Coefficient >
     
    ViewRow for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < ColIndex, RowIndex, Coefficient >
    
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{    
    type ViewMajorIntoIter  =   < Self::ViewMajor as IntoIterator >::IntoIter;   
    type ViewMajor          =   IterTwoType<
                                        Empty< ( ColIndex, Coefficient ) >,
                                        std::iter::Once< ( ColIndex, Coefficient ) >,
                                    > ;

    fn view_major( &self, index: Self::RowIndex ) -> Self::ViewMajor
    {
        match self.keymaj_to_ord( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ord ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ord_to_keymin( ord ),
                                                    self.ord_to_snzval( ord ),  
                                                )          
                            )
                    )
            }
        }
    }
}

impl < 'a, ColIndex, RowIndex, Coefficient >
     
    ViewRowAscend for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < ColIndex, RowIndex, Coefficient >
    
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{    
    type ViewMajorAscendIntoIter = < Self::ViewMajorAscend as IntoIterator >::IntoIter;
    type ViewMajorAscend =  IterTwoType<
                                    Empty< ( ColIndex, Coefficient ) >,
                                    std::iter::Once< ( ColIndex, Coefficient ) >,
                                > ;

    fn view_major_ascend( &self, index: Self::RowIndex ) -> Self::ViewMajorAscend
    {
        match self.keymaj_to_ord( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ord ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ord_to_keymin( ord ),
                                                    self.ord_to_snzval( ord ),  
                                                )          
                            )
                    )
            }
        }
    }
}

impl < 'a, ColIndex, RowIndex, Coefficient >
     
    ViewRowDescend for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < ColIndex, RowIndex, Coefficient >
    
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{    
    type ViewMajorDescendIntoIter   =     < Self::ViewMajorDescend as IntoIterator >::IntoIter;  
    type ViewMajorDescend           =  IterTwoType<
                                                Empty< ( ColIndex, Coefficient ) >,
                                                std::iter::Once< ( ColIndex, Coefficient ) >,
                                            > ;

    fn view_major_descend( &self, index: Self::RowIndex ) -> Self::ViewMajorDescend
    {
        match self.keymaj_to_ord( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ord ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ord_to_keymin( ord ),
                                                    self.ord_to_snzval( ord ),  
                                                )          
                            )
                    )
            }
        }
    }
}


impl < 'a, ColIndex, RowIndex, Coefficient >
     
    ViewCol for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < ColIndex, RowIndex, Coefficient >
    
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{    
    type ViewMinorIntoIter  =   < Self::ViewMinor as IntoIterator >::IntoIter;
    type ViewMinor          =   IterTwoType<
                                        Empty< ( RowIndex, Coefficient ) >,
                                        std::iter::Once< ( RowIndex, Coefficient ) >,
                                    > ;

    fn view_minor( &self, index: Self::ColIndex ) -> Self::ViewMinor
    {
        match self.keymin_to_ord( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ord ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ord_to_keymaj( ord ),
                                                    self.ord_to_snzval( ord ),  
                                                )          
                            )
                    )
            }
        }
    }
}

impl < 'a, ColIndex, RowIndex, Coefficient >
     
    ViewColAscend for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < ColIndex, RowIndex, Coefficient >
    
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{    
    type ViewMinorAscendIntoIter = < Self::ViewMinorAscend as IntoIterator >::IntoIter;
    type ViewMinorAscend =  IterTwoType<
                                    Empty< ( RowIndex, Coefficient ) >,
                                    std::iter::Once< ( RowIndex, Coefficient ) >,
                                > ;

    fn view_minor_ascend( &self, index: Self::ColIndex ) -> Self::ViewMinorAscend
    {
        match self.keymin_to_ord( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ord ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ord_to_keymaj( ord ),
                                                    self.ord_to_snzval( ord ),  
                                                )          
                            )
                    )
            }
        }
    }
}

impl < 'a, ColIndex, RowIndex, Coefficient >
     
    ViewColDescend for

    &'a GeneralizedMatchingArrayWithMajorOrdinals
        < ColIndex, RowIndex, Coefficient >
    
    where   ColIndex:     Clone + Hash + std::cmp::Eq,
            RowIndex:     Clone + Hash + std::cmp::Eq,
            Coefficient:     Clone,
{    
    type ViewMinorDescendIntoIter   =     < Self::ViewMinorDescend as IntoIterator >::IntoIter; 
    type ViewMinorDescend           =  IterTwoType<
                                                Empty< ( RowIndex, Coefficient ) >,
                                                std::iter::Once< ( RowIndex, Coefficient ) >,
                                            > ;

    fn view_minor_descend( &self, index: Self::ColIndex ) -> Self::ViewMinorDescend
    {
        match self.keymin_to_ord( &index ) {
            None => {
                IterTwoType::Iter1( std::iter::empty() )
            }
            Some( ord ) => {
                IterTwoType::Iter2( 
                        std::iter::once(        (   
                                                    self.ord_to_keymaj( ord ),
                                                    self.ord_to_snzval( ord ),  
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


/// A generalized matching matrix
#[derive(Debug, Clone)]
pub struct GeneralizedMatchingArray
                < ColIndex, RowIndex, EntryMaj, EntryMin, Coefficient >
    where   ColIndex:     Clone + Hash + std::cmp::Eq,  // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            RowIndex:     Clone + Hash + std::cmp::Eq,  // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            EntryMin:   Clone + KeyValGet< ColIndex, Coefficient > + KeyValNew< ColIndex, Coefficient >, // clone is ncessary in order for HashMap< RowIndex, EntryMin > to implement `EvaluateFunction< RowIndex, EntryMin >`
            EntryMaj:   Clone + KeyValGet< RowIndex, Coefficient > + KeyValNew< RowIndex, Coefficient >, // clone is ncessary in order for HashMap< ColIndex, EntryMaj > to implement `EvaluateFunction< ColIndex, EntryMaj >`
            Coefficient:     Clone,                       
{  
    keymin_to_entry:    HashMap< ColIndex, EntryMaj >,
    keymaj_to_entry:    HashMap< RowIndex, EntryMin >,
    phantom_snzval:     PhantomData< Coefficient >
}


impl < 'a, ColIndex, RowIndex, EntryMin, EntryMaj, Coefficient > 

    GeneralizedMatchingArray
        < ColIndex, RowIndex, EntryMaj, EntryMin, Coefficient >

    where   ColIndex:     Clone + Hash + std::cmp::Eq, // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            RowIndex:     Clone + Hash + std::cmp::Eq, // cloning seems unavoidable for push/insert operations, where we first need a major/minor key as a value, then, later, as a key 
            EntryMin:   Clone + KeyValGet< ColIndex, Coefficient > + KeyValNew< ColIndex, Coefficient >,
            EntryMaj:   Clone + KeyValGet< RowIndex, Coefficient > + KeyValNew< RowIndex, Coefficient >,   
            Coefficient:     Clone, 
{
    /// Create a new, empty matching array.
    /// 
    /// See [`GeneralizedMatchingArrayWithMajorOrdinals`] for usage.
    pub fn new() -> GeneralizedMatchingArray< ColIndex, RowIndex, EntryMaj, EntryMin, Coefficient > 
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
                    keymin:     ColIndex, 
                    keymaj:     RowIndex, 
                    snzval:     Coefficient 
                )
    {
        let a = self.keymin_to_entry.insert( keymin.clone(), EntryMaj::new( keymaj.clone(), snzval.clone() ) ); 
        let b = self.keymaj_to_entry.insert( keymaj, EntryMin::new( keymin, snzval ) );               
        if a.is_some()  { panic!("attempted to over-write the entry matched to a minor key in a generalized matching matrix")  }
        if b.is_some()  { panic!("attempted to over-write the entry matched to a major key in a generalized matching matrix")  }
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
                    ColIndex, &EntryMaj, RowIndex, 
                    & HashMap< ColIndex, EntryMaj >,                     
                    ExtractKey< RowIndex, Coefficient > 
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
                    RowIndex, & EntryMin, ColIndex, 
                    & HashMap< RowIndex, EntryMin >,                     
                    ExtractKey< ColIndex, Coefficient > 
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
                    ColIndex, & EntryMaj, Coefficient, 
                    & HashMap< ColIndex, EntryMaj >,                     
                    ExtractVal< RowIndex, Coefficient > 
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
                    &RowIndex, & EntryMin, Coefficient, 
                    &HashMap< RowIndex, EntryMin >,                     
                    ExtractVal< ColIndex, Coefficient > 
                > 
    {
        ComposeFunctions::new( & self.keymaj_to_entry, ExtractVal::new() )
    }    




    
    /// Returns an immutable reference `& HashMap< ColIndex, EntryMaj >` representing the partial bijection from 
    /// minor keys to pairs of form `(associated_major_key, associated_nonzero_scalar)`.
    pub fn hashmap_ref_keymin_to_entryymaj( &self ) -> & HashMap< ColIndex, EntryMaj > { & self.keymin_to_entry }

    /// Returns an immutable reference `& HashMap< RowIndex, EntryMaj >` representing the partial bijection from 
    /// major keys to pairs of form `(associated_minor_key, associated_nonzero_scalar)`.
    pub fn hashmap_ref_keymaj_to_entryymin( &self ) -> & HashMap< RowIndex, EntryMin > { & self.keymaj_to_entry }




    /// Returns the nonzero entry that corresponds to the major key, or returns `None` if the major key is unmatched.
    pub fn keymaj_to_entrymin_refopt( &self, keymaj: &RowIndex ) -> Option< & EntryMin > { 
        self.keymaj_to_entry.get( keymaj )
    }

    /// Returns the nonzero entry that corresponds to the minor key, or returns `None` if the minor key is unmatched.
    pub fn keymin_to_entrymaj_refopt( &self, keymin: &ColIndex ) -> Option< & EntryMaj > { 
        self.keymin_to_entry.get( keymin )
    }    
 
    /// The nonzero snzval associated with this `keymaj` index.
    pub fn keymaj_to_snzval_opt( &self, keymaj: & RowIndex ) -> Option< Coefficient >
        // where Coefficient: Clone    
    {
        self.keymaj_to_entry.get( keymaj ).map( |x| x.val() )
    }

    /// The nonzero snzval associated with this `keymin` index.
    pub fn keymin_to_snzval_opt( &self, keymin: & ColIndex ) -> Option< Coefficient >
        // where Coefficient: Clone    
    {
        self.keymin_to_entry.get( keymin ).map( |x| x.val() )
    }    

    /// Returns true iff the given minor key is matched.
    pub fn contains_keymin(  &self, keymin: &ColIndex ) -> bool {
        self.keymin_to_entry.contains_key( keymin )
    }

    /// Returns true iff the given major key is matched.
    pub fn contains_keymaj(  &self, keymaj: &RowIndex ) -> bool {
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
        matching_array.push( 1, 0, 3 );
        matching_array.push( 2, 1, 4 );
        matching_array.push( 0, 2, 5 );    
        
        assert_eq!( matching_array.vec_snzval_ref(), &vec![ 3, 4, 5 ] );

        assert_eq!( matching_array.keymin_to_keymaj( &0 ), Some(1) );
        assert_eq!( matching_array.keymaj_to_keymin( &0 ), Some(2) );     

        assert_eq!( matching_array.keymin_to_ord( &1 ), Some(1) );
        assert_eq!( matching_array.keymaj_to_ord( &1 ), Some(0) ); 

        assert_eq!( matching_array.ord_to_keymin( 0 ), 0 );
        assert_eq!( matching_array.ord_to_keymaj( 0 ), 1 );     

        assert_eq!( matching_array.keymin_to_snzval( &0 ), 3 );
        assert_eq!( matching_array.keymaj_to_snzval( &0 ), 5 );       
        
        assert_eq!( matching_array.ord_to_snzval( 0 ), 3 );
    }    
}



