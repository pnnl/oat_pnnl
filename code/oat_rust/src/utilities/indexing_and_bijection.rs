//! Bijective functions.

use std::iter::FromIterator;



//  ---------------------------------------------------------------------------
//  PERMUTATIONS
//  ---------------------------------------------------------------------------

/// Returns a permutation that sorts the elements of the vector.
pub fn  sort_perm< T: Ord >( vec: & Vec< T > ) -> Vec< usize > {
    let mut sortand     =   Vec::from_iter(
                                vec.iter().enumerate().map(|x| (x.1, x.0) )
                            );
    sortand.shrink_to_fit();
    sortand.sort();
    
    Vec::from_iter( sortand.iter().map(|x| x.1.clone()) )
}

/// Given a vector of length `n+1` representing a permutation on {0, .., n}, 
/// returns a vector that represents the inverse permutation.
pub fn  inverse_perm( vec: & Vec< usize > ) -> Vec< usize > {
    let mut inv_perm    =   Vec::from_iter( std::iter::repeat(0).take( vec.len()) );
    inv_perm.shrink_to_fit();
    for (ind_count, ind) in vec.iter().enumerate() {
        inv_perm[ *ind ] = ind_count;
    }
    inv_perm
}



//  ---------------------------------------------------------------------------
//  WORKING WITH VECTORS
//  ---------------------------------------------------------------------------

//  -----------
//  SUPER INDEX
//  -----------


pub trait SuperIndex < T >
{
    /// Return `self[index]` if `index < self.len()`; otherwise return `super_value`.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::indexing_and_bijection::SuperIndex;
    /// 
    /// let v = vec![0, 1, 2];
    /// assert_eq!( 1, v.sindex(1, 0) );
    /// assert_eq!( 0, v.sindex(5, 0) );
    /// ```
    fn sindex( &self, index: usize, super_value: T ) -> T;
}

impl < T > SuperIndex < T > for Vec< T >
    where T: Clone 
{
    fn sindex( &self, index: usize, super_value: T ) -> T { if self.len() > index { self[ index ].clone() } else { super_value } }
}

//  ----------
//  COMPOSE
//  ----------

pub fn compose_f_after_g< T: Clone > ( f: &Vec< T >, g: &Vec< usize > ) -> Vec< T > {
    Vec::from_iter( g.iter().map(|x| f[ *x ].clone() ) )
}

//  ----------
//  LAST INDEX
//  ----------


pub trait EndIndex< T > 
{

    // THIS FUNCTION IS OBVIATED BY https://doc.rust-lang.org/std/primitive.slice.html#method.last
    // /// Last value of a vector.
    // fn end_val( &self ) -> Option< T >;

    // THIS FUNCTION IS OBVIATED BY https://doc.rust-lang.org/std/primitive.slice.html#method.last_mut
    // /// Mutable reference to last value of a vector.
    // fn end_val_mut< 'a >( &'a mut self ) -> Option< &'a mut T >;

    /// Last ordinal for a vector
    fn end_index ( &self ) -> Option< usize > ;

}

impl < T > EndIndex< T > for Vec< T > 
{
    // /// Last value of a vector.
    // fn end_val( &self ) -> Option< T >  {
    //     println!("Deprecated: The same functionality here could be achieved with the `last_mut` method on slices; prefer that.");        
    //     match self.is_empty() { 
    //         true    =>  None, 
    //         false   =>  Some( self[ self.len() - 1].clone() ) 
    //     }
    // }

    // /// Mutable reference to last value of a vector.
    // fn end_val_mut< 'a >( &'a mut self ) -> Option< &'a mut T > {
    //     println!("Deprecated: The same functionality here could be achieved with the `last_mut` method on slices; prefer that.");
    //     match self.end_index() { 
    //         None      =>  None, 
    //         Some(i)   =>  Some( &mut self[i] ) 
    //     }
    // }

    /// Last ordinal for a vector
    fn end_index( &self ) -> Option< usize > { 
        match self.is_empty() { 
            true    =>  None, 
            false   =>  Some(self.len() - 1) 
        }
    }
}



//  ---------------------------------------------------------------------------
//  SUPER VECTORS
//  ---------------------------------------------------------------------------

/// Returns a constant value for all indices greater than the length of the 
/// internall stored vector.
#[derive(Clone, Debug, PartialEq)]
pub struct SuperVec< T > {
    pub vec: Vec< T >,
    pub val: T
}

impl < T > SuperVec < T > 
    where T : Clone + PartialEq
{
    pub fn val( &self, index: usize ) -> T {

        println!("PREFER USING SINDEX TRAIT ON VECTORS (DEFINED ABOVE) TO CREATING THIS STRUCT");

        if index < self.vec.len() { self.vec[ index ].clone() }
        else { self.val.clone() }
    }
}





#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_sort_perm()
    {
        // initialize vector
        let v           =   vec![1, 2, 3, 4, 2, 1, 2, 2, 1, 4, 3];
        
        // obtain permutations
        let new_to_old  =   sort_perm( & v );
        let old_to_new  =   inverse_perm( & new_to_old );

        // determine grond truth
        let mut v_sorted    =   v.clone();
        v_sorted.sort();

        let ascend      =   Vec::from_iter( 0..v.len() );

        println!("{:?}", compose_f_after_g(&v, &new_to_old));
        println!("{:?}", compose_f_after_g(&old_to_new, &new_to_old));
        println!("{:?}", compose_f_after_g(&new_to_old, &old_to_new));   
        
        assert_eq!(     &compose_f_after_g(&v, &new_to_old), 
                        &v_sorted                                   );

        assert_eq!(     &compose_f_after_g(&old_to_new, &new_to_old), 
                        &ascend                                     );                        

        assert_eq!(     &compose_f_after_g(&new_to_old, &old_to_new), 
                        &ascend                                     );                                                
        
        
    }     


}    
