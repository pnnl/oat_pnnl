//! Prime-order field operators.

use crate::rings::operator_traits::{Semiring, Ring, DivisionRing};





//  ---------------------------------------------------------
//  2   ELEMENT FIELD
//  =========================================================

/// Ring operator for the Boolean ring.  Elements of the ring are `bool`s, i.e. `true` and `false`.
#[derive(Debug, Clone, Copy)]
pub struct BooleanFieldOperator{}

impl BooleanFieldOperator {
    /// Create a new instance of `BooleanFieldOperator`.
    /// 
    /// The commands `BooleanFieldOperator::new()` and `BooleanFieldOperator{}` are equivalent.  The primary advantage of having this
    /// `new` function is that overlaps with the syntax of other ring objects that are harder to construct; 
    /// this makes it easier to interchange one ring object with another.
    pub fn new() -> BooleanFieldOperator { BooleanFieldOperator{} }
}

impl Semiring<bool> for BooleanFieldOperator 
{
    fn is_0( &self, x: bool ) -> bool { ! x         }
    fn is_1( &self, x: bool ) -> bool {   x.clone() }
    fn zero() -> bool { false }
    fn one()  -> bool { true  }

    fn add( &self, x : bool, y : bool ) -> bool { x ^ y }
    fn multiply( &self, x : bool, y: bool ) -> bool { x && y }
}

impl Ring<bool> for BooleanFieldOperator
{
    fn subtract( &self, x : bool, y: bool ) -> bool { x ^ y }
    fn negate( &self, x : bool ) -> bool { x }  // this one is tricky; you want to try logical negation, but you relly need to perform additive negation
}

impl DivisionRing<bool> for BooleanFieldOperator
{
    /// NOTE: THIS DIVISION IS UNSAFE; DESIGNERS MAY WANT TO ADD OPTION TO CHECK FOR DIVISION BY
    /// ZERO
    fn divide( &self, x : bool, y: bool ) -> bool { match y { true => {x.clone()}, false => { panic!("division by zero")}  } }
    
    /// NOTE: THIS DIVISION IS UNSAFE; DESIGNERS MAY WANT TO ADD OPTION TO CHECK FOR DIVISION BY
    /// ZERO
    fn invert( &self, x : bool ) -> bool { match x { true => {x.clone()}, false => { panic!("division by zero")}  } }
}



#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_gf2() {
        
        let ring                        =   BooleanFieldOperator{};

        assert!(    !   ring.is_0( true     ) );
        assert!(        ring.is_0( false    ) );
        assert!(        ring.is_1( true     ) );
        assert!(    !   ring.is_1( false    ) );        
        assert!(        ring.invert( true   ) );
        assert!(        ring.negate( true     ) );
        assert!(    !   ring.negate( false    ) );        
        assert!(    !   ring.add( false, false ) );        
        assert!(        ring.add( false, true  ) );        
        assert!(        ring.add( true,  false ) );                
        assert!(    !   ring.add( true,  true  ) );
        assert!(    !   ring.subtract( false, false ) );        
        assert!(        ring.subtract( false, true  ) );        
        assert!(        ring.subtract( true,  false ) );                
        assert!(    !   ring.subtract( true,  true  ) );  
        assert!(    !   ring.multiply( false, false ) );        
        assert!(    !   ring.multiply( false, true  ) );        
        assert!(    !   ring.multiply( true,  false ) );                
        assert!(        ring.multiply( true,  true  ) );                 
        assert!(    !   ring.divide( false, true  ) );            
        assert!(        ring.divide( true,  true  ) );                  

    }

}





//  ---------------------------------------------------------
//  2   PRIME ORDER FIELD
//  =========================================================



//  ---------------------------------------------------------------------------
//  MODULAR ARITHMETIC

/// Multiplicative inverse of `a` modulo `p`, where `a` and `p` are positive coprime integers.
/// 
/// Concretely, returns the unique integer `m` such that `a * m` is congruent to 1, modulo `p`.
/// 
/// This implementation of the extended Euclidean algorithm is adapted from 
/// pseudocode found here: <https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm>.
/// 
/// NOTE TO DEVELOPERS !!!: this implementation is based on the pseudocode for "computingn multiplicative
/// inverses in modular structures" at <https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm>.  There is
/// another algorithm which is proved to require the least number of iterations possible, described as
/// the method of least absolute remainders in <https://en.wikipedia.org/wiki/Euclidean_algorithm#:~:text=In%20mathematics%2C%20the%20Euclidean%20algorithm,them%20both%20without%20a%20remainder.&text=When%20that%20occurs%2C%20they%20are,of%20the%20original%20two%20numbers>.
/// 
/// If `a` lies outside [0,p), then it is first reduced mod `p`, then processed as normal.
/// 
/// Panics if 
///   * `a` is congruent to 0 mod `p`
///   * `a` has no inverse, modulo `p`
/// 
/// Does **not** panic if
///   * `p` is not prime.  (It is potentially expensive to check that `p` is prime;
///     the method used here checks divisibility by every odd integer less than the
///     square root)
/// 
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::rings::operator_structs::field_prime_order::invert_mod_p;
/// 
/// let primes = [2, 3, 5, 7, 313, 947];
/// for p in primes {
///     for a in 1 .. p {
///         let a_inv   = invert_mod_p( a, p );
///         assert_eq!( 1, (a * a_inv) % p );
///     }
/// }
/// ```
pub fn invert_mod_p( a:usize, p:usize ) -> usize {

    // preallocate (*all* of these values will be used)
    let mut t    = 0 as isize;       let mut r    = p as isize;
    let mut tnew = 1 as isize;       let mut rnew = a as isize;

    // ensure that a lies in [0,p)
    rnew    =   rnew % p as isize;
    
    // preallocate (*none* of these values will be used)
    let mut quotient; 
    let mut priorstep_t; 
    let mut priorstep_r; 
    let mut priorstep_tnew; 
    let mut priorstep_rnew; 

    // execute the algorithm
    while rnew != 0 {
        quotient = r / rnew;

        priorstep_t     =   t.clone();
        priorstep_tnew  =   tnew.clone();
        priorstep_r     =   r.clone();
        priorstep_rnew  =   rnew.clone();        

        t = priorstep_tnew;   tnew = priorstep_t - quotient * priorstep_tnew; 
        r = priorstep_rnew;   rnew = priorstep_r - quotient * priorstep_rnew; 
    }
    
    // format the output
    if      r > 1 { panic!("{:?} has no inverse modulo {:?}", a, p)  }
    else if t < 0 { return (t + (p as isize))  as usize }
    else          { return  t                  as usize }

 }



//  ---------------------------------------------------------------------------
//  PRIME-ORDER FIELDS

/// Prime-order field ring operations (add, substract, multiply, etc).  
/// Elements of the field are of type `usize`.
/// 
/// Note that any two fields of the same prime order are isomorphic.
/// 
/// The modulus is stored as a private attribute; it can be accessed via the [`PrimeOrderFieldOperator::modulus`]  function.
/// 
/// Semiring/Ring/DivisionRing operations panic if they receive an integer greater than `p` as input.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    /// use oat_rust::rings::operator_traits::{Semiring};
/// 
/// // Constructor panics if passed a modulus that is not prime
/// let field = PrimeOrderFieldOperator::new(7); // no panic
/// let result = std::panic::catch_unwind(|| PrimeOrderFieldOperator::new(4) ); // panic!
/// assert!(result.is_err());
/// 
/// // Operations panic on arguments >= `p`
/// let sum = field.add( 0, 2 ); // no panic
/// assert_eq!( sum, 2 );
/// let result = std::panic::catch_unwind(|| field.add( 0, 9 ) ); // panic!
/// assert!(result.is_err());
/// 
/// ```
#[derive(Debug, Clone, Copy)]
pub struct PrimeOrderFieldOperator{ modulus: usize }


impl PrimeOrderFieldOperator {

    /// Returns the modulus of the field.
    /// `
    pub fn modulus( &self ) -> usize { self.modulus.clone() }
    
    
    /// Create a new operator for a field of order `p`.  Panics if `p` is not prime.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    /// 
    /// let field = PrimeOrderFieldOperator::new(3); // no panic
    /// let result = std::panic::catch_unwind(|| PrimeOrderFieldOperator::new(4) ); // panic!
    /// assert!(result.is_err());
    /// 
    /// ```
    pub fn new( modulus: usize ) -> PrimeOrderFieldOperator { 
                match primes::is_prime( modulus as u64 ) {
                    true  =>  { return PrimeOrderFieldOperator{ modulus: modulus } },
                    false =>  { panic!() },
                }
            }
}

impl Semiring<usize> for PrimeOrderFieldOperator 
{
    fn is_0( &self, x: usize ) -> bool { x == 0         }
    fn is_1( &self, x: usize ) -> bool { x == 1         }
    fn zero() -> usize { 0 }
    fn one()  -> usize { 1 }

    fn add( &self, x : usize, y : usize ) -> usize { 
        if x >= self.modulus { panic!( "{:?} lies outside the field of order {:?}", x, self.modulus() )}
        else 
        if y >= self.modulus { panic!( "{:?} lies outside the field of order {:?}", y, self.modulus() )}        
        
        (x + y) % self.modulus 
    }

    fn multiply( &self, x : usize, y: usize ) -> usize { 
        if x >= self.modulus { panic!( "{:?} lies outside the field of order {:?}", x, self.modulus() )}
        else 
        if y >= self.modulus { panic!( "{:?} lies outside the field of order {:?}", y, self.modulus() )}

        (x * y) % self.modulus     
    }
}

impl Ring<usize> for PrimeOrderFieldOperator
{

    /// # Examples
    /// 
    /// ```
    /// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    /// use oat_rust::rings::operator_traits::{Semiring, Ring, DivisionRing};
    /// 
    /// let primes = [2, 3, 5, 7, 13];
    /// for p in primes {
    ///     let field = PrimeOrderFieldOperator::new(p);
    ///     for a in 0 .. p {
    ///         for b in 0 .. p {
    ///         assert_eq!( b, field.subtract( field.add(a, b), a)  );
    ///         }
    ///     }
    /// }    
    /// ```
    fn subtract( &self, x : usize, y: usize ) -> usize { 
        if x >= self.modulus { panic!( "{:?} lies outside the field of order {:?}", x, self.modulus() )}
        else 
        if y >= self.modulus { panic!( "{:?} lies outside the field of order {:?}", y, self.modulus() )}
        match x >= y {
            true  => { x - y },
            false => { x + self.modulus - y },
        }   
    }

    /// # Examples
    /// 
    /// ```
    /// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    /// use oat_rust::rings::operator_traits::{Semiring, Ring, DivisionRing};
    /// 
    /// let primes = [2, 3, 5, 7, 313];
    /// for p in primes {
    ///     let field = PrimeOrderFieldOperator::new(p);
    ///     for a in 0 .. p {
    ///         assert_eq!( 0, field.add(a, field.negate(a)));
    ///     }
    /// }
    /// ```    
    fn negate( &self, x : usize ) -> usize { 
        if x > 0 {
            match x < self.modulus {
                true  => { return self.modulus - x },
                false => { panic!( "{:?} lies outside the field of order {:?}", x, self.modulus() )}
            }
        }
        else { return 0 }        
    } 
}

impl DivisionRing<usize> for PrimeOrderFieldOperator
{
    /// NOTE: THIS DIVISION IS UNSAFE; DESIGNERS MAY WANT TO ADD OPTION TO CHECK FOR DIVISION BY
    /// ZERO
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    /// use oat_rust::rings::operator_traits::{Semiring, Ring, DivisionRing};
    /// 
    /// let primes = [2, 3, 5, 7, 313, 947];
    /// for p in primes {
    ///     let field = PrimeOrderFieldOperator::new(p);
    ///     for a in 1 .. p {
    ///         let a_inv   = field.invert(a);
    ///         assert_eq!( 1, field.multiply(a, a_inv));
    ///     }
    /// }
    /// ```
    fn divide( &self, x : usize, y: usize ) -> usize {  self.multiply(x, self.invert(y)) }
    
    /// NOTE: THIS DIVISION IS UNSAFE; DESIGNERS MAY WANT TO ADD OPTION TO CHECK FOR DIVISION BY
    /// ZERO
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::rings::operator_structs::field_prime_order::PrimeOrderFieldOperator;
    /// use oat_rust::rings::operator_traits::{Semiring, Ring, DivisionRing};
    /// 
    /// let primes = [2, 3, 5, 7, 313, 947];
    /// for p in primes {
    ///     let field = PrimeOrderFieldOperator::new(p);
    ///     for a in 1 .. p {
    ///         let a_inv   = field.invert(a);
    ///         assert_eq!( 1, field.multiply(a, a_inv));
    ///     }
    /// }
    /// ```
    fn invert( &self, x : usize ) -> usize {  invert_mod_p( x, self.modulus() ) }
}





#[cfg(test)]
mod tests_prime_order_field {

    // Modular arithmetic
    // ====================================================================================================================
    
    // * SEE DOC TESTS


    // Prime order fields
    // ====================================================================================================================
    
    // * SEE DOC TESTS

}
