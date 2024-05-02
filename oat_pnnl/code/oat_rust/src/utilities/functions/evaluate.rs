//! Evaluate functions.


use std::collections::HashMap;
use std::hash::Hash;
use std::marker::PhantomData;
 




/// Evaluate a function.
/// 
/// This trait is similar in spirit to `Fn`, `FnOnce` and `FnMut`.  A key problem with `Fn`, `FnOnce` and `FnMut` 
/// is that manual implementation is experimental [E0183](https://doc.rust-lang.org/error-index.html#E0183).  The
/// `EvaluateFunction` trait can act as a stand-in for one of these traits when manual implementation is required.
pub trait EvaluateFunction< Input, Output >{
    fn evaluate_function( &self, input: Input ) -> Output;
}


//  IDENTITY FUNCTION
//  -------------------------------------------------------------------------------------------------

/// A struct that implements the identify function via the `EvaluateFunction` trait
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::IdentityFunction;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
/// 
/// let identity_function = IdentityFunction{};
/// assert_eq!( 1, identity_function.evaluate_function(1) );
/// ```
pub struct IdentityFunction{}

impl IdentityFunction{
    pub fn new() -> IdentityFunction { IdentityFunction {} }
}

impl < Input > 

    EvaluateFunction
        < Input, Input > for 
    
    IdentityFunction
{
    fn evaluate_function( &self, input: Input ) -> Input { input }
}




//  IMPLEMENT ON WRAPPERS -- OF FNMUT
//  -------------------------------------------------------------------------------------------------


/// Wrapper for objects that implement `Fn( Input ) -> Output`; the wrapper implements `EvaluateFunction< Input, Output >`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionFnMutWrapper};
/// 
/// let fun = | x: usize | -> usize { 2 * x };
/// let wrapped_fun = EvaluateFunctionFnMutWrapper::new( fun );
/// assert_eq!( wrapped_fun.evaluate_function(1), 2 );
/// ```
/// 
/// # Design notes
/// 
/// We found that implementing `EvaluateFunction` on this wrapper helped the compiler to produce helpful error messages when code breaks;
/// however this took place at early stages of development.  It may be worth reviewing whether these wrappers are truly necessary.
pub struct EvaluateFunctionFnMutWrapper< Function > { unwrapped_function: Function }

impl < Function > EvaluateFunctionFnMutWrapper< Function > { 
    /// Create new function wrapper.
    pub fn new( unwrapped_function: Function ) -> EvaluateFunctionFnMutWrapper< Function > { EvaluateFunctionFnMutWrapper {unwrapped_function} } 
}

impl < Function, Input, Output >

    EvaluateFunction
        < Input, Output > for

    EvaluateFunctionFnMutWrapper
        < Function >

    where 
        Function:   Fn( Input ) -> Output
{
    fn evaluate_function( & self, input: Input ) -> Output {
        (self.unwrapped_function)( input )
    }
}       

impl < Function > 

    Clone for

    EvaluateFunctionFnMutWrapper
        < Function >   
        
    where
        Function:   Clone
{
    fn clone( &self ) -> EvaluateFunctionFnMutWrapper<Function>
    { EvaluateFunctionFnMutWrapper::new( self.unwrapped_function.clone() ) }
}


//  IMPLEMENT ON WRAPPERS -- OF IMPLEMENTORS
//  -------------------------------------------------------------------------------------------------


/// Wrapper for types that implement the [`EvaluateFunction`] trait in multiple different ways.
/// 
/// Some times, such as `HashMap` implement the [`EvaluateFunction`] trait in multiple different ways.  
/// This can lead to errors where the compiler is uncertain which function the object represents.
/// By wrapping the object in an `EvaluateFunctionDisambiguator` struct, we produce an object that
/// only implements one instance of the [`EvaluateFunction`] trait, thus removing any ambiguity.
pub struct EvaluateFunctionDisambiguator< Function: EvaluateFunction< Input, Output >, Input, Output > 
{ 
    unwrapped_function: Function,
    phantom_input:      PhantomData< Input >,
    phantom_output:     PhantomData< Output >,
}

impl < Input, Output, Function: EvaluateFunction< Input, Output > > 

    EvaluateFunctionDisambiguator
        < Function, Input, Output, >
{
    pub fn new( unwrapped_function: Function ) -> Self { EvaluateFunctionDisambiguator{ unwrapped_function, phantom_input: PhantomData, phantom_output: PhantomData } }
}        

impl < Input, Output, Function: EvaluateFunction< Input, Output > > 

    EvaluateFunction
        < Input, Output > for

    EvaluateFunctionDisambiguator
        < Function, Input, Output, >
{
    fn evaluate_function( &self, input: Input ) -> Output { self.unwrapped_function.evaluate_function(input) }
}    





//  IMPLEMENT ON VECTORS
//  -------------------------------------------------------------------------------------------------


//  UNWRAPPED VALUES
//  --------------------


impl < Output: Clone > 

    EvaluateFunction
        < usize, Output > for

    Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> Output {
        return self[ input ].clone()
    }    
}



impl < 'a, Output: Clone > 

    EvaluateFunction
        < usize, Output > for

    &'a Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> Output {
        return self[ input ].clone()
    }    
}



impl < 'a, Output > 

    EvaluateFunction
        < usize, &'a Output > for

    &'a Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> &'a Output {
        return & self[ input ]
    }    
}


//  WRAPPED IN OPTIONS
//  --------------------


impl < Output: Clone > 

    EvaluateFunction
        < usize, Option< Output > > for

    Vec< Output >
{
    fn evaluate_function( & self, input: usize ) -> Option< Output > {
        match input < self.len() {
            true    =>  { return Some( self[ input ].clone() ) },
            false   =>  { return None }
        }        
    }    
}



impl < 'a, Output: Clone > 

    EvaluateFunction
        < usize, Option< Output > > for

    &'a Vec< Output >
{
    /// Return `self[input].clone()`.
    fn evaluate_function( & self, input: usize ) -> Option< Output > {
        match input < self.len() {
            true    =>  { return Some( self[ input ].clone() ) },
            false   =>  { return None }
        }        
    }    
}



impl < 'a, Output > 

    EvaluateFunction
        < usize, Option< &'a Output > > for

    &'a Vec< Output >
{
    /// Return `Some( &self[ input ] )` if `input < self.len()`, and `None` otherwise.
    fn evaluate_function( & self, input: usize ) -> Option< &'a Output > {
        match input < self.len() {
            true    =>  { return Some( & self[ input ] ) },
            false   =>  { return None }
        }        
    }     
}



//  IMPLEMENT ON HASHMAPS
//  -------------------------------------------------------------------------------------------------


// 0 - Implement the trait for < Input, Output > on HashMap
impl < Input, Output > 

    EvaluateFunction
        < Input, Output > for 
        
    HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, i32 >
    ///     = EvaluateFunctionDisambiguator::new( hash );     
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1usize ), -1  );  
    /// ```
    /// 
    /// # Design note
    /// 
    /// This implementation clones the output of the hashmap by default; to do otherwise 
    /// may ultimately require the involvement of lifetime parameters, and a more complicated implementation.
    fn evaluate_function( & self, input: Input ) -> Output {
        return self.get( &input ).unwrap().clone()
    }
}

// 1 - Implement the trait for < &Input, Output > on HashMap
impl < Input, Output > 

    EvaluateFunction
        < &Input, Option< Output > > for 
        
    HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// assert_eq!( hash.evaluate_function( &1 ), Some( -1 )  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// This implementation clones the output of the hashmap by default; to do otherwise 
    /// may ultimately require the involvement of lifetime parameters, and a more complicated implementation.
    fn evaluate_function( & self, input: &Input ) -> Option< Output > {
        return self.get( &input ).cloned()
    }
}



// 2 - Implement the trait for < Input, Output > on HashMap
impl < Input, Output > 

    EvaluateFunction
        < Input, Option< Output > > for 
        
    HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;    
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, Option< i32 > >
    ///     = EvaluateFunctionDisambiguator::new( hash );
    ///
    /// assert_eq!( disambiguated.evaluate_function( 1usize ), Some( -1i32 )  );
    /// ```
    fn evaluate_function( & self, input: Input ) -> Option< Output > {
        return self.get( &input ).cloned()
    }
}



// 3 - Implement the trait for < &Input, Option< &'a Output > > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < &Input, Option< &'a Output > > for 
        
    &'a HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, &usize, Option<&i32>, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );
    /// 
    /// assert_eq!( disambiguated.evaluate_function( &1 ), Some( & -1 )  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: &Input ) -> Option< &'a Output > {
        return self.get( &input )
    }
}



// 3.1 - Implement the trait for < &Input, Option< &'a Output > > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < &Input, &'a Output > for 
        
    &'a HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, &usize, &i32, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );
    /// 
    /// assert_eq!( disambiguated.evaluate_function( &1 ), & -1  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: &Input ) -> &'a Output {
        return self.get( &input ).unwrap()
    }
}



// 3.2 - Implement the trait for < &Input, Option< &'a Output > > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, &'a Output > for 
        
    &'a HashMap
        < Input, Output > 

    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output        
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, &i32, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1 ), & -1  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: Input ) -> &'a Output {
        return self.get( &input ).unwrap()
    }
}






// 4 - Implement the trait for < Input, &'a Output > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, Option< &'a Output > > for 
        
    &'a HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<&i32>, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );        
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1 ), Some( &-1 )  );
    /// ```
    /// 
    /// # Design note
    /// 
    /// It would be difficult to implement a version of `EvaluateFunction<_, _>` where `evaulate_function` returns a reference
    /// on `HashMap` directly, because a lifetime parameter is likely needed on the output.  Hence we implement `EvaluateFunction<_, _>`
    /// on `& HashMap`.
    fn evaluate_function( & self, input: Input ) -> Option< &'a Output > {
        return self.get( &input )
    }
}



// 5 - Implement the trait for < Input, &'a Output > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, Option< Output > > for 
        
    &'a HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;                
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<i32>, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );  
    /// 
    /// assert_eq!( disambiguated.evaluate_function(1), Some( -1 )  );
    /// ```
    fn evaluate_function( & self, input: Input ) -> Option< Output > {
        return self.get( &input ).cloned()
    }
}



// 6 - Implement the trait for < Input, Output > on &'a HashMap
impl < 'a, Input, Output > 

    EvaluateFunction
        < Input, Output > for 
        
    &'a HashMap
        < Input, Output > 
    
    where
        Input:  std::cmp::Eq + Hash,
        Output: Clone,
{    
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
    /// use std::collections::HashMap;        
    /// 
    /// let mut hash = HashMap::new();
    /// hash.insert( 1usize, -1); 
    /// 
    /// // we must sometimes disambiguate the types of the input and output
    /// let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, i32, > 
    ///     = EvaluateFunctionDisambiguator::new( &hash );      
    /// 
    /// assert_eq!( disambiguated.evaluate_function( 1 ), -1  );
    /// ```
    fn evaluate_function( & self, input: Input ) -> Output {
        return self.get( &input ).unwrap().clone()
    }
}








//  =========================================================================================================
//  TESTS
//  =========================================================================================================


//  ---------------------------------------------------------------------
//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {

    // 0
    // EvaluateFunction
    //     < Input, Output > for 
        
    // HashMap
    //     < Input, Output > 
    #[test]
    fn test_0() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        
        // we must sometimes disambiguate the types of the input and output        
        let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, i32 >
            = EvaluateFunctionDisambiguator::new( hash );     
        
        assert_eq!( disambiguated.evaluate_function( 1usize ), -1  );   
    } 
    

    // 1
    // EvaluateFunction
    //     < &Input, Option< Output > > for 
        
    // HashMap
    //     < Input, Output >         
    #[test]
    fn test_1() {
        use crate::utilities::functions::evaluate::EvaluateFunction;
        use std::collections::HashMap;
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        
        assert_eq!( hash.evaluate_function( &1 ), Some( -1 )  );
    }
        
    // 2
    // EvaluateFunction
    //     < Input, Option< Output > > for 
        
    // HashMap
    //     < Input, Output > 
    #[test]
    fn test_2() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;    

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output        
        let disambiguated: EvaluateFunctionDisambiguator::< HashMap< usize, i32 >, usize, Option< i32 > >
            = EvaluateFunctionDisambiguator::new( hash );
        
        assert_eq!( disambiguated.evaluate_function( 1usize ), Some( -1i32 )  );
    }
        

    // 3
    // EvaluateFunction
    //     < &Input, Option< &'a Output > > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_3() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output        
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, &usize, Option<&i32>, > 
            = EvaluateFunctionDisambiguator::new( &hash );
        
        assert_eq!( disambiguated.evaluate_function( &1 ), Some( & -1 )  );
    }


    // 4
    // EvaluateFunction
    //     < Input, Option< &'a Output > > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_4() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<&i32>, > 
            = EvaluateFunctionDisambiguator::new( &hash );        
        
        assert_eq!( disambiguated.evaluate_function( 1 ), Some( &-1 )  );
    }
        

    // 5
    // EvaluateFunction
    //     < Input, Option< Output > > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_5() {
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;                
        
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 

        // we must sometimes disambiguate the types of the input and output
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, Option<i32>, > 
            = EvaluateFunctionDisambiguator::new( &hash );  
        
        assert_eq!( disambiguated.evaluate_function(1), Some( -1 )  );
    }
            
    // 6
    // EvaluateFunction
    //     < Input, Output > for 
        
    // &'a HashMap
    //     < Input, Output > 
    #[test]
    fn test_6() {        

        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionDisambiguator};
        use std::collections::HashMap;        

        let mut hash = HashMap::new();
        hash.insert( 1usize, -1); 
        
        // we must sometimes disambiguate the types of the input and output
        let disambiguated: EvaluateFunctionDisambiguator< &HashMap<usize, i32>, usize, i32, > 
            = EvaluateFunctionDisambiguator::new( &hash );      

        assert_eq!( disambiguated.evaluate_function( 1 ), -1  );
    
    }    
}