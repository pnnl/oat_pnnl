//! Structs that implement `EvaluateFunction`


use super::evaluate::EvaluateFunction;


//  IDENTITY FUNCTION
//  -------------------------------------------------------------------------------------------------

/// Implements the `EvaluateFunction` trait, for the function that consumes a vector and returns the reversed vector.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::utilities::functions::misc_functions::ReverseVector;
/// use oat_rust::utilities::functions::evaluate::EvaluateFunction;
/// 
/// let reverser = ReverseVector::new();
/// assert_eq!( vec![2,1], reverser.evaluate_function(vec![1,2]) );
/// ```
#[derive(Clone)]
pub struct ReverseVector{}

impl ReverseVector{
    pub fn new() -> ReverseVector { ReverseVector {} }
}

impl < T > 

    EvaluateFunction
        < Vec<T>, Vec<T> > for 
    
    ReverseVector
{
    fn evaluate_function( &self, input: Vec<T> ) -> Vec<T> { let mut a = input; a.reverse(); a }
}