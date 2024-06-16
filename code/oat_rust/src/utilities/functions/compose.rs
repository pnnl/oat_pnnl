//! Compose functions.

use std::marker::PhantomData;

use super::evaluate::EvaluateFunction;






/// Represents the composition of two functions.
pub struct ComposeFunctions< Input1, Output1, Output2, Function1, Function2 > 
    where 
        Function1:  EvaluateFunction< Input1, Output1 >,
        Function2:  EvaluateFunction< Output1, Output2 >,
{
    function1:          Function1,
    function2:          Function2,
    phantom_input1:     PhantomData< Input1 >,
    phantom_output1:    PhantomData< Output1 >,
    phantom_output2:    PhantomData< Output2 >,        
}

impl < Input1, Output1, Output2, Function1, Function2 > 

    ComposeFunctions
        < Input1, Output1, Output2, Function1, Function2 > 

    where 
        Function1:  EvaluateFunction< Input1, Output1 >,
        Function2:  EvaluateFunction< Output1, Output2 >,
{
    /// Creates a new [`ComposeFunctions`] struct representing the composition of the two 
    /// input functions.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::utilities::functions::compose::ComposeFunctions;
    /// use oat_rust::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionFnMutWrapper};
    /// 
    /// let function1 = | x: usize | -> usize { 2 * x }; // double
    /// let function2 = | x: usize | -> usize { x.pow(2) }; // square
    /// let function1_wrapper = EvaluateFunctionFnMutWrapper::new( function1 );
    /// let function2_wrapper = EvaluateFunctionFnMutWrapper::new( function2 );        
    /// 
    /// let composition1 = ComposeFunctions::new( function1_wrapper.clone(), function2_wrapper.clone() ); // double then square
    /// let composition2 = ComposeFunctions::new( function2_wrapper.clone(), function1_wrapper.clone() ); // square then double
    /// 
    /// assert_eq!( 4, composition1.evaluate_function(1) );
    /// assert_eq!( 2, composition2.evaluate_function(1) );
    /// ```
    pub fn new( function1: Function1, function2: Function2 ) ->  ComposeFunctions< Input1, Output1, Output2, Function1, Function2 >  {
        ComposeFunctions{ function1, function2, phantom_input1: PhantomData, phantom_output1: PhantomData, phantom_output2: PhantomData }
    }
}     

impl < Input1, Output1, Output2, Function1, Function2 > 

    EvaluateFunction
        < Input1, Output2 > for

    ComposeFunctions
        < Input1, Output1, Output2, Function1, Function2 > 

    where 
        Function1:  EvaluateFunction< Input1, Output1 >,
        Function2:  EvaluateFunction< Output1, Output2 >,        
{
    fn evaluate_function( &self, input: Input1 ) -> Output2 {
        self.function2.evaluate_function(
                self.function1.evaluate_function( input )
            )
    }
}
















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
    fn test_compose_functions() {
        use crate::utilities::functions::compose::ComposeFunctions;
        use crate::utilities::functions::evaluate::{EvaluateFunction, EvaluateFunctionFnMutWrapper};
        
        let function1 = | x: usize | -> usize { 2 * x }; // double
        let function2 = | x: usize | -> usize { x.pow(2) }; // square
        let function1_wrapper = EvaluateFunctionFnMutWrapper::new( function1 );
        let function2_wrapper = EvaluateFunctionFnMutWrapper::new( function2 );        
        
        let composition1 = ComposeFunctions::new( function1_wrapper.clone(), function2_wrapper.clone() ); // double then square
        let composition2 = ComposeFunctions::new( function2_wrapper.clone(), function1_wrapper.clone() ); // square then double
        
        assert_eq!( 4, composition1.evaluate_function(1) );
        assert_eq!( 2, composition2.evaluate_function(1) );
    }    
}