//! Simplices that have associated filtration values

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
// use pyo3_log;
use pyo3::types::IntoPyDict;

// use polars::prelude::*;
// use polars::df;




use oat_rust::algebra::chains::barcode::Bar;
use oat_rust::utilities::optimization::minimize_l1::minimize_l1;
use oat_rust::topology::simplicial::simplices::filtered::SimplexFiltered;

use itertools::Itertools;
use num::rational::Ratio;
use ordered_float::OrderedFloat;

use std::fmt::Debug;
use std::hash::Hash;



//  =========================================
//  A FILTERED SIMPLEX OBJECT
//  =========================================


#[pyclass]
#[derive(Clone)]
pub struct SimplexFilteredPy{ simplex: SimplexFiltered<OrderedFloat<f64>> }

#[pymethods]
impl SimplexFilteredPy {
    pub fn filtration( &self ) -> f64 { self.simplex.filtration().clone().into_inner() }    
    pub fn vertices( &self ) -> Vec<u16> { self.simplex.vertices().clone() }    
    pub fn dimension( &self ) -> usize { self.simplex.dimension() }
}

impl SimplexFilteredPy {
    pub fn simplex( &self ) -> &SimplexFiltered<OrderedFloat<f64>> { &self.simplex }    
    pub fn new( simplex: SimplexFiltered<OrderedFloat<f64>> ) -> SimplexFilteredPy 
    { SimplexFilteredPy{ simplex }  }
}

//  =========================================
//  TYPE CONVERSION FOR CHAINS INDEXED BY FILTERED SIMPLICES (A CONVENIENCE FUNCTION)
//  =========================================

pub fn convert_chain_to_python( chain: & Vec< (SimplexFiltered< OrderedFloat<f64> >, Ratio<isize> ) > ) 
        -> 
        Vec< ( SimplexFilteredPy, (isize, isize) ) >  {
        chain.iter()
            .map(|x| 
                    ( 
                        SimplexFilteredPy::new( x.0.clone() ), 
                        (x.1.numer().clone(), x.1.denom().clone()) 
                    ) 
                )
                .collect_vec()        
}  


//  =========================================
//  BARS IN A BARCODE INDEXED BY FILTERED SIMPLICES
//  =========================================

#[pyclass]
#[derive(Clone,Debug)]
pub struct BarPySimplexFilteredRational{
                    bar:    Bar<
                                    SimplexFiltered< OrderedFloat< f64 > >, 
                                    ( SimplexFiltered< OrderedFloat< f64 > >, Ratio<isize> )
                                >
}
#[pymethods]
impl BarPySimplexFilteredRational{
    /// id_number(&self, /)
    /// --
    ///
    /// Displays an integer used to uniquely identify the bar.
    pub fn id_number(&self) -> usize { self.bar.id_number().clone() }
    #[pyo3(text_signature = "($self)")]    
    pub fn dimension(&self) -> isize { self.bar.dimension() }
    #[pyo3(text_signature = "($self)")]    
    pub fn birth(&self) -> f64 { self.bar.birth_f64() }
    #[pyo3(text_signature = "($self)")]    
    pub fn death(&self) -> f64 { self.bar.death_f64() }    
    #[pyo3(text_signature = "($self)")]   
    /// birth_column(&self, /)
    /// --
    ///
    /// The "positive simplex" associated with this bar
    pub fn birth_column(&self) -> SimplexFilteredPy { SimplexFilteredPy::new( self.bar.birth_column().clone() ) }
    #[pyo3(text_signature = "($self)")]    
    /// birth_column(&self, /)
    /// --
    ///
    /// The "negative simplex" associated with this bar (or None, if there is no negative simplex)
    pub fn death_column(&self) -> Option< SimplexFilteredPy > { self.bar.death_column().clone().map(|x| SimplexFilteredPy::new(x) ) }
    #[pyo3(text_signature = "($self)")]    
    pub fn cycle_representative(&self) -> Option< Vec< ( SimplexFilteredPy, (isize, isize) ) > > {
                self.bar.cycle_representative().clone().map(|x| convert_chain_to_python(&x) ) 
            }
    #[pyo3(text_signature = "($self)")]            
    pub fn bounding_chain(&self) -> Option< Vec< ( SimplexFilteredPy, (isize, isize) ) > > {
                self.bar.bounding_chain().clone().map(|x| convert_chain_to_python(&x) ) 
            }       
}


impl BarPySimplexFilteredRational{
    /// Return a reference to the internally stored bar
    pub fn peek(&self) -> &Bar<
                                    SimplexFiltered< OrderedFloat< f64 > >, 
                                    ( SimplexFiltered< OrderedFloat< f64 > >, Ratio<isize> )
                                > 
        { & self.bar }

    /// Return the internally stored bar
    pub fn disolve( self ) 
        -> 
        Bar<
                SimplexFiltered< OrderedFloat< f64 > >, 
                ( SimplexFiltered< OrderedFloat< f64 > >, Ratio<isize> )
            >
    { self.bar }
}

//  =========================================
//  BARCODES INDEXED BY FILTERED SIMPLICES
//  =========================================

/// The barcode of the homological persistence module of a filtered simplicial complex.
#[pyclass]
#[derive(Clone)]
pub struct BarcodePySimplexFilteredRational{
                    barcode: oat_rust::algebra::chains::barcode::Barcode<
                                    SimplexFiltered< OrderedFloat< f64 > >, 
                                    ( SimplexFiltered< OrderedFloat< f64 > >, Ratio<isize> )
                                >
                }
// python-free methods
impl BarcodePySimplexFilteredRational{

    pub fn new( barcode: oat_rust::algebra::chains::barcode::Barcode<
                                    SimplexFiltered< OrderedFloat< f64 > >, 
                                    ( SimplexFiltered< OrderedFloat< f64 > >, Ratio<isize> )
                                > 
            ) -> Self {
        BarcodePySimplexFilteredRational{ barcode }
    }    
}

// python-bound methods
#[pymethods]
impl BarcodePySimplexFilteredRational{

    /// Construct a barcode from a list of bars
    /// 
    /// If you already have a collection of bars, then wrapping them in a `BarcodePySimplexFilteredRational`
    /// will allow you to use all the associated convenience functions (e.g. calculating the maximum
    /// finite endpoint).
    #[new]
    pub fn py_new( list: Vec< BarPySimplexFilteredRational > ) -> Self {
        let barcode 
            = oat_rust::algebra::chains::barcode::Barcode::new( 
                    list.into_iter().map(|x| x.disolve() ) 
                );
        BarcodePySimplexFilteredRational{ barcode }
    }

    /// bar( bar_id_number: usize, /)
    /// --
    ///
    /// Retreive a deep copy of an internally stored [`Bar`]
    // #[args(bar_id_number = "0",)]
    #[pyo3( signature = ( bar_id_number, / )) ]
    pub fn bar(&self, bar_id_number: usize ) -> PyResult< BarPySimplexFilteredRational > { 
        Ok( BarPySimplexFilteredRational{ bar: self.barcode.bar(bar_id_number).clone() } )
    }    

    /// bars($self, /)
    /// --
    ///
    /// A list containing deep copies of every bar in the barcode
    pub fn bars( &self ) -> Vec< BarPySimplexFilteredRational > { 
        self.barcode.iter().map(|x| BarPySimplexFilteredRational{ bar: x.clone() } ).collect() 
    }    

    /// bars_in_dim( &self, dim: isize /)
    /// --
    ///
    /// A list containing deep copies of every *dimension-`\dim` bar in the barcode
    #[args( bar_id_number = "0", )]
    pub fn bars_in_dim( &self, dim: isize ) -> Vec< BarPySimplexFilteredRational > { 
        self.barcode.iter().filter(|x| x.dimension()==dim ).map(|x| BarPySimplexFilteredRational{ bar: x.clone() } ).collect() 
    }        


    /// Returns a vector of triples `(birth, death, id)`, where `id` is the uniue id of the bar.
    pub fn intervals( &self, dim: isize) -> Vec< (f64, f64, usize ) > {
        self.barcode.intervals_f64( dim )
    }

    /// A vector of tuples `(t, betti_t)` where `t` is an endpoint of an
    /// interval in the barcode, and `betti_t` is the dimension `dim` betti
    /// number at filtration parameter `t`.
    pub fn betti_curve( &self, dim: isize ) -> Vec< ( f64, usize ) > {
        self.barcode.betti_curve(dim).into_iter().map(|x| (x.0.into_inner(), x.1)).collect_vec()
    } 

    ///  Return a sorted list of all endpoints of intervals in the barcode.
    pub fn endpoints( &self ) -> Vec< f64 >{ 
        self.barcode.endpoints_ordf64().into_iter().map(|x|x.into_inner()).collect_vec() 
    }

    /// The maximum finite endpoint of any interval (or `None`, if no such endpoint exists)
    pub fn max_finite_endpoint( &self ) -> Option<f64> { 
        self.barcode.max_finite_endpoint().map(|x| x.into_inner()) 
    }

    /// Either returns the maximum value of any finite endpoint of any interval in the barcode
    /// OR, if no such endpoint exists, the return the user-supplied `default`
    pub fn max_finite_endpoint_or( &self, default: f64 ) -> f64 { 
        self.barcode.max_finite_endpoint().map(|x| x.into_inner()).unwrap_or( default )
    }    

}    





//  =========================================
//  OPTIMIZED CYCLE REPRESENTATIVES
//  =========================================


// #[pyclass]
// #[derive(Clone)]
// pub struct MinimalCyclePySimplexFilteredRational {
//     minimal_cycle:  MinimalCycle< SimplexFilteredPy  >
// }

// #[pymethods]
// impl MinimalCyclePySimplexFilteredRational {
//     /// The initial cycle, prior to optimization
//     pub fn cycle_initial(&self) -> Vec< (SimplexFilteredPy, f64) > { self.minimal_cycle.cycle_initial().clone() }
//     /// The optimized cycle
//     pub fn cycle_optimal(&self) -> Vec< (SimplexFilteredPy, f64) > { self.minimal_cycle.cycle_optimal().clone() }    
//     /// The objective value of the initial cycle
//     /// 
//     /// Equal to the L1 norm of the linear combination of simplices, where each entry is scaled by the diameter of the corresponding simplex.
//     pub fn objective_initial(&self) -> f64 { self.minimal_cycle.objective_initial().clone() }
//     /// The objective value of the optimal cycle
//     /// 
//     /// Equal to the L1 norm of the linear combination of simplices, where each entry is scaled by the diameter of the corresponding simplex.
//     pub fn objective_optimal(&self) -> f64 { self.minimal_cycle.objective_optimal().clone() }    
//     /// A chain whose boundary equals the difference between the initial cycle and the optimal cycle
//     pub fn bounding_difference(&self) -> Vec< (SimplexFilteredPy, f64) > { self.minimal_cycle.bounding_difference().clone() }
// }

// impl MinimalCyclePySimplexFilteredRational {
//     pub fn new( minimal_cycle: MinimalCycle< SimplexFilteredPy  >   ) -> Self { MinimalCyclePySimplexFilteredRational{ minimal_cycle} }
// }