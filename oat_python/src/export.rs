//! Export data to Python

use itertools::Itertools;
use num::rational::Ratio;
use ordered_float::OrderedFloat;

use pyo3::prelude::*;
use pyo3::types::PyDict;

use oat_rust::topology::simplicial::simplices::filtered::SimplexFiltered;
use sprs::CsMatBase;


//  =========================================================
//  THE TRAIT
//  =========================================================


/// Generic wrapper used to export objects to Python
/// 
/// It's often necessary to wrap an object `t: T` in a wrapper struct
/// `Wrapper` in order to implement the PyO3 `IntoPy` or `IntoPyObject` traits.
/// This struct offers a convenient way to do so: for given `T`:
/// 1. Implement `IntoPy` or `IntoPyObject` on `ForExport<T>`
/// 2. Then convert any object `t: T` to an object `wrapper: ForExport< T >`
/// via the `Export` trait; concretely `wrapper = t.export()`.
#[derive(Copy,Clone,Debug,Eq,PartialEq,Ord,PartialOrd)]
pub struct ForExport< T > {
    pub data: T,
}

/// Provides a method `export` for any type `T`
pub trait Export where Self: Sized {
    fn export( self ) -> ForExport< Self >;   
}

impl < T > Export for T {
    fn export( self ) -> ForExport< Self > { ForExport{ data: self } }    
}

//  =========================================================
//  DELETEABLE
//  =========================================================

pub fn export_ratio( r: & Ratio< isize >) -> PyResult<Py<PyAny>> {
    Python::with_gil(|py| {
        let frac = py.import("fractions")?;
        frac.call_method("Fraction", ( r.numer().clone(), r.denom().clone() ), None)
            .map(Into::into)
    })
}

#[pyfunction]
pub fn export_fraction(numer: isize, denom: isize) -> PyResult<Py<PyAny>> {
    Python::with_gil(|py| {
        let frac = py.import("fractions")?;
        frac.call_method("Fraction", (numer, denom), None)
            .map(Into::into)
    })
}

#[pyfunction]
pub fn export_fraction_with_current_gil<'py>(py: Python<'py>, numer: isize, denom: isize) -> PyResult<&'py PyAny> {
    let frac = py.import("fractions")?;
    frac.call_method("Fraction", (numer, denom), None)
}

#[pyfunction]
pub fn my_dict<'py>(py: Python<'py>) -> &'py PyDict {
    let d = PyDict::new( py );
    d.set_item( 1, 4 ).ok();
    d 
}


#[pyfunction]
pub fn return_frame<'py>(py: Python<'py>) -> ForExport< Vec< ( SimplexFiltered< OrderedFloat<f64> >, Ratio< isize > ) > >  { // PyResult<Py<PyAny>> { // 
    let v 
        = vec![ 
                ( 
                    SimplexFiltered{ filtration: OrderedFloat(0.0), vertices: vec![0] },
                    Ratio::new( 1, 1 ),
                ) 
            ];
    return v.export()
    // let v = v.export();
    // return v.into_py(py)
}


//  =========================================================
//  Ratio< i64 >
//  =========================================================


impl IntoPy< PyResult<Py<PyAny>> > for ForExport< Ratio< isize > > {
    fn into_py(self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let frac = py.import("fractions")?;
        let numer = self.data.numer().clone();
        let denom = self.data.denom().clone();        
        frac.call_method("Fraction", ( numer, denom ), None)
            .map(Into::into)
    }
}

impl ToPyObject for ForExport< Ratio< isize > > {
    fn to_object( &self, py: Python<'_>) -> Py<PyAny> {
        export_ratio( & self.data ).ok().unwrap()
    }
}

//  =========================================================
//  OrderedFloat< f64 >
//  =========================================================

impl IntoPy< f64 > for ForExport< OrderedFloat< f64 > > {
    fn into_py(self, _py: Python<'_>) -> f64 {
        self.data.into_inner()
    }
}

impl ToPyObject for ForExport< OrderedFloat< f64 > > {
    fn to_object( &self, py: Python<'_>) -> PyObject {
        self.data.clone().into_inner().to_object(py)
    }
}


//  =========================================================
//  SimplexFiltered
//  =========================================================


impl IntoPy
        < Py<PyAny> > for 
        
    ForExport
        < 
                SimplexFiltered< OrderedFloat<f64> > 
        > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn into_py(self, py: Python<'_>) -> Py<PyAny> {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.vertices() ).ok().unwrap();
        dict.set_item( "filtration",  self.data.filtration().into_inner() ).ok().unwrap();
        return dict.into()
    }
}

impl ToPyObject for 
        
    ForExport
        < 
                SimplexFiltered< OrderedFloat<f64> > 
        > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn to_object(&self, py: Python<'_>) -> PyObject {
        return self.clone().into_py(py);
    }
}



//  =========================================================
//  Vec< (SimplexFiltered, Ratio) >
//  =========================================================


impl IntoPy
        < Py<PyAny> > for 
        
    ForExport
        < Vec< ( SimplexFiltered< OrderedFloat<f64> >, Ratio< isize > ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn into_py(self, py: Python<'_>) -> Py<PyAny> {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.vertices() ).collect_vec() ).ok().unwrap();
        dict.set_item( "filtration",  self.data.iter().map( |(s,_)| s.filtration().into_inner() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone().export() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}

impl ToPyObject for 
        
    ForExport
        < Vec< ( SimplexFiltered< OrderedFloat<f64> >, Ratio< isize > ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn to_object(&self, py: Python<'_>) -> PyObject {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.vertices() ).collect_vec() ).ok().unwrap();
        dict.set_item( "filtration",  self.data.iter().map( |(s,_)| s.filtration().into_inner() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone().export() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}


//  =========================================================
//  Vec< SimplexFiltered >
//  =========================================================


impl IntoPy
        < Py<PyAny> > for 
        
    ForExport
        < 
            Vec< 
                SimplexFiltered< OrderedFloat<f64> > 
            > 
        > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn into_py(self, py: Python<'_>) -> Py<PyAny> {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |s| s.vertices() ).collect_vec() ).ok().unwrap();
        dict.set_item( "filtration",  self.data.iter().map( |s| s.filtration().into_inner() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}

impl ToPyObject for 
        
    ForExport
        < 
            Vec< 
                SimplexFiltered< OrderedFloat<f64> > 
            > 
        > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn to_object(&self, py: Python<'_>) -> PyObject {
        return self.clone().into_py(py);
    }
}


//  =========================================================
//  Vec< (SimplexFiltered, f64) >
//  =========================================================


impl IntoPy
        < Py<PyAny> > for 
        
    ForExport
        < Vec< ( SimplexFiltered< OrderedFloat<f64> >, f64 ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn into_py(self, py: Python<'_>) -> Py<PyAny> {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.vertices() ).collect_vec() ).ok().unwrap();
        dict.set_item( "filtration",  self.data.iter().map( |(s,_)| s.filtration().into_inner() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}

impl ToPyObject for 
        
    ForExport
        < Vec< ( SimplexFiltered< OrderedFloat<f64> >, f64 ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn to_object(&self, py: Python<'_>) -> PyObject {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.vertices() ).collect_vec() ).ok().unwrap();
        dict.set_item( "filtration",  self.data.iter().map( |(s,_)| s.filtration().into_inner() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}


//  =========================================================
//  Vec< (Vec<isize>, Ratio) >
//  =========================================================


impl IntoPy
        < Py<PyAny> > for 
        
    ForExport
        < Vec< ( Vec<isize>, Ratio< isize > ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn into_py(self, py: Python<'_>) -> Py<PyAny> {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.clone() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone().export() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}

impl ToPyObject for 
        
    ForExport
        < Vec< ( Vec<isize>, Ratio< isize > ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn to_object(&self, py: Python<'_>) -> PyObject {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.clone() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone().export() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}


//  =========================================================
//  Vec< (Vec<isize>, f64) >
//  =========================================================


impl IntoPy
        < Py<PyAny> > for 
        
    ForExport
        < Vec< ( Vec<isize>, f64 ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn into_py(self, py: Python<'_>) -> Py<PyAny> {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.clone() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}

impl ToPyObject for 
        
    ForExport
        < Vec< ( Vec<isize>, f64 ) > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn to_object(&self, py: Python<'_>) -> PyObject {
        let dict = PyDict::new(py);
        dict.set_item( "simplex",     self.data.iter().map( |(s,_)| s.clone() ).collect_vec() ).ok().unwrap();
        dict.set_item( "coefficient", self.data.iter().map( |(_,z)| z.clone() ).collect_vec() ).ok().unwrap();
        let pandas = py.import("pandas").ok().unwrap();       
        pandas.call_method("DataFrame", ( dict, ), None)
            .map(Into::into).ok().unwrap()
    }
}



//  =========================================================
//  CSMAT
//  =========================================================


impl IntoPy
        < Py<PyAny> > for 
        
    ForExport
        < CsMatBase< Ratio<isize>, usize, Vec<usize>, Vec<usize>, Vec<Ratio<isize>> > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn into_py(self, py: Python<'_>) -> Py<PyAny> {

        let shape = self.data.shape();
        let (indptr, indices, data) = self.data.into_raw_storage();
        let data =  data.into_iter().map(|x| x.export().into_py(py).ok().unwrap() ).collect_vec();

        let sparse = py.import("scipy.sparse").ok().unwrap();
        return sparse.call_method("csr_matrix", 
            (
                (
                    data,
                    indices,
                    indptr,
                ),
                shape
            ), 
            None
        ).map(Into::into).ok().unwrap()
    }
}

impl ToPyObject for 
        
    ForExport
        < CsMatBase< Ratio<isize>, usize, Vec<usize>, Vec<usize>, Vec<Ratio<isize>> > > 
{
    /// Returns a Padas data frame with columns `simplex`, `filtration`, and `coefficient`
    fn to_object(&self, py: Python<'_>) -> PyObject {
        return self.clone().into_py(py)
    }
}