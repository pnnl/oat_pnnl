
use crate::algebra::matrices::query::{   ViewRow,
                                        ViewRowAscend,
                                        ViewRowDescend,
                                        ViewCol, 
                                        ViewColAscend,
                                        ViewColDescend,
                                        WhichMajor,
                                        MajorDimension};
use std::iter;


pub struct VecCsv< ColIndex, Coefficient >
{
    major_dimension: MajorDimension, 
    min_ind: Vec< Vec< ColIndex > > ,
    snz_val: Vec< Vec< Coefficient > >
}


impl    < ColIndex, Coefficient >
        VecOfVec 
        < ColIndex, Coefficient > 
{
    // Make new (empty) VecOfVec. 
    pub fn new( major_dimension: MajorDimension ) -> Self  
    {
        atrixOracle { scalar: scalar,
                             major_dimension: major_dimension,
                             phantom: PhantomData 
                            }
    }
}



