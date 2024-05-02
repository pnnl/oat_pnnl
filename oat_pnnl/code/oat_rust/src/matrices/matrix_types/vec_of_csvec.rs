
use crate::matrices::matrix_oracle_traits::{   OracleMajor,
                                        OracleMajorAscend,
                                        OracleMajorDescend,
                                        OracleMinor, 
                                        OracleMinorAscend,
                                        OracleMinorDescend,
                                        WhichMajor,
                                        MajorDimension};
use std::iter;


pub struct VecCsv< KeyMin, SnzVal >
{
    major_dimension: MajorDimension, 
    min_ind: Vec< Vec< KeyMin > > ,
    snz_val: Vec< Vec< SnzVal > >
}


impl    < KeyMin, SnzVal >
        VecOfVec 
        < KeyMin, SnzVal > 
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



