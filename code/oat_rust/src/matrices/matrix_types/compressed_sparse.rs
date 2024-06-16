

//  NOTES
//
//  Why have this data type?
//  * this overall structure for data is pretty common/standard
//  * why have this instead of just using standard csc/csr
//      * sometimes there is a natural bijection between *one* set of indices and {0, ..., n}, but
//      not the other set of indices.  This occurs naturally in the persistent cohomology
//      algorithm, where rows are naturally indexed by integers, but assigning integers to column
//      indices can be quite laborious.
