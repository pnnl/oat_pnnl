
//! General advice when working with matrices
//! 
//! # Many matrix objects have extra functionality
//! 
//! *You don't always have to use the oracle trait*  The oracle traits are powerful tools for building code that
//! can handle any type of sparse matrix data structure.  However there are many situations where only a single
//! data structure will ever be used.  That data structure may have a much wider range of methods and data than
//! can be accessed via the oracle traits.  It can be a great advantage to make use of this.