//! Tips for debuggin Rust
//! 
//! 
//! 
//! # cargo clean
//! 
//! Every once in a while the Rust compiler may generate a bug.
//! Often when this happens you'll receive an error message saying there is a problem with the compiler.
//! Moreover, the problem may persist even when you recompile.
//! In these cases, it may help to run `cargo clean`.  This will clear all compiled code.
//! 
//! # determining the type of an object
//! 
//! Sometimes it can be difficult to determine the type of a variable; this [stack exchange article](https://stackoverflow.com/questions/21747136/how-do-i-print-the-type-of-a-variable) suggests several solutions.