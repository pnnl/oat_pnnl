//! Gentle introduction to traits
//! 
//! One way to sort a list `L = [1,2,3]` in Python is to call `L.sort()`.
//! The function `.sort()` is an example of a *method*, which many
//! programming languages support.
//! 
//! Most languages treat the methods defined on two different objects as
//! different, even if they have the same name.  For example the `.sort()`
//! method used on lists is technically different from the `.sort()` method 
//! used on Numpy arrays.  However, sometimes it can be useful to write code
//! that will work with *any* object that has a `.sort()` method.
//! 
//! This is where traits come in.  A `trait` in Rust essentially consists of
//! - a function name, e.g. `my_method`, and 
//! - a function signature, which specifies the inputs that the function `my_method` should take
//! and the outputs that it should return.
//!
//! After defining a trait, you can define a method with that trait name
//! on multiple different objects -- as long as the method has the correct
//! function signature.  *Then you can write a function `F` that takes any
//! input `I` such that `I` has a method called `my_method`.*
//! 
//! This becomes important in OAT because the most basic operations on matrices, vectors, and
//! vector entries are all governed by trait methods.  For example, looking up a 
//! row or column of matrix is a trait method.
//! 
//! For more information on traits, see the [introduction to traits in the Rust book](https://doc.rust-lang.org/book/ch10-02-traits.html).