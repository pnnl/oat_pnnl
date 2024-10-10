//! Three step process for writing and running your own Rust code.
//! 
//! Developing code in Rust typically follows the following general workflow.
//! 
//! 1. Create a new crate.  
//! 2. Add program files and folders to the crate.
//! 3. Compile and run your program.
//! 
//! There's a lot to learn about this process when you're first starting out: you'll want to use additional packages, split programs into separate files,
//! compile code to executable, etc.  
//! 
//! The `intro_to_rust` crate, which ships with the OAT suite, can help get you started.
//! - The `README.md` file for `intro_to_rust` has many useful tips
//! - Feel free to copy parts of the demo for your own project!
//! 



// //! # Overview
// //! 
// //! Developing code in Rust typically follows the following general workflow.
// //! 
// //! 1. Create a new "crate."  In concrete terms, a crate is just a folder on your computer that contains
// //! some Rust files.  To start a new crate, open a shell, cd into the folder where you want your
// //! crate to be located, and run `cargo new my_crate`.  This will generate a new folder
// //! and fill it with a few files to help you get started.
// //! 
// //! 2. Add new Rust files to the crate.  Use these files to build the programs you want to
// //! execute.
// //! 
// //! 3. Once you've written your program(s), executing them is a two-step process.  First
// //! you "compile" the code; this means running a command that will read your program
// //! files and generate a new file (or set of files) called "binary".  Second, you run
// //! the binary.  



// //! Rust offers a number of commands for these two steps, depending on your particular needs.
// //! Some important points to be aware of include
// //! 
// //! - debug compilation: the default for several commands is to compile a "debug" version of
// //! your program.  The debug version offers very helpful options for finding and understanding
// //! errors in your code, but it runs *much* slower than the alternative, "release" version.
// //! 
// //! - library versus binary: there are two types of Rust crates: libraries and binaries.
// //! A single library crate can be used to create several different binaries, and you
// //! use different commands to compile, depending on which type you have.
//! 
//! # Detailed explanation
//! 
//! For a detailed explanation of code development in Rust, the best place to start is the first example in the online
//! [Rust Book](https://doc.rust-lang.org/book/ch01-00-getting-started.html)
