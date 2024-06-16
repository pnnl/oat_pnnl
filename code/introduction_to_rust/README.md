
# Introduction to Rust libraries

This is a gental introduction-by-example to Rust libraries.  The repository is itself Rust library, which you can

- download
- run
- explore
- modify
- learn from
- copy!

## How to read this document

Blocks of text formated as follows

```bash
here is some text
```

should be interepreted as commands that you can run in a terminal.

## Create a new library

There are two types of crates: binary and library.  This is a library crate.  To create a new library, first make sure that [Rust is installed](https://www.rust-lang.org/tools/install).  Then open a shell, cd to the directory where you want the crate to be located, and run `cargo new my_crate --lib` (replace `my_crate` with whatever you'd like the crate to be called).

## Build the documentaiton

Cd into the crate (i.e., into this repository) and run

```bash
cargo doc --no-deps --open
```

This will open the documentation home page in a new browser window.  Explore the files inside this crate to see how comments in different sections of the source code get turned into text in the documentation.

## Run a program file

Each file and (and each folder) inside `src/bin` can be compiled and run.  For example, `src/bin` contains a file called `introduce_flowers.rs`.  To run this program, cd into the crate, and run

```bash
cargo run --bin introduce_flowers
```

It should print

```bash
I am a daisy.
I am a violet.
I am a cherry blossom.
```

You can substitute `introduce_flowers` with `introduce_trees` or `introduce_numbers` in order to run the files `introduce_trees` or `introduce_numbers`.


## Run a program file *faster*

The programs in this library are short and fast.  But more complicated programs might take longer.  Rust allows you to compile versions of a program that run significantly faster, using the `--release` option.  For example, to run a release version of `introduce_flowers.rs`, use

```bash
cargo run --bin introduce_flowers --release
```

## Make an executable file

To compile a program without running it, replace `run` with `build`, e.g.

```bash
cargo build --bin introduce_flowers --release
```

The executable file can then be found under `src/target/release/introduce_flowers`.  If instead we ran

```bash
cargo build --bin introduce_flowers
```

then executable would be saved to `src/target/debug/introduce_flowers`.  The files in the debug folder run much slower.

## Organize programs into files and folders

The file structure of this crate is

```bash
.
|____Cargo.toml
|____Cargo.lock
|____README.md
|____src
| |____flower_functions
| | |____daisy.rs
| | |____violet.rs
| | |____blossom_functions
| | | |____cherry.rs
| | | |____mod.rs
| | |____mod.rs
| |____bin
| | |____introduce_flowers.rs
| | |____introduce_trees.rs
| | |____introduce_numbers.rs
| |____lib.rs
| |____oak.rs
| |____willow.rs

```

Here are some (oversimplified) guidelines for the organization of files and folders:

- files that you wish to turn into executables go in `src/bin/`
- other program files go in `src/`; each file shoud end in `.rs`
- you can group files in `src/` by placing them into folders; each folder should contain a `mod.rs` file (c.f. the `blossoms` folder)
- the `mod.rs` file lists other files within the same folder, which you want to make "visible" to other files outside the folder
- the `lib.rs` is essentially the same as a `mod.rs` file; it just gets a special name because it sits in `src/` and not a subfolder of `src/`.

See also this page on [project layout](https://doc.rust-lang.org/cargo/guide/project-layout.html).  Note that library crates do not have a `main.rs` file, unlike binary crates.

## Use programs from other crates

Suppose we wanted to use the [abs](https://docs.rs/num/0.4.0/num/fn.abs.html) function from the [num](https://docs.rs/num/0.4.0/num/index.html) crate in our `introduce_numbers.rs` file.  To do so, we need to take two steps

1. Modify the `Cargo.toml` file by adding a new line under `[dependencies]`, e.g.

    ```bash
    [dependencies]
    num = "0.4.0"
    ```

   This method works for adding crates that are available from the online registry [crates.io](https://crates.io/).  If you want to import a crate from your local computer, the syntax is

   ```bash
   [dependencies]
   other_crate = { path = "file/path/to/other_crate" }
   ```

   See the [Rust Book](https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html) for details.  

2. Place a `use` statement in any file where you wish to use features from the imported packaged.  See for example `src/bin/introduce_fractions.rs`

## Delete all compiled files, as well as documentation

```bash
cargo clean
```

This can be helpful if you get a compilation error (which is rare but does happen), made a mistake compiling the documentation, or want to reduce the size of the folder by deleting everything but source code.

## Additional resources

- [difference between library crates and binary crates](https://stackoverflow.com/questions/62365876/whats-the-difference-between-binary-and-library-in-rust)
- [project layout](https://doc.rust-lang.org/cargo/guide/project-layout.html)
- [organizing code](https://rust-classes.com/chapter_4_3.html)
- [running a binary file in a library crate](https://www.reddit.com/r/rust/comments/tiaor0/is_there_any_way_to_compile_and_run_a_single_rust/). 
