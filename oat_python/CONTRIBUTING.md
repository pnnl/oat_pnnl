# Open Applied Topology


Welcome to the OAT community!  We're glad you're here.

Users can find a great sense of satisfaction and accomplishment in helping fellow users and/or modifying open source software; that includes

- adding, subtracting, or changing the code
- catching typos and clarifying explanations
- joining discussions
- reporting problems
- and more!

Here's some information to get started.

- [Code of conduct](#code-of-conduct)
- [Get help](#get-help)  
- [Project overview](#project-overview)
- [How-to](#general-tips)
  - [Report an issue or request a feature](#report-an-issue-or-requiest-a-feature)
  - [Contribute new code](#contribute-new-code)  
- [Introduction to Rust](#general-tips)
- [Introduction to PyO3](#general-tips)



# Code of conduct

Your safety and well being are the top priorities of the OAT community.  The same holds for every one of our users.  See the [code of conduct](./CODE_OF_CONDUCT.md) for what that means.

# Get help

If you're stuck or don't know where to begin, then you're in good company -- we've all been there!  We're here to help, and we'd love to hear from you:

- open a issue report on Github
- email us at <gregory.roek@pnnl.gov>

# Project overview

- OAT
  - A high-performance, low-level software package written in pure Rust
  - Registered on `crates.io` as `oat_rust`
- OAT-Python
  - Powerful tools for interactive visualization and analysis.  Provides Python bindings for `oat_rust` using [pyO3](https://pyo3.rs/).
  - Registered on `PyPi` as `oat_python`
  - Registered on `crates.io` as `oat_python`
  - This package has 
    - a Rust component, stored in `oat_python/src`, and 
    - a Python component, stored in `oat_python/oat_python`. 
- Tutorials
  - Jupyter notebook tutorials, available on colab

These components have the following folder structure

```
tutorials       <-- Jupyter notebook tutorials

oat_rust        
├── src         <-- OAT source code, in Rust
└── developer   <-- documents and resources for devlopers

oat_python
├── oat_python       <-- OAT-Python source code, in Python
└── src         <-- OAT-Python source code, in Rust
```

# How to

The world of open source is wide; it can be a lot to take in!  If this is your first time working with an open source project, then welcome!  If you're an experienced contributor, then welcome back!  Either way, this online resource might help you [get oriented, or bursh up](https://opensource.guide/how-to-contribute/) on the process.

## Report an issue or request a feature

Here are the [steps to creating an issue on github](https://docs.github.com/en/issues/tracking-your-work-with-issues/quickstart)

- search for related issues on Github. You might be able to get answer without the hassle of creating an issue
- describe the current behavior and explain which behavior you expected to see instead and why. At this point you can also tell which alternatives do not work for you.  
  - (if applicable) provide error messages
  - (if applicable) provide a step by step description of the problem; if possible include code that others can use to reproduce it
  - You may want to [include screenshots and animated GIFs](https://www.cockos.com/licecap/) which help you demonstrate the steps or point out the part which the suggestion is related to. You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast) or [this tool](https://github.com/GNOME/byzanz) on Linux.
  - provide clear, specific title
  - include details on your setup (operating system, python version, etc.)
- use the most recent version of this library and the source language (e.g. Rust); that fixes a lot of problems  
- here are [more details on getting the most out of issue reporting!](https://marker.io/blog/how-to-write-bug-report)

## Contribute new code

Here is a [step-by-step guide to writing new code, and submiting it to the project](https://docs.github.com/en/get-started/quickstart/contributing-to-projects)

The more you know about a software library, the easier it is to get started writing code.  The best way to learn about this project is its the documentation!  See `README.md` to get started.





# Introduction to Rust

Rust is a low-level programming language with powerful features for scientific computing, e.g. memory safety and helpful error messages.  It has been voted [most-loved language](https://insights.stackoverflow.com/survey/2021) by the worldwide developer community since 2015.

* **Installation**  The  [Rust website](https://www.rust-lang.org/learn/get-started) has directions!

* **Search for what you need in the documentation** All Rust documenation has a consistent format, which you can search in a web browser.  
  - Lists of objects, functions, etc., appear in two places: either the bottom of a page, or in the menu bar on the left.
  - You can also use the search bar at the top to pull up a list of related terms.  
  - The question mark button to the right of the bar gives other helpful tips (for example, you can search for functions based on their type signature).

* **Save up to 90% coding time, with VS Code** If you are new to Rust, we strongly recommend the [VS Code editor](https://code.visualstudio.com/docs/languages/rust), with the [`rust-analyzer`](https://marketplace.visualstudio.com/items?itemName=rust-lang.rust-analyzer) extension.  *Everyone we know uses VS code for Rust development*.  It has great features for

  - highlighting
  - smart debugging messages
  - command completion
  - bracket coloring; since type parameters are defined with `<>` brackets, it may be useful to add these to your `configuration.json` file, as per this [discussion](https://github.com/microsoft/vscode/issues/132476) and [example](https://github.com/microsoft/vscode/blob/997228d528d687fd17cbf3ba117a0d4f668c1393/extensions/javascript/tags-language-configuration.json#L11)

  

* **Debugging** Rust is very good about providing helpful error messages, as a rule.
It is often possible to find help just by copying these messages into a web search.
The OAT developers have also collected a short list of [debugging tips and tricks](crate::developer::rust_debugging), which you may find useful.

* **Tips**

  - long type definitions: because OAT is highly modular, its type can become quite verbose.  The Rust compiler offers a lot of helpful information about types when running your code, but it abbreviates some types when printing to terminal; note however, that whenever it does this, it also writes a full type description in a file, and prints the file path in the command shell.
  - the [`debugit`](https://docs.rs/debugit/latest/debugit/) crate is useful for debugging


# Introduction to PyO3

PyO3 is a software package to link Rust code with Python.  There are three good ways to learn about PyO3 for this project:

- [PyO3 homepage](https://pyo3.rs/)
- [PyO3 API](https://docs.rs/pyo3) - this contains slightly different information than the homepage
- the [OAT-Python git respository](https://pyo3.rs/) - many examples that aren't covered in the previous two sources can be found and copy/pasted from the source code!