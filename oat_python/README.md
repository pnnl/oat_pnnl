# Open Applied Topology in Python

Open Applied Topology in Python (OAT-Python) is a package for fast, user-friendly algebra and topology.  It provides Python wrappers for [oat_rust](https://github.com/OpenAppliedTopology/oat_rust), the Rust library Open Applied Topology.

## Documentation and tutorials

Documentation for OAT-Python is currently under development. 

**Python users** The best resources currently available are
1. The Jupyter notebook tutorials published on OAT's [main git repository](https://github.com/OpenAppliedTopology/OAT/oat), and 
2. Python's `help()` function.

**Rust Developers** If you are interested in modifying the Rust code in OAT-Python, see the API documenation for `oat_python` available at [Crates.io](https://crates.io).

**Python Developers** Documentation for the Python API is not available at this time (we are working to address this). In the meantime we have done our best to document the code, and encourage you to check out the source code directly.

## Python Installation

1. Download and install the most recent version of [Rust](https://www.rust-lang.org/).  Make sure your installation is up to date by running `rustup update` in a command shell.

2. Obtain a copy of OAT suite, which contains folders for OAT-Rust and OAT-Python.  Don't move the contents out of this folder.

3. Create a virtual Python environment, e.g. using Anaconda, or open one that you already have.  In this example, let's assume the environment name is `myenv`.  Activate `myenv`, and run

    ```bash
    pip install maturin
    ```

    A number of warning messages may appear; this is normal. 

    **If you already have maturin installed, make sure it is up to date!**

4. With `myenv` activated, CD into the `oat_python` folder, and run

    ```bash
    maturin develop --release
    ```
    
5. OAT-Python should now be installed.  Try running the Jupyter notebooks with `myenv`!


## Contributing

For information on **contributing**, see [`CONTRIBUTING`](https://github.com/OpenAppliedTopology/oat_python/blob/main/CONTRIBUTING).

## License

For inforamtion on copyright and licensing, see [`LICENSE.md`](https://github.com/OpenAppliedTopology/oat_python/blob/main/LICENSE.md) and [`DISCLAIMER.md`](https://github.com/OpenAppliedTopology/oat_python/blob/main/DISCLAIMER.md).

## Attributions

OAT is an extension of the ExHACT library. See [`ATTRIBUTIONS.md`](https://github.com/OpenAppliedTopology/oat_python/blob/main/ATTRIBUTIONS.md) for details.
