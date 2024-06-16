# oat_python

This package is a Python wrapper for the library Open Applied Topology in Rust (oat_rust), which provides performant algorithms and data structures for applied algebraic topology.  It is part of the Open Applied Topology (OAT) project.

## Installation

1. Download and install the most recent version of Rust.
2. Obtain copies of oat_rust and oat_python.
3. Open `oat_python/Cargo.toml` and update

    ```bash
    oat_rust = { path = "file/path/to/oat_rust" }
    ```

    with the path to your local oat_rust folder.

4. Create a virtual python environment `myenv`, e.g. using Anaconda.  Activate the environment, and run

    ```bash
    pip install maturin
    ```

    A number of warning messages may appear; this is normal.

5. With `myenv` activated, CD into the `oat_python` folder, and run

    ```bash
    maturin develop --release
    ```

6. oat_python should now be installed.  Try running the following to get the unique elements in a list of lists:

    ```python
    >>> import oat_python as oat
    >>> oat.unique_rows( [ [0,0], [0,0], [0,0] ] )
    [[0, 0]]
    ```
