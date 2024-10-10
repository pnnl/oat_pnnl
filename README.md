# Welcome to Open Applied Topology at PNNL!

Welcome the Pacific Northwest National Laboratory (PNNL) central repository for Open Applied Topology (OAT)!

OAT is a collection of open-source software packages for applied topology.

- [OAT-Python](oat_python) Provides a user-friendly Python module for persistent homology, cycle representatives, and more.
- [OAT-Rust](oat_rust) Provides a high-performance backend for OAT-Python. Also a standalone library which allows users to compile programs as executable.
- [OAT-Jupyter](oat_jupyter) Provides a variety of Jupyter notebook tutorials.
- [Introduction to Rust](introduction_to_rust_pnnl) Provides a gentle introduction to programming in Rust.


## Installation and documentation

Steps to install each package and access its documentation can be found in the corresponding README.md.

See each package's README for installation instructions. 

## Tutorials

Jupyter notebook tutorials can be found in the [oat_jupyter](oat_jupyter) folder. An [OAT-Python](oat_python) installation is needed to run these notebooks; see the README.md file in the [oat_python](oat_python) for installation instructions.

## Updates

This repository contains the latest contributions from Pacific Northwest National Laboratory, but not the OAT community as a whole. **For the latest code, or to submit an issue report, please visit OAT's [homepage](https://openappliedtopology.github.io).** or visit one of the corresponding repositories

# Developers

OAT is developed by a multi-institution collaboration including researchers from twenty universities, governement, and industries. The lead contributor is currently Pacific Northwest National Laboratory, which maintains the following public repository for code and documentation (link TBD).

# Funding and attribution

Funding for this has been generously provided by the National Science Foundation through awards DMS-1854748, DMS-1854683. and DMS-1854703, and by Pacific Northwest National Laboratory.


OAT is an extension of ExHACT and SOLAR libraries available [here](https://github.com/ExHACT).  See `ATTRIBUTION.mdS` for details.

# Legacy  code

Earlier iterations of the OAT project separated the library into two halves, Sparse Oracle Linear Algebra in Rust (SOLAR) and Sparse Oracle Homological Algebra in Rust (SOHAR). Later iterations merged the Python and Rust repositories for SOHAR into the corresponding repositories for SOLAR, and renamed these repositories OAT-Rust and OAT-Python. The following repositories contain the legacy code from SOHAR, which is now deprecated.

- SOHAR-Rust [sohar_rust_pnnl](https://github.com/pnnl/sohar_rust_pnnl)
- SOHAR-Python [sohar_python_pnnl](https://github.com/pnnl/sohar_python_pnnl)