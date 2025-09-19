# Welcome to Open Applied Topology at PNNL!

Welcome the Pacific Northwest National Laboratory (PNNL) central repository for Open Applied Topology (OAT)!

This repository contains the latest contributions from Pacific Northwest National Laboratory, but not the latest version of OAT. **For the latest code, or to submit an issue report, please visit OAT's [homepage](https://openappliedtopology.github.io).** or visit one of the corresponding repositories

# Obtaining the code

The oat_pnnl code base has two major components, each with its own dedicated git repository

- [oat_rust_pnnl](https://github.com/pnnl/oat_rust_pnnl) (a fork of [oat_rust](https://github.com/OpenAppliedTopology/oat_rust))
- [oat_python_pnnl](https://github.com/pnnl/oat_python_pnnl) (a fork of [oat_python](https://github.com/OpenAppliedTopology/oat_python))

You can obtain these components directly from their github repositories, or from this repository. The folders labeled [oat_rust_pnnl](oat_rust_pnnl) and [oat_python_pnnl](oat_python_pnnl) inside [oat_pnnl](https://github.com/pnnl/oat_pnnl) are actually git [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) which have been set to track the corresponding standalone repositories.  Accessing the code in these submodules requires a few additional steps beyond `git clone`, but don't worry! The process is simple and well documented in the git [literature](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

# Developers

OAT is developed by a multi-institution collaboration including researchers from twenty universities, governement, and industries. The lead contributor is currently Pacific Northwest National Laboratory.

# Funding and attribution

Funding for this has been generously provided by the National Science Foundation through awards DMS-1854748, DMS-1854683. and DMS-1854703, and by Pacific Northwest National Laboratory.

OAT is an extension of ExHACT and SOLAR libraries available [here](https://github.com/ExHACT).  See `ATTRIBUTIONS.md` for details.

# Legacy  code

Earlier iterations of the OAT project separated the library into two halves, Sparse Oracle Linear Algebra in Rust (SOLAR) and Sparse Oracle Homological Algebra in Rust (SOHAR). Later iterations merged the Python and Rust repositories for SOHAR into the corresponding repositories for SOLAR, and renamed these repositories OAT-Rust and OAT-Python. The following repositories contain the legacy code from SOHAR, which is now deprecated.

- SOHAR-Rust [sohar_rust_pnnl](https://github.com/pnnl/sohar_rust_pnnl)
- SOHAR-Python [sohar_python_pnnl](https://github.com/pnnl/sohar_python_pnnl)


This respository also contains deprecated code from prior versions of the project in the [deprecated](deprecated) folder. This includes

- [introduction_to_rust](introduction_to_rust) A small repository illustrating basic design of a Rust project. This code has been deprecated, as large language models now provide a better resource.
- [oat_jupyter](oat_jupyter) A collection of jupyter notebooks illustrating basic functionality of the OAT library, published circa 2023. These notebooks have now been absorbed into the Sphinx documentation for oat_python under *Tutorials*.