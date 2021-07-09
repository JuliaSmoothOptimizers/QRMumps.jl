# A [Julia](http://julialang.org) Interface to [qr_mumps](http://buttari.perso.enseeiht.fr/qr_mumps/)

| **Documentation** | **Linux/macOS/Windows/FreeBSD** | **Coverage** | **DOI** |
|:-----------------:|:-------------------------------:|:------------:|:-------:|
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] | [![build-gh][build-gh-img]][build-gh-url] [![build-cirrus][build-cirrus-img]][build-cirrus-url] | [![codecov][codecov-img]][codecov-url] | [![doi][doi-img]][doi-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaSmoothOptimizers.github.io/QRMumps.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://JuliaSmoothOptimizers.github.io/QRMumps.jl/dev
[build-gh-img]: https://github.com/JuliaSmoothOptimizers/QRMumps.jl/workflows/CI/badge.svg?branch=main
[build-gh-url]: https://github.com/JuliaSmoothOptimizers/QRMumps.jl/actions
[build-cirrus-img]: https://img.shields.io/cirrus/github/JuliaSmoothOptimizers/QRMumps.jl?logo=Cirrus%20CI
[build-cirrus-url]: https://cirrus-ci.com/github/JuliaSmoothOptimizers/QRMumps.jl
[codecov-img]: https://codecov.io/gh/JuliaSmoothOptimizers/QRMumps.jl/branch/main/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/JuliaSmoothOptimizers/QRMumps.jl
[doi-img]: https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4716114-blue.svg
[doi-url]: https://doi.org/10.5281/zenodo.4716114

## How to Cite

If you use QRMumps.jl in your work, please cite using the format given in [`CITATION.bib`](https://github.com/JuliaSmoothOptimizers/QRMumps.jl/blob/main/CITATION.bib).

## How to install

```julia
julia> ]
pkg> add QRMumps
pkg> test QRMumps
```

## Content

[qr_mumps](http://buttari.perso.enseeiht.fr/qr_mumps/) is a software package for the solution of sparse, linear systems on multicore computers.
It implements a direct solution method based on the QR or Cholesky factorization of the input matrix. 
Therefore, it is suited to solving sparse least-squares problems, to computing the minimum-norm solution of sparse, underdetermined problems and to solving symmetric, positive-definite sparse linear systems.
It can obviously be used for solving square unsymmetric problems in which case the stability provided by the use of orthogonal transformations comes at the cost of a higher operation count with respect to solvers based on, e.g., the LU factorization such as [MUMPS](http://mumps.enseeiht.fr/).
It supports real and complex, single or double precision arithmetic.

## Custom Installation

**Note: qr_mumps is already precompiled with Yggdrasil for all platforms.**

To use your custom qr_mumps, set the environmental variables `JULIA_QR_MUMPS_LIBRARY_PATH`
to point the shared library. Note that **qr\_mumps** version 3.0.2 is needed.

**Very important note: you must set these environment variables before
calling `using QRMumps` in every Julia session.**

For example:
```julia
ENV["JULIA_QR_MUMPS_LIBRARY_PATH"] = "~/Applications/qr_mumps-3.0.2/build/lib"
using QRMumps
```

Alternatively, you can set these permanently through your operating system.
