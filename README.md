# A [Julia](http://julialang.org) Interface to [qr_mumps](http://buttari.perso.enseeiht.fr/qr_mumps/)

| **Documentation** | **Linux/macOS/Windows/FreeBSD** | **Coverage** |
|:-----------------:|:-------------------------------:|:------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSmoothOptimizers.github.io/QRMumps.jl/dev) | ![CI](https://github.com/JuliaSmoothOptimizers/QRMumps.jl/workflows/CI/badge.svg?branch=master) [![Build Status](https://img.shields.io/cirrus/github/JuliaSmoothOptimizers/QRMumps.jl?logo=Cirrus%20CI)](https://cirrus-ci.com/github/JuliaSmoothOptimizers/QRMumps.jl) | [![codecov.io](https://codecov.io/github/JuliaSmoothOptimizers/QRMumps.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaSmoothOptimizers/QRMumps.jl?branch=master) |

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
