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

To use your custom qr_mumps, set the environment variable `JULIA_QRMUMPS_LIBRARY_PATH`
to point to the folder holding the qr_mumps shared libraries before `using QRMumps`.
Note that the same version of qr_mumps as used by the `qr_mumps_jll` artifact is needed.
To initialize qr_mumps with a custom installation, the function [`qrm_init`](https://juliasmoothoptimizers.github.io/QRMumps.jl/stable/api/#QRMumps.qrm_init) must be called prior to any other functions of QRMumps.jl.

For example:
```bash
brew tap dpo/mumps-jl
brew install qr_mumps
export JULIA_QRMUMPS_LIBRARY_PATH=$(brew --prefix)/opt/qr_mumps/lib
```

Apple Silicon users should remember to use `arch x86_64 brew` to refer to Intel binaries run through Rosetta, as we do not (yet) ship Silicon binaries of qr_mumps via Homebrew.

The `JULIA_QRMUMPS_LIBRARY_PATH` environment variable may be set permanently in the shell's startup file, or in `$HOME/.julia/config/startup.jl`.

## References

> Emmanuel Agullo, Alfredo Buttari, Abdou Guermouche, and Florent Lopez (2016).
> Implementing multifrontal sparse solvers for multicore architectures with sequential task flow runtime systems.
> ACM Trans. Math. Softw., 43(2):13:1–13:22.
> [10.1145/2898348](http://dx.doi.org/10.1145/2898348)

> Emmanuel Agullo, Alfredo Buttari, Abdou Guermouche, and Florent Lopez (2015).
> Task-based multifrontal QR solver for GPU-accelerated multicore architectures.
> In HiPC, pages 54–63. IEEE Computer Society.
> [10.1109/HiPC.2015.27](http://dx.doi.org/10.1109/HiPC.2015.27)

## Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/QRMumps.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers) organization, so questions about any of our packages are welcome.
