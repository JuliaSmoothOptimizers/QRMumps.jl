module QRMumps

using Libdl, SparseArrays, LinearAlgebra

import Base: \
import LinearAlgebra: mul!

if haskey(ENV, "JULIA_QRMUMPS_LIBRARY_PATH")
  const libsqrm = joinpath(ENV["JULIA_QRMUMPS_LIBRARY_PATH"], "libsqrm.$dlext")
  const libdqrm = joinpath(ENV["JULIA_QRMUMPS_LIBRARY_PATH"], "libdqrm.$dlext")
  const libcqrm = joinpath(ENV["JULIA_QRMUMPS_LIBRARY_PATH"], "libcqrm.$dlext")
  const libzqrm = joinpath(ENV["JULIA_QRMUMPS_LIBRARY_PATH"], "libzqrm.$dlext")
  const libqrm_common = joinpath(ENV["JULIA_QRMUMPS_LIBRARY_PATH"], "libqrm_common.$dlext")
  const QRMUMPS_INSTALLATION = "CUSTOM"
else
  using OpenBLAS32_jll
  using qr_mumps_jll
  const QRMUMPS_INSTALLATION = "YGGDRASIL"
end

function __init__()
  if QRMUMPS_INSTALLATION == "YGGDRASIL"
    qrm_init()
    if VERSION ≥ v"1.10"
      config = LinearAlgebra.BLAS.lbt_get_config()
      if !any(lib -> lib.interface == :lp64, config.loaded_libs)
        LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
      end
    end
  end
end

include("wrapper/qr_mumps_common.jl")
include("wrapper/libqrmumps.jl")
include("wrapper/qr_mumps_api.jl")

include("utils.jl")

export qrm_spmat, qrm_spfct,
    qrm_init, qrm_finalize,
    qrm_spmat_init!, qrm_spmat_init, qrm_spmat_destroy!,
    qrm_spfct_init!, qrm_spfct_init, qrm_spfct_destroy!,
    qrm_analyse!, qrm_analyse,
    qrm_update!, qrm_factorize!,
    qrm_solve!, qrm_solve,
    qrm_apply!, qrm_apply,
    qrm_spfct_get_rp, qrm_spfct_get_cp, qrm_spfct_get_r,
    qrm_spmat_mv!, mul!, qrm_spmat_nrm,
    qrm_vecnrm!, qrm_vecnrm,
    qrm_spbackslash!, qrm_spbackslash, \,
    qrm_spposv!, qrm_spposv,
    qrm_least_squares!, qrm_least_squares,
    qrm_least_squares_semi_normal!, qrm_least_squares_semi_normal,
    qrm_min_norm!, qrm_min_norm,
    qrm_min_norm_semi_normal!, qrm_min_norm_semi_normal,
    qrm_residual_norm!, qrm_residual_norm,
    qrm_residual_orth!, qrm_residual_orth,
    qrm_refine!, qrm_refine, qrm_set, qrm_get,
    qrm_user_permutation!

@doc raw"""
    qrm_init(ncpu, ngpu)

This routine initializes qr\_mumps and should be called prior to any other qr\_mumps routine.
This function is automatically called if you use qr\_mumps precompiled with Yggdrasil.

    qrm_init()

`ncpu` and `ngpu` are optional arguments and their default value are, respectively, `1` and `0`.

#### Input Arguments :

* `ncpu`: number of working threads on CPU.
* `ngpu`: number of working threads on GPU.
"""
function qrm_init end

@doc raw"""
    qrm_finalize()

This routine finalizes qr\_mumps and no other qr\_mumps routine should be called afterwards.
"""
function qrm_finalize end

@doc raw"""
    qrm_spmat_init!(spmat, A; sym=false)
    spmat = qrm_spmat_init!(spmat, m, n, rows, cols, vals; sym=false)

This routine initializes a **qrm_spmat** type data structure from a **sparseMatrixCSC**.

#### Input Arguments :

In the first form,

* `spmat`: the **qrm_spmat** sparse matrix to be initialized.
* `A` : a Julia sparse matrix stored in either **SparseMatrixCSC** or **SparseMatrixCOO** format (see SparseMatricesCOO.jl for the second case).
* `sym` : a boolean to specify if the matrix is symmetric / hermitian (true) or unsymmetric (false).

In the second form, the matrix `A` is specified using

* `m`: the number of rows.
* `n`: the number of columns.
* `rows`: the array of row indices of nonzero elements.
* `cols`: the array of column indices of nonzero elements.
* `vals`: the array of values of nonzero elements.
"""
function qrm_spmat_init! end

@doc raw"""
    spmat = qrm_spmat_init(A; sym=false)
    spmat = qrm_spmat_init(m, n, rows, cols, vals; sym=false)
"""
function qrm_spmat_init end

@doc raw"""
    qrm_spmat_destroy!(spmat)

This routine cleans up a **qrm_spmat** type data structure.

#### Input Argument :

* `spfct`: the sparse matrix to be destroyed.
"""
function qrm_spmat_destroy! end

@doc raw"""
    qrm_spfct_init!(spmat, spfct)

This routine initializes a **qrm_spfct** type data structure.
This amounts to setting all the control parameters to the default values.

#### Input Arguments :

* `spmat`: the input matrix of type `qrm_spmat`.
* `spfct`: the sparse factorization object to be initialized.
"""
function qrm_spfct_init! end

@doc raw"""
    spfct = qrm_spfct_init(spmat)
"""
function qrm_spfct_init end

@doc raw"""
    qrm_spfct_destroy!(spfct)

This routine cleans up a **qrm_spfct** type data structure by deleting the result of a sparse factorization.

#### Input Argument :

* `spfct`: the sparse factorization object to be destroyed.
"""
function qrm_spfct_destroy! end

@doc raw"""
    qrm_update!(spmat, A)
    qrm_update!(spmat, vals)

This routine updates a **qrm_spmat** type data structure from a **sparseMatrixCSC**.
`spmat` and `A` must have the same sparsity pattern.

#### Input Arguments :

In the first form,

* `spmat`: the **qrm_spmat** sparse matrix to be updated.
* `A` : a Julia sparse matrix stored in **SparseMatrixCSC** format.

In the second form,

* `vals`: the array of values of nonzero elements of `A`.
"""
function qrm_update! end

@doc raw"""
    qrm_analyse!(spmat, spfct; transp='n')

This routine performs the analysis phase and updates spfct.

#### Input Arguments :

* `spmat`: the input matrix of type `qrm_spmat`.
* `spfct`: the sparse factorization object of type `qrm_spfct`.
* `transp`: whether the input matrix should be transposed or not. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_analyse! end

@doc raw"""
    spfct = qrm_analyse(spmat; transp='n')
"""
function qrm_analyse end

@doc raw"""
    qrm_factorize!(spmat, spfct; transp='n')

This routine performs the factorization phase. It can only be executed if the analysis is already done.

#### Input Arguments :

* `spmat`: the input matrix of type `qrm_spmat`.
* `spfct`: the sparse factorization object of type `qrm_spfct`.
* `transp`: whether the input matrix should be transposed or not. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_factorize! end

@doc raw"""
    qrm_solve!(spfct, b, x; transp='n')

This routine solves the triangular system `Rx = b` or `Rᵀx = b`. It can only be executed once the factorization is done.

#### Input Arguments :

* `spfct`: the sparse factorization object resulting from the qrm_factorize! function.
* `b`: the right-hand side(s).
* `x`: the solution vector(s).
* `transp`: whether to solve for R, Rᵀ or Rᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_solve! end

@doc raw"""
    x = qrm_solve(spfct, b; transp='n')
"""
function qrm_solve end

@doc raw"""
    qrm_apply!(spfct, b; transp='n')

This routine computes `z = Qb`, `z = Qᵀb` or `z = Qᴴb` in place and overwrites b. It can only be executed once the factorization is done.

#### Input Arguments :

* `spfct`: the sparse factorization object resulting from the qrm_factorize! function.
* `b`: the vector(s) to which Q, Qᵀ or Qᴴ is applied.
* `transp`: whether to apply Q, Qᵀ or Qᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_apply! end

@doc raw"""
    z = qrm_apply(spfct, b; transp='n')
"""
function qrm_apply end

@doc raw"""
    qrm_spmat_mv!(spmat, alpha, x, beta, y; transp='n')

This subroutine performs a matrix-vector product of the type `y = αAx + βy`, `y = αAᵀx + βy` or `y = αAᴴx + βy`.

#### Input Arguments :

* `spmat`: the input matrix.
* `alpha`, `beta` : the α and β scalars
* `x`: the x vector(s).
* `y`: the y vector(s).
* `transp`: whether to multiply by A, Aᵀ or Aᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_spmat_mv! end

@doc raw"""
    qrm_spmat_nrm(spmat; ntype='f')

This routine computes the one-norm ``\|A\|_1``, the infinity-norm ``\|x\|_{\infty}`` or the two-norm ``\|x\|_2`` of a matrix.

#### Input Arguments :

* `spmat`: the input matrix.
* `ntype`: the type of norm to be computed. It can be either `'i'`, `'1'` or `'f'` for the infinity, one and Frobenius norms, respectively.
"""
function qrm_spmat_nrm end

@doc raw"""
    qrm_vecnrm!(x, nrm; ntype='2')

This routine computes the one-norm ``\|x\|_1``, the infinity-norm ``\|x\|_{\infty}`` or the two-norm ``\|x\|_2`` of a vector.

#### Input Arguments :

* `x`: the x vector(s).
* `nrm`: the computed norm(s). If x is a matrix this argument has to be a vector and each of its elements will contain the norm of the corresponding column of x.
* `ntype`: the type of norm to be computed. It can be either `'i'`, `'1'` or `'2'` for the infinity, one and two norms, respectively.
"""
function qrm_vecnrm! end

@doc raw"""
    nrm = qrm_vecnrm(x; ntype='2')
"""
function qrm_vecnrm end

@doc raw"""
    qrm_spbackslash!(spmat, b, x; transp='n')
"""
function qrm_spbackslash! end

@doc raw"""
    x = qrm_spbackslash(spmat, b; transp='n')
"""
function qrm_spbackslash end

@doc raw"""
    qrm_spposv!(spmat, b, x)

This function can be used to solve a linear symmetric, positive definite problem.
It is a shortcut for the sequence

    x = b
    qrm_analyse!(spmat, spfct; transp='n')
    qrm_factorize!(spmat, spfct; transp='n')
    qrm_solve!(spfct, x, x; transp='t')
    qrm_solve!(spfct, x, x; transp='t')

#### Input Arguments :

* `spmat`: the input matrix.
* `b`: the right-hand side(s).
* `x`: the solution vector(s).
"""
function qrm_spposv! end

@doc raw"""
    x = qrm_spposv(spmat, b)
"""
function qrm_spposv end

@doc raw"""
    qrm_least_squares!(spmat, b, x; transp='n')

This function can be used to solve a linear least squares problem

```math
\min \|Ax − b\|_2
```

in the case where the input matrix is square or overdetermined.
It is a shortcut for the sequence

    qrm_analyse!(spmat, spfct; transp='n')
    qrm_factorize!(spmat, spfct; transp='n')
    qrm_apply!(spfct, b; transp='t')
    qrm_solve!(spfct, b, x; transp='n')

#### Input Arguments :

* `spmat`: the input matrix.
* `b`: the ight-hand side(s).
* `x`: the solution vector(s).
* `transp`: whether to use A, Aᵀ or Aᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_least_squares! end

@doc raw"""
    x = qrm_least_squares(spmat, b; transp='n')
"""
function qrm_least_squares end

@doc raw"""
qrm_least_squares_semi_normal!(spmat, b, x, z, Δx, y; transp='n')

This function can be used to solve a linear least squares problem

```math
\min \|Ax − b\|_2
```

in the case where A is square or overdetermined.
Contrary to `qrm_least_squares!`, this function allows to solve the problem without storing the Q-factor of the QR factorization of A.

It is a shortcut for the sequence

    qrm_analyse!(spmat, spfct, transp  = 'n')
    qrm_factorize!(spmat, spfct, transp = 'n')
    qrm_spmat_mv!(spmat, T(1), b, T(0), z, transp = 't')
    qrm_solve!(spfct, z, y, transp = 't')
    qrm_solve!(spfct, y, x, transp = 'n')
    qrm_refine!(spmat, spfct, x, z, Δx, y)

Note that the Q-factor is not used in this sequence; only A and R. 

#### Input Arguments

* `spmat`: the input matrix.
* `spfct`: a sparse factorization object of type `qrm_spfct`.
* `b`: the right-hand side(s).
* `x`: the solution vector(s).
* `Δx`: an auxiliary vector (or matrix if b and x are matrices) used to compute the solution, the size of this vector (resp. matrix) is the same as x.
* `z`: an auxiliary vector (or matrix if b and x are matrices) used to store the value z = Aᵀb, the size of this vector (resp. matrix) is the same as x and can be unitialized when the function is called.
* `y`: an auxiliary vector (or matrix if b and x are matrices) used to compute the solution, the size of this vector (resp. matrix) is the same as b.
* `transp`: whether to use A, Aᵀ or Aᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_least_squares_semi_normal! end

@doc raw"""
    x = qrm_least_squares_semi_normal(spmat, b)
"""
function qrm_least_squares_semi_normal end

@doc raw"""
    qrm_min_norm!(spmat, b, x; transp='n')

This function can be used to solve a linear minimum norm problem

```math
\min \|x\|_2 \quad s.t. \quad Ax = b
```

in the case where the input matrix is square or underdetermined.
It is a shortcut for the sequence

    qrm_analyse!(spmat, spfct; transp='t')
    qrm_factorize!(spmat, spfct; transp='t')
    qrm_solve!(spfct, b, x; transp='t')
    qrm_apply!(spfct, x; transp='n')

#### Input Arguments :

* `spmat`: the input matrix.
* `b`: the right-hand side(s).
* `x`: the solution vector(s).
* `transp`: whether to use A, Aᵀ or Aᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_min_norm! end

"""
    x = qrm_min_norm(spmat, b; transp='n')
"""
function qrm_min_norm end

@doc raw"""
    qrm_min_norm_semi_normal!(spmat, spfct, b, x, Δx, y; transp='n')

This function can be used to solve a linear minimum-norm problem

```math
\min \|x\|_2 \quad s.t. \quad Ax = b
```
in the case where A is square or underdetermined.
Contrary to `qrm_min_norm!`, this function allows to solve the problem without storing the Q-factor in the QR factorization of Aᵀ.
It is a shortcut for the sequence

    qrm_analyse!(spmat, spfct, transp = 't')
    qrm_factorize!(spmat, spfct, transp = 't')
    qrm_solve!(spfct, b, Δx, transp = 't')
    qrm_solve!(spfct, Δx, y, transp = 'n')
    qrm_spmat_mv!(spmat, T(1),  y, T(0), x, transp = 't')

Remark that the Q-factor is not used in this sequence but rather A and R. 

#### Input Arguments :

* `spmat`: the input matrix.
* `spfct`: a sparse factorization object of type `qrm_spfct`.
* `b`: the right-hand side(s).
* `x`: the solution vector(s).
* `Δx`: an auxiliary vector (or matrix if x and b are matrices) used to compute the solution, the size of this vector (resp. matrix) is the same as x.
* `y`: an auxiliary vector (or matrix if x and b are matrices) used to compute the solution, the size of this vector (resp. matrix) is the same as b.
* `transp`: whether to use A, Aᵀ or Aᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_min_norm_semi_normal! end

@doc raw"""
    x = qrm_min_norm_semi_normal(spmat, b)
"""
function qrm_min_norm_semi_normal end

@doc raw"""
    qrm_residual_norm!(spmat, b, x, nrm; transp='n')

This function computes the scaled norm of the residual ``\frac{\|b - Ax\|_{\infty}}{\|b\|_{\infty} + \|x\|_{\infty} \|A\|_{\infty}}``, i.e., the normwise backward error.

#### Input Arguments :

* `spmat`: the input matrix.
* `b`: the right-hand side(s).
* `x`: the solution vector(s).
* `nrm`: the computed norm(s).
* `transp`: whether to use A, Aᵀ or Aᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_residual_norm! end

@doc raw"""
    nrm = qrm_residual_norm(spmat, b, x; transp='n')
"""
function qrm_residual_norm end

@doc raw"""
    qrm_residual_orth!(spmat, r, nrm; transp='n')

Computes the quantity ``\frac{\|A^T r\|_2}{\|r\|_2}`` which can be used to evaluate the quality of the solution of a least squares problem.

#### Input Arguments :

* `spmat`: the input matrix.
* `r`: the residual(s).
* `nrm`: the computed norm(s).
* `transp`: whether to use A, Aᵀ or Aᴴ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_residual_orth! end

@doc raw"""
    nrm = qrm_residual_orth(spmat, r; transp='n')
"""
function qrm_residual_orth end

@doc raw"""
    qrm_refine!(spmat, spfct, x, z, Δx, y)

Given an approximate solution x of the linear system RᵀRx ≈ z where R is the R-factor of some QR factorization of size (m, n), compute a refined solution.

### Input Arguments :

* `spmat`: the input matrix.
* `spfct`: a sparse factorization object of type `qrm_spfct`.
* `x`: the approximate solution vector(s), the size of this vector is n (or n×k if there are multiple solutions).
* `z`: the RHS vector(s) of the linear system, the size of this vector is n (or n×k if there are multiple RHS).
* `Δx`: an auxiliary vector (or matrix if x and z are matrices) used to compute the refinement, the size of this vector (resp. matrix) is n (resp. n×k).
* `y`: an auxiliary vector (or matrix if x and z are matrices) used to compute the refinement, the size of this vector (resp. matrix) is m (resp. m×k).
"""
function qrm_refine! end

@doc raw"""
    x_refined = qrm_refine(spmat, spfct, x, z)
"""
function qrm_refine end

@doc raw"""
    qrm_set(str, val)
    qrm_set(spfct, str, val)

Set control parameters that define the behavior of `qr_mumps`.

#### Input Arguments :

* `spfct`: a sparse factorization object of type `qrm_spfct`.
* `str`: a string describing the parameter to set.
* `val`: the parameter value.
"""
function qrm_set end

@doc raw"""
    val = qrm_get(str)
    val = qrm_get(spfct, str)

Returns the value of a control parameter or an information parameter.

#### Input Arguments :

* `spfct`: a sparse factorization object of type `qrm_spfct`.
* `str`: a string describing the parameter to get.
"""
function qrm_get end

@doc raw"""
    rp = qrm_spfct_get_rp(spfct)

Returns the row permutation.

#### Input Arguments :

* `spfct`: a sparse factorization object of type `qrm_spfct`.
"""
function qrm_spfct_get_rp end

@doc raw"""
    cp = qrm_spfct_get_cp(spfct)

Returns the column permutation.

#### Input Arguments :

* `spfct`: a sparse factorization object of type `qrm_spfct`.
"""
function qrm_spfct_get_cp end

@doc raw"""
    R = qrm_spfct_get_r(spfct)

Returns the R factor as a **SparseMatrixCSC** matrix.

#### Input Arguments :

* `spfct`: a sparse factorization object of type `qrm_spfct`.
"""
function qrm_spfct_get_r end

end # module
