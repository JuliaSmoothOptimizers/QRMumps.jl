module qr_mumps

using Libdl, SparseArrays, LinearAlgebra

import Base: \
    
if haskey(ENV, "JULIA_QR_MUMPS_LIBRARY_PATH")
  println("Custom Installation")
  const libsqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libsqrm.$dlext")
  const libdqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libdqrm.$dlext")
  const libcqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libcqrm.$dlext")
  const libzqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libzqrm.$dlext")
  const libqrm_common = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libqrm_common.$dlext")
else
  println("Yggdrasil Installation")
  using qr_mumps_jll
end

include("wrapper/qr_mumps_common.jl")
include("wrapper/qr_mumps_api.jl")

export qrm_spmat, qrm_spfct,
    qrm_init, qrm_finalize,
    qrm_spmat_init!, qrm_spmat_init, qrm_spmat_destroy!,
    qrm_spfct_init!, qrm_spfct_init, qrm_spfct_destroy!,
    qrm_analyse!, qrm_analyse,
    qrm_update!, qrm_factorize!,
    qrm_solve!, qrm_solve,
    qrm_apply!, qrm_apply,
    qrm_spmat_mv!, qrm_spmat_nrm,
    qrm_vecnrm!, qrm_vecnrm,
    qrm_spbackslash!, qrm_spbackslash, \,
    qrm_spposv!, qrm_spposv,
    qrm_least_squares!, qrm_least_squares,
    qrm_min_norm!, qrm_min_norm,
    qrm_residual_norm!, qrm_residual_norm,
    qrm_residual_orth!, qrm_residual_orth,
    qrm_set, qrm_get

@doc raw"""
    qrm_init(ncpu, ngpu)

This routine initializes qr\_mumps and should be called prior to any other qr\_mumps routine.

    qrm_init()

`ncpu` and `ngpu` are optional arguments and their default value are, respectively, `1` and `0`.

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

This routine initializes a **qrm_spmat** type data structure from a **sparseMatrixCSC**.

#### Input Arguments :

* `spmat`: the **qrm_spmat** sparse matrix to be initialized.
* `A` : a Julia sparse matrix stored in **SparseMatrixCSC** format.
* `sym` : a boolean to specify if the matrix is symmetric / hermitian (true) or unsymmetric (false).
"""
function qrm_spmat_init! end

@doc raw"""
    spmat = qrm_spmat_init(A; sym=false)
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

This routine updates a **qrm_spmat** type data structure from a **sparseMatrixCSC**.
`spmat` and `A` must have the same sparsity pattern.

#### Input Arguments :

* `spmat`: the **qrm_spmat** sparse matrix to be updated.
* `A` : a Julia sparse matrix stored in **SparseMatrixCSC** format.
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
* `transp`: whether to solve for R or Rᵀ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_solve! end

@doc raw"""
    x = qrm_solve(spfct, b; transp='n')
"""
function qrm_solve end

@doc raw"""
    qrm_apply!(spfct, b; transp='n')

This routine computes `z = Qb` or `z = Qᵀb` in place and overwrites b. It can only be executed once the factorization is done.

#### Input Arguments :

* `spfct`: the sparse factorization object resulting from the qrm_factorize! function.
* `b`: the vector(s) to which Q or Qᵀ is applied.
* `transp`: whether to apply Q or Qᵀ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_apply! end

@doc raw"""
    z = qrm_apply(spfct, b; transp='n')
"""
function qrm_apply end

@doc raw"""
    qrm_spmat_mv!(spmat, alpha, x, beta, y; transp='n')

This subroutine performs a matrix-vector product of the type y = αAx + βy or y = αAᵀx + βy.

#### Input Arguments :

* `spmat`: the input matrix.
* `alpha`, `beta` : the α and β scalars
* `x`: the x vector(s).
* `y`: the y vector(s).
* `transp`: whether to multiply by A or Aᵀ. Can be either `'t'`, `'c'` or `'n'`.
"""
function qrm_spmat_mv! end

@doc raw"""
    qrm_spmat_nrm(spmat, nrm; ntype='f')

This routine computes the one-norm ``\|A\|_1``, the infinity-norm ``\|x\|_{\infty}`` or the two-norm ``\|x\|_2`` of a matrix.

#### Input Arguments :

* `spmat`: the input matrix.
* `nrm`: the computed norm.
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
    qrm_spbackslash!(spmat, b, x)

TO DO !
"""
function qrm_spbackslash! end

@doc raw"""
    x = qrm_spbackslash(spmat, b)
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
    qrm_least_squares!(spmat, b, x)

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
"""
function qrm_least_squares! end

@doc raw"""
    x = qrm_least_squares(spmat, b)
"""
function qrm_least_squares end

@doc raw"""
    qrm_min_norm!(spmat, b, x)

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
"""
function qrm_min_norm! end

"""
    x = qrm_min_norm(spmat, b)
"""
function qrm_min_norm end

@doc raw"""
    qrm_residual_norm!(spmat, b, x, nrm)

This function computes the scaled norm of the residual ``\frac{\|b - Ax\|_{\infty}}{\|b\|_{\infty} + \|x\|_{\infty} \|A\|_{\infty}}``, i.e., the normwise backward error.

#### Input Arguments :

* `spmat`: the input matrix.
* `b`: the right-hand side(s).
* `x`: the solution vector(s).
* `nrm`: the computed norm(s).
"""
function qrm_residual_norm! end

@doc raw"""
    nrm = qrm_residual_norm(spmat, b, x)
"""
function qrm_residual_norm end

@doc raw"""
    qrm_residual_orth!(spmat, r, nrm)

Computes the quantity ``\frac{\|A^T r\|_2}{\|r\|_2}`` which can be used to evaluate the quality of the solution of a least squares problem.

#### Input Arguments :

* `spmat`: the input matrix.
* `r`: the residual(s).
* `nrm`: the computed norm(s).
"""
function qrm_residual_orth! end

@doc raw"""
    nrm = qrm_residual_orth(spmat, r)
"""
function qrm_residual_orth end

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

end # module
