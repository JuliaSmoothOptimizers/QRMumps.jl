`QRMumps.jl` supports different column orderings to reduce fill-in in the factors.
By default, an ordering is chosen automatically by `qr_mumps`, but it is also possible to set the ordering manually.
This ordering must always be set before the analysis phase.

```@example ordering1
using QRMumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val, 7, 5)
b = [22.0, 5.0, 13.0, 8.0, 25.0, 5.0, 9.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# Use the natural ordering (no permutation is used)
qrm_set(spfct, "qrm_ordering", 1)

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)
z = qrm_apply(spfct, b, transp='t')
x = qrm_solve(spfct, z)

error_norm = norm(x - x_star)
r = A * x - b
optimality_residual_norm = norm(A' * r)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Optimality residual norm ‖Aᵀr‖ = %10.5e\n", optimality_residual_norm)
```

```@example ordering2
using QRMumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0+im, 2.0-im, 3.0+im, 1.0-im, 1.0+im, 2.0-im, 4.0+im, 1.0-im, 5.0+im, 1.0-im, 3.0+im, 6.0-im, 1.0+im]

A = sparse(irn, jcn, val, 7, 5)
b = [1.0+im, 2.0+im, 3.0+im, 4.0+im, 5.0+im, 6.0+im, 7.0+im]
z = copy(b)
x = zeros(ComplexF64, 5)

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# Provide your own column permutation
permutation = Cint[i for i = 5:-1:1]
qrm_user_permutation!(spfct, permutation)
qrm_set(spfct, "qrm_ordering", 2)

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)
qrm_apply!(spfct, z, transp='c')
qrm_solve!(spfct, z, x)

r = A * x - b
optimality_residual_norm = norm(A' * r)

@printf("Optimality residual norm ‖Aᵀr‖ = %10.5e\n", optimality_residual_norm)
```

```@example ordering3
using QRMumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5]
jcn = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5]
val = [4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0]

A = sparse(irn, jcn, val, 5, 5)
A_U = triu(A)
b = [5.0, 6.0, 6.0, 6.0, 5.0]
x_star = [1.0, 1.0, 1.0, 1.0, 1.0]

qrm_init()

spmat = qrm_spmat_init(A_U, sym=true)
spfct = qrm_spfct_init(spmat)

# Compute a column permutation with COLAMD
qrm_set(spfct, "qrm_ordering", 3)

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)
z = qrm_solve(spfct, b, transp='t')
x = qrm_solve(spfct, z)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example ordering4
using QRMumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 2, 3]
jcn = [1, 2, 3, 2, 3, 3]
val = [7.0, im, -5im, 8.0, 5.0, 10.0]

A = Hermitian(sparse(irn, jcn, val, 3, 3), :U)
b = [11.0-6im, 32.0+12im, 35.0+20im]
x_star = [1.0+im, 2.0+im, 3.0+im]
x = copy(b)

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# Compute a column permutation with METIS
qrm_set(spfct, "qrm_ordering", 4)

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)

qrm_solve!(spfct, x, x, transp='c')
qrm_solve!(spfct, x, x)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```
