```@example ls1
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val, 7, 5)
b = [22.0, 5.0, 13.0, 8.0, 25.0, 5.0, 9.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

qrm_init()

spmat = qrm_spmat_init(A)
x = qrm_least_squares(spmat, b)

error_norm = norm(x - x_star)
r = A * x - b
optimality_residual_norm = norm(A' * r)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Optimality residual norm ‖Aᵀr‖ = %10.5e\n", optimality_residual_norm)
```

```@example ls2
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val, 7, 5)
b = [22.0, 5.0, 13.0, 8.0, 25.0, 5.0, 9.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

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

```@example ls3
using qr_mumps, LinearAlgebra, SparseArrays, Printf

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

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)
qrm_apply!(spfct, z, transp='c')
qrm_solve!(spfct, z, x)

r = A * x - b
optimality_residual_norm = norm(A' * r)

@printf("Optimality residual norm ‖Aᵀr‖ = %10.5e\n", optimality_residual_norm)
```
