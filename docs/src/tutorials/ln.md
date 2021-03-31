```@example ln1
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5]
jcn = [3, 5, 7, 1, 4, 6, 2, 6, 5, 6, 3, 4, 7]
val = [2.0, 3.0, 5.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0, 2.0]

A = sparse(irn, jcn, val, 5, 7)
b = [56.0, 21.0, 16.0, 22.0, 25.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

qrm_init()

spmat = qrm_spmat_init(A)
x = qrm_min_norm(spmat, b)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example ln2
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5]
jcn = [3, 5, 7, 1, 4, 6, 2, 6, 5, 6, 3, 4, 7]
val = [2.0, 3.0, 5.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0, 2.0]

A = sparse(irn, jcn, val, 5, 7)
b = [56.0, 21.0, 16.0, 22.0, 25.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

qrm_analyse!(spmat, spfct, transp='t')
qrm_factorize!(spmat, spfct, transp='t')
z = qrm_solve(spfct, b, transp='t')
x = qrm_apply(spfct, z)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example ln3
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5]
jcn = [3, 5, 7, 1, 4, 6, 2, 6, 5, 6, 3, 4, 7]
val = [1.0-im, 2.0+im, 3.0-im, 1.0+im, 1.0-im, 2.0+im, 4.0-im, 1.0+im, 5.0-im, 1.0+im, 3.0-im, 6.0+im, 1.0-im]

A = sparse(irn, jcn, val, 5, 7)
b = [1.0+im, 2.0+im, 3.0+im, 4.0+im, 5.0+im]
x = zeros(ComplexF64, 7)

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

qrm_analyse!(spmat, spfct, transp='c')
qrm_factorize!(spmat, spfct, transp='c')
qrm_solve!(spfct, b, x, transp='c')
qrm_apply!(spfct, x)

residual_norm = norm(A * x - b)

@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```
