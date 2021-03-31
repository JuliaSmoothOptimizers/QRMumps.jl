```@example spd1
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 1, 2, 3, 3, 4, 4, 5]
jcn = [1, 3, 4, 5, 2, 3, 5, 4, 5, 5]
val = [53.0, 8.0, 4.0, 3.0, 10.0, 6.0, 8.0, 26.0, 5.0, 14.0]

A = Symmetric(sparse(irn, jcn, val, 5, 5), :U)
b = [108.0, 20.0, 66.0, 133.0, 117.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

qrm_init()

spmat = qrm_spmat_init(A)
x = qrm_spposv(spmat, b)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example spd2
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5]
jcn = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5]
val = [4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0]

A = sparse(irn, jcn, val, 5, 5)
A_L = tril(A)
b = [5.0, 6.0, 6.0, 6.0, 5.0]
x_star = [1.0, 1.0, 1.0, 1.0, 1.0]

qrm_init()

spmat = qrm_spmat_init(A_L, sym=true)
x = qrm_spposv(spmat, b)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example spd3
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 3, 4, 5, 2, 3, 5, 4, 5, 5]
jcn = [1, 1, 1, 1, 2, 3, 3, 4, 4, 5]
val = [53.0, 8.0, 4.0, 3.0, 10.0, 6.0, 8.0, 26.0, 5.0, 14.0]

A = Symmetric(sparse(irn, jcn, val, 5, 5), :L)
b = [108.0, 20.0, 66.0, 133.0, 117.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)
z = qrm_solve(spfct, b, transp='t')
x = qrm_solve(spfct, z)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example spd4
using qr_mumps, LinearAlgebra, SparseArrays, Printf

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

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)
z = qrm_solve(spfct, b, transp='t')
x = qrm_solve(spfct, z)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example spd5
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 2, 2, 3, 3, 3]
jcn = [1, 1, 2, 1, 2, 3]
val = [7.0, -im, 8.0, 5im, 5.0, 10.0]

A = Hermitian(sparse(irn, jcn, val, 3, 3), :L)
b = [11.0-6im, 32.0+12im, 35.0+20im]
x_star = [1.0+im, 2.0+im, 3.0+im]
x = copy(b)

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)

qrm_solve!(spfct, x, x, transp='c')
qrm_solve!(spfct, x, x)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```

```@example spd6
using qr_mumps, LinearAlgebra, SparseArrays, Printf

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

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)

qrm_solve!(spfct, x, x, transp='c')
qrm_solve!(spfct, x, x)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)
```
