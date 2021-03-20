# Tutorials

## Symmetric and positive definite linear systems

### Basic

```@example spd1
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 1, 2, 3, 3, 4, 4, 5]
jcn = [1, 3, 4, 5, 2, 3, 5, 4, 5, 5]
val = [53.0, 8.0, 4.0, 3.0, 10.0, 6.0, 8.0, 26.0, 5.0, 14.0]

A = Symmetric(sparse(irn, jcn, val), :U)
b = [108.0, 20.0, 66.0, 133.0, 117.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

qrm_init()

spmat = qrm_spmat_init(A)
x = qrm_spposv(spmat, b)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)

qrm_finalize()
```

### Full

```@example spd2
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 1, 2, 3, 3, 4, 4, 5]
jcn = [1, 3, 4, 5, 2, 3, 5, 4, 5, 5]
val = [53.0, 8.0, 4.0, 3.0, 10.0, 6.0, 8.0, 26.0, 5.0, 14.0]

A = Symmetric(sparse(irn, jcn, val), :U)
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

qrm_finalize()
```

## Least-squares problems

### Basic

```@example ls1
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val)
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

qrm_finalize()
```

### Full

```@example ls2
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val)
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

qrm_finalize()
```

## Least-norm problems

### Basic

```@example ln1
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5]
jcn = [3, 5, 7, 1, 4, 6, 2, 6, 5, 6, 3, 4, 7]
val = [2.0, 3.0, 5.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0, 2.0]

A = sparse(irn, jcn, val)
b = [56.0, 21.0, 16.0, 22.0, 25.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

qrm_init()

spmat = qrm_spmat_init(A)
x = qrm_min_norm(spmat, b)

error_norm = norm(x - x_star)
residual_norm = norm(A * x - b)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖Ax - b‖ = %10.5e\n", residual_norm)

qrm_finalize()
```

### Full

```@example ln2
using qr_mumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5]
jcn = [3, 5, 7, 1, 4, 6, 2, 6, 5, 6, 3, 4, 7]
val = [2.0, 3.0, 5.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0, 2.0]

A = sparse(irn, jcn, val)
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

qrm_finalize()
```
