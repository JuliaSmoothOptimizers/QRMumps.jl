```@example qless1
using QRMumps, LinearAlgebra, SparseArrays, Printf

m, n = 5, 7
irn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
jcn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val, m, n)
b = [40.0, 10.0, 44.0, 98.0, 87.0]
y_star = [0.0, 1.0, 2.0, 3.0, 4.0]
y₁ = zeros(n)
y₂ = zeros(m)

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

qrm_set(spfct, "qrm_keeph", 0)
qrm_analyse!(spmat, spfct; transp='t')
qrm_factorize!(spmat, spfct, transp='t')

qrm_solve!(spfct, b, y₁; transp='t')
qrm_solve!(spfct, y₁, y₂; transp='n')

error_norm = norm(y₂ - y_star)

@printf("Error norm ‖y* - y‖ = %10.5e\n", error_norm)
```

```@example qless2
using QRMumps, LinearAlgebra, SparseArrays, Printf

m, n = 7, 5
irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 3.0, 1.0, 3.0, 2.0, 1.0]

A = sparse(irn, jcn, val, m, n)
b = [22.0, 2.0, 13.0, 8.0, 17.0, 6.0, 5.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

qrm_set(spfct, "qrm_keeph", 0)

qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)

z = A'*b
x₁ = qrm_solve(spfct, z; transp = 't')
x₂ = qrm_solve(spfct, x₁; transp = 'n')

error_norm = norm(x₂ - x_star)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
```