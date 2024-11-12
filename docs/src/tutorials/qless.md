```@example qless1
# The Q-less QR factorization may be used to solve the least-norm problem
#
#   minimize ‖x‖  subject to  Ax=b
#
# while saving storage because Q is not formed.
# Thus it is appropriate for large problems where storage is at a premium.
# The normal equations of the second kind AAᵀy = b are the optimality conditions of the least-norm problems, where x = Aᵀy.
# If Aᵀ = QR, they can be equivalently written RᵀRy = b.
#
# The stability of this procedure is comparable to the method that uses Q---see
#
# C. C. Paige, An error analysis of a method for solving matrix equations,
# Mathematics of Computations, 27, pp. 355-359, 1973, DOI 10.2307/2005623.

using LinearAlgebra, Printf, SparseArrays
using QRMumps

# Initialize data
m, n = 5, 7
irn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
jcn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val, m, n)
b = [40.0, 10.0, 44.0, 98.0, 87.0]
x_star = [16.0, 1.0, 10.0, 3.0, 19.0, 3.0, 2.0]
y₁ = zeros(n)
y = zeros(m)
x = zeros(n)

# Initialize QRMumps
qrm_init()

# Initialize data structures
spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# Specify that we want the Q-less QR factorization
qrm_set(spfct, "qrm_keeph", 0)

# Perform symbolic analysis of Aᵀ and factorize Aᵀ = QR
qrm_analyse!(spmat, spfct; transp='t')
qrm_factorize!(spmat, spfct, transp='t')

# Solve RᵀR y = b in two steps:
# 1. Solve Rᵀy₁ = b  
qrm_solve!(spfct, b, y₁; transp='t')

# 2. Solve Ry = y₁
qrm_solve!(spfct, y₁, y; transp='n')


# Compute the least norm solution of Ax = b
x .= A'*y

# Compute error norm and residual norm
error_norm = norm(x - x_star)
residual_norm = norm(b - A*x)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖b - Ax‖ = %10.5e\n", residual_norm)
```

```@example qless2
# The Q-less QR factorization may be used to solve the least-square problem
#
#   minimize ‖Ax - b‖
#
# while saving storage because Q is not formed.
# Thus it is appropriate for large problems where storage is at a premium.
# The normal equations AᵀAx = Aᵀb are the optimality conditions of the least squares problems.
# If A = QR, they can be equivalently written RᵀRx = Aᵀb.
#
# This procedure is backward stable if we perform one step of iterative refinement---see
#
# Å. Björck, Stability analysis of the method of seminormal equations for linear least squares problems,
# Linear Algebra and its Applications, 88–89, pp. 31-48, 1987, DOI 10.1016/0024-3795(87)90101-7.

using LinearAlgebra, Printf, SparseArrays 
using QRMumps

# Initialize data
m, n = 7, 5
irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 3.0, 1.0, 3.0, 2.0, 1.0]

A = sparse(irn, jcn, val, m, n)
b = [22.0, 2.0, 13.0, 8.0, 17.0, 6.0, 5.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

z = zeros(n)
x₁ = zeros(m)
x = zeros(n)

y = zeros(m)
r = zeros(n)
Δx₁ = zeros(m)
Δx = zeros(n)

# Initialize QRMumps
qrm_init()

# Initialize data structures
spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# Specify that we want the Q-less QR factorization
qrm_set(spfct, "qrm_keeph", 0)

# Perform symbolic analysis of A and factorize A = QR
qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)

# Compute the RHS of the semi-normal equations
mul!(z, A', b)

# Solve RᵀR x = z = Aᵀb in two steps:
# 1. Solve Rᵀx₁ = z  
qrm_solve!(spfct, z, x₁; transp = 't')

# 2. Solve Rx = x₁
qrm_solve!(spfct, x₁, x; transp = 'n')

error_norm = norm(x - x_star)
Aresidual_norm = norm(A'*(A*x - b))

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Normal equations residual norm ‖Aᵀ(Ax - b)‖= %10.5e\n", Aresidual_norm)

# As such, this method is not backward stable and we need to add an iterative refinement step:                                                          
# For this, we compute the least-squares solution Δx of min ‖r - AΔx‖, where r is the residual r = Aᵀb - AᵀA*x.
# We then update x := x + Δx

# Compute the residual in two steps to prevent allocating memory:
# 1. Compute y = b - Ax
mul!(y, A, x)
@. y = b - y

# 2. Compute r = Aᵀy = Aᵀ(b - A*x)
mul!(r, A', y)

# Solve the semi-normal equations as before
qrm_solve!(spfct, r, Δx₁; transp='t')
qrm_solve!(spfct, Δx₁, Δx; transp='n')

# Update the least squares solution
@. x = x + Δx

error_norm = norm(x - x_star)
Aresidual_norm = norm(A'*(A*x - b))

@printf("Error norm (iterative refinement step) ‖x* - x‖ = %10.5e\n", error_norm)
@printf("Normal equations residual norm (iterative refinement step) ‖Aᵀ(Ax - b)‖= %10.5e\n", Aresidual_norm)
```
