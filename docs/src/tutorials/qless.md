```@example qless1
using QRMumps, LinearAlgebra, SparseArrays, Printf

# Initialize data
m, n = 5, 7
irn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
jcn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 5.0, 1.0, 3.0, 6.0, 1.0]

A = sparse(irn, jcn, val, m, n)
b = [40.0, 10.0, 44.0, 98.0, 87.0]
y_star = [16.0, 1.0, 10.0, 3.0, 19.0, 3.0, 2.0]
y₁ = zeros(n)
y₂ = zeros(m)
y₃ = zeros(n)

Δy₁ = zeros(n)
Δy₂ = zeros(m)
Δy₃ = zeros(n)
r = zeros(m)

# Initialize QRMumps
qrm_init()

# Initialize data structures
spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# Specify not storing Q for the QR factorization
qrm_set(spfct, "qrm_keeph", 0)

# Perform symbolic analysis of Aᵀ and factorize Aᵀ = QR
qrm_analyse!(spmat, spfct; transp='t')
qrm_factorize!(spmat, spfct, transp='t')

# Solve Rᵀy₁ =  b  
qrm_solve!(spfct, b, y₁; transp='t')

# Solve Ry₂ = y₁
qrm_solve!(spfct, y₁, y₂; transp='n')

# Overall, RᵀRy₂ = b. Equivalently, RᵀQᵀQRy₂ = b or AAᵀy₂ = b

# Compute least norm solution of min ‖b - Ay‖
y₃ .= A'*y₂

# Compute error norm and residual norm
error_norm = norm(y₃ - y_star)
residual_norm = norm(b - A*y₃)

@printf("Error norm ‖y* - y‖ = %10.5e\n", error_norm)
@printf("Residual norm ‖b - Ay‖ = %10.5e\n", residual_norm)

# We can improve this solution with iterative refinement: see Björck 1967 - Iterative refinement of linear least squares solutions I
#                                                             in BIT Numerical Mathematics.
# For this, we compute the least norm solution Δy₃ of min ‖r - AΔy‖, where r is the residual r = b - A*y₃.
# We then update y₃ := y₃ + Δy₃
r .= b - A*y₃

# Solve RᵀΔy₁ =  r 
qrm_solve!(spfct, r, Δy₁; transp='t')

# Solve RΔy₂ = Δy₁
qrm_solve!(spfct, Δy₁, Δy₂; transp='n')

# Overall, RᵀRΔy₂ = r. Equivalently, RᵀQᵀQRΔy₂ = r or AAᵀΔy₂ = r

# Compute least norm solution of min ‖r - Ay‖
Δy₃ .= A'*Δy₂

# Update the least norm solution
y₃ .= y₃ .+ Δy₃

error_norm = norm(y₃ - y_star)
residual_norm = norm(b - A*y₃)

@printf("Error norm - iterative refinement ‖y* - y‖ = %10.5e\n", error_norm)
@printf("Residual norm - iterative refinement ‖b - Ay‖ = %10.5e\n", residual_norm)
```

```@example qless2
using QRMumps, LinearAlgebra, SparseArrays, Printf

# Initialize data
m, n = 7, 5
irn = [1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7]
jcn = [1, 3, 5, 2, 3, 5, 1, 4, 4, 5, 2, 1, 3]
val = [1.0, 2.0, 3.0, 1.0, 1.0, 2.0, 4.0, 1.0, 3.0, 1.0, 3.0, 2.0, 1.0]

A = sparse(irn, jcn, val, m, n)
b = [22.0, 2.0, 13.0, 8.0, 17.0, 6.0, 5.0]
x_star = [1.0, 2.0, 3.0, 4.0, 5.0]

# Initialize QRMumps
qrm_init()

# Initialize data structures
spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# Specify not storing Q for the QR factorization
qrm_set(spfct, "qrm_keeph", 0)

# Perform symbolic analysis of A and factorize A = QR
qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)

z = A'*b

# Solve Rᵀx₁ = z = Aᵀb
x₁ = qrm_solve(spfct, z; transp = 't')

# Solve Rx₂ = x₁ 
x₂ = qrm_solve(spfct, x₁; transp = 'n')

# Overall, RᵀRx₂ = Aᵀb. Equivalently, RᵀQᵀQRx₂ = Aᵀb or AᵀAx₂ = Aᵀb
error_norm = norm(x₂ - x_star)

@printf("Error norm ‖x* - x‖ = %10.5e\n", error_norm)
```