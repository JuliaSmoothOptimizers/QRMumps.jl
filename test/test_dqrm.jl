m = 5
n = 3
A = sprand(m, n, .7)
b = rand(m)

# To avoid error in running finalizer: ErrorException("stream is closed or unusable")
blas_set_num_threads(1)

qrm_mat = QrmType(A)
transp = m < n ? 't' : 'n'
qrm_analyze(qrm_mat, transp=transp)
println(qrm_mat.gstats)

qrm_factorize(qrm_mat, transp=transp)
if m ≥ n
  # Solve the least-squares problem min ‖Rx - Q'b‖.
  qrm_apply(qrm_mat, b, transp='t')
  x = qrm_solve(qrm_mat, b, transp='n')
else
  # Solve the least-norm problem min ‖x‖ s.t. R'Q'x=b.
  x = qrm_solve(qrm_mat, b, transp='t')
  qrm_apply(qrm_mat, x, transp='n')
end

# x = qrm_least_squares(A, b)

println(x)
r = b - A * x
@printf("KKT residual: %7.1e\n", norm(A' * r))

x_ex = full(A) \ b  # ???
println(x_ex)

r_ex = b - A * x_ex
@printf("KKT residual (backslash): %7.1e\n", norm(A' * r_ex))
@printf("Relative error: %7.1e\n", norm(x - x_ex) / norm(x_ex))
