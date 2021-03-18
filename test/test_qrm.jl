qrm_init()
m = 200
n = 100
p = 5

for T in (Float32, Float64)
  tol = (T == Float32) ? 1e-4 : 1e-12
  for I in (Int32 , Int64)
    A = sprand(m, n, 0.3)
    A = convert(SparseMatrixCSC{T,I}, A)
    b = rand(T, m)
    B = rand(T, m, p)

    spmat = qrm_spmat_init(A)
    spfct = qrm_analyse(spmat)
    qrm_factorize!(spmat, spfct)

    z = qrm_apply(spfct, b, transp='t')
    x = qrm_solve(spfct, z, transp='n')
    r = b - A * x
    @test norm(A' * r) ≤ tol

    Z = qrm_apply(spfct, B, transp='t')
    X = qrm_solve(spfct, Z, transp='n')
    R = B - A * X
    @test norm(A' * R) ≤ tol
  end
end

for T in (ComplexF32, ComplexF64)
  tol = (T == ComplexF32) ? 1e-4 : 1e-12
  for I in (Int32 , Int64)
    A = sprand(m, n, 0.3)
    A = convert(SparseMatrixCSC{T,I}, A)
    b = rand(T, m)
    B = rand(T, m, p)

    spmat = qrm_spmat_init(A)
    spfct = qrm_analyse(spmat)
    qrm_factorize!(spmat, spfct)

    z = qrm_apply(spfct, b, transp='c')
    x = qrm_solve(spfct, z, transp='n')
    r = b - A * x
    @test norm(A' * r) ≤ tol

    Z = qrm_apply(spfct, B, transp='c')
    X = qrm_solve(spfct, Z, transp='n')
    R = B - A * X
    @test norm(A' * R) ≤ tol
  end
end

qrm_finalize()
