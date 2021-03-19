qrm_init()
m = 200
n = 100
p = 5

for T in (Float32, Float64, ComplexF32, ComplexF64)
  tol = (real(T) == Float32) ? 1e-4 : 1e-12
  for I in (Int32 , Int64)
    A = sprand(m, n, 0.3)
    A = convert(SparseMatrixCSC{T,I}, A)
    b = rand(T, m)
    B = rand(T, m, p)

    spmat = qrm_spmat_init(A)
    spfct = qrm_analyse(spmat)
    qrm_factorize!(spmat, spfct)

    transp = (T <: Real) ? 't' : 'c'
    z = qrm_apply(spfct, b, transp=transp)
    x = qrm_solve(spfct, z, transp='n')
    r = b - A * x
    @test norm(A' * r) ≤ tol

    Z = qrm_apply(spfct, B, transp=transp)
    X = qrm_solve(spfct, Z, transp='n')
    R = B - A * X
    @test norm(A' * R) ≤ tol

    x = spmat \ b
    r = b - A * x
    @test norm(A' * r) ≤ tol

    X = spmat \ B
    R = B - A * X
    @test norm(A' * R) ≤ tol

    x = spfct \ b
    r = b - A * x
    @test norm(A' * r) ≤ tol

    X = spfct \ B
    R = B - A * X
    @test norm(A' * R) ≤ tol

    spmat = qrm_spmat_init(T)
    qrm_spmat_init!(spmat, A)

    spfct = qrm_spfct_init(spmat)
    qrm_analyse!(spmat, spfct)

    qrm_factorize!(spmat, spfct)

    transp = (T <: Real) ? 't' : 'c'
    z = copy(b)
    qrm_apply!(spfct, z, transp=transp)
    qrm_solve!(spfct, z, x, transp='n')
    r = b - A * x
    @test norm(A' * r) ≤ tol

    Z = copy(B)
    qrm_apply!(spfct, Z, transp=transp)
    qrm_solve!(spfct, Z, X, transp='n')
    R = B - A * X
    @test norm(A' * R) ≤ tol
  end
end

qrm_finalize()
