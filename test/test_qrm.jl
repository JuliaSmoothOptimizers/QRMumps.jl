qrm_init(1, 0)
n = 100

for T in (Float32, Float64, ComplexF32, ComplexF64)
  for I in (Int32 , Int64)
    A = sprand(n, n, 0.1)
    A = convert(SparseMatrixCSC{T,I}, A)
    b = rand(T, n)
    x = rand(T, n)
    B = rand(T, n, n)
    X = rand(T, n, n)
    spmat = qrm_spmat_init(T)
    qrm_spmat_init!(spmat)
    spmat = qrm_spmat_init(A)
    qrm_spmat_init!(spmat, A)
    qrm_spmat_init(Symmetric(A))
    qrm_spmat_init(Hermitian(A))
    qrm_spmat_init!(spmat, Symmetric(A))
    qrm_spmat_init!(spmat, Hermitian(A))

    spfct = qrm_spfct_init(spmat)
    qrm_spfct_init!(spfct, spmat)

    # for transp ∈ ('n','t','c')
    #   qrm_analyse!(spmat, spfct, transp=transp)
    #   qrm_analyse(spmat, transp=transp)

    #   qrm_factorize!(spmat, spfct, transp=transp)

    #   qrm_solve(spfct, b, x, transp=transp)
    #   qrm_solve(spfct, B, X, transp=transp)
    #   qrm_solve!(spfct, b, x, transp=transp)
    #   qrm_solve!(spfct, B, X, transp=transp)
    # end
  end
end

qrm_finalize()
