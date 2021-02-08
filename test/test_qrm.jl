irn = Int32.([1, 1, 1, 1, 2, 3, 3, 4, 4, 5])
jcn = Int32.([1, 3, 4, 5, 2, 3, 5, 4, 5, 5])
val = [53.0, 8.0, 4.0, 3.0, 10.0, 6.0, 8.0, 26.0, 5.0, 14.0]
n   = Int32(5)
m   = Int32(5)
nz  = Int32(10)
sym = Int32(1)

for T in (Float32, Float64, ComplexF32, ComplexF64)
  qrm_init(1, 0)

  spmat = qrm_spmat{T}()
  qrm_spmat_init(spmat)
  @test spmat.h ≠ C_NULL
  spmat.irn = pointer(irn)
  spmat.jcn = pointer(jcn)
  spmat.val = pointer(T.(val))
  spmat.n   = n
  spmat.m   = m
  spmat.nz  = nz
  spmat.sym = sym

  spfct = qrm_spfct{T}()
  qrm_spfct_init(spfct, spmat)

  # A = sparse(Int64.(irn), Int64.(jcn), val, m, n)

  # for transp ∈ ('n', 't', 'c')
  #   qrm_analyse(spmat, spfct, transp=transp)
  #   qrm_factorize(spmat, spfct, transp=transp)

  #   b1 = rand(T, 5)
  #   x1 = zeros(T, 5)
  #   qrm_solve(spfct, b1, x1, transp=transp)
  #   qrm_apply(spfct, b1, transp=transp)
  #   show(x1)
  #   println()
  #   println(norm(b1 - A * x1))

  #   b2 = rand(T, 5, 3)
  #   x2 = zeros(T, 5, 3)
  #   qrm_solve(spfct, b2, x2, transp=transp)
  #   qrm_apply(spfct, b2, transp=transp)
  #   show(x2)
  #   println()
  #   println(norm(b2 - A * x2))
  # end

  # for ntype ∈ ('i', '1', 'f')
  #   nrm = one(real(T(3)))
  #   qrm_spmat_nrm(spmat, ntype, nrm)
  # end

  qrm_finalize()

end
