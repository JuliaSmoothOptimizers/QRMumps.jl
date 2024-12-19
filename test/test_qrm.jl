Random.seed!(1234)
m = 200
n = 100
p = 5
d = Dict(0 => "auto", 1 => "natural", 2 => "given", 3 => "colamd", 4 => "metis", 5 => "scotch")

@testset "least-squares problems" begin
  @testset "$T" for T in (Float32, Float64, ComplexF32, ComplexF64)

    tol = (real(T) == Float32) ? 1e-3 : 1e-12
    transp = (T <: Real) ? 't' : 'c'

    @testset "$INT" for INT in (Int32, Int64)
      A = sprand(T, m, n, 0.3)
      A = convert(SparseMatrixCSC{T,I}, A)

      A_transp = sprand(T, n, m, 0.3)
      A_transp = convert(SparseMatrixCSC{T,I}, A_transp)

      b = rand(T, m)
      B = rand(T, m, p)

      spmat = qrm_spmat_init(A)

      spfct = qrm_analyse(spmat)
      qrm_factorize!(spmat, spfct)
      spfct2 = (T <: Real) ? Transpose(spfct) : Adjoint(spfct)

      z = qrm_apply(spfct, b, transp=transp)
      x = qrm_solve(spfct, z, transp='n')
      r = b - A * x
      @test norm(A' * r) ≤ tol

      Z = qrm_apply(spfct, B, transp=transp)
      X = qrm_solve(spfct, Z, transp='n')
      R = B - A * X
      @test norm(A' * R) ≤ tol

      z = qrm_apply(spfct2, b)
      x = qrm_solve(spfct, z, transp='n')
      r = b - A * x
      @test norm(A' * r) ≤ tol

      Z = qrm_apply(spfct2, B)
      X = qrm_solve(spfct, Z, transp='n')
      R = B - A * X
      @test norm(A' * R) ≤ tol

      @testset "$(d[ordering])" for ordering in (0, 1, 2, 3, 4, 5)
        (ordering == 5) && continue  # SCOTCH v7.0.6 is not working with qr_mumps
        spmat = qrm_spmat_init(T)
        qrm_spmat_init!(spmat, A)

      spmat_transp = qrm_spmat_init(T)
      qrm_spmat_init!(spmat_transp, A_transp)

      spfct = qrm_spfct_init(spmat)
      qrm_analyse!(spmat, spfct)
      qrm_factorize!(spmat, spfct)
      spfct2 = (T <: Real) ? Transpose(spfct) : Adjoint(spfct)

        qrm_analyse!(spmat, spfct)
        qrm_factorize!(spmat, spfct)
        spfct2 = (T <: Real) ? Transpose(spfct) : Adjoint(spfct)

        if QRMumps.QRMUMPS_INSTALLATION == "YGGDRASIL"
          rp = qrm_spfct_get_rp(spfct)
          cp = qrm_spfct_get_cp(spfct)
          R = qrm_spfct_get_r(spfct)
          @test norm(A[rp,cp]' * A[rp,cp] - R[rp,cp]' * R[rp,cp]) ≤ tol
        end

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

        z = copy(b)
        qrm_apply!(spfct2, z)
        qrm_solve!(spfct, z, x, transp='n')
        r = b - A * x
        @test norm(A' * r) ≤ tol

        Z = copy(B)
        qrm_apply!(spfct2, Z)
        qrm_solve!(spfct, Z, X, transp='n')
        R = B - A * X
        @test norm(A' * R) ≤ tol
      end

      x = qrm_least_squares(spmat, b)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      x = qrm_least_squares_semi_normal(spmat, b)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      x = qrm_least_squares_semi_normal(spmat_transp, b, transp = transp)
      r = b - A_transp' * x
      @test norm(A_transp * r) ≤ tol

      X = qrm_least_squares(spmat, B)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      X = qrm_least_squares_semi_normal(spmat, B)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      X = qrm_least_squares_semi_normal(spmat_transp, B, transp = transp)
      R = B - A_transp' * X
      @test norm(A_transp * r) ≤ tol

      bc = copy(b)
      qrm_least_squares!(spmat, bc, x)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      z = similar(x)
      Δx = similar(x)
      y = similar(b)
      x = qrm_least_squares_semi_normal!(spmat, spfct, b, x, z, Δx, y)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      Bc = copy(B)
      qrm_least_squares!(spmat, Bc, X)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      Z = similar(X)
      ΔX = similar(X)
      Y = similar(B)
      qrm_least_squares_semi_normal!(spmat, spfct, B, X, Z, ΔX, Y)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      x = qrm_spbackslash(spmat, b)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      X = qrm_spbackslash(spmat, B)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      bc = copy(b)
      qrm_spbackslash!(spmat, bc, x)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      Bc = copy(B)
      qrm_spbackslash!(spmat, Bc, X)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      x = spmat \ b
      r = b - A * x
      @test norm(A' * r) ≤ tol

      X = spmat \ B
      R = B - A * X
      @test norm(A' * R) ≤ tol

      x = qrm_spbackslash(spfct, b)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      X = qrm_spbackslash(spfct, B)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      bc = copy(b)
      qrm_spbackslash!(spfct, bc, x)
      r = b - A * x
      @test norm(A' * r) ≤ tol

      Bc = copy(B)
      qrm_spbackslash!(spfct, Bc, X)
      R = B - A * X
      @test norm(A' * R) ≤ tol

      x = spfct \ b
      r = b - A * x
      @test norm(A' * r) ≤ tol

      X = spfct \ B
      R = B - A * X
      @test norm(A' * R) ≤ tol
    end
  end
end

@testset "least-norm problems" begin
  @testset "$T" for T in (Float32, Float64, ComplexF32, ComplexF64)

    tol = (real(T) == Float32) ? 1e-3 : 1e-12
    transp = (T <: Real) ? 't' : 'c'

    @testset "$INT" for INT in (Int32, Int64)
      A = sprand(T, n, m, 0.3)
      A = convert(SparseMatrixCSC{T,I}, A)

      A_transp = sprand(T, m, n, 0.3)
      A_transp = convert(SparseMatrixCSC{T,I}, A_transp)

      b = rand(T, n)
      B = rand(T, n, p)

      spmat = qrm_spmat_init(A)
      spfct = qrm_analyse(spmat, transp=transp)

      spmat_transp = qrm_spmat_init(T)
      qrm_spmat_init!(spmat_transp, A_transp)

      qrm_factorize!(spmat, spfct, transp=transp)
      spfct2 = (T <: Real) ? Transpose(spfct) : Adjoint(spfct)

      z = qrm_solve(spfct, b, transp=transp)
      x = qrm_apply(spfct, z, transp='n')
      r = b - A * x
      @test norm(r) ≤ tol

      Z = qrm_solve(spfct, B, transp=transp)
      X = qrm_apply(spfct, Z, transp='n')
      R = B - A * X
      @test norm(R) ≤ tol

      z = qrm_solve(spfct2, b)
      x = qrm_apply(spfct, z, transp='n')
      r = b - A * x
      @test norm(r) ≤ tol

      Z = qrm_solve(spfct2, B)
      X = qrm_apply(spfct, Z, transp='n')
      R = B - A * X
      @test norm(R) ≤ tol

      @testset "$(d[ordering])" for ordering in (0, 1, 2, 3, 4, 5)
        (ordering == 5) && continue  # SCOTCH v7.0.6 is not working with qr_mumps
        spmat = qrm_spmat_init(T)
        qrm_spmat_init!(spmat, A)

        spfct = qrm_spfct_init(spmat)
        qrm_set(spfct, "qrm_ordering", ordering)
        if ordering == 2
          permutation = Cint[i for i = n:-1:1]
          qrm_user_permutation!(spfct, permutation)
        end

        qrm_analyse!(spmat, spfct, transp=transp)
        qrm_factorize!(spmat, spfct, transp=transp)
        spfct2 = (T <: Real) ? Transpose(spfct) : Adjoint(spfct)

        if QRMumps.QRMUMPS_INSTALLATION == "YGGDRASIL"
          rp = qrm_spfct_get_rp(spfct)
          cp = qrm_spfct_get_cp(spfct)
          R = qrm_spfct_get_r(spfct)
          @test norm(A[cp,rp] * A[cp,rp]' - R[rp,cp]' * R[rp,cp]) ≤ tol
        end

        qrm_solve!(spfct, b, x, transp=transp)
        qrm_apply!(spfct, x, transp='n')
        r = b - A * x
        @test norm(r) ≤ tol

        qrm_solve!(spfct, B, X, transp=transp)
        qrm_apply!(spfct, X, transp='n')
        R = B - A * X
        @test norm(R) ≤ tol

        qrm_solve!(spfct2, b, x)
        qrm_apply!(spfct, x, transp='n')
        r = b - A * x
        @test norm(r) ≤ tol

        qrm_solve!(spfct2, B, X)
        qrm_apply!(spfct, X, transp='n')
        R = B - A * X
        @test norm(R) ≤ tol
      end

      x = qrm_min_norm(spmat, b)
      r = b - A * x
      @test norm(r) ≤ tol

      x = qrm_min_norm_semi_normal(spmat, b)
      r = b - A * x
      @test norm(r) ≤ tol

      x = qrm_min_norm_semi_normal(spmat_transp, b, transp = transp)
      r = b - A_transp' * x
      @test norm(r) ≤ tol

      X = qrm_min_norm(spmat, B)
      R = B - A * X
      @test norm(R) ≤ tol
      
      X = qrm_min_norm_semi_normal(spmat, B)
      R = B - A * X
      @test norm(R) ≤ tol

      X = qrm_min_norm_semi_normal(spmat_transp, B, transp = transp)
      R = B - A_transp' * X
      @test norm(R) ≤ tol

      qrm_min_norm!(spmat, b, x)
      r = b - A * x
      @test norm(r) ≤ tol
      
      Δx = similar(x)
      y = similar(b)
      qrm_min_norm_semi_normal!(spmat, spfct, b, x, Δx, y)
      r = b - A * x
      @test norm(r) ≤ tol

      qrm_min_norm!(spmat, B, X)
      R = B - A * X
      @test norm(R) ≤ tol

      ΔX = similar(X)
      Y = similar(B)
      qrm_min_norm_semi_normal!(spmat, spfct, B, X, ΔX, Y)
      R = B - A * X
      @test norm(R) ≤ tol

      x = qrm_spbackslash(spmat, b)
      r = b - A * x
      @test norm(r) ≤ tol

      X = qrm_spbackslash(spmat, B)
      R = B - A * X
      @test norm(R) ≤ tol

      qrm_spbackslash!(spmat, b, x)
      r = b - A * x
      @test norm(r) ≤ tol

      qrm_spbackslash!(spmat, B, X)
      R = B - A * X
      @test norm(R) ≤ tol

      x = spmat \ b
      r = b - A * x
      @test norm(r) ≤ tol

      X = spmat \ B
      R = B - A * X
      @test norm(R) ≤ tol

      x = qrm_spbackslash(spfct, b, transp=transp)
      r = b - A * x
      @test norm(r) ≤ tol

      X = qrm_spbackslash(spfct, B, transp=transp)
      R = B - A * X
      @test norm(R) ≤ tol

      qrm_spbackslash!(spfct, b, x, transp=transp)
      r = b - A * x
      @test norm(r) ≤ tol

      qrm_spbackslash!(spfct, B, X, transp=transp)
      R = B - A * X
      @test norm(R) ≤ tol

      x = spfct2 \ b
      r = b - A * x
      @test norm(r) ≤ tol

      X = spfct2 \ B
      R = B - A * X
      @test norm(R) ≤ tol
    end
  end
end

@testset "Symmetric and positive definite linear systems" begin
  @testset "$T" for T in (Float32, Float64, ComplexF32, ComplexF64)

    tol = (real(T) == Float32) ? 1e-3 : 1e-12
    transp = (T <: Real) ? 't' : 'c'

    @testset "$INT" for INT in (Int32, Int64)
      A = sprand(T, n, n, 0.01)
      A = convert(SparseMatrixCSC{T,INT}, A)
      A = A * A' + I
      A = (T <: Real) ? Symmetric(tril(A), :L) : Hermitian(tril(A), :L)
      b = rand(T, n)
      B = rand(T, n, p)

      spmat = qrm_spmat_init(A)
      spfct = qrm_analyse(spmat)
      qrm_factorize!(spmat, spfct)
      spfct2 = (T <: Real) ? Transpose(spfct) : Adjoint(spfct)

      z = qrm_solve(spfct, b, transp=transp)
      x = qrm_solve(spfct, z)
      r = b - A * x
      @test norm(r) ≤ tol

      Z = qrm_solve(spfct, B, transp=transp)
      X = qrm_solve(spfct, Z)
      R = B - A * X
      @test norm(R) ≤ tol

      z = qrm_solve(spfct2, b)
      x = qrm_solve(spfct, z)
      r = b - A * x
      @test norm(r) ≤ tol

      Z = qrm_solve(spfct2, B)
      X = qrm_solve(spfct, Z)
      R = B - A * X
      @test norm(R) ≤ tol

      A = sprand(T, n, n, 0.01)
      A = convert(SparseMatrixCSC{T,INT}, A)
      A = A * A' + I
      A = (T <: Real) ? Symmetric(triu(A), :U) : Hermitian(triu(A), :U)
      x = zeros(T, n)
      X = zeros(T, n, p)
      b = rand(T, n)
      B = rand(T, n, p)

      @testset "$(d[ordering])" for ordering in (0, 1, 2, 3, 4, 5)
        (ordering == 5) && continue  # SCOTCH v7.0.6 is not working with qr_mumps
        spmat = qrm_spmat_init(T)
        qrm_spmat_init!(spmat, A)

        spfct = qrm_spfct_init(spmat)
        qrm_set(spfct, "qrm_ordering", ordering)
        if ordering == 2
          permutation = Cint[i for i = n:-1:1]
          qrm_user_permutation!(spfct, permutation)
        end

        qrm_analyse!(spmat, spfct)
        qrm_factorize!(spmat, spfct)
        spfct2 = (T <: Real) ? Transpose(spfct) : Adjoint(spfct)

        if QRMumps.QRMUMPS_INSTALLATION == "YGGDRASIL"
          rp = qrm_spfct_get_rp(spfct)
          cp = qrm_spfct_get_cp(spfct)
          R = qrm_spfct_get_r(spfct)
          @test norm(A[rp,cp] - R[rp,cp]' * R[rp,cp]) ≤ tol
        end

        qrm_solve!(spfct, copy(b), x, transp=transp)
        qrm_solve!(spfct, x, x, transp='n')
        r = b - A * x
        @test norm(r) ≤ tol

        qrm_solve!(spfct, copy(B), X, transp=transp)
        qrm_solve!(spfct, X, X, transp='n')
        R = B - A * X
        @test norm(R) ≤ tol

        qrm_solve!(spfct2, copy(b), x)
        qrm_solve!(spfct, x, x, transp='n')
        r = b - A * x
        @test norm(r) ≤ tol

        qrm_solve!(spfct2, copy(B), X)
        qrm_solve!(spfct, X, X, transp='n')
        R = B - A * X
        @test norm(R) ≤ tol
      end

      x = qrm_spposv(spmat, b)
      r = b - A * x
      @test norm(r) ≤ tol

      X = qrm_spposv(spmat, B)
      R = B - A * X
      @test norm(R) ≤ tol

      bc = copy(b)
      qrm_spposv!(spmat, bc, x)
      r = b - A * x
      @test norm(r) ≤ tol

      Bc = copy(B)
      qrm_spposv!(spmat, Bc, X)
      R = B - A * X
      @test norm(R) ≤ tol

      x = qrm_spbackslash(spmat, b)
      r = b - A * x
      @test norm(r) ≤ tol

      X = qrm_spbackslash(spmat, B)
      R = B - A * X
      @test norm(R) ≤ tol

      bc = copy(b)
      qrm_spbackslash!(spmat, bc, x)
      r = b - A * x
      @test norm(r) ≤ tol

      Bc = copy(B)
      qrm_spbackslash!(spmat, Bc, X)
      R = B - A * X
      @test norm(R) ≤ tol

      x = spmat \ b
      r = b - A * x
      @test norm(r) ≤ tol

      X = spmat \ B
      R = B - A * X
      @test norm(R) ≤ tol

      x = qrm_spbackslash(spfct, b)
      r = b - A * x
      @test norm(r) ≤ tol

      X = qrm_spbackslash(spfct, B)
      R = B - A * X
      @test norm(R) ≤ tol

      bc = copy(b)
      qrm_spbackslash!(spfct, bc, x)
      r = b - A * x
      @test norm(r) ≤ tol

      Bc = copy(B)
      qrm_spbackslash!(spfct, Bc, X)
      R = B - A * X
      @test norm(R) ≤ tol

      x = spfct \ b
      r = b - A * x
      @test norm(r) ≤ tol

      X = spfct \ B
      R = B - A * X
      @test norm(R) ≤ tol
    end
  end
end

@testset "Auxiliary functions" begin
  for str ∈ QRMumps.GICNTL ∪ QRMumps.PICNTL ∪ QRMumps.RCNTL
    qrm_get(str)
    qrm_set(str, 1)
  end

  @testset "precision = $T" for T in (Float32, Float64, ComplexF32, ComplexF64)
    transp = (T <: Real) ? 't' : 'c'
    A = sprand(T, n, n, 0.3)
    spmat = qrm_spmat_init(A)
    spfct = qrm_analyse(spmat)
    qrm_factorize!(spmat, spfct)

    for str ∈ QRMumps.PICNTL ∪ QRMumps.RCNTL
      qrm_get(spfct, str)
      qrm_set(spfct, str, 1)
    end

    for str ∈ QRMumps.STATS
      qrm_get(spfct, str)
    end

    A = 2 * A
    qrm_update!(spmat, A)

    for ntype ∈ ('i', '1', 'f')
      qrm_spmat_nrm(spmat, ntype=ntype)
    end

    for ntype ∈ ('i', '1', '2')
      x = rand(T, 10)
      qrm_vecnrm(x, ntype=ntype)

      X = rand(T, 10, 5)
      nrm = qrm_vecnrm(X, ntype=ntype)

      X = 2 * X
      qrm_vecnrm!(X, nrm, ntype=ntype)
    end

    b = rand(T, n)
    x = rand(T, n)
    r = b - A * x
    qrm_residual_orth(spmat, r)
    qrm_residual_orth(spmat, r, transp=transp)
    qrm_residual_norm(spmat, b, x)
    qrm_residual_norm(spmat, b, x, transp=transp)

    B = rand(T, n, 5)
    X = rand(T, n, 5)
    R = B - A * X
    nrm = qrm_residual_orth(spmat, R)
    qrm_residual_orth!(spmat, R, nrm)
    qrm_residual_orth(spmat, R, transp=transp)
    nrm = qrm_residual_norm(spmat, B, X)
    qrm_residual_norm!(spmat, B, X, nrm)
    qrm_residual_norm(spmat, B, X, transp=transp)

    y = rand(T, n)
    y2 = copy(y)
    x = rand(T, n)
    Y = rand(T, n, 5)
    Y2 = copy(Y)
    X = rand(T, n, 5)
    α = rand(T)
    β = rand(T)

    mul!(y, spmat, x, α, β)
    mul!(y2, A, x, α, β)
    @test y ≈ y2
    mul!(Y, spmat, X, α, β)
    mul!(Y2, A, X, α, β)
    @test Y ≈ Y2

    mul!(y, Transpose(spmat), x, α, β)
    mul!(y2, Transpose(A), x, α, β)
    @test y ≈ y2
    mul!(Y, Transpose(spmat), X, α, β)
    mul!(Y2, Transpose(A), X, α, β)
    @test Y ≈ Y2

    mul!(y, Adjoint(spmat), x, α, β)
    mul!(y2, Adjoint(A), x, α, β)
    @test y ≈ y2
    mul!(Y, Adjoint(spmat), X, α, β)
    mul!(Y2, Adjoint(A), X, α, β)
    @test Y ≈ Y2
    
    A = sprand(T, m, n, 0.3)
    b = rand(T, n)
    spmat = qrm_spmat_init(A)
    spfct = qrm_analyse(spmat)
    qrm_set(spfct, "qrm_keeph", 0)
    qrm_factorize!(spmat, spfct)
    x₁ = qrm_solve(spfct,b, transp = transp)
    x = qrm_solve(spfct, x₁, transp = 'n')
    x_refined = qrm_refine(spmat, spfct, x, b)
    @test norm(b - A'*(A*x)) ≥ norm(b - A'*(A*x_refined))

  end
end

@testset "allocations" begin
  for T in (Float32, Float64, ComplexF32, ComplexF64)
    tol = (real(T) == Float32) ? 1e-3 : 1e-12
    transp = (T <: Real) ? 't' : 'c'

    @testset "$INT" for INT in (Int32, Int64)
      A = sprand(T, m, n, 0.3)
      A = convert(SparseMatrixCSC{T,INT}, A)

      spmat = qrm_spmat_init(T)
      qrm_spmat_init!(spmat, A)

      spfct = qrm_spfct_init(spmat)

      qrm_set(spfct, "qrm_rd_eps", tol)
      qrm_analyse!(spmat, spfct)
      qrm_factorize!(spmat, spfct)
      nbits = @allocated qrm_get(spfct, "qrm_rd_num")
      @test nbits == 0
    end
  end

end
