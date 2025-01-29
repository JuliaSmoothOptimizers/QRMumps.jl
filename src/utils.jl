function qrm_refine(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, x :: AbstractVector{T}, z :: AbstractVector{T}) where T
  Δx = similar(x)
  y = similar(x, spfct.fct.m)
  x_refined = copy(x)

  qrm_refine!(spmat, spfct, x_refined, z, Δx, y)
  return x_refined
end

function qrm_refine(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, x :: AbstractMatrix{T}, z :: AbstractMatrix{T}) where T
  Δx = similar(x)
  y = similar(x, (spfct.fct.m, size(x,2)))
  x_refined = copy(x)

  qrm_refine!(spmat, spfct, x_refined, z, Δx, y)
  return x_refined
end

function qrm_refine!(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, x :: AbstractVector{T}, z :: AbstractVector{T}, Δx :: AbstractVector{T}, y :: AbstractVector{T}) where T
  @assert length(x) == spfct.fct.n
  @assert length(z) == spfct.fct.n
  @assert length(Δx) == spfct.fct.n
  @assert length(y) == spfct.fct.m

  t = T <: Real ? 't' : 'c'
  transp = spfct.fct.n == spmat.mat.n ? 'n' : t # Check whether it is A^T = QR or A = QR.
  ntransp = transp == 't' || transp == 'c' ? 'n' : t

  qrm_spmat_mv!(spmat, T(1), x, T(0), y, transp = transp)
  qrm_spmat_mv!(spmat, T(1), y, T(0), Δx, transp = ntransp)
  @. Δx = z - Δx 

  qrm_solve!(spfct, Δx, y, transp = t)
  qrm_solve!(spfct, y, Δx, transp = 'n')
  @. x = x + Δx
end

function qrm_refine!(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, x :: AbstractMatrix{T}, z :: AbstractMatrix{T}, Δx :: AbstractMatrix{T}, y :: AbstractMatrix{T}) where T
  @assert size(x, 1) == spfct.fct.n
  @assert size(z, 1) == spfct.fct.n
  @assert size(Δx, 1) == spfct.fct.n
  @assert size(y, 1) == spfct.fct.m
  @assert size(x, 2) == size(z, 2)
  @assert size(Δx, 2) == size(z, 2)
  @assert size(y, 2) == size(z, 2)

  t = T <: Real ? 't' : 'c'
  transp = spfct.fct.n == spmat.mat.n ? 'n' : t # Check whether it is A^T = QR or A = QR.
  ntransp = transp == 't' || transp == 'c' ? 'n' : t

  qrm_spmat_mv!(spmat, T(1), x, T(0), y, transp = transp)
  qrm_spmat_mv!(spmat, T(1), y, T(0), Δx, transp = ntransp)
  @. Δx = z - Δx 

  qrm_solve!(spfct, Δx, y, transp = t)
  qrm_solve!(spfct, y, Δx, transp = 'n')
  @. x = x + Δx
end

function qrm_min_norm_semi_normal(spmat :: qrm_spmat{T}, b :: AbstractVecOrMat{T}; transp :: Char = 'n') where T

  n = transp == 'n' ? spmat.mat.n : spmat.mat.m

  spfct = qrm_spfct_init(spmat)

  if typeof(b) <: AbstractVector{T}
    x = similar(b, n)
  else 
    x = similar(b, (n, size(b, 2)))
  end
  Δx = similar(x)
  y = similar(b)

  qrm_min_norm_semi_normal!(spmat, spfct, b, x, Δx, y, transp = transp)
  return x
end

function qrm_min_norm_semi_normal!(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, b :: AbstractVector{T}, x :: AbstractVector{T}, Δx :: AbstractVector{T}, y :: AbstractVector{T}; transp :: Char = 'n') where T
  
  n = transp == 'n' ? spmat.mat.n : spmat.mat.m
  m = transp == 'n' ? spmat.mat.m : spmat.mat.n
  t = T <: Real ? 't' : 'c'
  ntransp = transp == 't' || transp == 'c' ? 'n' : t
  
  @assert length(x) == n
  @assert length(b) == m
  @assert length(Δx) == n
  @assert length(y) == m

  qrm_set(spfct, "qrm_keeph", 0)
  qrm_analyse!(spmat, spfct, transp = ntransp)
  qrm_factorize!(spmat, spfct, transp = ntransp)
  qrm_solve!(spfct, b, Δx, transp = t)
  qrm_solve!(spfct, Δx, y, transp = 'n')
  #x = A^T y
  qrm_spmat_mv!(spmat, T(1),  y, T(0), x, transp = ntransp)
end

function qrm_min_norm_semi_normal!(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, b :: AbstractMatrix{T}, x :: AbstractMatrix{T}, Δx :: AbstractMatrix{T}, y :: AbstractMatrix{T}; transp :: Char = 'n') where T
  
  n = transp == 'n' ? spmat.mat.n : spmat.mat.m
  m = transp == 'n' ? spmat.mat.m : spmat.mat.n
  t = T <: Real ? 't' : 'c'
  ntransp = transp == 't' || transp == 'c' ? 'n' : t
  
  @assert size(x, 1) == n
  @assert size(b, 1) == m
  @assert size(Δx, 1) == n
  @assert size(y, 1) == m
  @assert size(x, 2) == size(b, 2)
  @assert size(Δx, 2) == size(b, 2)
  @assert size(y, 2) == size(b, 2)

  qrm_set(spfct, "qrm_keeph", 0)
  qrm_analyse!(spmat, spfct, transp = ntransp)
  qrm_factorize!(spmat, spfct, transp = ntransp)
  qrm_solve!(spfct, b, Δx, transp = t)
  qrm_solve!(spfct, Δx, y, transp = 'n')
  #x = A^T y
  qrm_spmat_mv!(spmat, T(1),  y, T(0), x, transp = ntransp)
end

function qrm_least_squares_semi_normal(spmat :: qrm_spmat{T}, b :: AbstractVecOrMat{T}; transp :: Char = 'n') where T

  n = transp == 'n' ? spmat.mat.n : spmat.mat.m

  spfct = qrm_spfct_init(spmat)
  if typeof(b) <: AbstractVector{T}
    x = similar(b, n)
  else 
    x = similar(b, (n, size(b, 2)))
  end
  z = similar(x)
  Δx = similar(x)
  y = similar(b)

  qrm_least_squares_semi_normal!(spmat, spfct, b, x, z, Δx, y, transp = transp)
  return x
end

function qrm_least_squares_semi_normal!(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, b :: AbstractVector{T}, x :: AbstractVector{T}, z :: AbstractVector{T}, Δx :: AbstractVector{T}, y :: AbstractVector{T}; transp :: Char = 'n') where T
  
  n = transp == 'n' ? spmat.mat.n : spmat.mat.m
  m = transp == 'n' ? spmat.mat.m : spmat.mat.n
  t = T <: Real ? 't' : 'c'
  ntransp = transp == 't' || transp == 'c' ? 'n' : t

  @assert length(x) == n
  @assert length(b) == m
  @assert length(z) == n
  @assert length(Δx) == n
  @assert length(y) == m

  qrm_set(spfct, "qrm_keeph", 0)
  qrm_analyse!(spmat, spfct, transp  = transp)
  qrm_factorize!(spmat, spfct, transp = transp)
  # z = A^T b
  qrm_spmat_mv!(spmat, T(1),  b, T(0), z, transp = ntransp)
  #R^T y = z
  qrm_solve!(spfct, z, y, transp = t)
  qrm_solve!(spfct, y, x, transp = 'n')

  qrm_refine!(spmat, spfct, x, z, Δx, y)
end

function qrm_least_squares_semi_normal!(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, b :: AbstractMatrix{T}, x :: AbstractMatrix{T}, z :: AbstractMatrix{T}, Δx :: AbstractMatrix{T}, y :: AbstractMatrix{T}; transp :: Char = 'n') where T
  
  n = transp == 'n' ? spmat.mat.n : spmat.mat.m
  m = transp == 'n' ? spmat.mat.m : spmat.mat.n
  t = T <: Real ? 't' : 'c'
  ntransp = transp == 't' || transp == 'c' ? 'n' : t

  @assert size(x, 1) == n
  @assert size(b, 1) == m
  @assert size(z, 1) == n
  @assert size(Δx, 1) == n
  @assert size(y, 1) == m

  @assert size(x, 2) == size(b, 2)
  @assert size(z, 2) == size(b, 2)
  @assert size(Δx, 2) == size(b, 2)
  @assert size(y, 2) == size(b, 2)

  qrm_set(spfct, "qrm_keeph", 0)
  qrm_analyse!(spmat, spfct, transp  = transp)
  qrm_factorize!(spmat, spfct, transp = transp)
  # z = A^T b
  qrm_spmat_mv!(spmat, T(1),  b, T(0), z, transp = ntransp)
  qrm_solve!(spfct, z, y, transp = t)
  qrm_solve!(spfct, y, x, transp = 'n')

  qrm_refine!(spmat, spfct, x, z, Δx, y)
end

function qrm_shift_spmat(spmat :: qrm_spmat{T}, α :: T = T(0)) where T # Given a m×n matrix A, return the matrix Aₐ = [A √aI].
  @assert real(α) ≥ 0
  @assert imag(α) == 0

  m = spmat.mat.m
  n = spmat.mat.n
  irn = copy(spmat.irn)
  jcn = copy(spmat.jcn)
  val = copy(spmat.val)
  append!(irn, eltype(irn)(1):eltype(irn)(m))
  append!(jcn, eltype(irn)(1 + n):eltype(irn)(m + n))
  append!(val,  sqrt(α)*ones(T,m))
  shifted_spmat = qrm_spmat_init(m, n + m, irn, jcn, val)

  return qrm_shifted_spmat(shifted_spmat, α)
end

function qrm_update_shift_spmat!(shifted_spmat :: qrm_shifted_spmat{T}, α :: T) where T
  shifted_spmat.α = α
  shifted_spmat.spmat.val[shifted_spmat.spmat.mat.nz - shifted_spmat.spmat.mat.m + 1:end] .= sqrt(α)
end

function qrm_golub_riley(spmat :: qrm_spmat{T}, b :: AbstractVector{T}; α :: T = T(eps(real(T))), max_iter :: Int = 50, tol :: Real = eps(real(T)), transp :: Char = 'n') where T
  shifted_spmat = qrm_shift_spmat(spmat, α)
  spfct = qrm_spfct_init(shifted_spmat.spmat)
  n = shifted_spmat.spmat.mat.n
  m = shifted_spmat.spmat.mat.m

  x = similar(b, n)
  Δx = similar(b, n)
  y = similar(b)

  qrm_golub_riley!(
    shifted_spmat,
    spfct,
    x,
    b,
    Δx,
    y,
    α = α,
    max_iter = max_iter,
    tol = tol,
    transp = transp
  )
  return x[1:n-m]
end
#Given an underdetermined, rank defficient system Ax = b, compute the least-norm solution with Golub-Riley iteration
# TODO: add something for rank defficient least squares as well (i.e compute (Aᵀ)†z)
# TODO: add semi normal solves etc.
function qrm_golub_riley!(
  shifted_spmat :: qrm_shifted_spmat{T},
  spfct :: qrm_spfct{T},  
  x :: AbstractVector{T},
  b :: AbstractVector{T}, 
  Δx ::AbstractVector{T},
  y :: AbstractVector{T};  
  α :: T = T(eps(real(T))), 
  max_iter :: Int = 50, 
  tol :: Real = eps(real(T)),
  transp :: Char = 'n'
  ) where T
  
  spmat = shifted_spmat.spmat
  t = T <: Real ? 't' : 'c'
  ntransp = transp == 't' || transp == 'c' ? 'n' : t

  @assert real(α) ≥ 0
  @assert imag(α) == 0
  @assert length(b) == spmat.mat.m
  @assert length(x) == spmat.mat.n
  @assert length(y) == spmat.mat.m
  @assert length(Δx) == spmat.mat.n

  qrm_update_shift_spmat!(shifted_spmat, α)
  qrm_spfct_init!(spfct, spmat)
 
  qrm_set(spfct, "qrm_keeph", 0)
  qrm_analyse!(spmat, spfct, transp = transp)
  qrm_factorize!(spmat, spfct, transp = transp)

  k = 0
  solved = false
  x .= T(0)
  while k < max_iter && !solved

    qrm_spmat_mv!(spmat, T(1), x, T(0), y, transp = ntransp)
    @. y = b - y

    qrm_solve!(spfct, y, Δx, transp = t)
    qrm_solve!(spfct, Δx, y, transp = 'n')

    qrm_spmat_mv!(spmat, T(1), y, T(0), Δx, transp = transp)

    @. x = x + Δx
    
    solved = norm(Δx) ≤ tol*norm(x)
    k = k + 1
  end
end
