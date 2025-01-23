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
