# Given an approximate solution x of the system RᵀRx ≈ z, refine the solution
function qrm_refine(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, x :: AbstractVector{T}, z :: AbstractVector{T}) where T
  Δx = similar(x)
  y = similar(x, spfct.fct.m)
  x_refined = copy(x)

  qrm_refine!(spmat, spfct, x_refined, z, Δx, y)
  return x_refined
end

function qrm_refine!(spmat :: qrm_spmat{T}, spfct :: qrm_spfct{T}, x :: AbstractVector{T}, z :: AbstractVector{T}, Δx :: AbstractVector{T}, y :: AbstractVector{T}) where T
  @assert length(x) == spfct.fct.n
  @assert length(z) == spfct.fct.n
  @assert length(Δx) == spfct.fct.n
  @assert length(y) == spfct.fct.m

  transp = T <: Real ? 't' : 'c'

  qrm_spmat_mv!(spmat, T(1), x, T(0), y, transp = 'n')
  qrm_spmat_mv!(spmat, T(1), y, T(0), Δx, transp = transp)
  @. Δx = z - Δx 

  qrm_solve!(spfct, Δx, y, transp = transp)
  qrm_solve!(spfct, y, Δx, transp = 'n')
  @. x = x + Δx
end