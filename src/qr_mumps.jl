module qr_mumps

using CEnum, Libdl

if haskey(ENV, "JULIA_QR_MUMPS_LIBRARY_PATH")
  println("Custom Installation")
  const libsqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libsqrm.$dlext")
  const libdqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libdqrm.$dlext")
  const libcqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libcqrm.$dlext")
  const libzqrm = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libzqrm.$dlext")
  const libqrm_common = joinpath(ENV["JULIA_QR_MUMPS_LIBRARY_PATH"], "libqrm_common.$dlext")
else
  using qr_mumps_jll
end

include("wrapper/qr_mumps_common.jl")
include("wrapper/qr_mumps_api.jl")

export qrm_spmat, qrm_spfct,
       qrm_spmat_init, qrm_spmat_destroy,
       qrm_spfct_init, qrm_spfct_destroy,
       qrm_readmat, qrm_analyse, qrm_factorize, qrm_solve,
       qrm_apply, qrm_matmul, qrm_spmat_nrm, qrm_vecnrm,
       qrm_spbackslash, qrm_spposv,
       qrm_least_squares, qrm_min_norm,
       qrm_residual_norm, qrm_residual_orth,
       qrm_spfct_seti, qrm_spfct_geti, qrm_spfct_getii,
       qrm_swtime, qrm_gseti, qrm_ggeti, qrm_ggetii,
       qrm_init, qrm_finalize

end # module
