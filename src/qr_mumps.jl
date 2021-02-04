module qr_mumps

using qr_mumps_jll, CEnum

include("wrapper/qr_mumps_common.jl")
include("wrapper/qr_mumps_api.jl")

export qrm_spmat, qrm_spfct,
       qrm_spmat_init, qrm_spmat_destroy,
       qrm_spfct_init, qrm_spfct_destroy,
       qrm_readmat, qrm_analyse, qrm_factorize, qrm_solve,
       qrm_apply, qrm_matmul, qrm_spmat_nrm, qrm_vecnrm,
       qrm_spbackslash, qrm_spposv,
       qrm_least_squares, qrm_min_norm,
       qrm_residual_norm, qrm_residual_orth
       # qrm_spfct_seti, qrm_spfct_geti, qrm_spfct_getii,
       # qrm_swtime, qrm_gseti, qrm_ggeti, qrm_ggetii,
       # qrm_init, qrm_finalize

end # module
