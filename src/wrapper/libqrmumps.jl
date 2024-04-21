function sqrm_spmat_init_c(spmat)
    @ccall libsqrm.sqrm_spmat_init_c(spmat::Ref{c_spmat{Float32}})::Cint
end

function sqrm_spmat_destroy_c(spmat)
    @ccall libsqrm.sqrm_spmat_destroy_c(spmat::Ref{c_spmat{Float32}})::Cint
end

function sqrm_spfct_init_c(spfct, spmat)
    @ccall libsqrm.sqrm_spfct_init_c(spfct::Ref{c_spfct{Float32}},
                                     spmat::Ref{c_spmat{Float32}})::Cint
end

function sqrm_spfct_destroy_c(spfct)
    @ccall libsqrm.sqrm_spfct_destroy_c(spfct::Ref{c_spfct{Float32}})::Cint
end

function sqrm_analyse_c(spmat, spfct, transp)
    @ccall libsqrm.sqrm_analyse_c(spmat::Ref{c_spmat{Float32}},
                                  spfct::Ref{c_spfct{Float32}}, transp::Cchar)::Cint
end

function sqrm_factorize_c(spmat, spfct, transp)
    @ccall libsqrm.sqrm_factorize_c(spmat::Ref{c_spmat{Float32}},
                                    spfct::Ref{c_spfct{Float32}}, transp::Cchar)::Cint
end

function sqrm_solve_c(spfct, transp, b, x, nrhs)
    @ccall libsqrm.sqrm_solve_c(spfct::Ref{c_spfct{Float32}}, transp::Cchar, b::Ptr{Cfloat},
                                x::Ptr{Cfloat}, nrhs::Cint)::Cint
end

function sqrm_apply_c(spfct, transp, b, nrhs)
    @ccall libsqrm.sqrm_apply_c(spfct::Ref{c_spfct{Float32}}, transp::Cchar, b::Ptr{Cfloat},
                                nrhs::Cint)::Cint
end

function sqrm_spmat_mv_c(spmat, transp, alpha, x, beta, y, nrhs)
    @ccall libsqrm.sqrm_spmat_mv_c(spmat::Ref{c_spmat{Float32}}, transp::Cchar,
                                   alpha::Cfloat, x::Ptr{Cfloat}, beta::Cfloat,
                                   y::Ptr{Cfloat}, nrhs::Cint)::Cint
end

function sqrm_spmat_nrm_c(spmat, ntype, nrm)
    @ccall libsqrm.sqrm_spmat_nrm_c(spmat::Ref{c_spmat{Float32}}, ntype::Cchar,
                                    nrm::Ptr{Cfloat})::Cint
end

function sqrm_vecnrm_c(x, n, nrhs, ntype, nrm)
    @ccall libsqrm.sqrm_vecnrm_c(x::Ptr{Cfloat}, n::Cint, nrhs::Cint, ntype::Cchar,
                                 nrm::Ptr{Cfloat})::Cint
end

function sqrm_spbackslash_c(spmat, b, x, nrhs, transp)
    @ccall libsqrm.sqrm_spbackslash_c(spmat::Ref{c_spmat{Float32}}, b::Ptr{Cfloat},
                                      x::Ptr{Cfloat}, nrhs::Cint, transp::Cchar)::Cint
end

function sqrm_spfct_backslash_c(spfct, b, x, nrhs, transp)
    @ccall libsqrm.sqrm_spfct_backslash_c(spfct::Ref{c_spfct{Float32}}, b::Ptr{Cfloat},
                                          x::Ptr{Cfloat}, nrhs::Cint, transp::Cchar)::Cint
end

function sqrm_spposv_c(spmat, b, x, nrhs)
    @ccall libsqrm.sqrm_spposv_c(spmat::Ref{c_spmat{Float32}}, b::Ptr{Cfloat},
                                 x::Ptr{Cfloat}, nrhs::Cint)::Cint
end

function sqrm_least_squares_c(spmat, b, x, nrhs, transp)
    @ccall libsqrm.sqrm_least_squares_c(spmat::Ref{c_spmat{Float32}}, b::Ptr{Cfloat},
                                        x::Ptr{Cfloat}, nrhs::Cint, transp::Cchar)::Cint
end

function sqrm_min_norm_c(spmat, b, x, nrhs, transp)
    @ccall libsqrm.sqrm_min_norm_c(spmat::Ref{c_spmat{Float32}}, b::Ptr{Cfloat},
                                   x::Ptr{Cfloat}, nrhs::Cint, transp::Cchar)::Cint
end

function sqrm_residual_norm_c(spmat, b, x, nrhs, nrm, transp)
    @ccall libsqrm.sqrm_residual_norm_c(spmat::Ref{c_spmat{Float32}}, b::Ptr{Cfloat},
                                        x::Ptr{Cfloat}, nrhs::Cint, nrm::Ptr{Cfloat},
                                        transp::Cchar)::Cint
end

function sqrm_residual_orth_c(spmat, r, nrhs, nrm, transp)
    @ccall libsqrm.sqrm_residual_orth_c(spmat::Ref{c_spmat{Float32}}, r::Ptr{Cfloat},
                                        nrhs::Cint, nrm::Ptr{Cfloat}, transp::Cchar)::Cint
end

function sqrm_spfct_trsm_c(spfct, transp, b, x, nrhs)
    @ccall libsqrm.sqrm_spfct_trsm_c(spfct::Ref{c_spfct{Float32}}, transp::Cchar,
                                     b::Ptr{Cfloat}, x::Ptr{Cfloat}, nrhs::Cint)::Cint
end

function sqrm_spfct_sytrs_c(spfct, b, x, nrhs)
    @ccall libsqrm.sqrm_spfct_sytrs_c(spfct::Ref{c_spfct{Float32}}, b::Ptr{Cfloat},
                                      x::Ptr{Cfloat}, nrhs::Cint)::Cint
end

function sqrm_spfct_set_i4_c(spfct, string, val)
    @ccall libsqrm.sqrm_spfct_set_i4_c(spfct::Ref{c_spfct{Float32}}, string::Cstring,
                                       val::Cint)::Cint
end

function sqrm_spfct_set_r4_c(spfct, string, val)
    @ccall libsqrm.sqrm_spfct_set_r4_c(spfct::Ref{c_spfct{Float32}}, string::Cstring,
                                       val::Cfloat)::Cint
end

function sqrm_spfct_get_i4_c(spfct, string, val)
    @ccall libsqrm.sqrm_spfct_get_i4_c(spfct::Ref{c_spfct{Float32}}, string::Cstring,
                                       val::Ptr{Cint})::Cint
end

function sqrm_spfct_get_r4_c(spfct, string, val)
    @ccall libsqrm.sqrm_spfct_get_r4_c(spfct::Ref{c_spfct{Float32}}, string::Cstring,
                                       val::Ptr{Cfloat})::Cint
end

function sqrm_spfct_get_i8_c(spfct, string, val)
    @ccall libsqrm.sqrm_spfct_get_i8_c(spfct::Ref{c_spfct{Float32}}, string::Cstring,
                                       val::Ptr{Clonglong})::Cint
end

function sqrm_spfct_get_schur_c(spfct, s, i, j, m, n)
    @ccall libsqrm.sqrm_spfct_get_schur_c(spfct::Ref{c_spfct{Float32}}, s::Ptr{Cfloat},
                                          i::Cint, j::Cint, m::Cint, n::Cint)::Cint
end

function sqrm_spfct_get_r_c(spfct, spmat)
    @ccall libsqrm.sqrm_spfct_get_r_c(spfct::Ref{c_spfct{Float32}},
                                      spmat::Ref{c_spmat{Float32}})::Cint
end

function sqrm_spfct_get_cp_c(spfct, cp)
    @ccall libsqrm.sqrm_spfct_get_cp_c(spfct::Ref{c_spfct{Float32}},
                                       cp::Ptr{Ptr{Cint}})::Cint
end

function sqrm_spfct_get_rp_c(spfct, rp)
    @ccall libsqrm.sqrm_spfct_get_rp_c(spfct::Ref{c_spfct{Float32}},
                                       rp::Ptr{Ptr{Cint}})::Cint
end

function dqrm_spmat_init_c(spmat)
    @ccall libdqrm.dqrm_spmat_init_c(spmat::Ref{c_spmat{Float64}})::Cint
end

function dqrm_spmat_destroy_c(spmat)
    @ccall libdqrm.dqrm_spmat_destroy_c(spmat::Ref{c_spmat{Float64}})::Cint
end

function dqrm_spfct_init_c(spfct, spmat)
    @ccall libdqrm.dqrm_spfct_init_c(spfct::Ref{c_spfct{Float64}},
                                     spmat::Ref{c_spmat{Float64}})::Cint
end

function dqrm_spfct_destroy_c(spfct)
    @ccall libdqrm.dqrm_spfct_destroy_c(spfct::Ref{c_spfct{Float64}})::Cint
end

function dqrm_analyse_c(spmat, spfct, transp)
    @ccall libdqrm.dqrm_analyse_c(spmat::Ref{c_spmat{Float64}},
                                  spfct::Ref{c_spfct{Float64}}, transp::Cchar)::Cint
end

function dqrm_factorize_c(spmat, spfct, transp)
    @ccall libdqrm.dqrm_factorize_c(spmat::Ref{c_spmat{Float64}},
                                    spfct::Ref{c_spfct{Float64}}, transp::Cchar)::Cint
end

function dqrm_solve_c(spfct, transp, b, x, nrhs)
    @ccall libdqrm.dqrm_solve_c(spfct::Ref{c_spfct{Float64}}, transp::Cchar,
                                b::Ptr{Cdouble}, x::Ptr{Cdouble}, nrhs::Cint)::Cint
end

function dqrm_apply_c(spfct, transp, b, nrhs)
    @ccall libdqrm.dqrm_apply_c(spfct::Ref{c_spfct{Float64}}, transp::Cchar,
                                b::Ptr{Cdouble}, nrhs::Cint)::Cint
end

function dqrm_spmat_mv_c(spmat, transp, alpha, x, beta, y, nrhs)
    @ccall libdqrm.dqrm_spmat_mv_c(spmat::Ref{c_spmat{Float64}}, transp::Cchar,
                                   alpha::Cdouble, x::Ptr{Cdouble}, beta::Cdouble,
                                   y::Ptr{Cdouble}, nrhs::Cint)::Cint
end

function dqrm_spmat_nrm_c(spmat, ntype, nrm)
    @ccall libdqrm.dqrm_spmat_nrm_c(spmat::Ref{c_spmat{Float64}}, ntype::Cchar,
                                    nrm::Ptr{Cdouble})::Cint
end

function dqrm_vecnrm_c(x, n, nrhs, ntype, nrm)
    @ccall libdqrm.dqrm_vecnrm_c(x::Ptr{Cdouble}, n::Cint, nrhs::Cint, ntype::Cchar,
                                 nrm::Ptr{Cdouble})::Cint
end

function dqrm_spbackslash_c(spmat, b, x, nrhs, transp)
    @ccall libdqrm.dqrm_spbackslash_c(spmat::Ref{c_spmat{Float64}}, b::Ptr{Cdouble},
                                      x::Ptr{Cdouble}, nrhs::Cint, transp::Cchar)::Cint
end

function dqrm_spfct_backslash_c(spfct, b, x, nrhs, transp)
    @ccall libdqrm.dqrm_spfct_backslash_c(spfct::Ref{c_spfct{Float64}}, b::Ptr{Cdouble},
                                          x::Ptr{Cdouble}, nrhs::Cint, transp::Cchar)::Cint
end

function dqrm_spposv_c(spmat, b, x, nrhs)
    @ccall libdqrm.dqrm_spposv_c(spmat::Ref{c_spmat{Float64}}, b::Ptr{Cdouble},
                                 x::Ptr{Cdouble}, nrhs::Cint)::Cint
end

function dqrm_least_squares_c(spmat, b, x, nrhs, transp)
    @ccall libdqrm.dqrm_least_squares_c(spmat::Ref{c_spmat{Float64}}, b::Ptr{Cdouble},
                                        x::Ptr{Cdouble}, nrhs::Cint, transp::Cchar)::Cint
end

function dqrm_min_norm_c(spmat, b, x, nrhs, transp)
    @ccall libdqrm.dqrm_min_norm_c(spmat::Ref{c_spmat{Float64}}, b::Ptr{Cdouble},
                                   x::Ptr{Cdouble}, nrhs::Cint, transp::Cchar)::Cint
end

function dqrm_residual_norm_c(spmat, b, x, nrhs, nrm, transp)
    @ccall libdqrm.dqrm_residual_norm_c(spmat::Ref{c_spmat{Float64}}, b::Ptr{Cdouble},
                                        x::Ptr{Cdouble}, nrhs::Cint, nrm::Ptr{Cdouble},
                                        transp::Cchar)::Cint
end

function dqrm_residual_orth_c(spmat, r, nrhs, nrm, transp)
    @ccall libdqrm.dqrm_residual_orth_c(spmat::Ref{c_spmat{Float64}}, r::Ptr{Cdouble},
                                        nrhs::Cint, nrm::Ptr{Cdouble}, transp::Cchar)::Cint
end

function dqrm_spfct_trsm_c(spfct, transp, b, x, nrhs)
    @ccall libdqrm.dqrm_spfct_trsm_c(spfct::Ref{c_spfct{Float64}}, transp::Cchar,
                                     b::Ptr{Cdouble}, x::Ptr{Cdouble}, nrhs::Cint)::Cint
end

function dqrm_spfct_sytrs_c(spfct, b, x, nrhs)
    @ccall libdqrm.dqrm_spfct_sytrs_c(spfct::Ref{c_spfct{Float64}}, b::Ptr{Cdouble},
                                      x::Ptr{Cdouble}, nrhs::Cint)::Cint
end

function dqrm_spfct_set_i4_c(spfct, string, val)
    @ccall libdqrm.dqrm_spfct_set_i4_c(spfct::Ref{c_spfct{Float64}}, string::Cstring,
                                       val::Cint)::Cint
end

function dqrm_spfct_set_r4_c(spfct, string, val)
    @ccall libdqrm.dqrm_spfct_set_r4_c(spfct::Ref{c_spfct{Float64}}, string::Cstring,
                                       val::Cfloat)::Cint
end

function dqrm_spfct_get_i4_c(spfct, string, val)
    @ccall libdqrm.dqrm_spfct_get_i4_c(spfct::Ref{c_spfct{Float64}}, string::Cstring,
                                       val::Ptr{Cint})::Cint
end

function dqrm_spfct_get_r4_c(spfct, string, val)
    @ccall libdqrm.dqrm_spfct_get_r4_c(spfct::Ref{c_spfct{Float64}}, string::Cstring,
                                       val::Ptr{Cfloat})::Cint
end

function dqrm_spfct_get_i8_c(spfct, string, val)
    @ccall libdqrm.dqrm_spfct_get_i8_c(spfct::Ref{c_spfct{Float64}}, string::Cstring,
                                       val::Ptr{Clonglong})::Cint
end

function dqrm_spfct_get_schur_c(spfct, s, i, j, m, n)
    @ccall libdqrm.dqrm_spfct_get_schur_c(spfct::Ref{c_spfct{Float64}}, s::Ptr{Cdouble},
                                          i::Cint, j::Cint, m::Cint, n::Cint)::Cint
end

function dqrm_spfct_get_r_c(spfct, spmat)
    @ccall libdqrm.dqrm_spfct_get_r_c(spfct::Ref{c_spfct{Float64}},
                                      spmat::Ref{c_spmat{Float64}})::Cint
end

function dqrm_spfct_get_cp_c(spfct, cp)
    @ccall libdqrm.dqrm_spfct_get_cp_c(spfct::Ref{c_spfct{Float64}},
                                       cp::Ptr{Ptr{Cint}})::Cint
end

function dqrm_spfct_get_rp_c(spfct, rp)
    @ccall libdqrm.dqrm_spfct_get_rp_c(spfct::Ref{c_spfct{Float64}},
                                       rp::Ptr{Ptr{Cint}})::Cint
end

function cqrm_spmat_init_c(spmat)
    @ccall libcqrm.cqrm_spmat_init_c(spmat::Ref{c_spmat{ComplexF32}})::Cint
end

function cqrm_spmat_destroy_c(spmat)
    @ccall libcqrm.cqrm_spmat_destroy_c(spmat::Ref{c_spmat{ComplexF32}})::Cint
end

function cqrm_spfct_init_c(spfct, spmat)
    @ccall libcqrm.cqrm_spfct_init_c(spfct::Ref{c_spfct{ComplexF32}},
                                     spmat::Ref{c_spmat{ComplexF32}})::Cint
end

function cqrm_spfct_destroy_c(spfct)
    @ccall libcqrm.cqrm_spfct_destroy_c(spfct::Ref{c_spfct{ComplexF32}})::Cint
end

function cqrm_analyse_c(spmat, spfct, transp)
    @ccall libcqrm.cqrm_analyse_c(spmat::Ref{c_spmat{ComplexF32}},
                                  spfct::Ref{c_spfct{ComplexF32}}, transp::Cchar)::Cint
end

function cqrm_factorize_c(spmat, spfct, transp)
    @ccall libcqrm.cqrm_factorize_c(spmat::Ref{c_spmat{ComplexF32}},
                                    spfct::Ref{c_spfct{ComplexF32}}, transp::Cchar)::Cint
end

function cqrm_solve_c(spfct, transp, b, x, nrhs)
    @ccall libcqrm.cqrm_solve_c(spfct::Ref{c_spfct{ComplexF32}}, transp::Cchar,
                                b::Ptr{ComplexF32}, x::Ptr{ComplexF32}, nrhs::Cint)::Cint
end

function cqrm_apply_c(spfct, transp, b, nrhs)
    @ccall libcqrm.cqrm_apply_c(spfct::Ref{c_spfct{ComplexF32}}, transp::Cchar,
                                b::Ptr{ComplexF32}, nrhs::Cint)::Cint
end

function cqrm_spmat_mv_c(spmat, transp, alpha, x, beta, y, nrhs)
    @ccall libcqrm.cqrm_spmat_mv_c(spmat::Ref{c_spmat{ComplexF32}}, transp::Cchar,
                                   alpha::ComplexF32, x::Ptr{ComplexF32}, beta::ComplexF32,
                                   y::Ptr{ComplexF32}, nrhs::Cint)::Cint
end

function cqrm_spmat_nrm_c(spmat, ntype, nrm)
    @ccall libcqrm.cqrm_spmat_nrm_c(spmat::Ref{c_spmat{ComplexF32}}, ntype::Cchar,
                                    nrm::Ptr{Cfloat})::Cint
end

function cqrm_vecnrm_c(x, n, nrhs, ntype, nrm)
    @ccall libcqrm.cqrm_vecnrm_c(x::Ptr{ComplexF32}, n::Cint, nrhs::Cint, ntype::Cchar,
                                 nrm::Ptr{Cfloat})::Cint
end

function cqrm_spbackslash_c(spmat, b, x, nrhs, transp)
    @ccall libcqrm.cqrm_spbackslash_c(spmat::Ref{c_spmat{ComplexF32}}, b::Ptr{ComplexF32},
                                      x::Ptr{ComplexF32}, nrhs::Cint, transp::Cchar)::Cint
end

function cqrm_spfct_backslash_c(spfct, b, x, nrhs, transp)
    @ccall libcqrm.cqrm_spfct_backslash_c(spfct::Ref{c_spfct{ComplexF32}},
                                          b::Ptr{ComplexF32}, x::Ptr{ComplexF32},
                                          nrhs::Cint, transp::Cchar)::Cint
end

function cqrm_spposv_c(spmat, b, x, nrhs)
    @ccall libcqrm.cqrm_spposv_c(spmat::Ref{c_spmat{ComplexF32}}, b::Ptr{ComplexF32},
                                 x::Ptr{ComplexF32}, nrhs::Cint)::Cint
end

function cqrm_least_squares_c(spmat, b, x, nrhs, transp)
    @ccall libcqrm.cqrm_least_squares_c(spmat::Ref{c_spmat{ComplexF32}}, b::Ptr{ComplexF32},
                                        x::Ptr{ComplexF32}, nrhs::Cint, transp::Cchar)::Cint
end

function cqrm_min_norm_c(spmat, b, x, nrhs, transp)
    @ccall libcqrm.cqrm_min_norm_c(spmat::Ref{c_spmat{ComplexF32}}, b::Ptr{ComplexF32},
                                   x::Ptr{ComplexF32}, nrhs::Cint, transp::Cchar)::Cint
end

function cqrm_residual_norm_c(spmat, b, x, nrhs, nrm, transp)
    @ccall libcqrm.cqrm_residual_norm_c(spmat::Ref{c_spmat{ComplexF32}}, b::Ptr{ComplexF32},
                                        x::Ptr{ComplexF32}, nrhs::Cint, nrm::Ptr{Cfloat},
                                        transp::Cchar)::Cint
end

function cqrm_residual_orth_c(spmat, r, nrhs, nrm, transp)
    @ccall libcqrm.cqrm_residual_orth_c(spmat::Ref{c_spmat{ComplexF32}}, r::Ptr{ComplexF32},
                                        nrhs::Cint, nrm::Ptr{Cfloat}, transp::Cchar)::Cint
end

function cqrm_spfct_trsm_c(spfct, transp, b, x, nrhs)
    @ccall libcqrm.cqrm_spfct_trsm_c(spfct::Ref{c_spfct{ComplexF32}}, transp::Cchar,
                                     b::Ptr{ComplexF32}, x::Ptr{ComplexF32},
                                     nrhs::Cint)::Cint
end

function cqrm_spfct_sytrs_c(spfct, b, x, nrhs)
    @ccall libcqrm.cqrm_spfct_sytrs_c(spfct::Ref{c_spfct{ComplexF32}}, b::Ptr{ComplexF32},
                                      x::Ptr{ComplexF32}, nrhs::Cint)::Cint
end

function cqrm_spfct_set_i4_c(spfct, string, val)
    @ccall libcqrm.cqrm_spfct_set_i4_c(spfct::Ref{c_spfct{ComplexF32}}, string::Cstring,
                                       val::Cint)::Cint
end

function cqrm_spfct_set_r4_c(spfct, string, val)
    @ccall libcqrm.cqrm_spfct_set_r4_c(spfct::Ref{c_spfct{ComplexF32}}, string::Cstring,
                                       val::Cfloat)::Cint
end

function cqrm_spfct_get_i4_c(spfct, string, val)
    @ccall libcqrm.cqrm_spfct_get_i4_c(spfct::Ref{c_spfct{ComplexF32}}, string::Cstring,
                                       val::Ptr{Cint})::Cint
end

function cqrm_spfct_get_r4_c(spfct, string, val)
    @ccall libcqrm.cqrm_spfct_get_r4_c(spfct::Ref{c_spfct{ComplexF32}}, string::Cstring,
                                       val::Ptr{Cfloat})::Cint
end

function cqrm_spfct_get_i8_c(spfct, string, val)
    @ccall libcqrm.cqrm_spfct_get_i8_c(spfct::Ref{c_spfct{ComplexF32}}, string::Cstring,
                                       val::Ptr{Clonglong})::Cint
end

function cqrm_spfct_get_schur_c(spfct, s, i, j, m, n)
    @ccall libcqrm.cqrm_spfct_get_schur_c(spfct::Ref{c_spfct{ComplexF32}},
                                          s::Ptr{ComplexF32}, i::Cint, j::Cint, m::Cint,
                                          n::Cint)::Cint
end

function cqrm_spfct_get_r_c(spfct, spmat)
    @ccall libcqrm.cqrm_spfct_get_r_c(spfct::Ref{c_spfct{ComplexF32}},
                                      spmat::Ref{c_spmat{ComplexF32}})::Cint
end

function cqrm_spfct_get_cp_c(spfct, cp)
    @ccall libcqrm.cqrm_spfct_get_cp_c(spfct::Ref{c_spfct{ComplexF32}},
                                       cp::Ptr{Ptr{Cint}})::Cint
end

function cqrm_spfct_get_rp_c(spfct, rp)
    @ccall libcqrm.cqrm_spfct_get_rp_c(spfct::Ref{c_spfct{ComplexF32}},
                                       rp::Ptr{Ptr{Cint}})::Cint
end

function zqrm_spmat_init_c(spmat)
    @ccall libzqrm.zqrm_spmat_init_c(spmat::Ref{c_spmat{ComplexF64}})::Cint
end

function zqrm_spmat_destroy_c(spmat)
    @ccall libzqrm.zqrm_spmat_destroy_c(spmat::Ref{c_spmat{ComplexF64}})::Cint
end

function zqrm_spfct_init_c(spfct, spmat)
    @ccall libzqrm.zqrm_spfct_init_c(spfct::Ref{c_spfct{ComplexF64}},
                                     spmat::Ref{c_spmat{ComplexF64}})::Cint
end

function zqrm_spfct_destroy_c(spfct)
    @ccall libzqrm.zqrm_spfct_destroy_c(spfct::Ref{c_spfct{ComplexF64}})::Cint
end

function zqrm_analyse_c(spmat, spfct, transp)
    @ccall libzqrm.zqrm_analyse_c(spmat::Ref{c_spmat{ComplexF64}},
                                  spfct::Ref{c_spfct{ComplexF64}}, transp::Cchar)::Cint
end

function zqrm_factorize_c(spmat, spfct, transp)
    @ccall libzqrm.zqrm_factorize_c(spmat::Ref{c_spmat{ComplexF64}},
                                    spfct::Ref{c_spfct{ComplexF64}}, transp::Cchar)::Cint
end

function zqrm_solve_c(spfct, transp, b, x, nrhs)
    @ccall libzqrm.zqrm_solve_c(spfct::Ref{c_spfct{ComplexF64}}, transp::Cchar,
                                b::Ptr{ComplexF64}, x::Ptr{ComplexF64}, nrhs::Cint)::Cint
end

function zqrm_apply_c(spfct, transp, b, nrhs)
    @ccall libzqrm.zqrm_apply_c(spfct::Ref{c_spfct{ComplexF64}}, transp::Cchar,
                                b::Ptr{ComplexF64}, nrhs::Cint)::Cint
end

function zqrm_spmat_mv_c(spmat, transp, alpha, x, beta, y, nrhs)
    @ccall libzqrm.zqrm_spmat_mv_c(spmat::Ref{c_spmat{ComplexF64}}, transp::Cchar,
                                   alpha::ComplexF64, x::Ptr{ComplexF64}, beta::ComplexF64,
                                   y::Ptr{ComplexF64}, nrhs::Cint)::Cint
end

function zqrm_spmat_nrm_c(spmat, ntype, nrm)
    @ccall libzqrm.zqrm_spmat_nrm_c(spmat::Ref{c_spmat{ComplexF64}}, ntype::Cchar,
                                    nrm::Ptr{Cdouble})::Cint
end

function zqrm_vecnrm_c(x, n, nrhs, ntype, nrm)
    @ccall libzqrm.zqrm_vecnrm_c(x::Ptr{ComplexF64}, n::Cint, nrhs::Cint, ntype::Cchar,
                                 nrm::Ptr{Cdouble})::Cint
end

function zqrm_spbackslash_c(spmat, b, x, nrhs, transp)
    @ccall libzqrm.zqrm_spbackslash_c(spmat::Ref{c_spmat{ComplexF64}}, b::Ptr{ComplexF64},
                                      x::Ptr{ComplexF64}, nrhs::Cint, transp::Cchar)::Cint
end

function zqrm_spfct_backslash_c(spfct, b, x, nrhs, transp)
    @ccall libzqrm.zqrm_spfct_backslash_c(spfct::Ref{c_spfct{ComplexF64}},
                                          b::Ptr{ComplexF64}, x::Ptr{ComplexF64},
                                          nrhs::Cint, transp::Cchar)::Cint
end

function zqrm_spposv_c(spmat, b, x, nrhs)
    @ccall libzqrm.zqrm_spposv_c(spmat::Ref{c_spmat{ComplexF64}}, b::Ptr{ComplexF64},
                                 x::Ptr{ComplexF64}, nrhs::Cint)::Cint
end

function zqrm_least_squares_c(spmat, b, x, nrhs, transp)
    @ccall libzqrm.zqrm_least_squares_c(spmat::Ref{c_spmat{ComplexF64}}, b::Ptr{ComplexF64},
                                        x::Ptr{ComplexF64}, nrhs::Cint, transp::Cchar)::Cint
end

function zqrm_min_norm_c(spmat, b, x, nrhs, transp)
    @ccall libzqrm.zqrm_min_norm_c(spmat::Ref{c_spmat{ComplexF64}}, b::Ptr{ComplexF64},
                                   x::Ptr{ComplexF64}, nrhs::Cint, transp::Cchar)::Cint
end

function zqrm_residual_norm_c(spmat, b, x, nrhs, nrm, transp)
    @ccall libzqrm.zqrm_residual_norm_c(spmat::Ref{c_spmat{ComplexF64}}, b::Ptr{ComplexF64},
                                        x::Ptr{ComplexF64}, nrhs::Cint, nrm::Ptr{Cdouble},
                                        transp::Cchar)::Cint
end

function zqrm_residual_orth_c(spmat, r, nrhs, nrm, transp)
    @ccall libzqrm.zqrm_residual_orth_c(spmat::Ref{c_spmat{ComplexF64}}, r::Ptr{ComplexF64},
                                        nrhs::Cint, nrm::Ptr{Cdouble}, transp::Cchar)::Cint
end

function zqrm_spfct_trsm_c(spfct, transp, b, x, nrhs)
    @ccall libzqrm.zqrm_spfct_trsm_c(spfct::Ref{c_spfct{ComplexF64}}, transp::Cchar,
                                     b::Ptr{ComplexF64}, x::Ptr{ComplexF64},
                                     nrhs::Cint)::Cint
end

function zqrm_spfct_sytrs_c(spfct, b, x, nrhs)
    @ccall libzqrm.zqrm_spfct_sytrs_c(spfct::Ref{c_spfct{ComplexF64}}, b::Ptr{ComplexF64},
                                      x::Ptr{ComplexF64}, nrhs::Cint)::Cint
end

function zqrm_spfct_set_i4_c(spfct, string, val)
    @ccall libzqrm.zqrm_spfct_set_i4_c(spfct::Ref{c_spfct{ComplexF64}}, string::Cstring,
                                       val::Cint)::Cint
end

function zqrm_spfct_set_r4_c(spfct, string, val)
    @ccall libzqrm.zqrm_spfct_set_r4_c(spfct::Ref{c_spfct{ComplexF64}}, string::Cstring,
                                       val::Cfloat)::Cint
end

function zqrm_spfct_get_i4_c(spfct, string, val)
    @ccall libzqrm.zqrm_spfct_get_i4_c(spfct::Ref{c_spfct{ComplexF64}}, string::Cstring,
                                       val::Ptr{Cint})::Cint
end

function zqrm_spfct_get_r4_c(spfct, string, val)
    @ccall libzqrm.zqrm_spfct_get_r4_c(spfct::Ref{c_spfct{ComplexF64}}, string::Cstring,
                                       val::Ptr{Cfloat})::Cint
end

function zqrm_spfct_get_i8_c(spfct, string, val)
    @ccall libzqrm.zqrm_spfct_get_i8_c(spfct::Ref{c_spfct{ComplexF64}}, string::Cstring,
                                       val::Ptr{Clonglong})::Cint
end

function zqrm_spfct_get_schur_c(spfct, s, i, j, m, n)
    @ccall libzqrm.zqrm_spfct_get_schur_c(spfct::Ref{c_spfct{ComplexF64}},
                                          s::Ptr{ComplexF64}, i::Cint, j::Cint, m::Cint,
                                          n::Cint)::Cint
end

function zqrm_spfct_get_r_c(spfct, spmat)
    @ccall libzqrm.zqrm_spfct_get_r_c(spfct::Ref{c_spfct{ComplexF64}},
                                      spmat::Ref{c_spmat{ComplexF64}})::Cint
end

function zqrm_spfct_get_cp_c(spfct, cp)
    @ccall libzqrm.zqrm_spfct_get_cp_c(spfct::Ref{c_spfct{ComplexF64}},
                                       cp::Ptr{Ptr{Cint}})::Cint
end

function zqrm_spfct_get_rp_c(spfct, rp)
    @ccall libzqrm.zqrm_spfct_get_rp_c(spfct::Ref{c_spfct{ComplexF64}},
                                       rp::Ptr{Ptr{Cint}})::Cint
end

function qrm_swtime()
    @ccall libqrm_common.qrm_swtime()::Cdouble
end

function qrm_glob_set_i4_c(string, val)
    @ccall libqrm_common.qrm_glob_set_i4_c(string::Cstring, val::Cint)::Cint
end

function qrm_glob_set_r4_c(string, val)
    @ccall libqrm_common.qrm_glob_set_r4_c(string::Cstring, val::Cfloat)::Cint
end

function qrm_glob_get_i4_c(string, val)
    @ccall libqrm_common.qrm_glob_get_i4_c(string::Cstring, val::Ptr{Cint})::Cint
end

function qrm_glob_get_r4_c(string, val)
    @ccall libqrm_common.qrm_glob_get_r4_c(string::Cstring, val::Ptr{Cfloat})::Cint
end

function qrm_glob_get_i8_c(string, val)
    @ccall libqrm_common.qrm_glob_get_i8_c(string::Cstring, val::Ptr{Clonglong})::Cint
end

function qrm_init_c(ncpu, ngpu)
    @ccall libqrm_common.qrm_init_c(ncpu::Cint, ngpu::Cint)::Cint
end

function qrm_finalize_c()
    @ccall libqrm_common.qrm_finalize_c()::Cvoid
end
