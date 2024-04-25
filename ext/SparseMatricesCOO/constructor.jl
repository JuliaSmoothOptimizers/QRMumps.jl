function QRMumps.qrm_spmat_init(A :: SparseMatrixCOO{T,I}; sym :: Bool=false) where {T, I <: Integer}
    spmat = qrm_spmat{T}()
    qrm_spmat_init!(spmat, A; sym = sym)
    return spmat
end

function QRMumps.qrm_spmat_init!(spmat :: qrm_spmat{T}, A :: SparseMatrixCOO{T,I}; sym :: Bool=false) where {T, I <: Integer}
    qrm_spmat_init!(spmat, A.m, A.n, A.rows, A.cols, A.vals; sym = sym)
end