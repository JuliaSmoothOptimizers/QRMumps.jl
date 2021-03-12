@doc raw"""
This data type is used to store a sparse matrix in the COO (or coordinate) format through the irn, jcn
and val fields containing the row indices, column indices and values, respectively and the m, n and nz
containing the number of rows, columns and nonzeros, respectively. qr mumps uses a Fortran-style
1-based numbering and thus all row indices are expected to be between 1 and m and all the column
indices between 1 and n. Duplicate entries are summed during the factorization, out-of-bound entries
are ignored. The sym field is used to specify if the matrix is symmetric (> 0) or unsymmetric (= 0).
"""
mutable struct qrm_spmat{T} <: AbstractSparseMatrix{T, Cint}
  irn :: Ptr{Cint}
  jcn :: Ptr{Cint}
  val :: Ptr{T}
  m   :: Cint
  n   :: Cint
  nz  :: Cint
  sym :: Cint
  h   :: Ptr{Cvoid}

  function qrm_spmat{T}() where T
    spmat = new(C_NULL, C_NULL, C_NULL, 0, 0, 0, 0, C_NULL)
    return spmat
  end
end

Base.size(spmat :: qrm_spmat) = (spmat.m, spmat.n)
SparseArrays.nnz(spmat :: qrm_spmat) = spmat.nz

function Base.show(io :: IO, ::MIME"text/plain", spmat :: qrm_spmat)
  println(io, "Sparse matrix -- qrm_spmat of size ", size(spmat), " with ", nnz(spmat), " nonzeros.")
end

@doc raw"""
This type is used to set the parameters that control the behavior of a sparse factorization, to collect
information about its execution (number of flops, memory consumpnion etc) and store the result of 
the factorization, namely, the factors with all the symbolic information needed to use them in the
solve phase.
"""
mutable struct qrm_spfct{T} <: Factorization{T}
  m        :: Cint
  n        :: Cint
  nz       :: Cint
  sym      :: Cint
  cperm_in :: Ptr{Cint}
  icntl    :: NTuple{20, Cint}
  rcntl    :: NTuple{10, Cfloat}
  gstats   :: NTuple{10, Clonglong}
  h        :: Ptr{Cvoid}

  function qrm_spfct{T}() where T
    spfct = new(0, 0, 0, 0, C_NULL, ntuple(x -> Cint(0), 20), ntuple(x -> Cfloat(0), 10), ntuple(x -> Clonglong(0), 10), C_NULL)
    return spfct
  end
end

function Base.show(io :: IO, ::MIME"text/plain", spfct :: qrm_spfct)
  println(io, "Sparse factorization -- qrm_spfct")
end

const qrm_no_transp   = 'n'
const qrm_transp      = 't'
const qrm_conj_transp = 'c'

@cenum icntl::UInt32 begin
    qrm_ordering_ = 0
    qrm_sing_     = 1
    qrm_minamalg_ = 2
    qrm_mb_       = 3
    qrm_nb_       = 4
    qrm_ib_       = 5
    qrm_bh_       = 6
    qrm_keeph_    = 7
    qrm_rhsnb_    = 8
end

@cenum rcntl::UInt32 begin
    qrm_amalgthr_  = 0
    qrm_mem_relax_ = 1
    qrm_rd_eps_    = 2
end

@cenum ords::UInt32 begin
    qrm_auto     = 0
    qrm_natural_ = 1
    qrm_given_   = 2
    qrm_colamd_  = 3
    qrm_metis_   = 4
    qrm_scotch_  = 5
end

@cenum gstats::UInt32 begin
    qrm_e_facto_flops_   = 0
    qrm_e_nnz_r_         = 1
    qrm_e_nnz_h_         = 2
    qrm_facto_flops_     = 3
    qrm_nnz_r_           = 4
    qrm_nnz_h_           = 5
    qrm_e_facto_mempeak_ = 6
    qrm_rd_num_          = 7
end

@cenum yn::UInt32 begin
    qrm_no_  = 0
    qrm_yes_ = 1
end

function error_handling(err :: Cint)
  status = "Unknown qr_mumps error: $err"
  err == 1  && (status = "The provided sparse matrix format is not supported.")
  err == 3  && (status = "qrm_spfct.rcntl is invalid.")
  err == 4  && (status = "Trying to allocate an already allocated allocatable or pointer.")
  err == 5  && (status = "Memory allocation problem.")
  err == 6  && (status = "Memory allocation problem.")
  err == 8  && (status = "Input column permutation not provided or invalid.")
  err == 9  && (status = "The requested ordering method is unknown.")
  err == 10 && (status = "Internal error: insufficient size for array.")
  err == 11 && (status = "Internal error: Error in lapack routine.")
  err == 12 && (status = "Internal error: out of memory.")
  err == 13 && (status = "The analysis must be done before the factorization.")
  err == 14 && (status = "The factorization must be done before the solve.")
  err == 15 && (status = "This type of norm is not implemented.")
  err == 16 && (status = "Requested ordering method not available (i.e., has not been installed).")
  err == 17 && (status = "Internal error: error from call to subroutine...")
  err == 18 && (status = "An error has occured in a call to COLAMD.")
  err == 19 && (status = "An error has occured in a call to SCOTCH.")
  err == 20 && (status = "An error has occured in a call to Metis.")
  err == 23 && (status = "Incorrect argument to qrm_set or qrm_get.")
  err == 25 && (status = "Internal error: problem opening file.")
  err == 27 && (status = "Incompatible values in qrm_spfct.icntl.")
  err == 28 && (status = "Incorrect value for qrm_mb, qrm_nb or qrm_ib.")
  err == 29 && (status = "Incorrect value for qrm_spmat.m, qrm_spmat.n or qrm_spmat.nz.")
  err == 30 && (status = "qrm_apply cannot be called if the Q matrix is discarded.")
  err == 31 && (status = "StarPU initialization error.")
  err == 32 && (status = "Matrix is rank deficient.")
  return status
end
