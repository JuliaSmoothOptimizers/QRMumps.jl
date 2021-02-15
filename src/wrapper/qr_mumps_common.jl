mutable struct qrm_spmat{T}
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

mutable struct qrm_spfct{T}
  cperm_in :: Ptr{Cint}
  icntl    :: NTuple{20, Cint}
  rcntl    :: NTuple{10, Cfloat}
  gstats   :: NTuple{10, Clong}
  h        :: Ptr{Cvoid}

  function qrm_spfct{T}() where T
    spfct = new(C_NULL, ntuple(x -> Cint(0), 20), ntuple(x -> Cfloat(0), 10), ntuple(x -> Clong(0), 10), C_NULL)
    return spfct
  end
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
