module qr_mumps

export QrmType,
       qrm_analyse, qrm_analyze, qrm_factorize,
       qrm_apply, qrm_solve, qrm_solve!,
       qrm_least_squares, qrm_least_squares!,
       qrm_pseti

if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
  include("../deps/deps.jl")
else
  error("QR_MUMPS library not properly installed. Please run Pkg.build(\"qr_mumps\")")
end

"Exception type raised in case of error."
type QRMException <: Exception
  msg :: ASCIIString
end


type QrmType_Private
  irn :: Ptr{Cint}
  jcn :: Ptr{Cint}
  val :: Ptr{Cdouble}
  m   :: Cint
  n   :: Cint
  nz  :: Cint
  cperm_in :: Ptr{Cint}
  icntl :: Ptr{Cint}
  rcntl :: Ptr{Cdouble}
  gstats :: Ptr{Clong}
  h   :: Cint

  function QrmType_Private(irn :: Vector{Int}, jcn :: Vector{Int}, val :: Vector{Float64},
                           m :: Int, n :: Int, nz :: Int)
    qrm_priv = new(C_NULL, C_NULL, C_NULL, 0, 0, 0, C_NULL, C_NULL, C_NULL, C_NULL, 0)
    ccall((:dqrm_spmat_init_c, libdqrm), Void, (Ptr{QrmType_Private},), &qrm_priv)
    qrm_priv.irn = pointer(irn)
    qrm_priv.jcn = pointer(jcn)
    qrm_priv.val = pointer(val)
    qrm_priv.m = m
    qrm_priv.n = n
    qrm_priv.nz = nz
    # cperm_in can remain NULL for now.
    # icntl, rcntl and gstats are allocated inside QR_MUMPS.
    finalizer(qrm_priv, qrm_private_finalize)
    return qrm_priv
  end
end

type QrmType

  __qrm :: QrmType_Private
  irn :: Vector{Cint}  # to ensure references persist.
  jcn :: Vector{Cint}
  val :: Vector{Cdouble}

  function QrmType(irn :: Vector{Int}, jcn :: Vector{Int}, val :: Vector{Float64},
                   m :: Int, n :: Int, nz :: Int)

    qrm_priv = QrmType_Private(irn, jcn, val, m, n, nz)
    qrm = new(qrm_priv, irn, jcn, val)
    finalizer(qrm, qrm_finalize)
    return qrm
  end
end


# Constructor.
function QrmType(A :: SparseMatrixCSC{Float64,Int})
  m, n = size(A)
  nz = nnz(A)
  irn, jcn, val = findnz(A)
  return QrmType(irn, jcn, val, m, n, nz)
end

QrmType(A :: Array{Float64,2}) = QrmType(sparse(A))


# Destructors.
function qrm_private_finalize(qrm_priv :: QrmType_Private)
  # Don't print in the finalizer, or else:
  # jl_uv_writecb() ERROR: bad file descriptor EBADF
  ccall((:dqrm_spmat_destroy_c, libdqrm), Void, (Ptr{QrmType_Private},), &qrm_priv)
  qrm_priv.irn = C_NULL
  qrm_priv.jcn = C_NULL
  qrm_priv.val = C_NULL
  qrm_priv.cperm_in = C_NULL
  qrm_priv.icntl = C_NULL
  qrm_priv.rcntl = C_NULL
  qrm_priv.gstats = C_NULL
end

function qrm_finalize(qrm :: QrmType)
  qrm_private_finalize(qrm.__qrm)
end


# import Base.show, Base.print
# function show(io :: IO, qrm :: QrmType)
#   print(io, "QRM: m=$(qrm.m), n=$(qrm.n), nz=$(qrm.nz)");
# end
#
# function print(io :: IO, qrm :: QrmType)
#   print(io, "QRM: m=$(qrm.m), n=$(qrm.n), nz=$(qrm.nz)");
# end


function qrm_pseti(qrm :: QrmType, param :: ASCIIString, val :: Int)
  ccall((:dqrm_pseti_c, libdqrm), Void, (Ptr{QrmType}, Ptr{Uint8}, Cint), &qrm, param, val)
end


function qrm_analyse(qrm :: QrmType; transp :: Char='n')
  ccall((:dqrm_analyse_c, libdqrm), Void, (Ptr{QrmType}, Cchar), &qrm, transp)
end

# Zed's not dead!
qrm_analyze = qrm_analyse


function qrm_factorize(qrm :: QrmType; transp :: Char='n')
  ccall((:dqrm_factorize_c, libdqrm), Void, (Ptr{QrmType}, Cchar), &qrm, transp)
end


function qrm_apply(qrm :: QrmType, b :: Vector{Float64}; transp :: Char='t')
  # transp = 't' means we apply Q' to b.
  # This is consistent with transp = 'n' in the analysis and factorization.
  # It's qr_mumps that is inconsistent.
  ccall((:dqrm_apply_c, libdqrm), Void,
        (Ptr{QrmType}, Cchar,  Ptr{Cdouble}, Cint),
         &qrm,         transp, b,            1)
end


function qrm_solve(qrm :: QrmType, b :: Vector{Float64}; transp :: Char='n')
  x = zeros(Cdouble, transp == 'n' ? qrm.n : qrm.m)
  qrm_solve!(qrm, b, x, transp=transp)
  return x
end

function qrm_solve!(qrm :: QrmType, b :: Vector{Float64}, x :: Vector{Float64}; transp :: Char='n')
  length(x) ≥ (transp == 'n' ? qrm.n : qrm.m) || error("qr_mumps: size of x insufficient")
  ccall((:dqrm_solve_c, libdqrm), Void,
        (Ptr{QrmType}, Cchar,  Ptr{Cdouble}, Ptr{Cdouble}, Cint),
         &qrm,         transp, b,            x,            1)
end


function qrm_least_squares(qrm :: QrmType, b :: Vector{Float64})
  x = zeros(Cdouble, qrm.n)
  qrm_least_squares!(qrm, b, x)
  return x
end

function qrm_least_squares!(qrm :: QrmType, b :: Vector{Float64}, x :: Vector{Float64})
  length(x) ≥ qrm.n || error("qr_mumps: size of x insufficient")
  ccall((:dqrm_least_squares_c, libdqrm), Void,
        (Ptr{QrmType}, Ptr{Cdouble}, Ptr{Cdouble}, Cint),
         &qrm,         b,            x,            1)
end

function qrm_least_squares(A :: SparseMatrixCSC{Float64,Int}, b :: Vector{Float64})
  qrm = QrmType(A)
  return qrm_least_squares(qrm, b)
end


end # module
