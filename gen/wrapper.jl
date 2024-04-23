# Script to parse qr_mumps headers and generate Julia wrappers.
using qr_mumps_jll
using Clang
using Clang.Generators
using JuliaFormatter

function main()

  cd(@__DIR__)
  include_dir = joinpath(qr_mumps_jll.artifact_dir, "include")
  headers = map(header -> joinpath(include_dir, header), ["sqrm_c.h", "dqrm_c.h", "cqrm_c.h", "zqrm_c.h", "qrm_common_c.h"])

  options = load_options(joinpath(@__DIR__, "qrmumps.toml"))
  options["general"]["output_file_path"] = joinpath("..", "src", "wrapper", "libqrmumps.jl")
  options["general"]["output_ignorelist"] = ["[sdcz]qrm_spmat_type_c", "[sdcz]qrm_spfct_type_c", "icntl", "rcntl", "ords", "gstats", "yn"]

  args = get_default_args()
  push!(args, "-I$include_dir")
  
  ctx = create_context(headers, args, options)
  build!(ctx)

  path = options["general"]["output_file_path"]

  code = read(path, String)
  for pattern in ("::", ",", ")")
    qrm_spmat_c = "qrm_spmat_c" * pattern
    spmat = "spmat" * pattern
    code = replace(code, qrm_spmat_c => spmat)

    qrm_spfct_c = "qrm_spfct_c" * pattern
    spfct = "spfct" * pattern
    code = replace(code, qrm_spfct_c => spfct)
  end

  for (version, type) in [("s", "Float32"), ("d", "Float64"), ("c", "ComplexF32"), ("z", "ComplexF64")]
    type_spmat_c = "Ptr{" * version * "qrm_spmat_type_c}"
    type_spmat_julia = "Ref{c_spmat{" * type * "}}"
    code = replace(code, type_spmat_c => type_spmat_julia)

    type_spfct_c = "Ptr{" * version * "qrm_spfct_type_c}"
    type_spfct_julia = "Ref{c_spfct{" * type * "}}"
    code = replace(code, type_spfct_c => type_spfct_julia)
  end
  write(path, code)

  format_file(path, YASStyle())
  return nothing
end

# If we want to use the file as a script with `julia wrapper.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
