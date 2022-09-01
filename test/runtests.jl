using Test, QRMumps, LinearAlgebra, SparseArrays, Printf, Random

@info("QRMUMPS_INSTALLATION: $(QRMumps.QRMUMPS_INSTALLATION)")
(QRMumps.QRMUMPS_INSTALLATION == "CUSTOM") && qrm_init()

include("test_qrm.jl")
