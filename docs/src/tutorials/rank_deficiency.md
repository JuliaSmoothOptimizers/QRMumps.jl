```@example rd1
using QRMumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
jcn = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
val = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]

A = sparse(irn, jcn, val, 4, 3)

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# The control parameter `qrm_rd_eps` is a threshold to estimate the rank of the problem.
# If qrm_rd_eps > 0 the qrm_factorize routine will count the number of diagonal
# coefficients of the R factor whose absolute value is smaller than the provided value.
qrm_set(spfct, "qrm_rd_eps", 1e-12)

# Perform the analysis and factorization phases
qrm_analyse!(spmat, spfct)
qrm_factorize!(spmat, spfct)

qrm_get(spfct, "qrm_rd_eps")

# The information parameter `qrm_rd_num` contains the number of diagonal coefficients
# of the R factor whose absolute value is lower than `qrm_rd_eps`.
rank_deficiency = qrm_get(spfct, "qrm_rd_num")
```

```@example rd2
using QRMumps, LinearAlgebra, SparseArrays, Printf

irn = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
jcn = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]
val = [1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 6.0, 15.0, 7.0, 8.0, 9.0, 24.0]

A = sparse(irn, jcn, val, 3, 4)

qrm_init()

spmat = qrm_spmat_init(A)
spfct = qrm_spfct_init(spmat)

# The control parameter `qrm_rd_eps` is a threshold to estimate the rank of the problem.
# If qrm_rd_eps > 0 the qrm_factorize routine will count the number of diagonal
# coefficients of the R factor whose absolute value is smaller than the provided value.
qrm_set(spfct, "qrm_rd_eps", 1e-12)

# Perform the analysis and factorization phases
qrm_analyse!(spmat, spfct, transp='t')
qrm_factorize!(spmat, spfct, transp='t')

# The information parameter `qrm_rd_num` contains the number of diagonal coefficients
# of the R factor whose absolute value is lower than `qrm_rd_eps`.
rank_deficiency = qrm_get(spfct, "qrm_rd_num")
```
