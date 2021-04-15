# Control parameters

Control parameters define the behavior of **qr\_mumps** and can be set in two modes:

* **global mode**: in this mode it possible to either set generic parameters (e.g., the unit for output or error messages) or default parameter values (e.g., the ordering method to be used on the problem) that apply to all initialized **qrm\_spfct** factorizations.


* **problem mode**: these parameters control the behavior of **qr\_mumps** on a specific sparse factorization problem. Because the **qrm\_spfct\_init** routine sets the control parameters to their default values, these have to be modified after the sparse factorization object initialization.

All the control parameters can be set through the **qrm\_set** routine.

## Global parameters

* **qrm\_ncpu**: integer specifying the number of CPU cores to use for the subsequent **qr\_mumps** calls. It is an argument to the **qrm\_init** routine. Default is 1.


* **qrm\_ounit**: integer specifying the unit for output messages; if negative, output messages are suppressed. Default is 6 (stdout).


* **qrm\_eunit**: an integer specifying the unit for error messages; if negative, error messages are suppressed. Default is 0.

## Problem specific parameters

* **qrm\_ordering**: this parameter specifies what permutation to apply to the columns of the input matrix in order to reduce the fill-in and, consequently, the operation count of the factorization and solve phases. This parameter is used by **qr\_mumps** during the analysis phase and, therefore, has to be set before it starts. The following pre-defined values are accepted:
    * **qrm\_auto** (0) : the choice is automatically made by **qr\_mumps**. This is the default.
    * **qrm\_natural** (1) : no permutation is applied.
    * **qrm\_given** (2) : a column permutation is provided by the user through the **cperm\_in** attribute of a **qrm\_spfct** factorization.
    * **qrm\_colamd** (3) : the COLAMD software package (if installed) is used for computing the column permutation.
    * **qrm\_metis** (4) : the Metis software package (if installed) is used for computing the column permutation.
    * **qrm\_scotch** (5) : the SCOTCH software package (if installed) is used for computing the column permutation.


* **qrm\_keeph**: this parameter says whether the **Q** matrix should be kept for later use or discarded. This parameter is used by **qr\_mumps** during the factorization phase and, therefore, has to be set before it starts. Accepted value are:
    * **qrm\_yes** (1) : the **Q** matrix is kept. This is the default.
    * **qrm\_no** (0) : the **Q** matrix is discarded.


* **qrm\_mb** and **qrm\_nb**: These parameters define the block-size (rows and columns, respectively) for data partitioning and, thus, granularity of parallel tasks. Smaller values mean higher concurrence. This parameter, however, implicitly defines an upper bound for the granularity of call to BLAS and LAPACK routines (defined by the **qrm\_ib** parameter described below); therefore, excessively small values may result in poor performance. This parameter is used by **qr\_mumps** during the analysis and factorization phases and, therefore, has to be set before these start. The default value is 256 for both. Note that **qrm\_mb** has to be a multiple of **qrm\_nb**.


* **qrm\_ib**: this parameter defines the granularity of BLAS/LAPACK operations. Larger values mean better efficiency but imply more fill-in and thus more flops and memory consumption. The value of this parameter is upper-bounded by the **qrm\_nb** parameter described above. This parameter is used by **qr\_mumps** during the factorization phase and, therefore, has to be set before it starts. The default value is 32. It is strongly advised to choose, for this parameter, a submultiple of **qrm\_nb**.


* **qrm\_bh**: this parameter defines the type of algorithm for the communication-avoiding QR factorization of frontal matrices. Smaller values mean more concurrency but worse tasks efficiency; if lower or equal to zero the largest possible value is chosen for each front. Default value is -1.


* **qrm\_rhsnb**: in the case where multiple right-hand sides are passed to the **qrm\_apply** or the **qrm\_solve** routines, this parameter can be used to define a blocking of the right-hand sides. This parameter is used by **qr\_mumps** during the solve phase and, therefore, has to be set before it starts. By default, all the right-hand sides are treated in a single block.


* **qrm\_mem\_relax**: a value (≥ 1) that sets a relaxation parameter, with respect to the sequential peak, for the memory consumption in the factorization phase. If negative, the memory consumption is not bounded. Default value is −1.0.


* **qrm\_rd\_eps**: a value setting a threshold to estimate the rank of the problem. If > 0 the **qrm\_factorize** routine will count the number of diagonal coefficients of the **R** factor whose absolute value is smaller than the provided value. This number can be retrieved through the **qrm\_rd\_num** information parameter described in the next section.
