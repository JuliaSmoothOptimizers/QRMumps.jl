# Introduction

This package provides a Julia interface to [qr_mumps](http://qr_mumps.gitlab.io/), a software for solving sparse linear systems on multicore computers.
qr\_mumps implements a direct solution method based on the QR or Cholesky factorization of the input matrix. 
Therefore, it is suited to solving sparse least-squares problems, to computing the minimum-norm solution of sparse, underdetermined problems and to solving symmetric, positive-definite sparse linear systems. It can obviously be used for solving square unsymmetric problems in which case the stability provided by the use of orthogonal transformations comes at the cost of a higher operation count with respect to solvers based on, e.g., the LU factorization such as [MUMPS](https://mumps-solver.org/index.php). It supports real and complex, single or double precision arithmetic.

 As in all the sparse, direct solvers, the solution is achieved in three distinct phases:

* **Analysis**
    * In this phase an analysis of the structural properties of the input matrix is performed in preparation for the numerical factorization phase. This includes computing a column permutation which reduces the amount of fill-in coefficients (i.e., nonzeroes introduced by the factorization). This step does not perform any floating-point operation and is, thus, commonly much faster than the factorization and solve (depending on the number of right-hand sides) phases.

* **Factorization**
    * At this step, the actual QR or Cholesky factorization is computed. This step is the most computationally intense and, therefore, the most time consuming.

* **Solution**
    * Once the factorization is done, the factors can be used to compute the solution of the problem through two operations:
    * **Solve** : this operation computes the solution of the triangular system Rx=b or Rᵀx=b;
    * **Apply** : this operation applies the Q orthogonal matrix to a vector, i.e., y=Qx or y=Qᵀx. 

These three steps have to be done in order but each of them can be performed multiple times. If, for example, the problem has to be solved against multiple right-hand sides (not all available at once), the analysis and factorization can be done only once while the solution is repeated for each right-hand side. By the same token, if the coefficients of a matrix are updated but not its structure, the analysis can be performed only once for multiple factorization and solution steps.

**qr\_mumps** is based on the multifrontal factorization method. This method was first introduced by Duff and Reid as a method for the factorization of sparse, symmetric linear systems and, since then, has been the object of numerous studies and the method of choice for several, high-performance, software packages such as MUMPS and UMFPACK.
