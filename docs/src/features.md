# Features

## **Types of problems**

* **qr\_mumps** can handle unsymmetric and symmetric, positive definite problems. In the first case it will use a QR factorization whereas, in the second, it will use a Cholesky factorization. In order to choose one or the other method, **qr\_mumps** must be informed about the type of the problem through the **sym** argument of the **qrm\_spmat\_init** function: **false** means that the problem is unsymmetric and **true** means symmetric, positive definite. Note that in the second case, only half of the matrix must be provided, i.e., if the coefficient (i, j) is provided (j, i) must not be given.

## **Memory consumption control**

* **qr\_mumps** allows for controlling the amount of memory used in the parallel factorization stage. In the multifrontal method, the memory consumption varies greatly throughout the sequential factorization reaching a maximum value which is referred to as the sequential peak (_sp_). Parallelism can considerably increase this peak because, in order to feed the working threads, more data is allocated at the same time which results in higher concurrency. In **qr\_mumps** it is possible to bound the memory consumption of the factorization phase through the **qrm\_mem\_relax** parameter. If this parameter is set to a real value x ≥ 1, the memory consumption will be bounded by _x × sp_. Clearly, the tighter is this upper bound, the slower the factorization will proceed. Note that _sp_ only includes the memory consumed by the factorization operation; moreover, although in practice it is possible to precisely pre-compute this value in the analysis phase, this may be expensive and thus **qrm\_analyse** only computes a (hopefully) slight overestimation. The value of _sp_ is available upon completion of the analysis phase through the **qrm\_e\_facto\_mempeak** information parameter.

## **Fill-reducing permutations**

* **qr\_mumps** supports multiple methods for reducing the factorization fill-in through matrix column permutations. The choice is controlled through the **qrm\_ordering** control parameter. Nested-dissection based methods are available through the packages Metis and SCOTCH packages as well as average minimum degree through the COLAMD one. Nested-dissection based methods usually lead to lower fill-in which ultimately results in faster and less memory consuming factorization. COLAMD, instead, typically leads to a faster execution of the analysis phase although is not as effective in reducing the fill-in which may result in a slower and more memory consuming factorization. Because the overall execution time is commonly dominated by the factorization, nested-dissection methods are usually more effective especially for large size problems. **qr\_mumps** also allows the user to provide their own permutation.

# Optional features
The following features of the **qr\_mumps** software are currently unavailable in the Julia interface if the package is installed through Yggdrasil. If a custom **qr\_mumps** install is used which has support for the [StarPU](https://starpu.gitlabpages.inria.fr/) runtime, these features can be accessed, although they have not been thoroughly tested within Julia.

## **Multithreading**

* **qr\_mumps** is a parallel, multithreaded software based on the StarPU runtime system and it currently supports multicore or, more generally, shared memory multiprocessor computers. **qr\_mumps** does not run on distributed memory (e.g. clusters) parallel computers. Parallelism is achieved through a decomposition of the workload into fine-grained computational tasks which basically correspond to the execution of a BLAS or LAPACK operation on a blocks. It is strongly recommended to use sequential BLAS and LAPACK libraries and let **qr\_mumps** have full control of the parallelism. The granularity of the tasks is controlled by the **qrm\_mb** and **qrm\_nb** parameters which set the block size for partitioning internal data. Smaller values mean more parallelism; however, because this blocking factor is an upper-bound for the granularity of operations (or, more precisely for the granularity of calls to BLAS and LAPACK routines), it is recommended to choose reasonably large values in order to achieve high efficiency.

## **GPU acceleration**

* **qr\_mumps** can leverage the computing power of Nvidia GPU, commonly available on modern super-computing systems, to accelerate the solution of linear systems, especially large size ones. The use of GPUs is achieved through the StarPU runtime.
<!-- Note that it is possible to use multiple streams per GPU; this can be controlled through the StarPU **STARPU\_NWORKER\_PER\_CUDA** environment variable. -->

