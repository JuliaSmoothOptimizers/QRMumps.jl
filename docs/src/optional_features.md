# Optional features

The following features of the **qr\_mumps** software are currently unavailable in the Julia interface if the package is installed through Yggdrasil. If a custom **qr\_mumps** install is used which has support for the [StarPU](https://starpu.gitlabpages.inria.fr/) runtime, these features can be accessed, although they have not been thoroughly tested within Julia.

## **Multithreading**

* **qr\_mumps** is a parallel, multithreaded software based on the StarPU runtime system and it currently supports multicore or, more generally, shared memory multiprocessor computers. **qr\_mumps** does not run on distributed memory (e.g. clusters) parallel computers. Parallelism is achieved through a decomposition of the workload into fine-grained computational tasks which basically correspond to the execution of a BLAS or LAPACK operation on a blocks. It is strongly recommended to use sequential BLAS and LAPACK libraries and let **qr\_mumps** have full control of the parallelism. The granularity of the tasks is controlled by the **qrm\_mb** and **qrm\_nb** parameters which set the block size for partitioning internal data. Smaller values mean more parallelism; however, because this blocking factor is an upper-bound for the granularity of operations (or, more precisely for the granularity of calls to BLAS and LAPACK routines), it is recommended to choose reasonably large values in order to achieve high efficiency.

## **GPU acceleration**

* **qr\_mumps** can leverage the computing power of Nvidia GPU, commonly available on modern super-computing systems, to accelerate the solution of linear systems, especially large size ones. The use of GPUs is achieved through the StarPU runtime.
