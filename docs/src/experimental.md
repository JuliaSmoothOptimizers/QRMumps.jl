* **qrm\_ngpu**: integer specifying the number of GPUs to use for the subsequent **qr\_mumps** calls. It is an argument to the **qrm\_init** routine. Default is 0.

* **qrm\_pinth**: an integer value to control memory pinning when GPUs are used: all frontal matrices whose size (min(rows,cols)) is bigger than this value will be pinned.

Note that it is possible to use multiple streams per GPU; this can be controlled through the StarPU **STARPU\_NWORKER\_PER\_CUDA** environment variable.

## **GPU streams**

* When GPUs are used, it can be helpful (and it usually is) to use multiple streams per GPU to allow a single GPU to execute multiple tasks concurrently. Using multiple GPU streams is especially beneficial to achieve high GPU occupancy when a relatively small block size _mb_ is chosen to prevent CPU starvation. This can be controlled through the **STARPU\_NWORKER\_PER\_CUDA** StarPU environment variable. By default one stream is active per GPU device and higher performance can be commonly achieved with values of 2 up to 20.
