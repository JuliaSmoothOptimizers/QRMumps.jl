# Information parameters

Information parameters return information about the behavior of **qr\_mumps** and can be either global or problem specific.
All the information parameters can be gotten through the **qrm\_get** routine; problem specific control parameters can also be retrieved by manually reading the **gstats** attribute of a **qrm\_spfct** factorization.
The **qrm\_get** routine can also be used to retrieve the values of all the control parameters described in the previous section with the obvious usage.
The type of all information parameters is **Int64**.

## Global parameters

* **qrm\_max\_mem**: this parameter returns the maximum amount of memory allocated by **qr\_mumps** during its execution.


* **qrm\_tot\_mem**: this parameter returns the total amount of memory allocated by **qr\_mumps** at the moment when the **qrm\_get** routine is called.

## Problem specific parameters

* **qrm\_e\_facto\_flops**: this parameter returns an estimate, computed during the analysis phase, of the number of floating point operations performed during the factorization phase. This value is only available after the **qrm\_analyse** routine is executed.


* **qrm\_e\_nnz\_r**: this parameter returns an estimate, computed during the analysis phase, of the number of nonzero coefficients in the **R** factor. This value is only available after the **qrm\_analyse** routine is executed.


* **qrm\_e\_nnz\_h**: this parameter returns an estimate, computed during the analysis phase, of the number of nonzero coefficients in the **Q** matrix. This value is only available after the **qrm\_analyse** routine is executed.


* **qrm\_facto\_flops**: this parameter returns the number of floating point operations performed during the factorization phase. This value is only available after the **qrm\_analyse** routine is executed.


* **qrm\_nnz\_r**: this parameter returns the actual number of the nonzero coefficients in the **R** factor after the factorization is done. This value is only available after **the qrm\_factorize** routine is executed.


* **qrm\_nnz\_h**: this parameter returns the actual number of the nonzero coefficients in the **Q** matrix after the factorization is done. This value is only available after the **qrm\_factorize** routine is executed.


* **qrm\_e\_facto\_mempeak**: this parameter returns an estimate of the peak memory consumption of the factorization operation.


* **qrm\_rd\_num**: this information parameter returns the number of diagonal coefficients of the **R** factor whose absolute value is lower than **qrm\_rd\_eps** if this control parameter was set to a value greater than 0.
