This repository contains the extended encoding function acceleration of Classic McEliece Key Encapsulation Mechanism (CM KEM). 
It is based on the 4rth Round submission of the CM KEM to the NIST Post-Quantum Standardization contest. 
Specifically we use the manually-vectorized across 64-bits VEC version of the CM NIST submission.

The eEnc's accelerator's HLS source code is located at the "mceliece{LEVEL}/kernel/" branches. 
The directories ending with (e.g. "encryption_kernel_x_inj") are the PK preprocessor design variants, with an external parallelization of p=x, while the ones lacking any suffix (e.g. "encryption_kernel_x") are the FIFO resizer design variants having a respective external parallelization.

The accelerator can be independently synthesized and exported as a standalone Xilinx IP, and used as such in the context of a Vivado project.
The designer just has to import the specific kernel of interest in the Vitis tool and move on with the desired workflow.
The accelerator can also be integrated with the Vitis HLS Xilinx suite, using the Zynq UltraScale+ (zcu102) development board, or any of the Alveo acceleration cards.
We provide specifically tailored testbenches (see mceliece{LEVEL}/testench/) along with inputs for the encoding subroutine (mceliece{LEVEL}/testench/io_values) that are taken directly from the official NIST submission of CM.
These testbenches are writen in C and are platform independent.


The different directories "mceliece{LEVEL}", where level=(340864, 460896, 6688128, 6960119, 8192128), are hardware implementations for the different security levels implemented in CM KEM.
The "common/" directory contains libraries needed by the baseline software version and also the Xilinx IP of AES256 module that we used.
Do not hesitate to contact us for any doubt or help to reproduce our results!
This project is a work in progress and we will keep on updating both this README file as well as the respective source code.






