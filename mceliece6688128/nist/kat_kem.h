//HLS OPENCL specific libraries
#include <CL/opencl.h>
#include "params.h"

extern cl_int err;

cl_mem buffer_v_in;
cl_mem buffer_key_in;
cl_mem buffer_pk_in[PARALLELIZATION_FACTOR]; 
cl_mem buffer_e_out; 
cl_mem buffer_s_out[PARALLELIZATION_FACTOR];

cl_mem buffer_list_enc_in[PARALLELIZATION_FACTOR+2];
cl_mem buffer_list_enc_out[PARALLELIZATION_FACTOR+1];

unsigned char *ptr_v_in;
unsigned char *ptr_key_in;

unsigned char *ptr_pk_in[PARALLELIZATION_FACTOR];
unsigned char *ptr_s_out[PARALLELIZATION_FACTOR];
unsigned char *ptr_e_out;
extern cl_command_queue commands;
extern cl_kernel encryption_kernel;
