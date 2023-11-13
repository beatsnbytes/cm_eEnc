//HLS OPENCL specific libraries
#include <CL/opencl.h>
#include "params.h"

extern cl_int err;
extern cl_context context;

extern cl_kernel decryption_kernel;
extern cl_kernel encryption_kernel;

extern cl_mem buffer_v_in;
extern cl_mem buffer_key_in;
extern cl_mem buffer_pk_in[PARALLELIZATION_FACTOR]; 
extern cl_mem buffer_e_out; 
extern cl_mem buffer_s_out[PARALLELIZATION_FACTOR];

extern cl_mem buffer_list_enc_in[PARALLELIZATION_FACTOR+2];
extern cl_mem buffer_list_enc_out[PARALLELIZATION_FACTOR+1];

extern unsigned char *ptr_v_in;
extern unsigned char *ptr_key_in;

extern unsigned char *ptr_pk_in[PARALLELIZATION_FACTOR];
extern unsigned char *ptr_s_out[PARALLELIZATION_FACTOR];
extern unsigned char *ptr_e_out;
extern cl_command_queue commands;

extern cl_mem buffer_sk_irr_in;
extern cl_mem buffer_sk_rev_in;
extern cl_mem buffer_sk_fwd_in;

extern cl_mem buffer_c_in;
extern cl_mem buffer_check_out;
extern cl_mem buffer_e_out_dec;

extern cl_mem buffer_list_dec_in[4];
extern cl_mem buffer_list_dec_out[2];

extern unsigned char *ptr_sk_irr_in;
extern unsigned char *ptr_sk_rev_in;
extern unsigned char *ptr_sk_fwd_in;

extern unsigned char *ptr_c_in;
extern int *ptr_check_out;
extern unsigned char *ptr_e_out_dec;

