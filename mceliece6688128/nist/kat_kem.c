/*
   PQCgenKAT_kem.c
   Created by Bassham, Lawrence E (Fed) on 8/29/17.
   Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
   + mods from djb: see KATNOTES
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "rng.h"
#include "crypto_kem.h"
#include "operations.h"
#include "../common/custom_util.h"
#include "opencl_lib.h"
#include "hosts.h"

//Time measurement library
#include <sys/time.h>

//Includes the buffers, pointers sizes
//TODO move to common in source file?
#include "params.h"

//HLS OPENCL specific libraries
#include <CL/opencl.h>
#include <CL/cl_ext.h>


#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_CRYPTO_FAILURE  -4

void	fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L);

unsigned char entropy_input[48];
unsigned char seed[KATNUM][48];

double sum_keygen, sum_enc, sum_dec;
int times_keygen, times_enc, times_dec;

//GLOBAL variables needed by host functions
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

cl_command_queue commands;
cl_kernel encryption_kernel;

int main(int argc, char* argv[])
{
//Code specific for openCL and HLS functionality

    char *xclbin = argv[1];
    cl_context context;
    cl_int err;


    cl_program program;
    int i;


    //Initialization code
    platform_init(&program, &context, &commands, xclbin);

    //Kernel creation code

                                                    //ENCRYPTION//

    #ifdef ENCRYPTION_KERNEL
	{

        const char *kernel_name = "encryption_kernel";
        encryption_kernel = clCreateKernel(program, kernel_name, NULL);


		#ifdef ALVEO
        {
    		posix_memalign((unsigned char **)&ptr_v_in, 4096, sizeof(unsigned char)*16);
    		posix_memalign((unsigned char **)&ptr_key_in, 4096, sizeof(unsigned char)*32);
    		posix_memalign((unsigned char **)&ptr_e_out, 4096, sizeof(unsigned char)*(SYS_N/8));

    		for(i=0; i<PARALLELIZATION_FACTOR; i++){
    			posix_memalign((unsigned char **)&ptr_pk_in[i], 4096, sizeof(unsigned char)*(PK_NROWS*(PK_ROW_BYTES+12))/PARALLELIZATION_FACTOR);
    			posix_memalign((unsigned char **)&ptr_s_out[i], 4096, sizeof(unsigned char)*(SYND_BYTES)/PARALLELIZATION_FACTOR);
    		}

    	    buffer_v_in = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(unsigned char)*16, ptr_v_in, NULL);
    	    buffer_key_in = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(unsigned char)*32, ptr_key_in, NULL);

    	    for(i=0; i<PARALLELIZATION_FACTOR; i++){
    	        //PK
    	        buffer_pk_in[i] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(unsigned char)*(PK_NROWS*(PK_ROW_BYTES+12))/PARALLELIZATION_FACTOR, ptr_pk_in[i], NULL);
    	        //S
    	        buffer_s_out[i] = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, sizeof(unsigned char)*(SYND_BYTES)/PARALLELIZATION_FACTOR, ptr_s_out[i], NULL);
    	    }

    	    buffer_e_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, sizeof(unsigned char)*(SYS_N/8), ptr_e_out, NULL);


        }
		#endif


		#ifdef ZCU
        {
			buffer_v_in = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*16, NULL, &err);
			buffer_key_in = clCreateBuffer(context,  CL_MEM_READ_ONLY, sizeof(unsigned char)*32, NULL, &err);

			for(i=0; i<PARALLELIZATION_FACTOR; i++){
				//PK
				buffer_pk_in[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*(PK_NROWS*(PK_ROW_BYTES+12))/PARALLELIZATION_FACTOR, NULL, &err);
				//S
				buffer_s_out[i] = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char)*(SYND_BYTES)/PARALLELIZATION_FACTOR, NULL, &err);
			}

			buffer_e_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char)*(SYS_N/8), NULL, &err);
			ptr_v_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_v_in, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*16, 0, NULL, NULL, &err);
			ptr_key_in = (unsigned char *) clEnqueueMapBuffer(commands, buffer_key_in, true, CL_MAP_WRITE, 0, sizeof(unsigned char)*32, 0, NULL, NULL, &err);

			for(i=0; i<PARALLELIZATION_FACTOR; i++){
				//PK
				ptr_pk_in[i] = (unsigned char *) clEnqueueMapBuffer(commands, buffer_pk_in[i], true, CL_MAP_WRITE, 0, sizeof(unsigned char)*(PK_NROWS*(PK_ROW_BYTES+12))/PARALLELIZATION_FACTOR, 0, NULL, NULL, &err);
				//S_BUFF
				ptr_s_out[i] = (unsigned char *) clEnqueueMapBuffer(commands, buffer_s_out[i], true, CL_MAP_READ, 0, sizeof(unsigned char)*(SYND_BYTES)/PARALLELIZATION_FACTOR, 0, NULL, NULL, &err);
			}

			ptr_e_out = (unsigned char *) clEnqueueMapBuffer(commands, buffer_e_out, true, CL_MAP_READ, 0, sizeof(unsigned char)*(SYS_N/8), 0, NULL, NULL, &err);

        }
		#endif

        for(i=0; i<PARALLELIZATION_FACTOR; i++){    
            err = clSetKernelArg(encryption_kernel, i, sizeof(cl_mem), &buffer_s_out[i]);
            err = clSetKernelArg(encryption_kernel, i+(PARALLELIZATION_FACTOR+1), sizeof(cl_mem), &buffer_pk_in[i]);        
        }

        err = clSetKernelArg(encryption_kernel, PARALLELIZATION_FACTOR, sizeof(cl_mem), &buffer_e_out);

        err = clSetKernelArg(encryption_kernel, PARALLELIZATION_FACTOR*2+1, sizeof(cl_mem), &buffer_v_in);
        err = clSetKernelArg(encryption_kernel, PARALLELIZATION_FACTOR*2+2, sizeof(cl_mem), &buffer_key_in);

        for(i=0; i<PARALLELIZATION_FACTOR; i++){    
            buffer_list_enc_in[i] = buffer_pk_in[i]; 
        }
        buffer_list_enc_in[PARALLELIZATION_FACTOR] = buffer_v_in; 
        buffer_list_enc_in[PARALLELIZATION_FACTOR+1] = buffer_key_in;     


        for(i=0; i<PARALLELIZATION_FACTOR; i++){    
            buffer_list_enc_out[i] = buffer_s_out[i];  
        }
        buffer_list_enc_out[PARALLELIZATION_FACTOR] = buffer_e_out;

    }
    #endif



    //Original kat_kem.c source code of CME NIST submission

    FILE                *fp_req, *fp_rsp;
    int                 ret_val;
    unsigned char *ct = 0;
    unsigned char *ss = 0;
    unsigned char *ss1 = 0;
    unsigned char *pk = 0;
    unsigned char *sk = 0;

    struct timeval start_keygen, end_keygen, start_enc, end_enc, start_dec, end_dec;

    for (i=0; i<48; i++)
        entropy_input[i] = i;
    randombytes_init(entropy_input, NULL, 256);

    for (i=0; i<KATNUM; i++)
        randombytes(seed[i], 48);

    fp_req = fopen("kat_kem.req", "w");
    if (!fp_req)
        return KAT_FILE_OPEN_ERROR;

    for (i=0; i<KATNUM; i++) {
        fprintf(fp_req, "count = %d\n", i);
        fprintBstr(fp_req, "seed = ", seed[i], 48);
        fprintf(fp_req, "pk =\n");
        fprintf(fp_req, "sk =\n");
        fprintf(fp_req, "ct =\n");
        fprintf(fp_req, "ss =\n\n");
    }

    fp_rsp = fopen("kat_kem.rsp", "w");
    if (!fp_rsp)
        return KAT_FILE_OPEN_ERROR;

    fprintf(fp_rsp, "# kem/%s\n\n", crypto_kem_PRIMITIVE);

    for (i=0; i<KATNUM; i++) {
        if (!ct) ct = malloc(crypto_kem_CIPHERTEXTBYTES);
        if (!ct) abort();
        if (!ss) ss = malloc(crypto_kem_BYTES);
        if (!ss) abort();
        if (!ss1) ss1 = malloc(crypto_kem_BYTES);
        if (!ss1) abort();
        if (!pk) pk = malloc(crypto_kem_PUBLICKEYBYTES);
        if (!pk) abort();
        if (!sk) sk = malloc(crypto_kem_SECRETKEYBYTES);
        if (!sk) abort();

        randombytes_init(seed[i], NULL, 256);

        fprintf(fp_rsp, "count = %d\n", i);
        fprintBstr(fp_rsp, "seed = ", seed[i], 48);
       
        gettimeofday(&start_keygen, NULL);

        ret_val = crypto_kem_keypair(pk, sk);

        gettimeofday(&end_keygen, NULL);
        get_event_time(&start_keygen, &end_keygen, &sum_keygen, &times_keygen);


        if (ret_val != 0) {
            fprintf(stderr, "crypto_kem_keypair returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        fprintBstr(fp_rsp, "pk = ", pk, crypto_kem_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk = ", sk, crypto_kem_SECRETKEYBYTES);        

        gettimeofday(&start_enc, NULL);

        ret_val = crypto_kem_enc(ct, ss, pk);

        gettimeofday(&end_enc, NULL);
        get_event_time(&start_enc, &end_enc, &sum_enc, &times_enc);


        if (ret_val != 0) {
            fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        fprintBstr(fp_rsp, "ct = ", ct, crypto_kem_CIPHERTEXTBYTES);
        fprintBstr(fp_rsp, "ss = ", ss, crypto_kem_BYTES);
        
        fprintf(fp_rsp, "\n");

        gettimeofday(&start_dec, NULL);

        ret_val =  crypto_kem_dec(ss1, ct, sk);

        gettimeofday(&end_dec, NULL);
        get_event_time(&start_dec, &end_dec, &sum_dec, &times_dec);

        if (ret_val != 0) {
            fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        
        if ( memcmp(ss, ss1, crypto_kem_BYTES) ) {
            fprintf(stderr, "crypto_kem_dec returned bad 'ss' value\n");
            return KAT_CRYPTO_FAILURE;
        }
    }
	
    printf("\n\t**********TIMING RESULTS**********\t\n");    
	printf("Key Generation Part ");
	print_event_execution_time(&sum_keygen, &times_keygen);
	printf("Encapsulation Part ");
	print_event_execution_time(&sum_enc, &times_enc);
	printf("Decapsulation Part ");
	print_event_execution_time(&sum_dec, &times_dec);
	printf("Encryption Part ");
	print_event_execution_time(&sum_encrypt, &times_encrypt);
	printf("Decryption Part ");
	print_event_execution_time(&sum_decrypt, &times_decrypt);

	// printf("Encrypt Function");
	// print_event_execution_time(&sum_encryption_func, &times_encryption_func);
    // printf("Data to kernel\n");
    // print_kernel_execution_time(sum_list_event_migr_tokern, &times_event_migr_tokern, 1);
    printf("Kernel Enqueue\n");
    print_kernel_execution_time(sum_list_event_enq, &times_event_enq, 1);
    // printf("Data to host\n");
    // print_kernel_execution_time(sum_list_event_migr_tohost, &times_event_migr_tohost, 1);

    return KAT_SUCCESS;
}

void
fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L)
{
	unsigned long long i;

	fprintf(fp, "%s", S);

	for ( i=0; i<L; i++ )
		fprintf(fp, "%02X", A[i]);

	if ( L == 0 )
		fprintf(fp, "00");

	fprintf(fp, "\n");
}
