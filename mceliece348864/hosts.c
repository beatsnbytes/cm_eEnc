#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
//CUSTOM OPENCL libraries
#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "nist/kat_kem.h"

//TODO can I give direct path if I include the path of common with -I
#include "../common/custom_util.h"
#include "params.h"

int times_event_migr_tokern=0;
double sum_list_event_migr_tokern[1];

// int times_event_enq=0;
// double sum_list_event_enq[1];

int times_event_migr_tohost=0;
double sum_list_event_migr_tohost[1];

double sum_encrypt;
double times_encrypt;


int encrypt_host(unsigned char *s, unsigned char *e){
	cl_event event_enq, event_migr_tokern, event_migr_tohost;

	// clEnqueueMigrateMemObjects(commands, 3, buffer_list_enc_in, 0, 0, NULL, &event_migr_tokern);
	
	// struct timeval start_encrypt, end_encrypt;
	// gettimeofday(&start_encrypt, NULL);

	// clEnqueueTask(commands, encryption_kernel, 1, &event_migr_tokern, &event_enq);

	// gettimeofday(&end_encrypt, NULL);
	// get_event_time(&start_encrypt, &end_encrypt, &sum_encrypt, &times_encrypt);

	// clEnqueueMigrateMemObjects(commands, 2, buffer_list_enc_out, CL_MIGRATE_MEM_OBJECT_HOST, 1, &event_enq, &event_migr_tohost);

	// clWaitForEvents(1, &event_migr_tohost);

	// memcpy(s, ptr_s_out, sizeof(unsigned char)*(SYND_BYTES));
	// memcpy(e, ptr_e_out, sizeof(unsigned char)*(SYS_N/8));
}
