#include <gmp.h>
#define __gmp_const const

#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include<string.h>
#include"ap_int.h"

#include "scaling_kernel_testbench.h"


int main(){
    //Establish an initial return value. 0 = success
    int ret;
    int iter;
    int i, j, k;

    //Helper variables
    sk_irr_whole sk;
    vec recv;


    //Variables that will store the IO data
    //Inputs
    sk_irr_whole *sk_reference_input = (sk_irr_whole *) malloc(sizeof(sk)*(128));
    vec *recv_reference_input = (vec *) malloc(sizeof(vec)*(64));


    //Outputs
    vec *inv_golden_output = (vec *) malloc(sizeof(vec)* 64*(GFBITS));
    vec *out_golden_output = (vec *) malloc(sizeof(vec)* 64*(GFBITS));
    vec *inv_kernel_output = (vec *) malloc(sizeof(vec)* 64*(GFBITS));
    vec *out_kernel_output = (vec *) malloc(sizeof(vec)* 64*(GFBITS));


    // Call any preliminary functions required to prepare input for the test.
    FILE *sk_fptr;
    FILE *recv_fptr;
    FILE *inv_golden_fptr;
    FILE *out_golden_fptr;

	ret=0;

	sk_fptr = fopen("./io_values/scaling_1/sk_in.dat","r");
	recv_fptr = fopen("./io_values/scaling_1/recv_in.dat","r");
	inv_golden_fptr = fopen("./io_values/scaling_1/inv_golden_out.dat","r");
	out_golden_fptr = fopen("./io_values/scaling_1/out_golden_out.dat","r");

	//Read from the pk values file into the pk matrix
	for(i=0; i<128; i++){
		fscanf(sk_fptr, "%d\n", &sk);
		sk_reference_input[i] |= sk<<(i*8);
	}
	fclose(sk_fptr);



	for(i=0; i<64; i++){
		fscanf(recv_fptr, "%lu\n", &recv);
		recv_reference_input[i] = recv;
	}
	fclose(recv_fptr);


	scaling_kernel(out_kernel_output, inv_kernel_output, sk_reference_input, recv_reference_input);

	// Compare the results of the function against expected results

	for(i=0; i<64; i++){
		for(j=0; j<GFBITS; j++){
		if(*(out_kernel_output+(i*64)+j)!=*(out_golden_output+(i*64)+j)){
			printf("ERROR! Expected out[%d]=%lX but got %lX\n", i, *(out_golden_output+(i*64)+j), *(out_kernel_output+(i*64)+j));
			ret=1;
		}
	}

		for(i=0; i<64; i++){
			for(j=0; j<GFBITS; j++){
			if(*(inv_kernel_output+(i*64)+j)!=*(inv_golden_output+(i*64)+j)){
				printf("ERROR! Expected inv[%d]=%lX but got %lX\n", i, *(inv_golden_output+(i*64)+j), *(inv_kernel_output+(i*64)+j));
				ret=1;
			}
		}


	if (ret != 0) {
			printf("Iter %d Test failed  !!!\n", iter);
			ret=1;
	} else {
			printf("Iter %d Test passed !\n", iter);
	}

    
    return ret;
}

