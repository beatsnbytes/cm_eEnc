#include <gmp.h>
#define __gmp_const const

#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include<string.h>
#include"ap_int.h"

#include "benes_kernel_testbench.h"


int main () { 
    //Establish an initial return value. 0 = success
    int ret;
    int iter;
    int i, j, k;

    //Helper variables
    unsigned char sk;
    unsigned char s;
    vec r;


    //Variables that will store the IO data
    //Inputs
    unsigned char *sk_reference_input = (unsigned char *) malloc(sizeof(sk)*(6492));
    unsigned char *s_reference_input = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES));


    //Outputs
    vec *r_golden_output = (vec *) malloc(sizeof(vec)* 64);
    vec *r_kernel_output = (vec *) malloc(sizeof(vec)* 64);


    // Call any preliminary functions required to prepare input for the test.
    FILE *sk_fptr;
    FILE *s_fptr;
    FILE *r_kernel_fptr;
    FILE *r_golden_fptr;

    //TODO change for the file in the local machine where the TB will be executed at first

	ret=0;

	sk_fptr = fopen("./io_values/benes_1/sk_in.dat","r");
	s_fptr = fopen("./io_values/benes_1/s_in.dat","r");
	r_golden_fptr = fopen("./io_values/benes_1/r_out_golden.dat","r");

	//Read from the pk values file into the pk matrix
	for(i=0; i<6492; i++){
		fscanf(sk_fptr, "%d\n", &sk);
		sk_reference_input[i] = sk;
	}
	fclose(sk_fptr);

	for(i=0; i<SYND_BYTES; i++){
		fscanf(s_fptr, "%d\n", &s);
		s_reference_input[i] = s;
	}
	fclose(s_fptr);


	//Read from the r values file into the s matrix
	for(i=0; i<64; i++){
		fscanf(r_golden_fptr, "%ld\n", &r);
		r_golden_output[i] = r;
	}
	fclose(r_golden_fptr);

	benes_kernel(r_kernel_output, sk_reference_input, s_reference_input);

	// Compare the results of the function against expected results

	for(i=0; i<64; i++){
		if(r_kernel_output[i]!=r_golden_output[i]){
			printf("ERROR! Expected r[%d]=%lX but got %lX\n", i, r_golden_output[i].to_uint64(), r_kernel_output[i].to_uint64());
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

