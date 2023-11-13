#include <gmp.h>
#define __gmp_const const

#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include<string.h>
#include"ap_int.h"

#include "decryption_kernel_testbench.h"


int main () { 
    //Establish an initial return value. 0 = success
    int ret;
    int iter;
    int i, j, k;

    //Helper variables
    unsigned char sk;
    unsigned char c;
    unsigned char e;
    int check;


    //Variables that will store the IO data
    //Inputs
    unsigned char *sk_reference_input = (unsigned char *) malloc(sizeof(unsigned char)*(6492));
    unsigned char *c_reference_input = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES));


    //Outputs
    unsigned char *e_golden_output = (unsigned char *) malloc(sizeof(unsigned char)* (SYS_N/8));
    unsigned char *e_kernel_output = (unsigned char *) malloc(sizeof(unsigned char)* (SYS_N/8));
    int *check_golden_output = (int *) malloc(sizeof(int));
    int *check_kernel_output = (int *) malloc(sizeof(int));


    // Call any preliminary functions required to prepare input for the test.
    FILE *sk_fptr;
    FILE *c_fptr;
    FILE *e_golden_fptr;

    //TODO change for the file in the local machine where the TB will be executed at first

	ret=0;

	sk_fptr = fopen("./io_values/dec_1/sk_in.dat","r");
	c_fptr = fopen("./io_values/dec_1/c_in.dat","r");
	e_golden_fptr = fopen("./io_values/dec_1/e_out_golden.dat","r");

	//Read from the pk values file into the pk matrix
	for(i=0; i<6492; i++){
		fscanf(sk_fptr, "%d\n", &sk);
		sk_reference_input[i] = sk;
	}
	fclose(sk_fptr);

	for(i=0; i<SYND_BYTES; i++){
		fscanf(c_fptr, "%d\n", &c);
		c_reference_input[i] = c;
	}
	fclose(c_fptr);

    for(i=0; i<(SYS_N/8); i++){
        fscanf(e_golden_fptr, "%d\n", (e_golden_output + i));
    }
    fclose(e_golden_fptr);


	decryption_kernel(check_kernel_output, e_kernel_output, sk_reference_input, c_reference_input);

	// Compare the results of the function against expected results

	for(i=0; i<(SYS_N/8); i++){
		if(e_kernel_output[i]!=e_golden_output[i]){
			printf("ERROR! Expected e[%d]=%d but got %d\n", i, e_golden_output[i], e_kernel_output[i]);
			ret=1;
		}
	}

	if(*check_kernel_output!=*check_golden_output){
		printf("ERROR! Expected check=%d but got %d\n", *check_golden_output, *check_kernel_output);
		ret=1;
	}


	if (ret != 0) {
			printf("Iter %d Test failed  !!!\n", iter);
			ret=1;
	} else {
			printf("Iter %d Test passed !\n", iter);
	}

    
    return ret;
}

