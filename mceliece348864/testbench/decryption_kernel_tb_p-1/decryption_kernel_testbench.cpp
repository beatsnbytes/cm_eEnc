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
    int iter=0;
    int i, j, k;

    //Helper variables
    sk_2k_vec sk_2k_wide;
    sk_1k_vec sk_1k_wide;
    unsigned char c;
    unsigned char e;
    int check;


    //Variables that will store the IO data
    //Inputs
    unsigned char *c_reference_input = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES));
    sk_2k_vec *sk_reference_rev_input = (sk_2k_vec *) malloc(sizeof(sk_2k_vec)*(23));
    sk_2k_vec *sk_reference_fwd_input = (sk_2k_vec *) malloc(sizeof(sk_2k_vec)*(23));
    sk_1k_vec *sk_irr_reference_input = (sk_1k_vec *) malloc(sizeof(sk_1k_vec));


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
	for(iter=0; iter<4; iter++){

		if(iter==0){
			sk_fptr = fopen("./io_values/dec_1/sk_in.dat","r");
			c_fptr = fopen("./io_values/dec_1/c_in.dat","r");
			e_golden_fptr = fopen("./io_values/dec_1/e_out_golden.dat","r");
		}else if(iter==1){
			sk_fptr = fopen("./io_values/dec_2/sk_in.dat","r");
			c_fptr = fopen("./io_values/dec_2/c_in.dat","r");
			e_golden_fptr = fopen("./io_values/dec_2/e_out_golden.dat","r");
		}else if(iter==2){
			sk_fptr = fopen("./io_values/dec_3/sk_in.dat","r");
			c_fptr = fopen("./io_values/dec_3/c_in.dat","r");
			e_golden_fptr = fopen("./io_values/dec_3/e_out_golden.dat","r");
		}else{
			sk_fptr = fopen("./io_values/dec_4/sk_in.dat","r");
			c_fptr = fopen("./io_values/dec_4/c_in.dat","r");
			e_golden_fptr = fopen("./io_values/dec_4/e_out_golden.dat","r");
		}



		//SK_IRR
		//Read from the sk values file into the sk matrix

		if(iter==0){
			sk_fptr = fopen("./io_values/dec_1/sk_in.dat","r");
		}else if(iter==1){
			sk_fptr = fopen("./io_values/dec_2/sk_in.dat","r");
		}else if(iter==2){
			sk_fptr = fopen("./io_values/dec_3/sk_in.dat","r");
		}else{
			sk_fptr = fopen("./io_values/dec_4/sk_in.dat","r");
		}

		(*sk_irr_reference_input) = 0;
		for(i=0; i<128; i++){
			sk_1k_wide=0;
			fscanf(sk_fptr, "%d\n", &sk_1k_wide);
			(*sk_irr_reference_input) |= (sk_1k_wide)<<(i*8);
		}
		fclose(sk_fptr);


		//SK FORWARD & REVERSE
		if(iter==0){
			sk_fptr = fopen("./io_values/dec_1/sk_in.dat","r");
		}else if(iter==1){
			sk_fptr = fopen("./io_values/dec_2/sk_in.dat","r");
		}else if(iter==2){
			sk_fptr = fopen("./io_values/dec_3/sk_in.dat","r");
		}else{
			sk_fptr = fopen("./io_values/dec_4/sk_in.dat","r");
		}

		//Has to skip the 128 first lines (IRR_BYTES)
		for(i=0; i<128; i++){
			fscanf(sk_fptr, "%d\n", &sk_1k_wide);
		}

		//Read the remaining and save them in the reverse order
		for(i=22; i>=0; i--){
			sk_reference_rev_input[i] = 0;
			sk_reference_fwd_input[22-i] = 0;
			for(j=0; j<256; j++){
				sk_2k_wide=0;
				fscanf(sk_fptr, "%d\n", &sk_2k_wide);
				sk_reference_rev_input[i] |= (sk_2k_wide)<<(j*8);
				sk_reference_fwd_input[22-i] |= (sk_2k_wide)<<(j*8);
			}
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


		decryption_kernel_active(check_kernel_output, e_kernel_output, sk_irr_reference_input, sk_reference_rev_input,  sk_reference_fwd_input, c_reference_input);

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
	}
    
    return ret;
}

