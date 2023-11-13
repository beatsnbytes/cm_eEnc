#include <gmp.h>
#define __gmp_const const

#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include<string.h>
#include"ap_int.h"

#include "encryption_kernel_testbench.h"


int main () { 
    //Establish an initial return value. 0 = success
    int ret;
    int iter;
    int i, j, k;

    //Helper variables
    aligned_width pk[ALIGNED_BYTES];
    ap_uint<128> v;
    ap_uint<256> key;

    //Variables that will store the IO data
    //Inputs
    aligned_width *pk_reference_input_1 = (aligned_width *) malloc(sizeof(aligned_width)*((PK_NROWS*PK_ROW_BYTES)/(ALIGNED_BYTES*2)));
    aligned_width *pk_reference_input_2 = (aligned_width *) malloc(sizeof(aligned_width)*((PK_NROWS*PK_ROW_BYTES)/(ALIGNED_BYTES*2)));

    ap_uint<256> *key_reference_input = (ap_uint<256> *) malloc(sizeof(ap_uint<256>));
    ap_uint<128> *v_reference_input = (ap_uint<128> *) malloc(sizeof(ap_uint<128>));


    //Outputs
    unsigned char *s_golden_output = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES));
    unsigned char *s_kernel_output_1 = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES)/2);
    unsigned char *s_kernel_output_2 = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES)/2);

    unsigned char *e_golden_output = (unsigned char *) malloc(sizeof(unsigned char)*(SYS_N/8));
    unsigned char *e_kernel_output = (unsigned char *) malloc(sizeof(unsigned char)*(SYS_N/8));

    // Call any preliminary functions required to prepare input for the test.
    FILE *pk_fptr;
    FILE *v_fptr;
    FILE *key_fptr;

    FILE *s_golden_fptr;
    FILE *s_kernel_fptr;

    FILE *e_golden_fptr;
    FILE *e_kernel_fptr;

    //TODO change for the file in the local machine where the TB will be executed at first
for(iter=0; iter<4; iter++){

	ret=0;

	if(iter==0){
		pk_fptr = fopen("./io_values/enc_1/pk_in.dat","r");
		v_fptr = fopen("./io_values/enc_1/v_in.dat","r");
		key_fptr = fopen("./io_values/enc_1/key_in.dat","r");
		s_golden_fptr = fopen("./io_values/enc_1/s_out_golden.dat","r");
		s_kernel_fptr = fopen("./io_values/enc_1/s_out_kernel.dat","w+");
		e_golden_fptr = fopen("./io_values/enc_1/e_out_golden.dat","r");
		e_kernel_fptr = fopen("./io_values/enc_1/e_out_kernel.dat","w+");
	}else if(iter==1){
		pk_fptr = fopen("./io_values/enc_2/pk_in.dat","r");
		v_fptr = fopen("./io_values/enc_2/v_in.dat","r");
		key_fptr = fopen("./io_values/enc_2/key_in.dat","r");
		s_golden_fptr = fopen("./io_values/enc_2/s_out_golden.dat","r");
		s_kernel_fptr = fopen("./io_values/enc_2/s_out_kernel.dat","w+");
		e_golden_fptr = fopen("./io_values/enc_2/e_out_golden.dat","r");
		e_kernel_fptr = fopen("./io_values/enc_2/e_out_kernel.dat","w+");
	}else if(iter==2){
		pk_fptr = fopen("./io_values/enc_3/pk_in.dat","r");
		v_fptr = fopen("./io_values/enc_3/v_in.dat","r");
		key_fptr = fopen("./io_values/enc_3/key_in.dat","r");
		s_golden_fptr = fopen("./io_values/enc_3/s_out_golden.dat","r");
		s_kernel_fptr = fopen("./io_values/enc_3/s_out_kernel.dat","w+");
		e_golden_fptr = fopen("./io_values/enc_3/e_out_golden.dat","r");
		e_kernel_fptr = fopen("./io_values/enc_3/e_out_kernel.dat","w+");
	}else{
		pk_fptr = fopen("./io_values/enc_4/pk_in.dat","r");
		v_fptr = fopen("./io_values/enc_4/v_in.dat","r");
		key_fptr = fopen("./io_values/enc_4/key_in.dat","r");
		s_golden_fptr = fopen("./io_values/enc_4/s_out_golden.dat","r");
		s_kernel_fptr = fopen("./io_values/enc_4/s_out_kernel.dat","w+");
		e_golden_fptr = fopen("./io_values/enc_4/e_out_golden.dat","r");
		e_kernel_fptr = fopen("./io_values/enc_4/e_out_kernel.dat","w+");
	}
    //Read from the pk values file into the pk matrix
    for(i=0; i<((PK_NROWS*PK_ROW_BYTES)/(ALIGNED_BYTES*2)); i++){
    	(pk_reference_input_1[i]) = 0;
    	for(j=0; j<(ALIGNED_BYTES); j++){
    		pk[j] = 0;
    		fscanf(pk_fptr, "%d\n", (pk + j));
//    		if (i<2) printf("\npk_1 = %d\n", pk[j].to_uint());
    		(pk_reference_input_1[i]) |= pk[j]<<(j*8);
    	}
    }

    for(i=0; i<((PK_NROWS*PK_ROW_BYTES)/(ALIGNED_BYTES*2)); i++){
    	(pk_reference_input_2[i]) = 0;
    	for(j=0; j<(ALIGNED_BYTES); j++){
    		pk[j] = 0;
    		fscanf(pk_fptr, "%d\n", (pk + j));
//    		if (i<2) printf("\npk_2 = %d\n", pk[j].to_uint());
    		(pk_reference_input_2[i]) |= pk[j]<<(j*8);
    	}
    }

    fclose(pk_fptr);

    *v_reference_input=0;
	for(j=0; j<16; j++){
		v = 0;
		fscanf(v_fptr, "%d\n", &v);
		*v_reference_input |= v<<(j*8);
//		printf("v[%d]=%d\n", j, v.to_uint());
	}
    fclose(v_fptr);

    *key_reference_input=0;
	for(j=0; j<32; j++){
		key = 0;
		fscanf(key_fptr, "%d\n", &key);
//		printf("key = %d\n", key.to_uint());
		*key_reference_input |= key<<(j*8);
	}
    fclose(key_fptr);


    //Read from the s values file into the s matrix
    for(i=0; i<(SYND_BYTES); i++){
        fscanf(s_golden_fptr, "%d\n", (s_golden_output + i));
    }
    fclose(s_golden_fptr);

    for(i=0; i<(SYS_N/8); i++){
        fscanf(e_golden_fptr, "%d\n", (e_golden_output + i));
    }
    fclose(e_golden_fptr);

    encryption_kernel(s_kernel_output_1, s_kernel_output_2, e_kernel_output, pk_reference_input_1, pk_reference_input_2, v_reference_input, key_reference_input);

    // }

    // Capture the output results of the function, write to a file

    //Read from the s kernel output matrix and store to a file 
//    for(i=0; i<(SYND_BYTES/2); i++){
//        fprintf(s_kernel_fptr, "%d\n", *(s_kernel_output_1 + i));
//    }
    
    // Compare the results of the function against expected results
    
    for(i=0; i<(SYND_BYTES)/2; i++){
    	if(*(s_kernel_output_1 + i)!=*(s_golden_output + i)){
    		printf("ERROR! Expected s_1[%d]=%d but got %d\n", i, *(s_golden_output + i), *(s_kernel_output_1 + i));
    		ret=1;
    	}
    }

    for(i=0; i<(SYND_BYTES)/2; i++){
    	if(*(s_kernel_output_2 + i)!= *(s_golden_output + i + (SYND_BYTES)/2)){
    		printf("ERROR! Expected s_2[%d]=%d but got %d\n", i, *(s_golden_output + i + (SYND_BYTES)/2), *(s_kernel_output_2 + i));
    		ret=1;
    	}
    }


    for(i=0; i<(SYS_N/8); i++){
        fprintf(e_kernel_fptr, "%d\n", *(e_kernel_output + i));
    }

    // Compare the results of the function against expected results

    for(i=0; i<SYS_N/8; i++){
    	if(*(e_kernel_output + i)!=*(e_golden_output + i)){
    		printf("ERROR! Expected e[%d]=%d but got %d\n", i, *(e_golden_output + i), *(e_kernel_output + i));
    		ret=1;
    	}
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
