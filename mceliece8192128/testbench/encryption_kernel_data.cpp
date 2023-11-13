//#include "encryption_monolithic_dataflow.cpp"
#include "encryption_kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(){

  int ret_1=0;
  int ret_2=0;

  int i, j;
  int iter=1;

  // Call any preliminary functions required to prepare input for the test.
  char data_file_ext[] = ".dat";


  char iter_num[100];


FILE *fp;

unsigned char *dummy_pk = (unsigned char *) malloc(sizeof(unsigned char)*(crypto_kem_mceliece348864_ref_PUBLICKEYBYTES));
data_packed_pk pk_in[3072];

unsigned char *c_out = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES));
unsigned char *c_tb = (unsigned char *) malloc(sizeof(unsigned char)*(SYND_BYTES));

unsigned char *e_out = (unsigned char *) malloc( sizeof(unsigned char) * (SYS_N/8));
data_packed_e e_tb[4];
//unsigned char *e_tb = (unsigned char *) malloc( sizeof(unsigned char) * (SYS_N/8));

unsigned char *key_in = (unsigned char *) malloc( sizeof(unsigned char) * 32);
unsigned char *v_in = (unsigned char *) malloc( sizeof(unsigned char) * 16);


for(int iteration=0; iteration<iter; iteration++){

	//Esablish an initial return value. 0 = success
	int ret_1=0;
	int ret_2=0;

	char pk_data_path[] = "./matrix_data/pk_golden_in_";
	char v_data_path[] = "./matrix_data/v_in_";
	char key_data_path[] = "./matrix_data/key_in_";

	char e_out_data_path[] = "./matrix_data/e_out_";
	char e_tb_data_path[] = "./matrix_data/e_tb_";

	char c_out_data_path[] = "./matrix_data/c_out_";
	char c_tb_data_path[] = "./matrix_data/c_tb_";


    // convert iter to string
	sprintf(iter_num, "%d", iteration);
    //concat iter to file ext
    strcat(iter_num, data_file_ext);



	//pk_in
	strcat(pk_data_path, iter_num);
	fp = fopen(pk_data_path, "r");
	for(i=0; i<(crypto_kem_mceliece348864_ref_PUBLICKEYBYTES); i++){

		   fscanf(fp, "%d\n ", (dummy_pk + i));
	   }
	fclose(fp);

	for(i=0; i<768; i++){
		for(j=0; j<340; j++){
			pk_in[(i*4)+(j/85)].packed_values[j%85] = *(dummy_pk + i*340+j);
		}
	}




	//v_in
	strcat(v_data_path, iter_num);
	fp = fopen(v_data_path, "r");
	for(i=0; i<(16); i++){

		   fscanf(fp, "%d\n ", (v_in + i));
	   }
	fclose(fp);

	//key_in
	strcat(key_data_path, iter_num);
	fp = fopen(key_data_path, "r");
	for(i=0; i<(32); i++){

		   fscanf(fp, "%d\n ", (key_in + i));
	   }
	fclose(fp);

	//e_out
	strcat(e_out_data_path, iter_num);
	//printf("\n %s \n", e_out_data_path);
	fp = fopen(e_out_data_path, "r");
	for(i=0; i<(SYS_N/8); i++){

		   fscanf(fp, "%d\n ", (e_out + i));
	   }
	fclose(fp);

	//c_out
	strcat(c_out_data_path, iter_num);
	//printf("\n %s \n", c_out_data_path);
	fp = fopen(c_out_data_path, "r");
	for(i=0; i<(SYND_BYTES); i++){

		   fscanf(fp, "%d\n ", (c_out + i));
	   }
	fclose(fp);




	// Call the top-level function multiple times, passing input stimuli as needed.
	 //encryption_monolithic_dataflow(c_tb, e_tb, pk_in);
	 encryption_monolithic_packed(c_tb, e_tb, pk_in, v_in, key_in);

	//Capture the output results of the function, write to a file

	 //e_tb
	strcat(e_tb_data_path, iter_num);
	fp = fopen(e_tb_data_path, "w+");

	for(i=0; i<(SYS_N/8); i++){
		fprintf(fp, "%d\n", *(e_tb + i));
	   }
	fclose(fp);

	 //c_tb
	strcat(c_tb_data_path, iter_num);
	fp = fopen(c_tb_data_path, "w+");

	for(i=0; i<(SYND_BYTES); i++){
		fprintf(fp, "%d\n", *(c_tb + i));
	   }
	fclose(fp);


    int k;
    printf("encrypt e: positions");
    for (k = 0;k < SYS_N;++k)
      if (e_tb[(k/8)/109].packed_values[(k/8)%109] & (1 << (k&7))) //e_out[k/8]
        printf(" %d",k);
    printf("\n");


	for(i=0; i<SYS_N/8; i++){
		//if(*(e_out + i)!=*(e_tb + i)){
		if(*(e_out + i)!= e_tb[i/109].packed_values[i%109]){
			printf("ERROR in e at index %d. %d different than expected %d!\n", i, e_tb[i/109].packed_values[i%109], *(e_out + i));
			ret_1=1;
		}
	}
	printf("Test %d for e[] = %d\n", iteration, ret_1);

	for(i=0; i<SYND_BYTES; i++){
		if(*(c_out + i)!=*(c_tb + i)){
			printf("ERROR in c at index %d. %d different than expected %d!\n", i, *(c_tb + i), *(c_out + i));
			ret_2=1;
		}
	}
	printf("Test %d for c[] = %d\n", iteration, ret_2);


}

  if ((ret_1 | ret_2) != 0) {
        printf("Test failed  !!!\n");
  } else {
        printf("Tests passed !\n");
  }

  return (ret_1 | ret_2);
}

