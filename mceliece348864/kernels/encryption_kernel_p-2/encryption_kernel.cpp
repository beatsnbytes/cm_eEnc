#ifndef __SYNTHESIS__
#include <gmp.h>
#define __gmp_const const
#endif
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#define AP_INT_MAX_W 4096
#include"ap_int.h"
#include "xf_security/ecb.hpp"
#include "hls_stream.h"
#include "params.h"

#define ALIGNED_VALUES 85
#define BIT_WIDTH 680
#define ALIGNED_WIDTH 1024
#define ALIGNED_BYTES 128
#define PACK_FACTOR BIT_WIDTH/8
typedef ap_uint<BIT_WIDTH> custom_width;
typedef ap_uint<ALIGNED_WIDTH> aligned_width;




typedef ap_uint<128> v_vector;
typedef ap_uint<256> key_vector;

//TODO need to be changed
typedef ap_uint<12> gf_value;
typedef ap_uint<1> b_width;
typedef ap_uint<128> cipher_vector;


void dump_to_v_stream_aes(hls::stream<v_vector> &v_stream, hls::stream<v_vector> &v_stream_internal, v_vector v_local_vector){

	int i;
	//Issue here. when 0xFF not the same behavior as sw
	//Key and V do not change size so should I change this here?
    v_vector one_msb = ap_uint<128> ("0x01000000000000000000000000000000", 16);
    v_vector one_msb_alt = ap_uint<128> ("0x00010000000000000000000000000000", 16);
    v_vector change_dec = ap_uint<128> ("0xFF000000000000000000000000000000", 16);
    v_vector current_v=0;

    current_v = v_local_vector;

    for(int i=0; i<((2*SYS_T)/8); i++){
	#pragma HLS PIPELINE II=1
    	if(current_v.range(127, 120)==0xFF){
    		current_v += one_msb_alt;
    	}
        	current_v += one_msb;

    	v_stream.write(current_v);
    }

    v_stream_internal.write(current_v);

}

void dump_to_v_stream_update(hls::stream<v_vector> &v_stream_update, hls::stream<v_vector> &v_stream_internal){

	int i;
    v_vector one_msb = ap_uint<128> ("0x01000000000000000000000000000000", 16);

    v_vector current_v=0;
    current_v = v_stream_internal.read();

    for(int i=0; i<3; i++){
	#pragma HLS PIPELINE II=1
    	current_v += one_msb;
    	v_stream_update.write(current_v);
    }
}



void dump_to_key_stream(hls::stream<key_vector> &key_stream, hls::stream<key_vector> &key_stream_update, key_vector key_local_vector){

    key_vector current_key=0;

    current_key = key_local_vector;
    key_stream.write(current_key);
    key_stream_update.write(current_key);
}

void dump_to_update_stream(hls::stream<bool> &key_update_stream, bool key_update){

    key_update_stream.write(key_update);
}

void dump_to_output_ciphertext_stream( hls::stream<cipher_vector> &output_ciphertext_stream, hls::stream<cipher_vector> &ciphertext_stream){

	int i, j;
	cipher_vector ciphertext_block;

	for(i=0; i<((2*SYS_T)/8); i++){
	#pragma HLS PIPELINE II=1
		ciphertext_block = ciphertext_stream.read();
		output_ciphertext_stream.write(ciphertext_block);
	}


}

void dump_updated_key_v_values(v_vector *v_local_vector, key_vector *key_local_vector, hls::stream<v_vector> &ciphertext_stream){

	key_vector key_final_vector;
	v_vector v_final_vector;

	key_final_vector.range(127, 0) = ciphertext_stream.read();
	key_final_vector.range(255, 128) = ciphertext_stream.read();
	*key_local_vector = key_final_vector;

	v_final_vector = ciphertext_stream.read();
	*v_local_vector = v_final_vector;

}

void randombytes(hls::stream<cipher_vector> &output_ciphertext_stream,  v_vector *v_local_vector_out, key_vector *key_local_vector_out, v_vector v_local_vector_in, key_vector key_local_vector_in)
{

    static hls::stream<v_vector> v_stream_aes("v stream aes");
	#pragma HLS STREAM variable=v_stream_aes depth=16 //Contains 16 16B values that are incrementing
    static hls::stream<v_vector> v_stream_update("v stream update");
	#pragma HLS STREAM variable=v_stream_update depth=3 //Contains 3 16B values that are incrementing
    static hls::stream<v_vector> v_stream_internal("v stream internal");
	#pragma HLS STREAM variable=v_stream_internal depth=2 //It contains a single value to be passed to the update
    static hls::stream<key_vector> key_stream_aes("key stream aes");
	#pragma HLS STREAM variable=key_stream_aes depth=2 //Contains 1 key value
    static hls::stream<key_vector> key_stream_update("key stream update");
	#pragma HLS STREAM variable=key_stream_update depth=3 //Contains 1 key value
    static hls::stream<cipher_vector> ciphertext_stream("ciphertext stream");
	#pragma HLS STREAM variable=ciphertext_stream depth=16 //Contains 16 16B cipheretext values
    static hls::stream<cipher_vector> ciphertext_stream_update("ciphertext stream");
	#pragma HLS STREAM variable=ciphertext_stream_update depth=3 //Contains 3 16B cipheretext values

	#pragma HLS dataflow


    //Writes 16 values to the v_stream
    dump_to_v_stream_aes(v_stream_aes, v_stream_internal, v_local_vector_in); //reads v_local_vector
    dump_to_key_stream(key_stream_aes, key_stream_update, key_local_vector_in);  //reads key_local_vector

    //Reads 16 values from the v_stream and computes 16 consecutive ciphertext blocks
    xf::security::aes256EcbEncrypt(v_stream_aes, key_stream_aes, ciphertext_stream, ((SYS_T*2)/8));

    dump_to_output_ciphertext_stream(output_ciphertext_stream, ciphertext_stream);
    dump_to_v_stream_update(v_stream_update, v_stream_internal);  //reads v_local_vector

    //Reads 3 values from the v_stream and computes 3 consecutive ciphertext blocks (2 for the key update and 1 for the IV update)
    xf::security::aes256EcbEncrypt(v_stream_update, key_stream_update, ciphertext_stream_update, 3);

    dump_updated_key_v_values(v_local_vector_out, key_local_vector_out, ciphertext_stream_update); //writes both v_local_vector and key_local_vector
}

//TODO use sorter that can take streams?
void even_odd_sort(gf_value *a)
{
	unsigned int n = SYS_T;
	sort_loop:
	for (unsigned int i = 0; i < n / 2; i++) {
	#pragma HLS PIPELINE II=1
		for (unsigned int j = 0; j + 1 < n; j += 2) {
			if (a[j] > a[j + 1]) {
			gf_value t = a[j];
			a[j] = a[j + 1];
			a[j + 1] = t;
			}
		}
		for (unsigned int j = 1; j + 1 < n; j += 2) {
			if (a[j] > a[j + 1]) {
				gf_value t = a[j];
				a[j] = a[j + 1];
				a[j + 1] = t;
			}
		}
	}
}

//TODO split/merge streams?
//TODO correct the following loop getting the individual valuesin the correct order
void compose_gf_values(gf_value *gf_mat, hls::stream<cipher_vector> &aes_ciphertext_stream){

	int i, j;
	cipher_vector aes_block;
	gf_value gf_element;

	//TODO concat consecutive gf values together from here and continue down the line
	//TODO Every 16B I am throwing away 2....can I do it earlier and savew some overhead?
	LOOP_COMPOSE_GF_VALUES:
	for (i = 0; i <(SYS_T*2)/8; i++){
	#pragma HLS PIPELINE II=1
		aes_block = aes_ciphertext_stream.read(); //Each aes block is 128b=16B=8 load_gf() iterations
		gf_element=0;
		//TODO verify this change here
		for(j=0; j<16; j+=2){

			gf_element = aes_block.range(((j+2)*8)-(16-GFBITS+1), 0+j*8);
			gf_mat[i*8+(j/2)] = gf_element;
		}
	}
}

//TODO think of how can I pass indexes in stream...since the size will vsry there will be leftover elements
void index_counting(bool *approved, gf_value *indexes_mat, gf_value *gf_mat){

	//TODO concat gf values and check for higher value to compare??
	gf_value gf_element;
	int count = 0;
	int i;

	LOOP_COMPARE_TO_SYS_N:
	for (i = 0; i < (SYS_T*2); i++){
		if(gf_mat[i] < SYS_N){
			indexes_mat[count++] = gf_mat[i];
		}
	}

	if(count<SYS_T){
		*approved=0;
	}
}

//TODO implement stream version
void repetition_check(bool *approved, gf_value *indexes_mat){

	int i;
	LOOP_REPETITION_CHECK:
	for (i = 1; i < SYS_T; i++){
	#pragma HLS PIPELINE
		if (indexes_mat[i-1] == indexes_mat[i]){
			*approved = 0;
		}
	}
}

void dump_to_e_1_stream(hls::stream<unsigned char> &e_1_stream, hls::stream<unsigned char> &e_1_internal_stream){
	int i;

	for(i=0; i<(SYND_BYTES)/2; i++){
		e_1_stream.write(e_1_internal_stream.read());
	}
}

void dump_to_e_2_stream(hls::stream<custom_width> &e_2_stream_1, hls::stream<custom_width> &e_2_stream_2, unsigned char *e_2_mat){

	int i, j, k;
	custom_width e_2_stream_value=0;
	custom_width e_2_tmp=0;

	int cnt;
	for(k=0; k<PK_NROWS/2; k++){
	#pragma HLS PIPELINE
		cnt=0;
		for(i=(SYND_BYTES)+1; i<=SYS_N/8; i++){
			e_2_tmp = e_2_mat[i-1];

			e_2_stream_value |= e_2_tmp<<(((i-1-SYND_BYTES)%(PACK_FACTOR))*8);
			if(((i-SYND_BYTES)%(PACK_FACTOR))==0){
				e_2_stream_1.write(e_2_stream_value);
				e_2_stream_2.write(e_2_stream_value);
				e_2_stream_value = 0;
			}

		}
	}

}


void dump_to_e_out_stream(hls::stream<unsigned char> &e_out_stream, hls::stream<unsigned char> &e_out_internal_stream){
	int i;
	for(i=0; i<SYS_N/8; i++){
		e_out_stream.write(e_out_internal_stream.read());
	}
}

void write_to_input_streams(hls::stream<v_vector> &v_input_stream, hls::stream<key_vector> &key_input_stream, v_vector *v_local_vector, key_vector *key_local_vector){

	v_input_stream.write(*v_local_vector);
	key_input_stream.write(*key_local_vector);

}

void read_from_output_streams(v_vector *v_local_vector, key_vector *key_local_vector, hls::stream<v_vector> &v_output_stream, hls::stream<key_vector> &key_output_stream){

	*v_local_vector = v_output_stream.read();
	*key_local_vector = key_output_stream.read();
}


void load_v_key(v_vector *v_local_vector, key_vector *key_local_vector, v_vector *v, key_vector *key){

	*v_local_vector = *v;

	*key_local_vector = *key;

}

void dump_to_indexes_stream(hls::stream<gf_value> &indexes_stream, gf_value *indexes_mat){
	int i;

	for(i=0; i<SYS_T; i++){
		indexes_stream.write(indexes_mat[i]);
	}

}

void compute_indexes_mat(hls::stream<gf_value> &indexes_stream, v_vector *v, key_vector *key){

	static hls::stream<cipher_vector> aes_ciphertext_stream("aes cipherteext");
	#pragma HLS STREAM variable=aes_ciphertext_stream depth=16
	static hls::stream<gf_value> gf_values_stream("gf values");
	#pragma HLS STREAM variable=gf_values_stream depth=256 //There should be in total 128 values. In case of dataflow this wont be full probably.
    static hls::stream<v_vector> v_input_stream("v input stream gen_e");
	#pragma HLS STREAM variable=v_input_stream depth=1 //Stream with size 1 since we are performing 1 randomytes() at a time
    static hls::stream<v_vector> v_output_stream("v output stream gen_e");
	#pragma HLS STREAM variable=v_output_stream depth=1 //Stream with size 1 since we are performing 1 randomytes() at a time
    static hls::stream<key_vector> key_input_stream("key input stream gen_e");
	#pragma HLS STREAM variable=key_input_stream depth=1 //Stream with size 1 since we are performing 1 randomytes() at a time
    static hls::stream<key_vector> key_output_stream("key output stream gen_e");
	#pragma HLS STREAM variable=key_output_stream depth=1 //Stream with size 1 since we are performing 1 randomytes() at a time


    //TODO think of completely partitioning gf_mat and indexes_mat. Does it help? util without the partition?

    gf_value gf_mat[SYS_T*2];
#pragma HLS ARRAY_PARTITION variable=gf_mat factor=128


	v_vector v_local_vector_in, v_local_vector_out;
	key_vector key_local_vector_in, key_local_vector_out;

	gf_value indexes_mat[SYS_T*2];
#pragma HLS ARRAY_PARTITION variable=indexes_mat factor=128

	bool approved;

	load_v_key(&v_local_vector_out, &key_local_vector_out, v, key);


	while (1)
	{
		approved=1;

		key_local_vector_in = key_local_vector_out;
		v_local_vector_in = v_local_vector_out;

		randombytes(aes_ciphertext_stream, &v_local_vector_out, &key_local_vector_out, v_local_vector_in, key_local_vector_in);

		compose_gf_values(gf_mat, aes_ciphertext_stream);

		index_counting(&approved, indexes_mat, gf_mat);

		if (!approved){
			continue;
		}

		even_odd_sort(indexes_mat);

		// check for repetition
		repetition_check(&approved, indexes_mat);

		if (approved){
			break;
		}

	}

	dump_to_indexes_stream(indexes_stream, indexes_mat);

}

void compute_e_streams(hls::stream<unsigned char> &e_1_internal_stream_1, hls::stream<unsigned char> &e_1_internal_stream_2, unsigned char *e_2_mat, hls::stream<unsigned char> &e_out_internal_stream, hls::stream<gf_value> &indexes_stream){


	gf_value indexes_mat_msb;
	gf_value indexes_mat_lsb;
	ap_uint<32> mask, index_shift;
	unsigned char e_int, e_tmp;
	int i, j;

	gf_value indexes_mat[SYS_T];
	#pragma HLS ARRAY_PARTITION variable=indexes_mat factor=64



	LOOP_INDEXES:
	for(i=0; i<SYS_T; i++){
		indexes_mat[i]=indexes_stream.read();
	}
	int cnt=0;

	LOOP_COMPUTE_E_MAT:
	for (i = 0; i < SYS_N/8; i++)
	{
	#pragma HLS PIPELINE II=1
		e_int = 0;
		e_tmp = 0;

		for (j = 0; j < SYS_T; j++)
		{
			//TODO fix
			indexes_mat_msb = indexes_mat[j].range((GFBITS-1), 3);
			indexes_mat_lsb = indexes_mat[j].range(2, 0);
			//Same mask? Can it be done differently on HW?
			mask = i ^ (indexes_mat_msb);
			mask -= 1;
			mask >>= 31;
			mask = -mask;
			mask = mask & 0xFF;

			e_int |= (1 << (indexes_mat_lsb)) & mask;
		}
		e_out_internal_stream.write(e_int);
		if(i<SYND_BYTES/2){
			e_1_internal_stream_1.write(e_int);
		}else if((i>=SYND_BYTES/2) && (i<SYND_BYTES)){
			e_1_internal_stream_2.write(e_int);

		}else{
			e_2_mat[i] = e_int;
		}


	}


}


static void gen_e(hls::stream<unsigned char> &e_1_stream_1, hls::stream<unsigned char> &e_1_stream_2, hls::stream<custom_width> &e_2_stream_1, hls::stream<custom_width> &e_2_stream_2, hls::stream<unsigned char> &e_out_stream,  v_vector *v, key_vector *key)
{
	int i, j, eq, count;

    static hls::stream<unsigned char> e_out_internal_stream("e_out_internal stream");
	#pragma HLS STREAM variable=e_out_internal_stream depth=2

    static hls::stream<unsigned char> e_1_internal_stream_1("e_1_internal stream_1");
	#pragma HLS STREAM variable=e_1_internal_stream_1 depth=48
    static hls::stream<unsigned char> e_1_internal_stream_2("e_1_internal stream_2");
	#pragma HLS STREAM variable=e_1_internal_stream_2 depth=48

    static hls::stream<gf_value> indexes_stream("indexes_stream");
	#pragma HLS STREAM variable=indexes_stream depth=2

    unsigned char e_2_mat[SYS_N/8];
	#pragma HLS ARRAY_PARTITION variable=e_2_mat cyclic factor=64



#pragma HLS dataflow


	compute_indexes_mat(indexes_stream, v, key);
	compute_e_streams(e_1_internal_stream_1, e_1_internal_stream_2, e_2_mat, e_out_internal_stream, indexes_stream);

	dump_to_e_1_stream(e_1_stream_1, e_1_internal_stream_1);
	dump_to_e_1_stream(e_1_stream_2, e_1_internal_stream_2);

	dump_to_e_2_stream(e_2_stream_1, e_2_stream_2, e_2_mat);

	dump_to_e_out_stream(e_out_stream, e_out_internal_stream);
}

void store_e(unsigned char *e, hls::stream<unsigned char> &e_out_stream){
	int i;

	for(i=0; i<SYS_N/8; i++){
		*(e+i) = e_out_stream.read();
	}
}

//TODO break pk reading between 2 ports?
void axi_aligned_pk_reader(hls::stream<aligned_width> &pk_aligned_stream, aligned_width *pk){

	int i;
	LOOP_LOAD_PK:
	for(i=0; i<((PK_NROWS*PK_ROW_BYTES)/(ALIGNED_BYTES*2)); i++){
		#pragma HLS PIPELINE II=1
		pk_aligned_stream.write( *(pk+i) );
	}

}

void pk_stream_resizer(hls::stream<custom_width> &pk_stream, hls::stream<aligned_width> &pk_aligned_stream){

int i;
	int curr_idx = 0;
	int whole_idx = 0;
	int remnant_bytes=0;
	aligned_width aligned_value;
	custom_width whole=0;
	custom_width complement=0;
	custom_width remnant=0;

	custom_width current_pk;


	for(i=0; i<((PK_NROWS*PK_ROW_BYTES)/(ALIGNED_BYTES*2)); i++){
	#pragma HLS PIPELINE II=1
		aligned_value = pk_aligned_stream.read();
		curr_idx = BIT_WIDTH - remnant_bytes;
		complement = (aligned_value.range(curr_idx-1, 0));
		whole.range(remnant_bytes-1, 0) = remnant;
		whole.range(BIT_WIDTH-1, remnant_bytes) = complement;
		pk_stream.write(whole);
		if( curr_idx <= (ALIGNED_WIDTH -  BIT_WIDTH)){
			 whole_idx = (curr_idx + BIT_WIDTH);
			 whole =  aligned_value.range(whole_idx-1, curr_idx);
			 pk_stream.write(whole);
			 remnant_bytes = ALIGNED_WIDTH - whole_idx;
			 if(remnant_bytes==0){
				 remnant = aligned_value.bit(ALIGNED_WIDTH-1);
			 }else{
				remnant = aligned_value.range(ALIGNED_WIDTH-1, whole_idx);
			 }
		}else{
			 remnant_bytes = ALIGNED_WIDTH - curr_idx;
			 remnant = aligned_value.range(ALIGNED_WIDTH-1, curr_idx);
		}
	}
}

void pk_stream_resizer_alt(hls::stream<custom_width> &pk_stream, hls::stream<aligned_width> &pk_aligned_stream){

int i;
	int curr_idx = 0;
	int whole_idx = 0;
	int remnant_bytes=0;
	aligned_width aligned_value;
	custom_width whole=0;
	custom_width complement=0;
	custom_width current_pk;
	int remnant_bits=0;
	custom_width bundle=0;
	aligned_width remnant=0;

	for(i=0; i<((PK_NROWS*PK_ROW_BYTES)/(ALIGNED_VALUES*2)); i++){
	#pragma HLS PIPELINE II=1
		if(remnant_bits<BIT_WIDTH){
			aligned_value = pk_aligned_stream.read();

			curr_idx = BIT_WIDTH - remnant_bits;

			bundle.range(remnant_bits-1, 0) = remnant;
			bundle.range(BIT_WIDTH-1, remnant_bits) = aligned_value.range(curr_idx-1, 0);
			pk_stream.write(bundle);

			remnant = aligned_value.range(ALIGNED_WIDTH-1, curr_idx);

			remnant_bits = ALIGNED_WIDTH - curr_idx;
		}else{
			bundle = remnant.range(BIT_WIDTH-1, 0);
			pk_stream.write(bundle);

			remnant = remnant.range(remnant_bits-1, BIT_WIDTH);

			remnant_bits = remnant_bits - BIT_WIDTH;

		}

	}

}


void compute_n_store_s(unsigned char *s, hls::stream<unsigned char> &e_1_stream, hls::stream<b_width> &b_stream){

	int i, j;
	unsigned char s_mat_custom=0;
	//TODO check with 2 dimensional s matrix

	LOOP_S_COMPUTE_N_DUMP_TO_STREAM:
	for(i=0; i<(SYND_BYTES)/2; i++){
	#pragma HLS PIPELINE II=8
		s_mat_custom = e_1_stream.read();
		for(j=0; j<8; j++){
			s_mat_custom ^= ((b_stream.read()).bit(0) << (j%8));
		}
		*(s+i) =  (s_mat_custom);
	}

}

void compute_b_mat(hls::stream<b_width> &b_stream, hls::stream<custom_width> &pk_stream, hls::stream<custom_width> &e_2_stream)
{

	int i, j;
	custom_width b;
	custom_width b_1;

	custom_width pk_tmp, e_tmp;

	int idx=0;
	int odd_cnt=1;

	LOOP_PK_NROWS_ITERATE:
	for (i = 0; i < (PK_NROWS/2); i++)
	{
	#pragma HLS PIPELINE II=1
	
		b = 0;
		LOOP_PRODUCT:
		for (j = 0; j < (PK_NCOLS/(BIT_WIDTH)); j++){
			b ^= pk_stream.read() & e_2_stream.read();
		}

		b_stream.write(b.xor_reduce());

 
	}

}



void encryption_kernel(unsigned char *s_1, unsigned char *s_2, unsigned char *e, aligned_width *pk_1, aligned_width *pk_2, v_vector *v, key_vector *key)
{
	#pragma HLS INTERFACE m_axi     port=s_1     offset=slave depth=48  bundle=gmem
	#pragma HLS INTERFACE m_axi     port=s_2     offset=slave depth=48  bundle=gmem1
	#pragma HLS INTERFACE m_axi     port=e       offset=slave depth=436 bundle=gmem2
    #pragma HLS INTERFACE m_axi     port=pk_1    offset=slave depth=1020 bundle=gmem3
	#pragma HLS INTERFACE m_axi     port=pk_2    offset=slave depth=1020 bundle=gmem4
    #pragma HLS INTERFACE m_axi     port=v       offset=slave depth=1 bundle=gmem5
	#pragma HLS INTERFACE m_axi     port=key     offset=slave depth=1 bundle=gmem6
    #pragma HLS INTERFACE s_axilite port=s_1     bundle=control
	#pragma HLS INTERFACE s_axilite port=s_2     bundle=control
	#pragma HLS INTERFACE s_axilite port=e       bundle=control
	#pragma HLS INTERFACE s_axilite port=pk_1    bundle=control
	#pragma HLS INTERFACE s_axilite port=pk_2    bundle=control
    #pragma HLS INTERFACE s_axilite port=v       bundle=control
	#pragma HLS INTERFACE s_axilite port=key     bundle=control
    #pragma HLS INTERFACE s_axilite port=return  bundle=control

	static hls::stream<custom_width> pk_stream_1("pk stream 1");
	#pragma HLS STREAM variable=pk_stream_1 depth=2
	static hls::stream<custom_width> pk_stream_2("pk stream 2");
	#pragma HLS STREAM variable=pk_stream_2 depth=2

	static hls::stream<unsigned char> e_1_stream_1("e_1 stream 1");
	#pragma HLS STREAM variable=e_1_stream_1 depth=2
	static hls::stream<unsigned char> e_1_stream_2("e_1 stream 2");
	#pragma HLS STREAM variable=e_1_stream_2 depth=2

	static hls::stream<custom_width> e_2_stream_1("e_2 stream 1");
	#pragma HLS STREAM variable=e_2_stream_1 depth=2
	static hls::stream<custom_width> e_2_stream_2("e_2 stream 2");
	#pragma HLS STREAM variable=e_2_stream_2 depth=2

	static hls::stream<unsigned char> e_out_stream("e_out stream");
	#pragma HLS STREAM variable=e_out_stream depth=2

	static hls::stream<b_width> b_stream_1("b stream 1");
	#pragma HLS STREAM variable=b_stream_1 depth=2
	static hls::stream<b_width> b_stream_2("b stream 2");
	#pragma HLS STREAM variable=b_stream_2 depth=2

	static hls::stream<aligned_width> pk_aligned_stream_1("pk aligned stream 1");
	#pragma HLS STREAM variable=pk_aligned_stream_1 depth=2
	static hls::stream<aligned_width> pk_aligned_stream_2("pk aligned stream 2");
	#pragma HLS STREAM variable=pk_aligned_stream_2 depth=2


	#pragma HLS dataflow

	axi_aligned_pk_reader(pk_aligned_stream_1, pk_1);
	axi_aligned_pk_reader(pk_aligned_stream_2, pk_2);

	pk_stream_resizer(pk_stream_1, pk_aligned_stream_1);
	pk_stream_resizer(pk_stream_2, pk_aligned_stream_2);

	gen_e(e_1_stream_1, e_1_stream_2, e_2_stream_1, e_2_stream_2, e_out_stream, v, key);

	compute_b_mat(b_stream_1, pk_stream_1, e_2_stream_1);
	compute_b_mat(b_stream_2, pk_stream_2, e_2_stream_2);

	compute_n_store_s(s_1, e_1_stream_1, b_stream_1);
	compute_n_store_s(s_2, e_1_stream_2, b_stream_2);

	store_e(e, e_out_stream);

}

