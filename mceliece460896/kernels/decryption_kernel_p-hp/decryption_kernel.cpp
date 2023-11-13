#ifndef __SYNTHESIS__
#include <gmp.h>
#define __gmp_const const
#endif
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#define AP_INT_MAX_W 4096
#include"ap_int.h"
#include"hls_stream.h"
#include"../../params.h"

#define SK_SIZE 13608

typedef ap_uint<13> gf;
typedef ap_uint<64> vec;
typedef ap_uint<96> syst_vec;
typedef ap_uint<13> gf_vec;
typedef ap_uint<4096> sk_4k_vec;
typedef ap_uint<2048> sk_2k_vec;
typedef ap_uint<1024> sk_1k_vec;
typedef ap_uint<512> sk_05k_vec;
typedef ap_uint<1> bit;

typedef struct{
	vec mat[GFBITS];
}mat_struct;


//TODO maybe instead of 12 iterations play with mat[0] pointers
void stream_out(hls::stream<mat_struct> &out_stream, vec out_mat[64][GFBITS]){

	int k, s, l;
	int stream_idx=0;

	mat_struct tmp_struct;

	for (k = 0; k < 32; k++)
	{
		for(l=0; l<GFBITS; l++){
			tmp_struct.mat[l] = out_mat[k][l];
		}
		out_stream.write(tmp_struct);

		for(l=0; l<GFBITS; l++){
			tmp_struct.mat[l] = out_mat[k+32][l];
		}
		out_stream.write(tmp_struct);
	}

}



//TODO read in 512B=4096b from the CPU
//Change accordingly the below
//let for later
void load_sk_vec(hls::stream<sk_4k_vec> &sk_wide_mat_stream, sk_4k_vec *sk_values){

	int i;

	//Reading in512B=4096b
	LOOP_SK:
	for(i=0; i<25; i++){
		sk_wide_mat_stream.write(*(sk_values +i));
	}
}

/* bitsliced field squarings */
void vec_sq(vec * out, vec * in)
{
	int i;
	vec result[GFBITS], t;

	t = in[11] ^ in[12];

	result[0] = in[0] ^ in[11];
	result[1] = in[7] ^ t;
	result[2] = in[1] ^ in[7];
	result[3] = in[8] ^ t;
	result[4] = in[2] ^ in[7];
	result[4] = result[4] ^ in[8];
	result[4] = result[4] ^ t;
	result[5] = in[7] ^ in[9];
	result[6] = in[3] ^ in[8];
	result[6] = result[6] ^ in[9];
	result[6] = result[6] ^ in[12];
	result[7] = in[8] ^ in[10];
	result[8] = in[4] ^ in[9];
	result[8] = result[8] ^ in[10];
	result[9] = in[9] ^ in[11];
	result[10] = in[5] ^ in[10];
	result[10] = result[10] ^ in[11];
	result[11] = in[10] ^ in[12];
	result[12] = in[6] ^ t;

	for (i = 0; i < GFBITS; i++)
		out[i] = result[i];
}


/* field multiplication */
gf gf_mul(gf in0, gf in1)
{
	int i;

	vec tmp;
	vec t0;
	vec t1;
	vec t;

	t0 = in0;
	t1 = in1;

	tmp = t0 * (t1 & 1);

	for (i = 1; i < GFBITS; i++)
		tmp ^= (t0 * (t1 & (1 << i)));

	//

	t = tmp & 0x1FF0000;
	tmp ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);

	t = tmp & 0x000E000;
	tmp ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);

	return tmp & GFMASK;
}

  vec vec_set1_16b(uint16_t v)
{
	vec ret;

	ret = (vec)v;
	ret |= ret << 16;
	ret |= ret << 32;

	return ret;
}

void transpose_64x64_benes_out(hls::stream<vec> &out_mat_stream, vec *in_mat){

	int i, j;
	vec tmp;

	LOOP_TRANSPOSE:
	for(j=0; j<64; j++){
		#pragma HLS PIPELINE II=1
		for(i=0; i<64; i++){
			tmp.bit(i) = in_mat[i].bit(j);
		}
		out_mat_stream.write(tmp);
//		printf("out[%d] = %lX\n", j, tmp.to_uint64());
	}

}


void transpose_64x64_benes_in(vec *out_mat, hls::stream<vec> &in_mat_stream){

	int i, j;
	vec tmp;

	LOOP_TRANSPOSE:
	for(j=0; j<64; j++){
		#pragma HLS PIPELINE II=1
		for(i=0; i<64; i++){
			tmp = in_mat_stream.read();
			out_mat[j].bit(i) = tmp.bit(j);
		}
	}


}


void transpose_64x64(vec *out_mat, vec *in_mat){

	int i, j;

	LOOP_TRANSPOSE:
	for(j=0; j<64; j++){
		#pragma HLS PIPELINE II=1
		for(i=0; i<64; i++){
			out_mat[j].bit(i) = in_mat[i].bit(j);
		}
	}
}


  vec vec_setbits(vec b)
{
	vec ret = -b;

	return ret;
}

  vec vec_or_reduce(hls::stream<vec> &a)
{
	int i;
	vec ret;

	ret = a.read();
	for (i = 1; i < GFBITS; i++){
		ret |= a.read();
	}

	return ret;
}

//TODO can I make simples since I have vectors?
//what is it exactly trying to perform here? (popcount?)
  int vec_testz(vec a)
{
	a |= a >> 32;
	a |= a >> 16;
	a |= a >> 8;
	a |= a >> 4;
	a |= a >> 2;
	a |= a >> 1;

	return (a&1)^1;
}

  void vec_copy(vec * out, vec * in)
{
	int i;

	for (i = 0; i < GFBITS; i++)
		out[i] = in[i];
}


  void vec_copy128(vec * out, vec * in)
{
	int i;

	for (i = 0; i < 128; i++)
		out[i] = in[i];
}


void vec_copy64(vec * out, vec * in)
{
	int i;

	for (i = 0; i < 64; i++)
		out[i] = in[i];
}

void vec_mul_bm_out(hls::stream<vec> &out_stream, const vec * f, const vec * g)
{
	int i, j;
	vec buf[ 2*GFBITS-1 ];

	#pragma HLS array_partition variable=buf

	for (i = 0; i < 2*GFBITS-1; i++){
//	#pragma HLS pipeline
	#pragma HLS unroll
		buf[i] = 0;
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS PIPELINE
		for (j = 0; j < GFBITS; j++){
			buf[i+j] ^= f[i] & g[j];
		}
	}

	for (i = 2*GFBITS-2; i >= GFBITS; i--)
	#pragma HLS pipeline
	{
		buf[i-GFBITS+4] ^= buf[i];
		buf[i-GFBITS+3] ^= buf[i];
		buf[i-GFBITS+1] ^= buf[i];
		buf[i-GFBITS+0] ^= buf[i];
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS pipeline
		out_stream.write(buf[i]);
	}
}

void vec_mul_stream_in(vec * h, hls::stream<vec> &f_stream, const vec * g)
{
	int i, j;
	vec buf[ 2*GFBITS-1 ];
	#pragma HLS array_partition variable=buf

	vec f_tmp;


	for (i = 0; i < 2*GFBITS-1; i++){
//	#pragma HLS pipeline
	#pragma HLS unroll
		buf[i] = 0;
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS PIPELINE
		f_tmp = f_stream.read();
		for (j = 0; j < GFBITS; j++){
			buf[i+j] ^= f_tmp & g[j];
		}
	}

	for (i = 2*GFBITS-2; i >= GFBITS; i--)
	#pragma HLS pipeline
	{
		buf[i-GFBITS+4] ^= buf[i];
		buf[i-GFBITS+3] ^= buf[i];
		buf[i-GFBITS+1] ^= buf[i];
		buf[i-GFBITS+0] ^= buf[i];
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS unroll
		h[i] = buf[i];
	}
}

void vec_mul_stream_2_in(vec * h, hls::stream<vec> &f_stream, hls::stream<vec> &g_stream)
{
	int i, j;
	vec buf[ 2*GFBITS-1 ];
	#pragma HLS array_partition variable=buf

	vec f_tmp;


	for (i = 0; i < 2*GFBITS-1; i++){
//	#pragma HLS pipeline
	#pragma HLS unroll
		buf[i] = 0;
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS PIPELINE
		f_tmp = f_stream.read();
		for (j = 0; j < GFBITS; j++){
			buf[i+j] ^= f_tmp & g_stream.read();
		}
	}

	for (i = 2*GFBITS-2; i >= GFBITS; i--)
	#pragma HLS pipeline
	{
		buf[i-GFBITS+4] ^= buf[i];
		buf[i-GFBITS+3] ^= buf[i];
		buf[i-GFBITS+1] ^= buf[i];
		buf[i-GFBITS+0] ^= buf[i];
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS unroll
		h[i] = buf[i];
	}
}

void vec_mul_stream_inout(hls::stream<vec> &h_stream, hls::stream<vec> &f_stream, const vec * g)
{
	int i, j;
	vec buf[ 2*GFBITS-1 ];
	#pragma HLS array_partition variable=buf

	vec f_tmp;


	for (i = 0; i < 2*GFBITS-1; i++){
//	#pragma HLS pipeline
	#pragma HLS unroll
		buf[i] = 0;
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS PIPELINE
		f_tmp = f_stream.read();
		for (j = 0; j < GFBITS; j++){
			buf[i+j] ^= f_tmp & g[j];
		}
	}

	for (i = 2*GFBITS-2; i >= GFBITS; i--)
	#pragma HLS pipeline
	{
		buf[i-GFBITS+4] ^= buf[i];
		buf[i-GFBITS+3] ^= buf[i];
		buf[i-GFBITS+1] ^= buf[i];
		buf[i-GFBITS+0] ^= buf[i];
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS unroll
		h_stream.write(buf[i]);
	}
}


void vec_mul(vec * h, const vec * f, const vec * g)
{
	int i, j;
	vec buf[ 2*GFBITS-1 ];
	#pragma HLS array_partition variable=buf


	for (i = 0; i < 2*GFBITS-1; i++){
	#pragma HLS pipeline
		buf[i] = 0;
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS PIPELINE
		for (j = 0; j < GFBITS; j++){
			buf[i+j] ^= f[i] & g[j];
		}
	}

	for (i = 2*GFBITS-2; i >= GFBITS; i--)
	#pragma HLS pipeline
	{
		buf[i-GFBITS+4] ^= buf[i];
		buf[i-GFBITS+3] ^= buf[i];
		buf[i-GFBITS+1] ^= buf[i];
		buf[i-GFBITS+0] ^= buf[i];
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS PIPELINE
		h[i] = buf[i];
	}
}



//#define vec_add(z, x, y) for (b = 0; b < GFBITS; b++) { z[b] = x[b]^y[b]; }

//TODO reconsider argument passing
void vec_add(vec *z, vec *x, vec *y){

	int b;
	for (b = 0; b < GFBITS; b++) {
		z[b] = x[b]^y[b];
	}
}

void vec_add_stream_out(hls::stream<vec> &z_stream, vec *x, vec *y){

	int b;
	for (b = 0; b < GFBITS; b++) {
		z_stream.write(x[b]^y[b]);
	}
}

void loop_function_5(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &mask_stream_0, hls::stream<vec> &mask_stream_1){

	int i;
	vec tmp_in_0 = 0;
	vec tmp_in_1 = 0;

	//but in its output it has to pass 2 streams...one for the 0 and one for the 1
	for (i = 0; i < GFBITS; i++)
	{
		//It is reading the in[0] that can be streamed directly to the output since it is not overwritten
		tmp_in_0 = in_stream_0.read();
		tmp_in_1 = in_stream_1.read();

		out_stream_0.write(tmp_in_0);
		tmp_in_1 ^= (tmp_in_0 & mask_stream_0.read()) >> 32;

		tmp_in_1 ^= (tmp_in_1 & mask_stream_1.read()) << 32;
		out_stream_1.write(tmp_in_1);
	}

}


void loop_function_4(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &mask_stream_0, hls::stream<vec> &mask_stream_1, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1){

	int i, k;

	vec in_int[2][GFBITS];

	#pragma HLS array_partition variable=in_int

	//Is this ok in HLS?
//	vec *tmp_ptr_0 = &in_int[0][0];
//	vec *tmp_ptr_1 = &in_int[1][0];

	//Does this work?
	vec_mul_stream_2_in(&in_int[0][0], in_stream_0, s_stream_0); // scaling
	vec_mul_stream_2_in(&in_int[1][0], in_stream_1, s_stream_1); // scaling

	for (i = 0; i < GFBITS; i++){
		for (k = 4; k <= 4; k++)
		{
			in_int[0][i] ^= (in_int[0][i] & mask_stream_0.read()) << (1 << 4);
			in_int[0][i] ^= (in_int[0][i] & mask_stream_1.read()) << (1 << 4);

			in_int[1][i] ^= (in_int[1][i] & mask_stream_0.read()) << (1 << 4);
			in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << (1 << 4);
		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		in_int[1][i] ^= (in_int[0][i] & mask_stream_0.read()) >> 32;
		in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << 32;

		//Write the output streams for the next function
		out_stream_0.write(in_int[0][i]);
		out_stream_1.write(in_int[1][i]);
	}


}

void loop_function_3(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &mask_stream_0, hls::stream<vec> &mask_stream_1, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1){

	int i, k;
	vec in_int[2][GFBITS];

	#pragma HLS array_partition variable=in_int

	vec_mul_stream_2_in(&in_int[0][0], in_stream_0, s_stream_0); // scaling
	vec_mul_stream_2_in(&in_int[1][0], in_stream_1, s_stream_1); // scaling

	for (i = 0; i < GFBITS; i++){
		for (k = 3; k <= 4; k++)
		{
			in_int[0][i] ^= (in_int[0][i] & mask_stream_0.read()) << (1 << k);
			in_int[0][i] ^= (in_int[0][i] & mask_stream_1.read()) << (1 << k);

			in_int[1][i] ^= (in_int[1][i] & mask_stream_0.read()) << (1 << k);
			in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << (1 << k);
		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		in_int[1][i] ^= (in_int[0][i] & mask_stream_0.read()) >> 32;
		in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << 32;
		//Write the output streams for the next function
		out_stream_0.write(in_int[0][i]);
		out_stream_1.write(in_int[1][i]);
	}

}

void loop_function_2(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &mask_stream_0, hls::stream<vec> &mask_stream_1, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1){

	int i, k;

	vec in_int[2][GFBITS];
	#pragma HLS array_partition variable=in_int

	vec_mul_stream_2_in(&in_int[0][0], in_stream_0, s_stream_0); // scaling
	vec_mul_stream_2_in(&in_int[1][0], in_stream_1, s_stream_1); // scaling

	for (i = 0; i < GFBITS; i++){
		for (k = 2; k <= 4; k++)
		{
			in_int[0][i] ^= (in_int[0][i] & mask_stream_0.read()) << (1 << k);
			in_int[0][i] ^= (in_int[0][i] & mask_stream_1.read()) << (1 << k);

			in_int[1][i] ^= (in_int[1][i] & mask_stream_0.read()) << (1 << k);
			in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << (1 << k);
		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		in_int[1][i] ^= (in_int[0][i] & mask_stream_0.read()) >> 32;
		in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << 32;
		//Write the output streams for the next function
		out_stream_0.write(in_int[0][i]);
		out_stream_1.write(in_int[1][i]);
	}


}


void loop_function_1(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &mask_stream_0, hls::stream<vec> &mask_stream_1, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1){

	int i, k;

	vec in_int[2][GFBITS];
	#pragma HLS array_partition variable=in_int

	vec_mul_stream_2_in(&in_int[0][0], in_stream_0, s_stream_0); // scaling
	vec_mul_stream_2_in(&in_int[1][0], in_stream_1, s_stream_1); // scaling

	for (i = 0; i < GFBITS; i++){
		for (k = 1; k <= 4; k++)
		{
			in_int[0][i] ^= (in_int[0][i] & mask_stream_0.read()) << (1 << k);
			in_int[0][i] ^= (in_int[0][i] & mask_stream_1.read()) << (1 << k);

			in_int[1][i] ^= (in_int[1][i] & mask_stream_0.read()) << (1 << k);
			in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << (1 << k);
		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		in_int[1][i] ^= (in_int[0][i] & mask_stream_0.read()) >> 32;
		in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << 32;
		//Write the output streams for the next function
		out_stream_0.write(in_int[0][i]);
		out_stream_1.write(in_int[1][i]);
	}

}


void loop_function_0(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &mask_stream_0, hls::stream<vec> &mask_stream_1, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1){

	int i, k;
	vec in_int[2][GFBITS];
	#pragma HLS array_partition variable=in_int

	vec_mul_stream_2_in(&in_int[0][0], in_stream_0, s_stream_0); // scaling
	vec_mul_stream_2_in(&in_int[1][0], in_stream_1, s_stream_1); // scaling

	for (i = 0; i < GFBITS; i++){
		for (k = 0; k <= 4; k++)
		{
			in_int[0][i] ^= (in_int[0][i] & mask_stream_0.read()) << (1 << k);
			in_int[0][i] ^= (in_int[0][i] & mask_stream_1.read()) << (1 << k);

			in_int[1][i] ^= (in_int[1][i] & mask_stream_0.read()) << (1 << k);
			in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << (1 << k);
		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		in_int[1][i] ^= (in_int[0][i] & mask_stream_0.read()) >> 32;
		in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << 32;
		//Write the output streams for the next function
		out_stream_0.write(in_int[0][i]);
		out_stream_1.write(in_int[1][i]);
	}

}

void loop_function_n(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &mask_stream_0, hls::stream<vec> &mask_stream_1, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1, int iter){

	int i, k;
	vec in_int[2][GFBITS];

	vec_mul_stream_2_in(&in_int[0][0], in_stream_0, s_stream_0); // scaling
	vec_mul_stream_2_in(&in_int[1][0], in_stream_1, s_stream_1); // scaling

	for (i = 0; i < GFBITS; i++){
		for (k = iter; k <= 4; k++)
		{
			in_int[0][i] ^= (in_int[0][i] & mask_stream_0.read()) << (1 << k);
			in_int[0][i] ^= (in_int[0][i] & mask_stream_1.read()) << (1 << k);

			in_int[1][i] ^= (in_int[1][i] & mask_stream_0.read()) << (1 << k);
			in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << (1 << k);
		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		in_int[1][i] ^= (in_int[0][i] & mask_stream_0.read()) >> 32;
		in_int[1][i] ^= (in_int[1][i] & mask_stream_1.read()) << 32;
		//Write the output streams for the next function
		out_stream_0.write(in_int[0][i]);
		out_stream_1.write(in_int[1][i]);
	}

}



void s_streaming(hls::stream<vec> &s_stream_4_0, hls::stream<vec> &s_stream_4_1, hls::stream<vec> &s_stream_3_0, hls::stream<vec> &s_stream_3_1, hls::stream<vec> &s_stream_2_0, hls::stream<vec> &s_stream_2_1, hls::stream<vec> &s_stream_1_0, hls::stream<vec> &s_stream_1_1, hls::stream<vec> &s_stream_0_0, hls::stream<vec> &s_stream_0_1, vec s[5][2][GFBITS]){

	int i, j;

	for(i=0; i<GFBITS; i++){
		for(j=0; j<GFBITS; j++){
			s_stream_4_0.write(s[4][0][j]);
			s_stream_4_1.write(s[4][1][j]);

			s_stream_3_0.write(s[3][0][j]);
			s_stream_3_1.write(s[3][1][j]);

			s_stream_2_0.write(s[2][0][j]);
			s_stream_2_1.write(s[2][1][j]);

			s_stream_1_0.write(s[1][0][j]);
			s_stream_1_1.write(s[1][1][j]);

			s_stream_0_0.write(s[0][0][j]);
			s_stream_0_1.write(s[0][1][j]);
		}
	}

}

void mask_streaming(hls::stream<vec> &mask_stream_5_0, hls::stream<vec> &mask_stream_5_1, hls::stream<vec> &mask_stream_4_0, hls::stream<vec> &mask_stream_4_1, hls::stream<vec> &mask_stream_3_0, hls::stream<vec> &mask_stream_3_1, hls::stream<vec> &mask_stream_2_0, hls::stream<vec> &mask_stream_2_1, hls::stream<vec> &mask_stream_1_0, hls::stream<vec> &mask_stream_1_1, hls::stream<vec> &mask_stream_0_0, hls::stream<vec> &mask_stream_0_1, vec mask[6][2]){

	int i, j;

	//TODO find more effcient way to do it
	for(j=0; j<GFBITS; j++){

		mask_stream_5_0.write(mask[5][0]);
		mask_stream_5_1.write(mask[5][1]);

		mask_stream_4_0.write(mask[4][0]);
		mask_stream_4_1.write(mask[4][1]);
		mask_stream_4_0.write(mask[4][0]);
		mask_stream_4_1.write(mask[4][1]);


		mask_stream_3_0.write(mask[3][0]);
		mask_stream_3_1.write(mask[3][1]);
		mask_stream_3_0.write(mask[3][0]);
		mask_stream_3_1.write(mask[3][1]);
		mask_stream_3_0.write(mask[4][0]);
		mask_stream_3_1.write(mask[4][1]);
		mask_stream_3_0.write(mask[4][0]);
		mask_stream_3_1.write(mask[4][1]);


		mask_stream_2_0.write(mask[2][0]);
		mask_stream_2_1.write(mask[2][1]);
		mask_stream_2_0.write(mask[2][0]);
		mask_stream_2_1.write(mask[2][1]);
		mask_stream_2_0.write(mask[3][0]);
		mask_stream_2_1.write(mask[3][1]);
		mask_stream_2_0.write(mask[3][0]);
		mask_stream_2_1.write(mask[3][1]);
		mask_stream_2_0.write(mask[4][0]);
		mask_stream_2_1.write(mask[4][1]);
		mask_stream_2_0.write(mask[4][0]);
		mask_stream_2_1.write(mask[4][1]);


		mask_stream_1_0.write(mask[1][0]);
		mask_stream_1_1.write(mask[1][1]);
		mask_stream_1_0.write(mask[1][0]);
		mask_stream_1_1.write(mask[1][1]);
		mask_stream_1_0.write(mask[2][0]);
		mask_stream_1_1.write(mask[2][1]);
		mask_stream_1_0.write(mask[2][0]);
		mask_stream_1_1.write(mask[2][1]);
		mask_stream_1_0.write(mask[3][0]);
		mask_stream_1_1.write(mask[3][1]);
		mask_stream_1_0.write(mask[3][0]);
		mask_stream_1_1.write(mask[3][1]);
		mask_stream_1_0.write(mask[4][0]);
		mask_stream_1_1.write(mask[4][1]);
		mask_stream_1_0.write(mask[4][0]);
		mask_stream_1_1.write(mask[4][1]);


		mask_stream_0_0.write(mask[0][0]);
		mask_stream_0_1.write(mask[0][1]);
		mask_stream_0_0.write(mask[0][0]);
		mask_stream_0_1.write(mask[0][1]);
		mask_stream_0_0.write(mask[1][0]);
		mask_stream_0_1.write(mask[1][1]);
		mask_stream_0_0.write(mask[1][0]);
		mask_stream_0_1.write(mask[1][1]);
		mask_stream_0_0.write(mask[2][0]);
		mask_stream_0_1.write(mask[2][1]);
		mask_stream_0_0.write(mask[2][0]);
		mask_stream_0_1.write(mask[2][1]);
		mask_stream_0_0.write(mask[3][0]);
		mask_stream_0_1.write(mask[3][1]);
		mask_stream_0_0.write(mask[3][0]);
		mask_stream_0_1.write(mask[3][1]);
		mask_stream_0_0.write(mask[4][0]);
		mask_stream_0_1.write(mask[4][1]);
		mask_stream_0_0.write(mask[4][0]);
		mask_stream_0_1.write(mask[4][1]);



	}

	for(j=0; j<GFBITS; j++){
		mask_stream_4_0.write(mask[5][0]);
		mask_stream_4_1.write(mask[5][1]);

		mask_stream_3_0.write(mask[5][0]);
		mask_stream_3_1.write(mask[5][1]);

		mask_stream_2_0.write(mask[5][0]);
		mask_stream_2_1.write(mask[5][1]);

		mask_stream_1_0.write(mask[5][0]);
		mask_stream_1_1.write(mask[5][1]);

		mask_stream_0_0.write(mask[5][0]);
		mask_stream_0_1.write(mask[5][1]);

	}

}


//TODO has to read the wing streams and write accordingly to the output
void butterfly_tr_operation(mat_struct *op_a, mat_struct *op_b, const vec consts[128][GFBITS], unsigned int index){

	vec tmp[GFBITS];
	//TODO apply dataflow
	//Maybe pass as streams? This will lower the FF too..current is 67 cycles per layer
	vec_add(op_a->mat, op_a->mat, op_b->mat);
	vec_mul(tmp, op_a->mat, consts[index]);
	vec_add(op_b->mat, op_b->mat, tmp);


}


static void fft_tr_layer(hls::stream<mat_struct> &sink_0, hls::stream<mat_struct> &sink_1, hls::stream<mat_struct> &src_0, hls::stream<mat_struct> &src_1, int start, int reps)
{

	const vec consts[ 128 ][ GFBITS ] =
	{
#include "consts.data"
	};

	mat_struct wing1;
	mat_struct wing2;

	int iter=0;

	unsigned int index = start;

	for (unsigned int i = 0; i < 64; i++) {
		#pragma HLS pipeline
		wing1 = src_0.read();
		wing2 = src_1.read();

//		printf("Consts index=%d\n", index);
		//TODO apply dataflow inside the butterfly?
		butterfly_tr_operation(&wing1, &wing2, consts, index);

		if (i % 64 < 32) {
			sink_0.write(wing1);
			sink_0.write(wing2);

		}
		else {
			sink_1.write(wing1);
			sink_1.write(wing2);
		}


		// handle the consts matrix index
		if (++iter == reps) {
			iter = 0;
			if (++index == 2 * start + 1) {
				index = start;
			}
		}

	}

}

void split_struct_stream(hls::stream<mat_struct> &out_stream_0, hls::stream<mat_struct> &out_stream_1, hls::stream<mat_struct> &in_struct_stream){

	int i;
	mat_struct tmp;

	for(i=0; i<128; i++){
			tmp = in_struct_stream.read();
			if(i%2==0){
				out_stream_0.write(tmp);
			}else{
				out_stream_1.write(tmp);
			}
	}
}

void split_struct_stream_tr(hls::stream<mat_struct> &out_stream_0, hls::stream<mat_struct> &out_stream_1, hls::stream<mat_struct> &in_struct_stream){

	int i;
	mat_struct tmp;

	for(i=0; i<128; i++){
			tmp = in_struct_stream.read();
			if(i<64){
				out_stream_0.write(tmp);
			}else{
				out_stream_1.write(tmp);
			}
	}
}

void merge_tr_struct_stream(hls::stream<vec> &out_stream, hls::stream<mat_struct> &in_stream_0, hls::stream<mat_struct> &in_stream_1){

	int i, j;
	mat_struct tmp_0, tmp_1;

	int idx_0, idx_1;

	vec out_mat[128][GFBITS];

	#pragma HLS array_partition variable=out_mat dim=2 complete

	for(i=0; i<128; i+=2){
			tmp_0 = in_stream_0.read();
			tmp_1 = in_stream_1.read();
			for(j=0; j<GFBITS; j++){
				out_mat[i][j] = tmp_0.mat[j];
				out_mat[i+1][j] = tmp_1.mat[j];
			}
	}

	//TODO make more efficient
	//write according to the reversal matrix
	for(i=0; i<128; i++){
	#pragma HLS pipeline
		for(j=0; j<GFBITS; j++){
			out_stream.write(out_mat[i][j]);
//			out_stream_1.write(out_mat[j+64][i]);
		}
	}

}



//TODO has to read the wing streams and write accordingly to the output
void butterfly_operation(mat_struct *op_a, mat_struct *op_b, vec consts[128][GFBITS], int index){

	vec tmp[GFBITS];
	//TODO apply dataflow
	//Maybe pass as streams?
	vec_mul(tmp, op_b->mat, consts[index]);
	vec_add(op_a->mat, op_a->mat, tmp);
	vec_add(op_b->mat, op_b->mat, op_a->mat);



}


static void fft_layer(hls::stream<mat_struct> &sink_0, hls::stream<mat_struct> &sink_1, hls::stream<mat_struct> &src_0, hls::stream<mat_struct> &src_1, int start, int reps)
{

	vec consts[ 128 ][ GFBITS ] =
	{
#include "consts.data"
	};

	mat_struct wing1;
	mat_struct wing2;

	vec mat_a[GFBITS];
	vec mat_b[GFBITS];

	int idx_0;
	int idx_1;
	int iter=0;

	unsigned int index = start;

	for (unsigned int i = 0; i < 64; i++) {
		#pragma HLS pipeline

//		printf("Consts index=%d\n", index);

		if (i % 64 < 32) {
			wing1 = src_0.read();
			wing2 = src_0.read();

		}
		else {
			wing1 = src_1.read();
			wing2 = src_1.read();
		}

		//TODO apply dataflow inside the butterfly?
		butterfly_operation(&wing1, &wing2, consts, index);

		sink_0.write(wing1);
		sink_1.write(wing2);


		// handle the consts matrix index
		if (++iter == reps) {
			iter = 0;
			if (++index == 2 * start + 1) {
				index = start;
			}
		}


	}

}


void broadcast(hls::stream<mat_struct> &out_struct_stream,  hls::stream<vec> &in_stream){



	int i, j, k, s, b;

	vec tmp[ GFBITS ];
	vec pre[8][ GFBITS ];
	vec buf[128], buf_tr[128];

	vec consts_ptr = 2;

	vec consts[ 128 ][ GFBITS ] =
	{
#include "consts.data"
	};

	unsigned char reversal[128] =
	{
	  0, 64, 32, 96, 16, 80, 48, 112,
	  8, 72, 40, 104, 24, 88, 56, 120,
	  4, 68, 36, 100, 20, 84, 52, 116,
	  12, 76, 44, 108, 28, 92, 60, 124,
	  2, 66, 34, 98, 18, 82, 50, 114,
	  10, 74, 42, 106, 26, 90, 58, 122,
	  6, 70, 38, 102, 22, 86, 54, 118,
	  14, 78, 46, 110, 30, 94, 62, 126,
	  1, 65, 33, 97, 17, 81, 49, 113,
	  9, 73, 41, 105, 25, 89, 57, 121,
	  5, 69, 37, 101, 21, 85, 53, 117,
	  13, 77, 45, 109, 29, 93, 61, 125,
	  3, 67, 35, 99, 19, 83, 51, 115,
	  11, 75, 43, 107, 27, 91, 59, 123,
	  7, 71, 39, 103, 23, 87, 55, 119,
	  15, 79, 47, 111, 31, 95, 63, 127
	};

	ap_uint<16> beta[7] = {2522, 7827, 7801, 8035, 6897, 8167, 3476};

	vec in[2][GFBITS];
	vec out[128][GFBITS];


	//Read the matrix 0
	for (i = 0; i < GFBITS; i++)
	{
		in[0][i] = in_stream.read();
//		printf("in[0][%d] = %lX\n", i, in[0][i].to_uint64());
	}

	//Read the matrix 1
	for (i = 0; i < GFBITS; i++)
	{
		in[1][i] = in_stream.read();
//		printf("in[0][%d] = %lX\n", i, in[1][i].to_uint64());
	}

	for (i = 0; i < 7; i++)
	{
		for (j = 0; j < GFBITS; j++)
		{
			pre[i][j] = (beta[i] >> j) & 1;
			pre[i][j] = -pre[i][j];
		}

		vec_mul(pre[i], in[1], pre[i]);
	}

	for (i = 0; i < GFBITS; i++)
	{
		buf[0] = in[0][i];

		buf[1] = buf[0] ^ pre[0][i];      buf[32] = in[0][i] ^ pre[5][i];
		buf[3] = buf[1] ^ pre[1][i];      buf[96] = buf[32] ^ pre[6][i];
                                              buf[97] = buf[96] ^ pre[0][i];
		buf[2] = in[0][i] ^ pre[1][i];  buf[99] = buf[97] ^ pre[1][i];
		buf[6] = buf[2] ^ pre[2][i];      buf[98] = buf[99] ^ pre[0][i];
		buf[7] = buf[6] ^ pre[0][i];      buf[102] = buf[98] ^ pre[2][i];
		buf[5] = buf[7] ^ pre[1][i];      buf[103] = buf[102] ^ pre[0][i];
                                              buf[101] = buf[103] ^ pre[1][i];
		buf[4] = in[0][i] ^ pre[2][i];  buf[100] = buf[101] ^ pre[0][i];
		buf[12] = buf[4] ^ pre[3][i];     buf[108] = buf[100] ^ pre[3][i];
		buf[13] = buf[12] ^ pre[0][i];    buf[109] = buf[108] ^ pre[0][i];
		buf[15] = buf[13] ^ pre[1][i];    buf[111] = buf[109] ^ pre[1][i];
		buf[14] = buf[15] ^ pre[0][i];    buf[110] = buf[111] ^ pre[0][i];
		buf[10] = buf[14] ^ pre[2][i];    buf[106] = buf[110] ^ pre[2][i];
		buf[11] = buf[10] ^ pre[0][i];    buf[107] = buf[106] ^ pre[0][i];
		buf[9] = buf[11] ^ pre[1][i];     buf[105] = buf[107] ^ pre[1][i];
                                              buf[104] = buf[105] ^ pre[0][i];
		buf[8] = in[0][i] ^ pre[3][i];  buf[120] = buf[104] ^ pre[4][i];
		buf[24] = buf[8] ^ pre[4][i];     buf[121] = buf[120] ^ pre[0][i];
		buf[25] = buf[24] ^ pre[0][i];    buf[123] = buf[121] ^ pre[1][i];
		buf[27] = buf[25] ^ pre[1][i];    buf[122] = buf[123] ^ pre[0][i];
		buf[26] = buf[27] ^ pre[0][i];    buf[126] = buf[122] ^ pre[2][i];
		buf[30] = buf[26] ^ pre[2][i];    buf[127] = buf[126] ^ pre[0][i];
		buf[31] = buf[30] ^ pre[0][i];    buf[125] = buf[127] ^ pre[1][i];
		buf[29] = buf[31] ^ pre[1][i];    buf[124] = buf[125] ^ pre[0][i];
		buf[28] = buf[29] ^ pre[0][i];    buf[116] = buf[124] ^ pre[3][i];
		buf[20] = buf[28] ^ pre[3][i];    buf[117] = buf[116] ^ pre[0][i];
		buf[21] = buf[20] ^ pre[0][i];    buf[119] = buf[117] ^ pre[1][i];
		buf[23] = buf[21] ^ pre[1][i];    buf[118] = buf[119] ^ pre[0][i];
		buf[22] = buf[23] ^ pre[0][i];    buf[114] = buf[118] ^ pre[2][i];
		buf[18] = buf[22] ^ pre[2][i];    buf[115] = buf[114] ^ pre[0][i];
		buf[19] = buf[18] ^ pre[0][i];    buf[113] = buf[115] ^ pre[1][i];
		buf[17] = buf[19] ^ pre[1][i];    buf[112] = buf[113] ^ pre[0][i];
                                              buf[80] = buf[112] ^ pre[5][i];
		buf[16] = in[0][i] ^ pre[4][i]; buf[81] = buf[80] ^ pre[0][i];
		buf[48] = buf[16] ^ pre[5][i];    buf[83] = buf[81] ^ pre[1][i];
		buf[49] = buf[48] ^ pre[0][i];    buf[82] = buf[83] ^ pre[0][i];
		buf[51] = buf[49] ^ pre[1][i];    buf[86] = buf[82] ^ pre[2][i];
		buf[50] = buf[51] ^ pre[0][i];    buf[87] = buf[86] ^ pre[0][i];
		buf[54] = buf[50] ^ pre[2][i];    buf[85] = buf[87] ^ pre[1][i];
		buf[55] = buf[54] ^ pre[0][i];    buf[84] = buf[85] ^ pre[0][i];
		buf[53] = buf[55] ^ pre[1][i];    buf[92] = buf[84] ^ pre[3][i];
		buf[52] = buf[53] ^ pre[0][i];    buf[93] = buf[92] ^ pre[0][i];
		buf[60] = buf[52] ^ pre[3][i];    buf[95] = buf[93] ^ pre[1][i];
		buf[61] = buf[60] ^ pre[0][i];    buf[94] = buf[95] ^ pre[0][i];
		buf[63] = buf[61] ^ pre[1][i];    buf[90] = buf[94] ^ pre[2][i];
		buf[62] = buf[63] ^ pre[0][i];    buf[91] = buf[90] ^ pre[0][i];
		buf[58] = buf[62] ^ pre[2][i];    buf[89] = buf[91] ^ pre[1][i];
		buf[59] = buf[58] ^ pre[0][i];    buf[88] = buf[89] ^ pre[0][i];
		buf[57] = buf[59] ^ pre[1][i];    buf[72] = buf[88] ^ pre[4][i];
		buf[56] = buf[57] ^ pre[0][i];    buf[73] = buf[72] ^ pre[0][i];
		buf[40] = buf[56] ^ pre[4][i];    buf[75] = buf[73] ^ pre[1][i];
		buf[41] = buf[40] ^ pre[0][i];    buf[74] = buf[75] ^ pre[0][i];
		buf[43] = buf[41] ^ pre[1][i];    buf[78] = buf[74] ^ pre[2][i];
		buf[42] = buf[43] ^ pre[0][i];    buf[79] = buf[78] ^ pre[0][i];
		buf[46] = buf[42] ^ pre[2][i];    buf[77] = buf[79] ^ pre[1][i];
		buf[47] = buf[46] ^ pre[0][i];    buf[76] = buf[77] ^ pre[0][i];
		buf[45] = buf[47] ^ pre[1][i];    buf[68] = buf[76] ^ pre[3][i];
		buf[44] = buf[45] ^ pre[0][i];    buf[69] = buf[68] ^ pre[0][i];
		buf[36] = buf[44] ^ pre[3][i];    buf[71] = buf[69] ^ pre[1][i];
		buf[37] = buf[36] ^ pre[0][i];    buf[70] = buf[71] ^ pre[0][i];
		buf[39] = buf[37] ^ pre[1][i];    buf[66] = buf[70] ^ pre[2][i];
		buf[38] = buf[39] ^ pre[0][i];    buf[67] = buf[66] ^ pre[0][i];
		buf[34] = buf[38] ^ pre[2][i];    buf[65] = buf[67] ^ pre[1][i];
		buf[35] = buf[34] ^ pre[0][i];
		buf[33] = buf[35] ^ pre[1][i];    buf[64] = in[0][i] ^ pre[6][i];

		transpose_64x64(buf_tr +  0, buf +  0);
		transpose_64x64(buf_tr + 64, buf + 64);

		//This was added by me since transpse needs different iout and outout to function correctly
		vec_copy128(buf, buf_tr);


		for (j = 0; j < 128; j++){
			out[ reversal[j] ][i] = buf[j];
		}

	}

	for (i = 0; i < 128; i++){
		mat_struct tmp_struct;
		for (j = 0; j < GFBITS; j++){
			tmp_struct.mat[j] = out[i][j];
		}
		out_struct_stream.write(tmp_struct);
	}



}


void transpose_broadcast(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1,  hls::stream<vec> &out_stream_2,  hls::stream<vec> &out_stream_3, hls::stream<vec> &in_stream){

	int i, j, k, s, b;

		vec tmp[ GFBITS ];
		vec pre[6][2][ GFBITS ];
		vec buf[2][64];
		vec buf_tr[2][64];

		const vec consts[ 128 ][ GFBITS ] =
		{
	#include "consts.data"
		};

		vec consts_ptr = 128;

		const unsigned char reversal[128] =
		{
		  0, 64, 32, 96, 16, 80, 48, 112,
		  8, 72, 40, 104, 24, 88, 56, 120,
		  4, 68, 36, 100, 20, 84, 52, 116,
		  12, 76, 44, 108, 28, 92, 60, 124,
		  2, 66, 34, 98, 18, 82, 50, 114,
		  10, 74, 42, 106, 26, 90, 58, 122,
		  6, 70, 38, 102, 22, 86, 54, 118,
		  14, 78, 46, 110, 30, 94, 62, 126,
		  1, 65, 33, 97, 17, 81, 49, 113,
		  9, 73, 41, 105, 25, 89, 57, 121,
		  5, 69, 37, 101, 21, 85, 53, 117,
		  13, 77, 45, 109, 29, 93, 61, 125,
		  3, 67, 35, 99, 19, 83, 51, 115,
		  11, 75, 43, 107, 27, 91, 59, 123,
		  7, 71, 39, 103, 23, 87, 55, 119,
		  15, 79, 47, 111, 31, 95, 63, 127
		};

		const uint16_t beta[6] = {5246, 5306, 6039, 6685, 4905, 6755};

		mat_struct tmp_struct;

		vec in[ 128 ][ GFBITS ];
		vec out[ 4 ][ GFBITS ];

		for (i = 0; i < 128; i++){
			for (j = 0; j < GFBITS; j++){
				in[i][j] = in_stream.read();
			}
		}

	for (i = 0; i < GFBITS; i++)
	{
		for (k = 0; k < 128; k++)
			(&buf[0][0])[ k ] = in[ reversal[k] ][i];

		transpose_64x64(buf_tr[0], buf[0]);
		transpose_64x64(buf_tr[1], buf[1]);

		vec_copy64(buf[0], buf_tr[0]);
		vec_copy64(buf[1], buf_tr[1]);

		for (k = 0; k < 2; k++)
		{
			pre[0][k][i] = buf[k][32]; buf[k][33] ^= buf[k][32];
			pre[1][k][i] = buf[k][33]; buf[k][35] ^= buf[k][33];
			pre[0][k][i] ^= buf[k][35]; buf[k][34] ^= buf[k][35];
			pre[2][k][i] = buf[k][34]; buf[k][38] ^= buf[k][34];
			pre[0][k][i] ^= buf[k][38]; buf[k][39] ^= buf[k][38];
			pre[1][k][i] ^= buf[k][39]; buf[k][37] ^= buf[k][39];
			pre[0][k][i] ^= buf[k][37]; buf[k][36] ^= buf[k][37];
			pre[3][k][i] = buf[k][36]; buf[k][44] ^= buf[k][36];
			pre[0][k][i] ^= buf[k][44]; buf[k][45] ^= buf[k][44];
			pre[1][k][i] ^= buf[k][45]; buf[k][47] ^= buf[k][45];
			pre[0][k][i] ^= buf[k][47]; buf[k][46] ^= buf[k][47];
			pre[2][k][i] ^= buf[k][46]; buf[k][42] ^= buf[k][46];
			pre[0][k][i] ^= buf[k][42]; buf[k][43] ^= buf[k][42];
			pre[1][k][i] ^= buf[k][43]; buf[k][41] ^= buf[k][43];
			pre[0][k][i] ^= buf[k][41]; buf[k][40] ^= buf[k][41];
			pre[4][k][i] = buf[k][40]; buf[k][56] ^= buf[k][40];
			pre[0][k][i] ^= buf[k][56]; buf[k][57] ^= buf[k][56];
			pre[1][k][i] ^= buf[k][57]; buf[k][59] ^= buf[k][57];
			pre[0][k][i] ^= buf[k][59]; buf[k][58] ^= buf[k][59];
			pre[2][k][i] ^= buf[k][58]; buf[k][62] ^= buf[k][58];
			pre[0][k][i] ^= buf[k][62]; buf[k][63] ^= buf[k][62];
			pre[1][k][i] ^= buf[k][63]; buf[k][61] ^= buf[k][63];
			pre[0][k][i] ^= buf[k][61]; buf[k][60] ^= buf[k][61];
			pre[3][k][i] ^= buf[k][60]; buf[k][52] ^= buf[k][60];
			pre[0][k][i] ^= buf[k][52]; buf[k][53] ^= buf[k][52];
			pre[1][k][i] ^= buf[k][53]; buf[k][55] ^= buf[k][53];
			pre[0][k][i] ^= buf[k][55]; buf[k][54] ^= buf[k][55];
			pre[2][k][i] ^= buf[k][54]; buf[k][50] ^= buf[k][54];
			pre[0][k][i] ^= buf[k][50]; buf[k][51] ^= buf[k][50];
			pre[1][k][i] ^= buf[k][51]; buf[k][49] ^= buf[k][51];
			pre[0][k][i] ^= buf[k][49]; buf[k][48] ^= buf[k][49];
			pre[5][k][i] = buf[k][48]; buf[k][16] ^= buf[k][48];
			pre[0][k][i] ^= buf[k][16]; buf[k][17] ^= buf[k][16];
			pre[1][k][i] ^= buf[k][17]; buf[k][19] ^= buf[k][17];
			pre[0][k][i] ^= buf[k][19]; buf[k][18] ^= buf[k][19];
			pre[2][k][i] ^= buf[k][18]; buf[k][22] ^= buf[k][18];
			pre[0][k][i] ^= buf[k][22]; buf[k][23] ^= buf[k][22];
			pre[1][k][i] ^= buf[k][23]; buf[k][21] ^= buf[k][23];
			pre[0][k][i] ^= buf[k][21]; buf[k][20] ^= buf[k][21];
			pre[3][k][i] ^= buf[k][20]; buf[k][28] ^= buf[k][20];
			pre[0][k][i] ^= buf[k][28]; buf[k][29] ^= buf[k][28];
			pre[1][k][i] ^= buf[k][29]; buf[k][31] ^= buf[k][29];
			pre[0][k][i] ^= buf[k][31]; buf[k][30] ^= buf[k][31];
			pre[2][k][i] ^= buf[k][30]; buf[k][26] ^= buf[k][30];
			pre[0][k][i] ^= buf[k][26]; buf[k][27] ^= buf[k][26];
			pre[1][k][i] ^= buf[k][27]; buf[k][25] ^= buf[k][27];
			pre[0][k][i] ^= buf[k][25]; buf[k][24] ^= buf[k][25];
			pre[4][k][i] ^= buf[k][24]; buf[k][8] ^= buf[k][24];
			pre[0][k][i] ^= buf[k][8]; buf[k][9] ^= buf[k][8];
			pre[1][k][i] ^= buf[k][9]; buf[k][11] ^= buf[k][9];
			pre[0][k][i] ^= buf[k][11]; buf[k][10] ^= buf[k][11];
			pre[2][k][i] ^= buf[k][10]; buf[k][14] ^= buf[k][10];
			pre[0][k][i] ^= buf[k][14]; buf[k][15] ^= buf[k][14];
			pre[1][k][i] ^= buf[k][15]; buf[k][13] ^= buf[k][15];
			pre[0][k][i] ^= buf[k][13]; buf[k][12] ^= buf[k][13];
			pre[3][k][i] ^= buf[k][12]; buf[k][4] ^= buf[k][12];
			pre[0][k][i] ^= buf[k][4]; buf[k][5] ^= buf[k][4];
			pre[1][k][i] ^= buf[k][5]; buf[k][7] ^= buf[k][5];
			pre[0][k][i] ^= buf[k][7]; buf[k][6] ^= buf[k][7];
			pre[2][k][i] ^= buf[k][6]; buf[k][2] ^= buf[k][6];
			pre[0][k][i] ^= buf[k][2]; buf[k][3] ^= buf[k][2];
			pre[1][k][i] ^= buf[k][3]; buf[k][1] ^= buf[k][3];

			pre[0][k][i] ^= buf[k][1]; out[k][i] = buf[k][0] ^ buf[k][1];
		}
	}

	for (j = 0; j < GFBITS; j++){
		tmp[j] = vec_setbits((beta[0] >> j) & 1);
	}

	vec_mul(out[2], pre[0][0], tmp);
	vec_mul(out[3], pre[0][1], tmp);

	for (i = 1; i < 6; i++)
	{
		for (j = 0; j < GFBITS; j++){
			tmp[j] = vec_setbits((beta[i] >> j) & 1);
		}

		vec_mul(pre[i][0], pre[i][0], tmp);
		vec_mul(pre[i][1], pre[i][1], tmp);

		for (b = 0; b < GFBITS; b++)
		{
			out[2][b] ^= pre[i][0][b];
			out[3][b] ^= pre[i][1][b];
		}
	}

	//Write the 4 output streams
	for (i = 0; i < GFBITS; i++){
		out_stream_0.write(out[0][i]);
		out_stream_1.write(out[1][i]);
		out_stream_2.write(out[2][i]);
		out_stream_3.write(out[3][i]);
//		printf("out[0][%d] = %lX\n", i, out[0][i].to_uint64());
//		printf("out[1][%d] = %lX\n", i, out[1][i].to_uint64());
//		printf("out[2][%d] = %lX\n", i, out[2][i].to_uint64());
//		printf("out[3][%d] = %lX\n", i, out[3][i].to_uint64());

	}




}




static inline uint16_t mask_nonzero(gf a)
{
	uint32_t ret = a;

	ret -= 1;
	ret >>= 31;
	ret -= 1;

	return ret;
}

static inline uint16_t mask_leq(uint16_t a, uint16_t b)
{
	uint32_t a_tmp = a;
	uint32_t b_tmp = b;
	uint32_t ret = b_tmp - a_tmp;

	ret >>= 31;
	ret -= 1;

	return ret;
}

  void vec_cmov(vec *out, vec *in, uint16_t mask)
{
	int i;

	vec m0, m1;

	m0 = vec_set1_16b(mask);
	m1 = ~m0;

	for (i = 0; i < GFBITS; i++)
	{
		out[i] = (in[i] & m0) | (out[i] & m1);
		out[i] = (in[i] & m0) | (out[i] & m1);
	}
}

  void interleave(vec *in, int idx0, int idx1, vec *mask, int b)
{
	int s = 1 << b;

	vec x, y;

	x = (in[idx0] & mask[0]) | ((in[idx1] & mask[0]) << s);
	y = ((in[idx0] & mask[1]) >> s) | (in[idx1] & mask[1]);

	in[idx0] = x;
	in[idx1] = y;
}


  void write_buf_value(gf *out, vec buf[16], int idx){

  	int k;

  	for (k = 0; k <  4; k++){
  	#pragma HLS PIPELINE
  		out[ k*16 + idx] = buf[idx].range((k+1)*16-4, k*16);
  	}

  }


/* input: in, field elements in bitsliced form */
/* output: out, field elements in non-bitsliced form */
void get_coefs(gf *out, hls::stream<vec> &in_stream)
{
	int i, k;

	vec mask[4][2];
	vec buf[16];
	#pragma HLS array_partition variable=buf

  	for (i =  0; i < GFBITS; i++){
  		buf[i] = in_stream.read();
  	}

  	for (i = GFBITS; i < 16; i++){
  		buf[i] = 0;
  	}

  	mask[0][0] = 0x5555555555555555;

  	mask[0][1] = 0xAAAAAAAAAAAAAAAA;

  	mask[1][0] = 0x3333333333333333;

  	mask[1][1] = 0xCCCCCCCCCCCCCCCC;

  	mask[2][0] = 0x0F0F0F0F0F0F0F0F;

  	mask[2][1] = 0xF0F0F0F0F0F0F0F0;

  	mask[3][0] = 0x00FF00FF00FF00FF;

  	mask[3][1] = 0xFF00FF00FF00FF00;

  	interleave(buf,  0,  8, mask[3], 3);
  	interleave(buf,  1,  9, mask[3], 3);
  	interleave(buf,  2, 10, mask[3], 3);
  	interleave(buf,  3, 11, mask[3], 3);
  	interleave(buf,  4, 12, mask[3], 3);
  	interleave(buf,  5, 13, mask[3], 3);
  	interleave(buf,  6, 14, mask[3], 3);
  	interleave(buf,  7, 15, mask[3], 3);

  	interleave(buf,  0,  4, mask[2], 2);
  	interleave(buf,  1,  5, mask[2], 2);
  	interleave(buf,  2,  6, mask[2], 2);
  	interleave(buf,  3,  7, mask[2], 2);
  	interleave(buf,  8, 12, mask[2], 2);
  	interleave(buf,  9, 13, mask[2], 2);
  	interleave(buf, 10, 14, mask[2], 2);
  	interleave(buf, 11, 15, mask[2], 2);

  	interleave(buf,  0,  2, mask[1], 1);
  	interleave(buf,  1,  3, mask[1], 1);
  	interleave(buf,  4,  6, mask[1], 1);
  	interleave(buf,  5,  7, mask[1], 1);
  	interleave(buf,  8, 10, mask[1], 1);
  	interleave(buf,  9, 11, mask[1], 1);
  	interleave(buf, 12, 14, mask[1], 1);
  	interleave(buf, 13, 15, mask[1], 1);

  	interleave(buf,  0,  1, mask[0], 0);
  	write_buf_value(out, buf, 0);
  	write_buf_value(out, buf, 1);

  	interleave(buf,  2,  3, mask[0], 0);
  	write_buf_value(out, buf, 2);
  	write_buf_value(out, buf, 3);

  	interleave(buf,  4,  5, mask[0], 0);
  	write_buf_value(out, buf, 4);
  	write_buf_value(out, buf, 5);

  	interleave(buf,  6,  7, mask[0], 0);
  	write_buf_value(out, buf, 6);
  	write_buf_value(out, buf, 7);

  	interleave(buf,  8,  9, mask[0], 0);
  	write_buf_value(out, buf, 8);
  	write_buf_value(out, buf, 9);

  	interleave(buf, 10, 11, mask[0], 0);
  	write_buf_value(out, buf, 10);
  	write_buf_value(out, buf, 11);

  	interleave(buf, 12, 13, mask[0], 0);
  	write_buf_value(out, buf, 12);
  	write_buf_value(out, buf, 13);

  	interleave(buf, 14, 15, mask[0], 0);
  	write_buf_value(out, buf, 14);
  	write_buf_value(out, buf, 15);
}


//gf vec_reduce(vec in[2][GFBITS])
//{
//	int i;
//	vec tmp;
//	gf ret = 0;
//
////Give the d_vec immediately a value here?
//
//	for (i = GFBITS-1; i >= 0; i--)
//	#pragma HLS pipeline
//	{
//		tmp = in[0][i] ^ in[1][i];
//		ret.bit(i) = tmp.xor_reduce();
//	}
//
//	return ret;
//}

gf vec_reduce(vec in[2][GFBITS])
{
	int i;
	vec tmp;
	gf ret = 0;

	for (i = GFBITS-1; i >= 0; i--)
	{
		tmp = in[0][i] ^ in[1][i];
		ret.bit(i) = tmp.xor_reduce();

//		tmp ^= tmp >> 32;
//		tmp ^= tmp >> 16;
//		tmp ^= tmp >> 8;
//		tmp ^= tmp >> 4;
//		tmp ^= tmp >> 2;
//		tmp ^= tmp >> 1;
//
//		ret <<= 1;
//		ret |= tmp & 1;
	}

	return ret;
}

void update(vec in[2][GFBITS], const gf e)
{
	int i;
	vec tmp;

	for (i = 0; i < GFBITS; i++)
	{
//		tmp=0;
		tmp = (e >> i) & 1;
//		tmp.bit(0) = e.bit(i);

		in[0][i] = (in[0][i] >> 1) | (in[1][i] << 63);
		in[1][i] = (in[1][i] >> 1) | (tmp      << 63);
	}
}


void store_to_stream(hls::stream<vec> &out_stream, vec *tmp){

	int i, j;

		for (i = 0; i < GFBITS; i++){
			out_stream.write(tmp[i]);
//			printf("out[%d] = %lX\n", i, tmp[i].to_uint64());
		}

}

/* input: in, sequence of field elements */
/* output: out, minimal polynomial of in */
//void bm(vec out[][ GFBITS ], vec in[][ GFBITS ])
void bm(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &in_stream_2, hls::stream<vec> &in_stream_3)
{
	#pragma HLS Inline recursive

	int i;
	ap_uint<16> N, L;
	ap_uint<16> mask, mask_1, mask_2;
	vec one = 0x1;



	vec prod[2][GFBITS];
	vec interval[2][GFBITS];
	vec dd[2][GFBITS], bb[2][GFBITS];
	vec B[2][GFBITS], C[2][GFBITS];
	vec B_tmp[2][GFBITS], C_tmp[2][GFBITS];


	gf d, b;
	gf coefs[256];

	vec tmp_0, tmp_1;

	// initialization

#pragma HLS array_partition variable=prod
#pragma HLS array_partition variable=dd
#pragma HLS array_partition variable=bb
#pragma HLS array_partition variable=B
#pragma HLS array_partition variable=C
#pragma HLS array_partition variable=B_tmp
#pragma HLS array_partition variable=C_tmp
#pragma HLS array_partition variable=interval

	get_coefs(&coefs[  0], in_stream_0);
	get_coefs(&coefs[ 64], in_stream_1);
	get_coefs(&coefs[128], in_stream_2);
	get_coefs(&coefs[192], in_stream_3);


	C[0][0] = 0;
	C[1][0] = one << 63;
	B[0][0] = 0; 
	B[1][0] = one << 62;

	for(i = 1; i < GFBITS; i++){
	#pragma HLS pipeline
		C[0][i] = 0;
		C[1][i] = 0;
		B[0][i] = 0;
		B[1][i] = 0;
	}

	b = 1;
	L = 0;

	//

	for (i = 0; i < GFBITS; i++){
	#pragma HLS pipeline
		interval[0][i] = 0;
		interval[1][i] = 0;
	}

	for (N = 0; N < SYS_T*2; N++)
	{
	#pragma HLS pipeline
		update(interval, coefs[N]);

		vec_mul(prod[0], C[0], interval[0]);
		vec_mul(prod[1], C[1], interval[1]);

		d = vec_reduce(prod);
//		printf("d[%d]=%lX\n", N, d.to_uint64());

		mask_1 = mask_nonzero(d);
//		printf("mask1[%d]=%lX\n", N.to_uint64(), mask_1.to_uint64());
		mask_2 = mask_leq(L*2, N);
//		printf("mask2[%d]=%lX\n", N.to_uint64(), mask_2.to_uint64());

		mask = mask_1 & mask_2;
//		mask = mask_nonzero(d) & mask_leq(L*2, N);

		for (i = 0; i < GFBITS; i++) 
		{
			dd[0][i] = dd[1][i] = vec_setbits((d >> i) & 1);
			bb[0][i] = bb[1][i] = vec_setbits((b >> i) & 1);

		}
		
		vec_mul(B_tmp[0], dd[0], B[0]);
		vec_mul(B_tmp[1], dd[1], B[1]);
		vec_mul(C_tmp[0], bb[0], C[0]);
		vec_mul(C_tmp[1], bb[1], C[1]);

		vec_cmov(B[0], C[0], mask);
		vec_cmov(B[1], C[1], mask);

		update(B, 0);


		for (i = 0; i < GFBITS; i++)
		{ 
			C[0][i] = B_tmp[0][i] ^ C_tmp[0][i];
			C[1][i] = B_tmp[1][i] ^ C_tmp[1][i];

		}

		b = (d & mask) | (b & ~mask);
		L = ((N+1-L) & mask) | (L & ~mask);

	}


	for (i = 0; i < GFBITS; i++)
	#pragma HLS pipeline
	{
		tmp_0 = (C[0][i] >> 31) | (C[1][i] << 33);
		tmp_1 =  C[1][i] >> 31;

		out_stream_0.write(tmp_0);
		out_stream_1.write(tmp_1);

	}
}



/* input: in, polynomial in bitsliced form */
/* output: in, result of applying the radix conversions on in */
void radix_conversions_stream(hls::stream<vec> &out_stream, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1)
{
	int i, j, k, l;

	vec in[2][GFBITS];

	const vec mask[5][2] = 
	{
		{0x8888888888888888, 0x4444444444444444},
		{0xC0C0C0C0C0C0C0C0, 0x3030303030303030},
		{0xF000F000F000F000, 0x0F000F000F000F00},
		{0xFF000000FF000000, 0x00FF000000FF0000},
		{0xFFFF000000000000, 0x0000FFFF00000000}
	};

	const vec s[5][2][GFBITS] = 
	{
#include "scalars_2x.data"
	};

	
	for (j = 0; j <= 5; j++)
	{
		if(j==0){
			for (i = 0; i < GFBITS; i++)
			{
				in[1][i] = in_stream_1.read();
				in[0][i] = in_stream_0.read();
			}
		}
		for (i = 0; i < GFBITS; i++)
		{
			in[1][i] ^= in[1][i] >> 32;
			in[0][i] ^= in[1][i] << 32;
		}
	


		for (i = 0; i < GFBITS; i++)
		for (k = 4; k >= j; k--)
		{
			in[0][i] ^= (in[0][i] & mask[k][0]) >> (1 << k);
			in[0][i] ^= (in[0][i] & mask[k][1]) >> (1 << k);
			in[1][i] ^= (in[1][i] & mask[k][0]) >> (1 << k);
			in[1][i] ^= (in[1][i] & mask[k][1]) >> (1 << k);
		}

		if (j < 5)
		{
			vec_mul(in[0], in[0], s[j][0]);
			vec_mul(in[1], in[1], s[j][1]);
		} else {
			for (l = 0; l < GFBITS; l++)
			{
				out_stream.write(in[0][l]);
			}
			for (l = 0; l < GFBITS; l++)
			{
				out_stream.write(in[1][l]);
			}
		}

	}

}


void butterflies_stream(hls::stream<vec> &out_stream, hls::stream<vec> &in_stream){



	int i, j, k, s, b;

	vec tmp[ GFBITS ];
	vec pre[8][ GFBITS ];
	vec buf[128], buf_tr[128];

	vec consts_ptr = 2;

	vec consts[ 128 ][ GFBITS ] =
	{
#include "consts.data"
	};

	unsigned char reversal[128] =
	{ 
	  0, 64, 32, 96, 16, 80, 48, 112, 
	  8, 72, 40, 104, 24, 88, 56, 120, 
	  4, 68, 36, 100, 20, 84, 52, 116, 
	  12, 76, 44, 108, 28, 92, 60, 124, 
	  2, 66, 34, 98, 18, 82, 50, 114, 
	  10, 74, 42, 106, 26, 90, 58, 122, 
	  6, 70, 38, 102, 22, 86, 54, 118, 
	  14, 78, 46, 110, 30, 94, 62, 126, 
	  1, 65, 33, 97, 17, 81, 49, 113, 
	  9, 73, 41, 105, 25, 89, 57, 121, 
	  5, 69, 37, 101, 21, 85, 53, 117, 
	  13, 77, 45, 109, 29, 93, 61, 125, 
	  3, 67, 35, 99, 19, 83, 51, 115, 
	  11, 75, 43, 107, 27, 91, 59, 123, 
	  7, 71, 39, 103, 23, 87, 55, 119, 
	  15, 79, 47, 111, 31, 95, 63, 127
	};

	ap_uint<16> beta[7] = {2522, 7827, 7801, 8035, 6897, 8167, 3476};

	vec in[2][GFBITS];
	vec out[128][GFBITS];


	//Read the matrix 0
	for (i = 0; i < GFBITS; i++)
	{
		in[0][i] = in_stream.read();
//		printf("in[0][%d] = %lX\n", i, in[0][i].to_uint64());
	}

	//Read the matrix 1
	for (i = 0; i < GFBITS; i++)
	{
		in[1][i] = in_stream.read();
//		printf("in[0][%d] = %lX\n", i, in[1][i].to_uint64());
	}

	for (i = 0; i < 7; i++)
	{
		for (j = 0; j < GFBITS; j++)
		{
			pre[i][j] = (beta[i] >> j) & 1; 
			pre[i][j] = -pre[i][j];
		}

		vec_mul(pre[i], in[1], pre[i]);
	}

	for (i = 0; i < GFBITS; i++)
	{
		buf[0] = in[0][i];             

		buf[1] = buf[0] ^ pre[0][i];      buf[32] = in[0][i] ^ pre[5][i];
		buf[3] = buf[1] ^ pre[1][i];      buf[96] = buf[32] ^ pre[6][i];
                                              buf[97] = buf[96] ^ pre[0][i];
		buf[2] = in[0][i] ^ pre[1][i];  buf[99] = buf[97] ^ pre[1][i];
		buf[6] = buf[2] ^ pre[2][i];      buf[98] = buf[99] ^ pre[0][i];
		buf[7] = buf[6] ^ pre[0][i];      buf[102] = buf[98] ^ pre[2][i];
		buf[5] = buf[7] ^ pre[1][i];      buf[103] = buf[102] ^ pre[0][i];
                                              buf[101] = buf[103] ^ pre[1][i];
		buf[4] = in[0][i] ^ pre[2][i];  buf[100] = buf[101] ^ pre[0][i];
		buf[12] = buf[4] ^ pre[3][i];     buf[108] = buf[100] ^ pre[3][i];
		buf[13] = buf[12] ^ pre[0][i];    buf[109] = buf[108] ^ pre[0][i];
		buf[15] = buf[13] ^ pre[1][i];    buf[111] = buf[109] ^ pre[1][i];
		buf[14] = buf[15] ^ pre[0][i];    buf[110] = buf[111] ^ pre[0][i];
		buf[10] = buf[14] ^ pre[2][i];    buf[106] = buf[110] ^ pre[2][i];
		buf[11] = buf[10] ^ pre[0][i];    buf[107] = buf[106] ^ pre[0][i];
		buf[9] = buf[11] ^ pre[1][i];     buf[105] = buf[107] ^ pre[1][i];
                                              buf[104] = buf[105] ^ pre[0][i];
		buf[8] = in[0][i] ^ pre[3][i];  buf[120] = buf[104] ^ pre[4][i];
		buf[24] = buf[8] ^ pre[4][i];     buf[121] = buf[120] ^ pre[0][i];
		buf[25] = buf[24] ^ pre[0][i];    buf[123] = buf[121] ^ pre[1][i];
		buf[27] = buf[25] ^ pre[1][i];    buf[122] = buf[123] ^ pre[0][i];
		buf[26] = buf[27] ^ pre[0][i];    buf[126] = buf[122] ^ pre[2][i];
		buf[30] = buf[26] ^ pre[2][i];    buf[127] = buf[126] ^ pre[0][i];
		buf[31] = buf[30] ^ pre[0][i];    buf[125] = buf[127] ^ pre[1][i];
		buf[29] = buf[31] ^ pre[1][i];    buf[124] = buf[125] ^ pre[0][i];
		buf[28] = buf[29] ^ pre[0][i];    buf[116] = buf[124] ^ pre[3][i];
		buf[20] = buf[28] ^ pre[3][i];    buf[117] = buf[116] ^ pre[0][i];
		buf[21] = buf[20] ^ pre[0][i];    buf[119] = buf[117] ^ pre[1][i];
		buf[23] = buf[21] ^ pre[1][i];    buf[118] = buf[119] ^ pre[0][i];
		buf[22] = buf[23] ^ pre[0][i];    buf[114] = buf[118] ^ pre[2][i];
		buf[18] = buf[22] ^ pre[2][i];    buf[115] = buf[114] ^ pre[0][i];
		buf[19] = buf[18] ^ pre[0][i];    buf[113] = buf[115] ^ pre[1][i];
		buf[17] = buf[19] ^ pre[1][i];    buf[112] = buf[113] ^ pre[0][i];
                                              buf[80] = buf[112] ^ pre[5][i];
		buf[16] = in[0][i] ^ pre[4][i]; buf[81] = buf[80] ^ pre[0][i];
		buf[48] = buf[16] ^ pre[5][i];    buf[83] = buf[81] ^ pre[1][i];
		buf[49] = buf[48] ^ pre[0][i];    buf[82] = buf[83] ^ pre[0][i];
		buf[51] = buf[49] ^ pre[1][i];    buf[86] = buf[82] ^ pre[2][i];
		buf[50] = buf[51] ^ pre[0][i];    buf[87] = buf[86] ^ pre[0][i];
		buf[54] = buf[50] ^ pre[2][i];    buf[85] = buf[87] ^ pre[1][i];
		buf[55] = buf[54] ^ pre[0][i];    buf[84] = buf[85] ^ pre[0][i];
		buf[53] = buf[55] ^ pre[1][i];    buf[92] = buf[84] ^ pre[3][i];
		buf[52] = buf[53] ^ pre[0][i];    buf[93] = buf[92] ^ pre[0][i];
		buf[60] = buf[52] ^ pre[3][i];    buf[95] = buf[93] ^ pre[1][i];
		buf[61] = buf[60] ^ pre[0][i];    buf[94] = buf[95] ^ pre[0][i];
		buf[63] = buf[61] ^ pre[1][i];    buf[90] = buf[94] ^ pre[2][i];
		buf[62] = buf[63] ^ pre[0][i];    buf[91] = buf[90] ^ pre[0][i];
		buf[58] = buf[62] ^ pre[2][i];    buf[89] = buf[91] ^ pre[1][i];
		buf[59] = buf[58] ^ pre[0][i];    buf[88] = buf[89] ^ pre[0][i];
		buf[57] = buf[59] ^ pre[1][i];    buf[72] = buf[88] ^ pre[4][i];
		buf[56] = buf[57] ^ pre[0][i];    buf[73] = buf[72] ^ pre[0][i];
		buf[40] = buf[56] ^ pre[4][i];    buf[75] = buf[73] ^ pre[1][i];
		buf[41] = buf[40] ^ pre[0][i];    buf[74] = buf[75] ^ pre[0][i];
		buf[43] = buf[41] ^ pre[1][i];    buf[78] = buf[74] ^ pre[2][i];
		buf[42] = buf[43] ^ pre[0][i];    buf[79] = buf[78] ^ pre[0][i];
		buf[46] = buf[42] ^ pre[2][i];    buf[77] = buf[79] ^ pre[1][i];
		buf[47] = buf[46] ^ pre[0][i];    buf[76] = buf[77] ^ pre[0][i];
		buf[45] = buf[47] ^ pre[1][i];    buf[68] = buf[76] ^ pre[3][i];
		buf[44] = buf[45] ^ pre[0][i];    buf[69] = buf[68] ^ pre[0][i];
		buf[36] = buf[44] ^ pre[3][i];    buf[71] = buf[69] ^ pre[1][i];
		buf[37] = buf[36] ^ pre[0][i];    buf[70] = buf[71] ^ pre[0][i];
		buf[39] = buf[37] ^ pre[1][i];    buf[66] = buf[70] ^ pre[2][i];
		buf[38] = buf[39] ^ pre[0][i];    buf[67] = buf[66] ^ pre[0][i];
		buf[34] = buf[38] ^ pre[2][i];    buf[65] = buf[67] ^ pre[1][i];
		buf[35] = buf[34] ^ pre[0][i];                                     
		buf[33] = buf[35] ^ pre[1][i];    buf[64] = in[0][i] ^ pre[6][i];

		transpose_64x64(buf_tr +  0, buf +  0);
		transpose_64x64(buf_tr + 64, buf + 64);

		//This was added by me since transpse needs different iout and outout to function correctly
		vec_copy128(buf, buf_tr);


		for (j = 0; j < 128; j++){
			out[ reversal[j] ][i] = buf[j];
		}
	}


	//butterflies
	for (i = 1; i <= 6; i++)
	{
		s = 1 << i;

		for (j = 0; j < 128; j += 2*s){
			for (k = j; k < j+s; k++)
			{
				vec_mul(tmp, out[k+s], consts[ consts_ptr + (k-j) ]);
				vec_add(out[k], out[k], tmp);
				vec_add(out[k+s], out[k+s], out[k]);

			}
		}

		consts_ptr += (1 << i);
	}

	//Write to the output
	for (i = 0; i < 128; i++){
		for (j = 0; j < GFBITS; j++){
			out_stream.write(out[i][j]);
//			printf("out[%d][%d] = %lX\n", i, j, out[i][j].to_uint64());
		}
	}




}

void write_output(hls::stream<vec> &out_stream, hls::stream<mat_struct> &in_stream_0, hls::stream<mat_struct> &in_stream_1){

	int i, b;

	vec tmp_val_0, tmp_val_1;
	mat_struct tmp_struct;

	//Pass directly as structs?

	for (i = 0; i < 64; i++){
	#pragma HLS pipeline
		tmp_struct = in_stream_0.read();
		for (b = 0; b < GFBITS; b++){
			tmp_val_0 = tmp_struct.mat[b];
			out_stream.write(tmp_val_0);

		}

	}

	for (i = 64; i < 128; i++){
	#pragma HLS pipeline
		tmp_struct = in_stream_1.read();
		for (b = 0; b < GFBITS; b++){
			tmp_val_1 = tmp_struct.mat[b];
			out_stream.write(tmp_val_1);

		}

	}




}




////////////////////FFT//////////////////////


vec load8(const unsigned char * in)
{
	int i;
	vec ret = in[7];

	for (i = 6; i >= 0; i--)
	{
		ret <<= 8;
		ret |= in[i];
	}

	return ret;
}


/* middle layers of the benes network */
static void layer_in(vec data[2][64], vec * bits, int lgs)
{
	int i, j, s;

	vec d;

	s = 1 << lgs;

	for (i = 0; i < 64; i += s*2)
	for (j = i; j < i+s; j++)
	{

		d = (data[0][j+0] ^ data[0][j+s]);
		d &= (*bits++);
		data[0][j+0] ^= d;
		data[0][j+s] ^= d;

		d = (data[1][j+0] ^ data[1][j+s]);
		d &= (*bits++);
		data[1][j+0] ^= d;
		data[1][j+s] ^= d;
	}
}

/* first and last layers of the benes network */
static void layer_ex(vec * data, vec * bits, int lgs)
{
	int i, j, s;

	vec d;

	s = 1 << lgs;

	for (i = 0; i < 128; i += s*2)
	for (j = i; j < i+s; j++)
	{

		d = (data[j+0] ^ data[j+s]);
		d &= (*bits++);
		data[j+0] ^= d;
		data[j+s] ^= d;
	}
}


/* bitsliced field inverses */
void vec_inv(vec * out, vec * in)
{
	vec tmp_11[ GFBITS ];
	vec tmp_1111[ GFBITS ];

	vec_copy(out, in);

	vec_sq(out, out);
	vec_mul(tmp_11, out, in); // ^11

	vec_sq(out, tmp_11);
	vec_sq(out, out);
	vec_mul(tmp_1111, out, tmp_11); // ^1111

	vec_sq(out, tmp_1111);
	vec_sq(out, out);
	vec_sq(out, out);
	vec_sq(out, out);
	vec_mul(out, out, tmp_1111); // ^11111111

	vec_sq(out, out);
	vec_sq(out, out);
	vec_sq(out, out);
	vec_sq(out, out);
	vec_mul(out, out, tmp_1111); // ^111111111111

	vec_sq(out, out); // ^1111111111110
}


/* input: r, sequence of bits to be permuted */
/*        bits, condition bits of the Benes network */
/*        rev, 0 for normal application; !0 for inverse */
/* output: r, permuted bits */
void benes(hls::stream<vec> &r_out, hls::stream<vec> &r_in, hls::stream<sk_4k_vec> &vec_bits)
{
	int i, iter;	//


	const unsigned char *cond_ptr;
	int inc, low;

	vec cond[64], cond_tr[64],cond_tr_2[64];

	vec r1[64], r2[64], r3[64], r4[64];
	vec r_mat[64];
	sk_4k_vec tmp;
	int sk_vec_idx = 0;

	vec r_int_v[2][64];
	vec r_int_h[2][64];
	vec b_int_v[64];
	vec b_int_h[64];

#pragma HLS array_partition variable=r_int_v
#pragma HLS array_partition variable=r_int_h


#pragma HLS array_partition variable=r_mat
#pragma HLS array_partition variable=r1
#pragma HLS array_partition variable=r2
#pragma HLS array_partition variable=r3
#pragma HLS array_partition variable=cond
#pragma HLS array_partition variable=cond_tr
#pragma HLS array_partition variable=cond_tr_2
	

	//Load the inout stream and perform the matrix split
	for (i = 0; i < 64; i++)
	{
		r_int_v[0][i] = r_in.read();
		r_int_v[1][i] = r_in.read();
	}
	// load_r(r_mat, r_in);


	transpose_64x64(r_int_h[0], r_int_v[0]);
	transpose_64x64(r_int_h[1], r_int_v[1]);
	// transpose_64x64(r1, r_mat);


for (iter = 0; iter <= 6; iter++)
#pragma HLS pipeline
	{
		//Reads 4kb
		tmp = vec_bits.read();
		for (i = 0; i < 64; i++)
		{
			b_int_v[i].range(63,0) = tmp.range(((i+1)*64)-1, (i*64));
		}

		transpose_64x64(b_int_h, b_int_v);
		layer_ex(r_int_h[0], b_int_h, iter);
	}

	transpose_64x64(r_int_v[0], r_int_h[0]);
	transpose_64x64(r_int_v[1], r_int_h[1]);

	for (iter = 0; iter <= 5; iter++) 
	{
		//Reads 4kb
		tmp = vec_bits.read();
		for (i = 0; i < 64; i++)
		{
			b_int_v[i].range(63,0) = tmp.range(((i+1)*64)-1, (i*64));
		}

		layer_in(r_int_v, b_int_v, iter);
	}

	for (iter = 4; iter >= 0; iter--) 
	{
		//Reads 4kb
		tmp = vec_bits.read();
		for (i = 0; i < 64; i++)
		{
			b_int_v[i].range(63,0) = tmp.range(((i+1)*64)-1, (i*64));
		}

		layer_in(r_int_v, b_int_v, iter);
	}

	transpose_64x64(r_int_h[0], r_int_v[0]);
	transpose_64x64(r_int_h[1], r_int_v[1]);

	for (iter = 6; iter >= 0; iter--)
	{
		//Reads 4kb
		tmp = vec_bits.read();
		for (i = 0; i < 64; i++)
		{
			b_int_v[i].range(63,0) = tmp.range(((i+1)*64)-1, (i*64));
		}

		transpose_64x64(b_int_h, b_int_v);
		layer_ex(r_int_h[0], b_int_h, iter);
	}

	transpose_64x64(r_int_v[0], r_int_h[0]);
	transpose_64x64(r_int_v[1], r_int_h[1]);

	for (i = 0; i < 64; i++)
	{
		r_out.write(r_int_v[0][i]);
		r_out.write(r_int_v[1][i]);
	}

}


//void benes_test(hls::stream<vec> &r_out, hls::stream<vec> &r_in, hls::stream<sk_1k_vec> &vec_bits)
//{
//	int i;	//
//
//
//	const unsigned char *cond_ptr;
//	int inc, low;
//
//	vec cond[64], cond_tr[64],cond_tr_2[64];
//
//	vec r1[64], r2[64], r3[64], r4[64];
//
//#pragma HLS array_partition variable=r1
//#pragma HLS array_partition variable=r2
//#pragma HLS array_partition variable=r3
//
//
//
//#pragma HLS array_partition variable=cond_tr
//#pragma HLS array_partition variable=cond_tr
//#pragma HLS array_partition variable=cond_tr_2
//
//vec r_mat[64];
//sk_1k_vec tmp_1, tmp_2;
//
//	int sk_vec_idx = 0;
//
//	load_r(r_mat, r_in);
//
//	transpose_64x64(r1, r_mat);
//
//
//	for (low = 0; low <= 5; low++)
//#pragma HLS pipeline
//	{
//		//This loop read 256B, so one whole sk_wide_vec.
//		//Split accordingly!
//		tmp_1 = vec_bits.read();
//		for (i = 0; i < 32; i++){
//			cond[i].range(31,0) = tmp_1.range(((i+1)*32)-1, (i*32));
//			cond[i].range(63,32) = 0;
//		}
//
//		tmp_2 = vec_bits.read();
//		for (i = 0; i < 32; i++){
//			cond[i+32].range(31,0) = tmp_2.range(((i+1)*32)-1, (i*32));
//			cond[i+32].range(63,32) = 0;
//		}
//
//		transpose_64x64(cond_tr, cond);
//		layer(r1, cond_tr, low);
//	}
//
//	transpose_64x64(r2, r1);
//
//	for (low = 0; low <= 5; low++)
//	#pragma HLS pipeline
//	{
//		tmp_1 = vec_bits.read();
//		for (i = 0; i < 16; i++){
//			cond_tr[i] = tmp_1.range(((i+1)*64)-1, (i*64));
//		}
//
//		tmp_2 = vec_bits.read();
//		for (i = 0; i < 16; i++){
//			cond_tr[i+16] = tmp_2.range(((i+1)*64)-1, (i*64));
//		}
//
//		layer(r2, cond_tr, low);
//	}
//
//	for (low = 4; low >= 0; low--)
//	#pragma HLS pipeline
//	{
//		tmp_1 = vec_bits.read();
//		for (i = 0; i < 16; i++){
//			cond_tr[i] = tmp_1.range(((i+1)*64)-1, (i*64));
//
//		}
//
//		tmp_2 = vec_bits.read();
//		for (i = 0; i < 16; i++){
//			cond_tr[i+16] = tmp_2.range(((i+1)*64)-1, (i*64));
//
//		}
//
//		layer(r2, cond_tr, low);
//	}
//
//	transpose_64x64(r3, r2);
//
//	for (low = 5; low >= 0; low--)
//	#pragma HLS pipeline
//	{
//		tmp_1 = vec_bits.read();
//		for (i = 0; i < 32; i++){
//			cond_tr[i].range(31,0) = tmp_1.range(((i+1)*32)-1, (i*32));
//			cond_tr[i].range(63,32) = 0;
//		}
//
//		tmp_2 = vec_bits.read();
//		for (i = 0; i < 32; i++){
//			cond_tr[i+32].range(31,0) = tmp_2.range(((i+1)*32)-1, (i*32));
//			cond_tr[i+32].range(63,32) = 0;
//		}
//
//		transpose_64x64(cond_tr_2, cond_tr);
//		layer(r3, cond_tr_2, low);
//	}
//
//	transpose_64x64_benes_out(r_out, r3);
//
//}


void load_irr_vec(hls::stream<sk_05k_vec> &out, sk_05k_vec * sk){

	int i;
	
	for(i=0; i<3; i++){
		out.write(*(sk+i));		
	}
}



 void irr_load_stream(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<sk_05k_vec> &in)
{
	int i, j;
	vec irr[ SYS_T + 1 ];
	sk_05k_vec tmp;

	vec out[GFBITS];

	vec v[2];

#pragma HLS array_partition variable=irr complete
#pragma HLS array_partition variable=tmp complete

	for(i=0; i<3; i++){
		tmp = in.read();
		for(j=0; j<(SYS_T/3); j++){
			#pragma HLS pipeline
			irr[i*32+j].range(63, 13) = 0;
			irr[i*32+j].range(12, 0) = tmp.range(((j+1)*16)-4, (j*16));
		}
	}

	irr[ SYS_T ] = 1;

	// for (i = 0; i < GFBITS; i++)
	// 	out[i] = 0;

	//TODO fix here! The following loop with ap_uints was giving wrong results
	// for (i = SYS_T; i >= 0; i--)
	// #pragma HLS pipeline
	// for (j = 0; j < GFBITS; j++)
	// {
	// 	out[j] <<= 1;
	// 	out[j] |= (irr[i] >> j) & 1;

	// }

	//Since my vec is arbitrarily long (96b here) I do not have to split in a s2 matrix.
	// for (j = 0; j < GFBITS; j++)
	// {
	// 	out_stream.write(out[j]);
	// }

	//Solution with 2 matrices
	//Maybe modify like this to stream immediately for lvl 1 too
	for (i = 0; i < GFBITS; i++)
	{
		v[0] = v[1] = 0;

		for (j =    63; j >=  0; j--) {
			v[0] <<= 1; 
			v[0] |= (irr[j] >> i) & 1;
		}
		for (j = SYS_T; j >= 64; j--) {
			v[1] <<= 1; 
			v[1] |= (irr[j] >> i) & 1;
		}

		// out[0][i] = v[0];
		// out[1][i] = v[1];

		out_stream_0.write(v[0]);
		out_stream_1.write(v[1]);

//		printf("v[0][%d] = %lX\n", i, v[0].to_uint64());
//		printf("v[1][%d] = %lX\n", i, v[1].to_uint64());

	}




	//TODO Try also this solution
	// for (i = SYS_T; i >= 1; i--)
	// #pragma HLS pipeline
	// for (j = 0; j < GFBITS; j++)
	// {
	// 	out[j] <<= 1;
	// 	out[j] |= (irr[i] >> j) & 1;
	// }

	// for (j = 0; j < GFBITS; j++)
	// #pragma HLS pipeline
	// {
	// 	out[j] <<= 1;
	// 	out[j] |= (irr[0] >> j) & 1;
	// 	out_stream.write(out[j]);
	// }



}







void vec_sq_loop(hls::stream<vec> &mat_out_stream, hls::stream<vec> &mat_in_stream){

	int i, j;
	vec tmp_vec_in[GFBITS];
	vec tmp_vec_out[GFBITS];


	//TODO Make more elegant. Apply dataflow inside the vec_sq and make a vec_sq stream
	for(j=0; j<GFBITS; j++){
		tmp_vec_in[j] = mat_in_stream.read();
	}

	for (i = 0; i < 128; i++){
	#pragma HLS pipeline
		vec_sq(tmp_vec_out, tmp_vec_in);

		for(j=0; j<GFBITS; j++){
//			printf("out[%d] = %lX\n", j, tmp_vec_out[j].to_uint64());
			mat_out_stream.write(tmp_vec_out[j]);

			//It has to be read one time less since it is read outside the loop once
			if(i<127){
				tmp_vec_in[j] = mat_in_stream.read();
			}
		}
	}
}

//TODO find a solution to stream with structs and not have to read in loops the inner matrix. can it be?
void stream_mat_rev(hls::stream<vec> &mat_out_rev_stream, vec mat_in[64][GFBITS]){

	int i, j;

	for (i = 127; i >= 0; i--){
		for(j=0; j<GFBITS; j++){
		#pragma HLS pipeline
			mat_out_rev_stream.write(mat_in[i][j]);
		}
	}


}

void mul_loop_1(hls::stream<vec> &mat_out_stream, hls::stream<vec> &mat_in_stream){

	int i;
	vec mat_int[128][GFBITS];

	#pragma HLS array_partition variable=mat_int dim=2 complete

	for (i = 0; i < GFBITS; i++){
	#pragma HLS pipeline
		mat_int[0][i] = mat_in_stream.read();
	}

	for (i = 1; i < 128; i++){
	#pragma HLS pipeline
		vec_mul_stream_in(mat_int[i], mat_in_stream, mat_int[i-1]);
	}

	stream_mat_rev(mat_out_stream, mat_int);
}


void mul_loop_2(hls::stream<vec> &mat_out_stream, hls::stream<vec> &inv_in_stream, hls::stream<vec> &eval_in_stream){

	int i, j;
	vec tmp[ GFBITS ];
	vec tmp_vec_mat[ GFBITS ];
	vec tmp_vec_mat_out[128][ GFBITS ];

	#pragma HLS array_partition variable=tmp complete
	#pragma HLS array_partition variable=tmp_vec_mat_out complete

	for(j=0; j<GFBITS; j++){
	#pragma HLS pipeline
		tmp_vec_mat[j] = inv_in_stream.read();
	}

	vec_inv(tmp, tmp_vec_mat);

	//TODO apply dataflow here? Or partition the array ot foster pipelining
	for (i = 126; i >= 0; i--)
	#pragma HLS pipeline
	{
		//TODO change the output inv to be able to stream it to the next
		vec_mul_stream_in(tmp_vec_mat_out[i+1], inv_in_stream, tmp);
		vec_mul_stream_in(tmp, eval_in_stream, tmp);

	}


	for(j=0; j<GFBITS; j++){
	#pragma HLS pipeline
		tmp_vec_mat_out[0][j] = tmp[j];
	}

	//Stream the rest of the values
	for(i=0; i<128; i++){
		for(j=0; j<GFBITS; j++){
		#pragma HLS pipeline
			mat_out_stream.write(tmp_vec_mat_out[i][j]);
		}
	}

}

void compute_outputs(hls::stream<mat_struct> &output_struct, hls::stream<vec> &inv_output_stream, hls::stream<vec> &inv_input_stream, hls::stream<vec> &recv_in_stream){

	int i, j;
	vec out_mat[128][GFBITS];
	vec tmp_recv;
	vec tmp_inv;

	mat_struct tmp_struct;

	for (i = 0; i < 128; i++){
	#pragma HLS pipeline
		tmp_recv = recv_in_stream.read();
		for (j = 0; j < GFBITS; j++){
			tmp_inv = inv_input_stream.read();
			out_mat[i][j] = tmp_inv & tmp_recv;
			inv_output_stream.write(tmp_inv);
			tmp_struct.mat[j] = out_mat[i][j];

//			printf("out[%d][%d] = %lX\n", i, j, out_mat[i][j].to_uint64());
//			printf("inv[%d][%d] = %lX\n", i, j, tmp_inv.to_uint64());
		}
		output_struct.write(tmp_struct);
	}

}

void eval_stream_duplicate(hls::stream<vec> &out_stream_fwd, hls::stream<vec> &out_stream_rev, hls::stream<vec> &in_stream){

	//TODO the following loops can be executed in parallel if I provide with two same in streams.
	//And provided I put the dataflow directive here
	int i, j;
	vec tmp_vec[128][GFBITS];

	//Stream the forward version
	for (i = 0; i < 128; i++){
	#pragma HLS pipeline
		for(j=0; j<GFBITS; j++){
			tmp_vec[i][j] = in_stream.read();
			out_stream_fwd.write(tmp_vec[i][j]);
		}
	}

	//Stream the reverse version
	for (i = 127; i >= 1; i--){
	#pragma HLS pipeline
		for(j=0; j<GFBITS; j++){
			out_stream_rev.write(tmp_vec[i][j]);
		}
	}
}

void recv_read(hls::stream<vec> &recv_out, hls::stream<vec> &recv_in){

	int i;
	vec tmp;

	for(i=0; i<64; i++){
		tmp = recv_in.read();
		recv_out.write(tmp);
	}
}

 void preprocess(hls::stream<vec> &recv, hls::stream<unsigned char> &s)
{
	int i;
	unsigned char r[ 1024 ];
	vec tmp;

	for (i = 0; i < SYND_BYTES; i++){
		r[i] = s.read();
//		printf("r[%d]=%lX\n", i, r[i]);
	}

	for (i = SYND_BYTES; i < 1024; i++)
		r[i] = 0;

	//TODO some loops are only zero....simplify
	for (i = 0; i < 128; i++){
		tmp = load8(r + i*8);
		recv.write(tmp);
//		printf("recv[%d]=%lX\n", i, tmp.to_uint64());
	}
}

//TODO Have extra unneeded iterations here (512-436)
//TODO maybe merge loops?
 void postprocess(unsigned char *e, hls::stream<unsigned char> &e_stream_out, hls::stream<vec> &err_stream_in)
{
	int i, j;
	unsigned char error8[ (1 << GFBITS)/8 ];

#pragma HLS array_partition variable=error8 factor=8

	vec tmp;

	for (i = 0; i < 128; i++){
		tmp = err_stream_in.read();
		for(j=0; j<8; j++){
			error8[i*8+j] = tmp.range((j+1)*8-1, j*8);
		}

	}

	for (i = 0; i < SYS_N/8; i++){
		e[i] = error8[i];
		e_stream_out.write(e[i]);
	}

}

 void scaling_inv(hls::stream<mat_struct> &out_struct, hls::stream<vec> &inv, hls::stream<vec> &recv)
{
	int i, j;
	vec recv_tmp;

	vec out[128][GFBITS];

	for (i = 0; i < 128; i++){
		recv_tmp = recv.read();
		for (j = 0; j < GFBITS; j++){
			out[i][j] = inv.read() & recv_tmp;
//			printf("out[%d][%d] = %lX\n", i, j, out[i][j].to_uint64());
		}
	}

//	stream_out(out_struct, out);
	int k, l;

	mat_struct tmp_struct;

	for (k = 0; k < 128; k++)
	{
		for(l=0; l<GFBITS; l++){
			tmp_struct.mat[l] = out[k][l];
		}
		out_struct.write(tmp_struct);
	}


}

void weight_check(hls::stream<bit> &check_weight, hls::stream<unsigned char> &e_stream, hls::stream<vec> &error_stream)
{
	int i, j;

	gf w0 = 0;
	gf w1 = 0;
	gf check;
	bit check_bit;

	vec tmp_error;
	ap_uint<8> tmp_e;

	//TODO change unsigned char to ap_uint<8>?

	//TODO In fact here it just adds all bits of error together!
	//TODO change accordingly?

	for (i = 0; i < (1 << GFBITS)/64; i++){
		tmp_error = error_stream.read();
		for(j=0; j<64; j++){
			w0 += (tmp_error.bit(j));
		}
	}

	for (i = 0; i < SYS_N/8; i++){
		tmp_e = e_stream.read();
		for(j=0; j<8; j++){
			w1 += (tmp_e.bit(j));
		}
	}

	check = (w0 ^ SYS_T) | (w1 ^ SYS_T);
	check -= 1;
	check_bit = check.bit(12);

	check_weight.write(check_bit);

}

 void synd_cmp(hls::stream<bit> &check_out, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1, hls::stream<vec> &s_stream_2, hls::stream<vec> &s_stream_3, hls::stream<vec> &s_cmp_stream_0, hls::stream<vec> &s_cmp_stream_1, hls::stream<vec> &s_cmp_stream_2, hls::stream<vec> &s_cmp_stream_3)
{
	int i, j;
	vec diff = 0;
	vec diff_0 = 0;
	vec diff_1 = 0;
	vec diff_2 = 0;
	vec diff_3 = 0;

	vec out_bit;

	for (j = 0; j < GFBITS; j++){
	#pragma HLS pipeline
		diff_0 |= (s_stream_0.read() ^ s_cmp_stream_0.read());
		diff_1 |= (s_stream_1.read() ^ s_cmp_stream_1.read());
		diff_2 |= (s_stream_2.read() ^ s_cmp_stream_2.read());
		diff_3 |= (s_stream_3.read() ^ s_cmp_stream_3.read());
	}

	diff = diff_0 | diff_1 | diff_2 | diff_3;

	diff |= diff >> 32;
	diff |= diff >> 16;
	diff |= diff >> 8;
	diff |= diff >> 4;
	diff |= diff >> 2;
	diff |= diff >> 1;
	diff = (diff&1)^1;

	out_bit = diff.bit(0);
	check_out.write(out_bit);
}


void load_s(hls::stream<unsigned char> &s_mat_stream, unsigned char *s){
	
	int i;

	LOOP_S_READ:
	for (i = 0; i < SYND_BYTES; i++){
		s_mat_stream.write(*(s+i));
	}

}

void store_e(unsigned char *e, unsigned char *e_mat){
	
	int i;

	LOOP_E_STORE:
	for (i = 0; i < (SYS_N/8); i++){
		*(e+i) = *(e_mat+i);
	}
}

void store_check_value(int *check_value, int *check){

	*check_value = *check;
}

void error_compute(hls::stream<vec> &out_stream, hls::stream<vec> &in_stream){


	vec allone = 0xffffffffffffffff;
	vec error_mat[128];
	int i, j;

	for (i = 0; i < 128; i++)
	{
		error_mat[i] = vec_or_reduce(in_stream);
		error_mat[i] ^= allone;
		out_stream.write(error_mat[i]);
//		printf("error[%d] = %lX\n", i, error_mat[i].to_uint64());
	}
}

void priv_stream_duplicate(hls::stream<vec> &out_stream_0_0, hls::stream<vec> &out_stream_0_1, hls::stream<vec> &out_stream_0_2, hls::stream<vec> &out_stream_0_3, hls::stream<vec> &out_stream_1_0, hls::stream<vec> &out_stream_1_1, hls::stream<vec> &out_stream_1_2, hls::stream<vec> &out_stream_1_3, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &in_stream_2, hls::stream<vec> &in_stream_3){

	int i, j;
	vec tmp_0, tmp_1, tmp_2, tmp_3;

	for (j = 0; j < GFBITS; j++){
		tmp_0 = in_stream_0.read();
		tmp_1 = in_stream_1.read();
		tmp_2 = in_stream_2.read();
		tmp_3 = in_stream_3.read();
		out_stream_0_0.write(tmp_0);
		out_stream_1_0.write(tmp_0);
		out_stream_0_1.write(tmp_1);
		out_stream_1_1.write(tmp_1);
		out_stream_0_2.write(tmp_2);
		out_stream_1_2.write(tmp_2);
		out_stream_0_3.write(tmp_3);
		out_stream_1_3.write(tmp_3);
	}

}

void error_stream_duplicate(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream){

	int i;
	vec tmp;

	for(i=0; i<128; i++){
		tmp = in_stream.read();
		out_stream_0.write(tmp);
		out_stream_1.write(tmp);
	}
}


void check_value_compute(int *check, hls::stream<bit> &check_synd, hls::stream<bit> &check_weight){

	bit check_value=0;
	bit one=1;

	check_value =  one - (check_synd.read() & check_weight.read());
	*check = (int)check_value;

}

static inline void radix_conversions_tr(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &out_stream_2, hls::stream<vec> &out_stream_3, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1, hls::stream<vec> &in_stream_2, hls::stream<vec> &in_stream_3)
{
	int i, j, k, l;

	const vec mask[6][2] = 
	{
		{0x2222222222222222, 0x4444444444444444},
		{0x0C0C0C0C0C0C0C0C, 0x3030303030303030},
		{0x00F000F000F000F0, 0x0F000F000F000F00},
		{0x0000FF000000FF00, 0x00FF000000FF0000},
		{0x00000000FFFF0000, 0x0000FFFF00000000},
		{0xFFFFFFFF00000000, 0x00000000FFFFFFFF}
	};

	const vec s[6][4][GFBITS] = 
	{
#include "scalars_4x.data"
	};
	
	vec in[4][GFBITS];

	//

	for (j = 6; j >= 0; j--)
	{
		if (j < 6)
		{
			vec_mul(in[0], in[0], s[j][0]); // scaling
			vec_mul(in[1], in[1], s[j][1]); // scaling
			vec_mul(in[2], in[2], s[j][2]); // scaling
			vec_mul(in[3], in[3], s[j][3]); // scaling
		}else{
			//Reading the inputs
			//TODO Do in 1 stream
			for(l=0; l<GFBITS; l++){
				in[0][l] = in_stream_0.read();
				in[1][l] = in_stream_1.read();
				in[2][l] = in_stream_2.read();
				in[3][l] = in_stream_3.read();
			}
		}

		for (k = j; k <= 4; k++)
		for (i = 0; i < GFBITS; i++)
		{
			in[0][i] ^= (in[0][i] & mask[k][0]) << (1 << k);
			in[0][i] ^= (in[0][i] & mask[k][1]) << (1 << k);
			in[1][i] ^= (in[1][i] & mask[k][0]) << (1 << k);
			in[1][i] ^= (in[1][i] & mask[k][1]) << (1 << k);
			in[2][i] ^= (in[2][i] & mask[k][0]) << (1 << k);
			in[2][i] ^= (in[2][i] & mask[k][1]) << (1 << k);
			in[3][i] ^= (in[3][i] & mask[k][0]) << (1 << k);
			in[3][i] ^= (in[3][i] & mask[k][1]) << (1 << k);
		}

		if (j <= 5)
		for (i = 0; i < GFBITS; i++)
		{
			in[1][i] ^= in[0][i] >> 32;
			in[1][i] ^= in[1][i] << 32;

			in[3][i] ^= in[2][i] >> 32;
			in[3][i] ^= in[3][i] << 32;
		}

		for (i = 0; i < GFBITS; i++){
//			 in[3][i] ^= in[2][i] ^= in[1][i];
			//The same functionality?
			in[2][i] ^= in[1][i];
			in[3][i] ^= in[2][i];
		}
	}

	for (i = 0; i < GFBITS; i++){
		out_stream_0.write(in[0][i]);
		out_stream_1.write(in[1][i]);
		out_stream_2.write(in[2][i]);
		out_stream_3.write(0);
//		printf("in[0][%d] = %lX\n", i, in[0][i].to_uint64());
//		printf("in[1][%d] = %lX\n", i, in[1][i].to_uint64());
//		printf("in[2][%d] = %lX\n", i, in[2][i].to_uint64());
//		printf("in[3][%d] = %lX\n", i, 0);
	}

}


static void butterflies_tr(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &out_stream_2, hls::stream<vec> &out_stream_3, hls::stream<mat_struct> &in_stream)
{
	int i, j, k, s, b;

	vec tmp[ GFBITS ];
	vec pre[6][2][ GFBITS ];
	vec buf[2][64];
	vec buf_tr[2][64];

	const vec consts[ 128 ][ GFBITS ] =
	{
#include "consts.data"
	};

	vec consts_ptr = 128;

	const unsigned char reversal[128] = 
	{ 
	  0, 64, 32, 96, 16, 80, 48, 112, 
	  8, 72, 40, 104, 24, 88, 56, 120, 
	  4, 68, 36, 100, 20, 84, 52, 116, 
	  12, 76, 44, 108, 28, 92, 60, 124, 
	  2, 66, 34, 98, 18, 82, 50, 114, 
	  10, 74, 42, 106, 26, 90, 58, 122, 
	  6, 70, 38, 102, 22, 86, 54, 118, 
	  14, 78, 46, 110, 30, 94, 62, 126, 
	  1, 65, 33, 97, 17, 81, 49, 113, 
	  9, 73, 41, 105, 25, 89, 57, 121, 
	  5, 69, 37, 101, 21, 85, 53, 117, 
	  13, 77, 45, 109, 29, 93, 61, 125, 
	  3, 67, 35, 99, 19, 83, 51, 115, 
	  11, 75, 43, 107, 27, 91, 59, 123, 
	  7, 71, 39, 103, 23, 87, 55, 119, 
	  15, 79, 47, 111, 31, 95, 63, 127
	};

	const uint16_t beta[6] = {5246, 5306, 6039, 6685, 4905, 6755};

	mat_struct tmp_struct;

	vec in[ 128 ][ GFBITS ];
	vec out[ 4 ][ GFBITS ];

	//Reading the streaming input
	for (i = 0; i < 128; i++){
		tmp_struct = in_stream.read();
		for (j = 0; j < GFBITS; j++){
			in[i][j] = tmp_struct.mat[j];
		}
	}


//TODO implemente the butterflies after initial check

	for (i = 6; i >= 0; i--)
	{
		s = 1 << i;
		consts_ptr -= s;

		for (j = 0; j < 128; j += 2*s)
		for (k = j; k < j+s; k++)
		{
			for (b = 0; b < GFBITS; b++){ 
				in[k][b] ^= in[k+s][b];
			}

			vec_mul(tmp, in[k], consts[ consts_ptr + (k-j) ]);

			for (b = 0; b < GFBITS; b++){
				in[k+s][b] ^= tmp[b];
			}

		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		for (k = 0; k < 128; k++)
			(&buf[0][0])[ k ] = in[ reversal[k] ][i];

		transpose_64x64(buf_tr[0], buf[0]);
		transpose_64x64(buf_tr[1], buf[1]);

		vec_copy64(buf[0], buf_tr[0]);
		vec_copy64(buf[1], buf_tr[1]);

		for (k = 0; k < 2; k++)
		{
			pre[0][k][i] = buf[k][32]; buf[k][33] ^= buf[k][32];
			pre[1][k][i] = buf[k][33]; buf[k][35] ^= buf[k][33];
			pre[0][k][i] ^= buf[k][35]; buf[k][34] ^= buf[k][35];
			pre[2][k][i] = buf[k][34]; buf[k][38] ^= buf[k][34];
			pre[0][k][i] ^= buf[k][38]; buf[k][39] ^= buf[k][38];
			pre[1][k][i] ^= buf[k][39]; buf[k][37] ^= buf[k][39];
			pre[0][k][i] ^= buf[k][37]; buf[k][36] ^= buf[k][37];
			pre[3][k][i] = buf[k][36]; buf[k][44] ^= buf[k][36];
			pre[0][k][i] ^= buf[k][44]; buf[k][45] ^= buf[k][44];
			pre[1][k][i] ^= buf[k][45]; buf[k][47] ^= buf[k][45];
			pre[0][k][i] ^= buf[k][47]; buf[k][46] ^= buf[k][47];
			pre[2][k][i] ^= buf[k][46]; buf[k][42] ^= buf[k][46];
			pre[0][k][i] ^= buf[k][42]; buf[k][43] ^= buf[k][42];
			pre[1][k][i] ^= buf[k][43]; buf[k][41] ^= buf[k][43];
			pre[0][k][i] ^= buf[k][41]; buf[k][40] ^= buf[k][41];
			pre[4][k][i] = buf[k][40]; buf[k][56] ^= buf[k][40];
			pre[0][k][i] ^= buf[k][56]; buf[k][57] ^= buf[k][56];
			pre[1][k][i] ^= buf[k][57]; buf[k][59] ^= buf[k][57];
			pre[0][k][i] ^= buf[k][59]; buf[k][58] ^= buf[k][59];
			pre[2][k][i] ^= buf[k][58]; buf[k][62] ^= buf[k][58];
			pre[0][k][i] ^= buf[k][62]; buf[k][63] ^= buf[k][62];
			pre[1][k][i] ^= buf[k][63]; buf[k][61] ^= buf[k][63];
			pre[0][k][i] ^= buf[k][61]; buf[k][60] ^= buf[k][61];
			pre[3][k][i] ^= buf[k][60]; buf[k][52] ^= buf[k][60];
			pre[0][k][i] ^= buf[k][52]; buf[k][53] ^= buf[k][52];
			pre[1][k][i] ^= buf[k][53]; buf[k][55] ^= buf[k][53];
			pre[0][k][i] ^= buf[k][55]; buf[k][54] ^= buf[k][55];
			pre[2][k][i] ^= buf[k][54]; buf[k][50] ^= buf[k][54];
			pre[0][k][i] ^= buf[k][50]; buf[k][51] ^= buf[k][50];
			pre[1][k][i] ^= buf[k][51]; buf[k][49] ^= buf[k][51];
			pre[0][k][i] ^= buf[k][49]; buf[k][48] ^= buf[k][49];
			pre[5][k][i] = buf[k][48]; buf[k][16] ^= buf[k][48];
			pre[0][k][i] ^= buf[k][16]; buf[k][17] ^= buf[k][16];
			pre[1][k][i] ^= buf[k][17]; buf[k][19] ^= buf[k][17];
			pre[0][k][i] ^= buf[k][19]; buf[k][18] ^= buf[k][19];
			pre[2][k][i] ^= buf[k][18]; buf[k][22] ^= buf[k][18];
			pre[0][k][i] ^= buf[k][22]; buf[k][23] ^= buf[k][22];
			pre[1][k][i] ^= buf[k][23]; buf[k][21] ^= buf[k][23];
			pre[0][k][i] ^= buf[k][21]; buf[k][20] ^= buf[k][21];
			pre[3][k][i] ^= buf[k][20]; buf[k][28] ^= buf[k][20];
			pre[0][k][i] ^= buf[k][28]; buf[k][29] ^= buf[k][28];
			pre[1][k][i] ^= buf[k][29]; buf[k][31] ^= buf[k][29];
			pre[0][k][i] ^= buf[k][31]; buf[k][30] ^= buf[k][31];
			pre[2][k][i] ^= buf[k][30]; buf[k][26] ^= buf[k][30];
			pre[0][k][i] ^= buf[k][26]; buf[k][27] ^= buf[k][26];
			pre[1][k][i] ^= buf[k][27]; buf[k][25] ^= buf[k][27];
			pre[0][k][i] ^= buf[k][25]; buf[k][24] ^= buf[k][25];
			pre[4][k][i] ^= buf[k][24]; buf[k][8] ^= buf[k][24];
			pre[0][k][i] ^= buf[k][8]; buf[k][9] ^= buf[k][8];
			pre[1][k][i] ^= buf[k][9]; buf[k][11] ^= buf[k][9];
			pre[0][k][i] ^= buf[k][11]; buf[k][10] ^= buf[k][11];
			pre[2][k][i] ^= buf[k][10]; buf[k][14] ^= buf[k][10];
			pre[0][k][i] ^= buf[k][14]; buf[k][15] ^= buf[k][14];
			pre[1][k][i] ^= buf[k][15]; buf[k][13] ^= buf[k][15];
			pre[0][k][i] ^= buf[k][13]; buf[k][12] ^= buf[k][13];
			pre[3][k][i] ^= buf[k][12]; buf[k][4] ^= buf[k][12];
			pre[0][k][i] ^= buf[k][4]; buf[k][5] ^= buf[k][4];
			pre[1][k][i] ^= buf[k][5]; buf[k][7] ^= buf[k][5];
			pre[0][k][i] ^= buf[k][7]; buf[k][6] ^= buf[k][7];
			pre[2][k][i] ^= buf[k][6]; buf[k][2] ^= buf[k][6];
			pre[0][k][i] ^= buf[k][2]; buf[k][3] ^= buf[k][2];
			pre[1][k][i] ^= buf[k][3]; buf[k][1] ^= buf[k][3];
		
			pre[0][k][i] ^= buf[k][1]; out[k][i] = buf[k][0] ^ buf[k][1];
		}
	}	

	for (j = 0; j < GFBITS; j++){
		tmp[j] = vec_setbits((beta[0] >> j) & 1);
	}

	vec_mul(out[2], pre[0][0], tmp);
	vec_mul(out[3], pre[0][1], tmp);

	for (i = 1; i < 6; i++)
	{
		for (j = 0; j < GFBITS; j++){ 
			tmp[j] = vec_setbits((beta[i] >> j) & 1);
		}

		vec_mul(pre[i][0], pre[i][0], tmp);
		vec_mul(pre[i][1], pre[i][1], tmp);

		for (b = 0; b < GFBITS; b++) 
		{
			out[2][b] ^= pre[i][0][b];
			out[3][b] ^= pre[i][1][b];
		}
	}

	//Write the 4 output streams
	for (i = 0; i < GFBITS; i++){
		out_stream_0.write(out[0][i]);
		out_stream_1.write(out[1][i]);
		out_stream_2.write(out[2][i]);
		out_stream_3.write(out[3][i]);
//		printf("out[0][%d] = %lX\n", i, out[0][i].to_uint64());
//		printf("out[1][%d] = %lX\n", i, out[1][i].to_uint64());
//		printf("out[2][%d] = %lX\n", i, out[2][i].to_uint64());
//		printf("out[3][%d] = %lX\n", i, out[3][i].to_uint64());

	}


}


//Interface is OK
void decryption_kernel(int *check, unsigned char *e, sk_05k_vec *sk_irr, sk_4k_vec *sk_rev, sk_4k_vec *sk_fwd, unsigned char *s){

#pragma HLS INTERFACE m_axi     port=e       offset=slave depth=576 bundle=gmem1
#pragma HLS INTERFACE m_axi     port=sk_irr  offset=slave depth=3 bundle=gmem2
#pragma HLS INTERFACE m_axi     port=sk_rev  offset=slave depth=25 bundle=gmem3
#pragma HLS INTERFACE m_axi     port=sk_fwd  offset=slave depth=25 bundle=gmem4
#pragma HLS INTERFACE m_axi     port=s       offset=slave depth=156 bundle=gmem5
#pragma HLS INTERFACE m_axi     port=check   offset=slave depth=1 bundle=gmem6
#pragma HLS INTERFACE s_axilite port=e       bundle=control
#pragma HLS INTERFACE s_axilite port=sk_irr  bundle=control
#pragma HLS INTERFACE s_axilite port=sk_rev  bundle=control
#pragma HLS INTERFACE s_axilite port=sk_fwd  bundle=control
#pragma HLS INTERFACE s_axilite port=s       bundle=control
#pragma HLS INTERFACE s_axilite port=check   bundle=control
#pragma HLS INTERFACE s_axilite port=return  bundle=control



	//MAIN FUNCTION
	static hls::stream<unsigned char> s_mat_stream("s_mat_stream");
	#pragma HLS STREAM variable=s_mat_stream depth=156

	static hls::stream<sk_4k_vec> sk_vec_rev_mat_stream("sk_vec_rev_mat_stream");
	#pragma HLS STREAM variable=sk_vec_rev_mat_stream depth=25

	static hls::stream<sk_4k_vec> sk_vec_fwd_mat_stream("sk_vec_fwd_mat_stream");
	#pragma HLS STREAM variable=sk_vec_fwd_mat_stream depth=25

	static hls::stream<sk_05k_vec> sk_irr_vec_stream("sk_irr_vec_stream");
	#pragma HLS STREAM variable=sk_irr_vec_stream depth=3

	static hls::stream<vec> recv_stream_preprocess("recv_stream_preprocess");
	#pragma HLS STREAM variable=recv_stream_preprocess depth=128

	static hls::stream<vec> recv_stream_benes("recv_stream_benes");
	#pragma HLS STREAM variable=recv_stream_benes depth=128

	static hls::stream<vec> inv_stream("inv_stream");
	#pragma HLS STREAM variable=inv_stream depth=1664

	static hls::stream<mat_struct> scaled_struct_stream_0("scaled_struct_stream_0");
	#pragma HLS STREAM variable=scaled_struct_stream_0 depth=64//2

	static hls::stream<mat_struct> scaled_struct_stream_3("scaled_struct_stream_3");
	#pragma HLS STREAM variable=scaled_struct_stream_3 depth=64

	static hls::stream<vec> s_priv_stream_0("s_priv_stream_0");
	#pragma HLS STREAM variable=s_priv_stream_0 depth=13

	static hls::stream<vec> s_priv_stream_1("s_priv_stream_1");
	#pragma HLS STREAM variable=s_priv_stream_1 depth=13

	static hls::stream<vec> s_priv_stream_2("s_priv_stream_2");
	#pragma HLS STREAM variable=s_priv_stream_2 depth=13

	static hls::stream<vec> s_priv_stream_3("s_priv_stream_3");
	#pragma HLS STREAM variable=s_priv_stream_2 depth=13

	static hls::stream<vec> s_priv_stream_0_0("s_priv_stream_0_0");
	#pragma HLS STREAM variable=s_priv_stream_0_0 depth=13

	static hls::stream<vec> s_priv_stream_0_1("s_priv_stream_0_1");
	#pragma HLS STREAM variable=s_priv_stream_0_1 depth=13

	static hls::stream<vec> s_priv_stream_0_2("s_priv_stream_0_2");
	#pragma HLS STREAM variable=s_priv_stream_0_2 depth=13

	static hls::stream<vec> s_priv_stream_0_3("s_priv_stream_0_3");
	#pragma HLS STREAM variable=s_priv_stream_0_3 depth=13

	static hls::stream<vec> s_priv_stream_1_0("s_priv_stream_1_0");
	#pragma HLS STREAM variable=s_priv_stream_1_0 depth=13

	static hls::stream<vec> s_priv_stream_1_1("s_priv_stream_1_1");
	#pragma HLS STREAM variable=s_priv_stream_1_1 depth=13

	static hls::stream<vec> s_priv_stream_1_2("s_priv_stream_1_2");
	#pragma HLS STREAM variable=s_priv_stream_1_2 depth=13

	static hls::stream<vec> s_priv_stream_1_3("s_priv_stream_1_3");
	#pragma HLS STREAM variable=s_priv_stream_1_3 depth=13

	static hls::stream<vec> s_priv_stream_2_0("s_priv_stream_2_0");
	#pragma HLS STREAM variable=s_priv_stream_2_0 depth=13

	static hls::stream<vec> s_priv_stream_2_1("s_priv_stream_2_1");
	#pragma HLS STREAM variable=s_priv_stream_2_1 depth=13

	static hls::stream<vec> s_priv_stream_3_0("s_priv_stream_3_0");
	#pragma HLS STREAM variable=s_priv_stream_3_0 depth=13

	static hls::stream<vec> s_priv_stream_3_1("s_priv_stream_3_1");
	#pragma HLS STREAM variable=s_priv_stream_3_1 depth=13

	static hls::stream<vec> s_priv_cmp_stream_0("s_priv_cmp_stream_0");
	#pragma HLS STREAM variable=s_priv_cmp_stream_0 depth=13

	static hls::stream<vec> s_priv_cmp_stream_1("s_priv_cmp_stream_1");
	#pragma HLS STREAM variable=s_priv_cmp_stream_1 depth=13

	static hls::stream<vec> s_priv_cmp_stream_2("s_priv_cmp_stream_2");
	#pragma HLS STREAM variable=s_priv_cmp_stream_2 depth=13

	static hls::stream<vec> s_priv_cmp_stream_3("s_priv_cmp_stream_3");
	#pragma HLS STREAM variable=s_priv_cmp_stream_3 depth=13

	static hls::stream<vec> locator_stream_0("locator_stream_0");
	#pragma HLS STREAM variable=locator_stream_0 depth=13

	static hls::stream<vec> locator_stream_1("locator_stream_1");
	#pragma HLS STREAM variable=locator_stream_1 depth=13

	static hls::stream<vec> eval_stream("eval_stream");
	#pragma HLS STREAM variable=eval_stream depth=768

	static hls::stream<vec> error_stream("error_stream");
	#pragma HLS STREAM variable=error_stream depth=64

	static hls::stream<vec> error_stream_0("error_stream_0");
	#pragma HLS STREAM variable=error_stream_0 depth=64

	static hls::stream<vec> error_stream_1("error_stream_1");
	#pragma HLS STREAM variable=error_stream_1 depth=64

	static hls::stream<vec> error_stream_benes_out("error_stream_benes_out");
	#pragma HLS STREAM variable=error_stream_benes_out depth=64

	static hls::stream<vec> error_stream_benes_out_0("error_stream_benes_out_0");
	#pragma HLS STREAM variable=error_stream_benes_out_0 depth=64

	static hls::stream<vec> error_stream_benes_out_1("error_stream_benes_out_1");
	#pragma HLS STREAM variable=error_stream_benes_out_1 depth=64

	static hls::stream<unsigned char> e_stream_int("e_stream_int");
	#pragma HLS STREAM variable=e_stream_int depth=436

	//SCALING
	static hls::stream<vec> irr_int_stream_0("irr_int_stream_0");
	#pragma HLS STREAM variable=irr_int_stream_0 depth=13

	static hls::stream<vec> irr_int_stream_1("irr_int_stream_1");
	#pragma HLS STREAM variable=irr_int_stream_1 depth=13

	static hls::stream<vec> eval_in_stream("eval_in_stream");
	#pragma HLS STREAM variable=eval_in_stream depth=768

	static hls::stream<vec> eval_out_stream("eval_out_stream");
	#pragma HLS STREAM variable=eval_out_stream depth=768

	static hls::stream<vec> eval_out_stream_fwd("eval_out_stream_fwd");
	#pragma HLS STREAM variable=eval_out_stream_fwd depth=768

	static hls::stream<vec> eval_out_stream_rev("eval_out_stream_rev");
	#pragma HLS STREAM variable=eval_out_stream_rev depth=768

	static hls::stream<vec> inv_int_rev_stream("inv_int_rev_stream");
	#pragma HLS STREAM variable=inv_int_rev_stream depth=768

	static hls::stream<vec> inv_int_fwd_stream("inv_int_fwd_stream");
	#pragma HLS STREAM variable=inv_int_fwd_stream depth=768

	//RADIX_CONVERSIONS_TR #1
	static hls::stream<vec> out_stream_5_0_0("out_stream_5_0_0");
	#pragma HLS STREAM variable=out_stream_5_0_0 depth=13

	static hls::stream<vec> out_stream_5_1_0("out_stream_5_1_0");
	#pragma HLS STREAM variable=out_stream_5_1_0 depth=13

	static hls::stream<vec> out_stream_4_0_0("out_stream_4_0_0");
	#pragma HLS STREAM variable=out_stream_4_0_0 depth=13

	static hls::stream<vec> out_stream_4_1_0("out_stream_4_1_0");
	#pragma HLS STREAM variable=out_stream_4_1_0 depth=13

	static hls::stream<vec> out_stream_3_0_0("out_stream_3_0_0");
	#pragma HLS STREAM variable=out_stream_3_0_0 depth=13

	static hls::stream<vec> out_stream_3_1_0("out_stream_3_1_0");
	#pragma HLS STREAM variable=out_stream_3_1_0 depth=13

	static hls::stream<vec> out_stream_2_0_0("out_stream_2_0_0");
	#pragma HLS STREAM variable=out_stream_2_0_0 depth=13

	static hls::stream<vec> out_stream_2_1_0("out_stream_2_1_0");
	#pragma HLS STREAM variable=out_stream_2_1_0 depth=13

	static hls::stream<vec> out_stream_1_0_0("out_stream_1_0_0");
	#pragma HLS STREAM variable=out_stream_2_0_0 depth=13

	static hls::stream<vec> out_stream_1_1_0("out_stream_1_1_0");
	#pragma HLS STREAM variable=out_stream_1_1_0 depth=13

	static hls::stream<vec> s_stream_4_0_0("s_stream_4_0_0");
	#pragma HLS STREAM variable=s_stream_4_0_0 depth=144

	static hls::stream<vec> s_stream_4_1_0("s_stream_4_1_0");
	#pragma HLS STREAM variable=s_stream_4_1_0 depth=144

	static hls::stream<vec> s_stream_3_0_0("s_stream_3_0_0");
	#pragma HLS STREAM variable=s_stream_3_0_0 depth=144

	static hls::stream<vec> s_stream_3_1_0("s_stream_3_1_0");
	#pragma HLS STREAM variable=s_stream_3_1_0 depth=144

	static hls::stream<vec> s_stream_2_0_0("s_stream_2_0_0");
	#pragma HLS STREAM variable=s_stream_2_0_0 depth=144

	static hls::stream<vec> s_stream_2_1_0("s_stream_2_1_0");
	#pragma HLS STREAM variable=s_stream_2_1_0 depth=144

	static hls::stream<vec> s_stream_1_0_0("s_stream_1_0_0");
	#pragma HLS STREAM variable=s_stream_1_0_0 depth=144

	static hls::stream<vec> s_stream_1_1_0("s_stream_1_1_0");
	#pragma HLS STREAM variable=s_stream_1_1_0 depth=144

	static hls::stream<vec> s_stream_0_0_0("s_stream_0_0_0");
	#pragma HLS STREAM variable=s_stream_0_0_0 depth=144

	static hls::stream<vec> s_stream_0_1_0("s_stream_0_1_0");
	#pragma HLS STREAM variable=s_stream_0_1_0 depth=144

	static hls::stream<vec> mask_stream_5_0_0("mask_stream_5_0_0");
	#pragma HLS STREAM variable=mask_stream_5_0_0 depth=13

	static hls::stream<vec> mask_stream_5_1_0("mask_stream_5_1_0");
	#pragma HLS STREAM variable=mask_stream_5_1_0 depth=13

	static hls::stream<vec> mask_stream_4_0_0("mask_stream_4_0_0");
	#pragma HLS STREAM variable=mask_stream_4_0_0 depth=36

	static hls::stream<vec> mask_stream_4_1_0("mask_stream_4_1_0");
	#pragma HLS STREAM variable=mask_stream_4_1_0 depth=36

	static hls::stream<vec> mask_stream_3_0_0("mask_stream_3_0_0");
	#pragma HLS STREAM variable=mask_stream_3_0_0 depth=60

	static hls::stream<vec> mask_stream_3_1_0("mask_stream_3_1_0");
	#pragma HLS STREAM variable=mask_stream_3_1_0 depth=60

	static hls::stream<vec> mask_stream_2_0_0("mask_stream_2_0_0");
	#pragma HLS STREAM variable=mask_stream_2_0_0 depth=84

	static hls::stream<vec> mask_stream_2_1_0("mask_stream_2_1_0");
	#pragma HLS STREAM variable=mask_stream_2_1_0 depth=84

	static hls::stream<vec> mask_stream_1_0_0("mask_stream_1_0_0");
	#pragma HLS STREAM variable=mask_stream_1_0_0 depth=108

	static hls::stream<vec> mask_stream_1_1_0("mask_stream_1_1_0");
	#pragma HLS STREAM variable=mask_stream_1_1_0 depth=108

	static hls::stream<vec> mask_stream_0_0_0("mask_stream_0_0_0");
	#pragma HLS STREAM variable=mask_stream_0_0_0 depth=132

	static hls::stream<vec> mask_stream_0_1_0("mask_stream_0_1_0");
	#pragma HLS STREAM variable=mask_stream_0_1_0 depth=132

	static hls::stream<vec> int_stream_0_0("int_stream_0_0");
	#pragma HLS STREAM variable=int_stream_0_0 depth=13

	static hls::stream<vec> int_stream_1_0("int_stream_1_0");
	#pragma HLS STREAM variable=int_stream_1_0 depth=13

	static hls::stream<vec> int_stream_2_0("int_stream_2_0");
	#pragma HLS STREAM variable=int_stream_2_0 depth=13

	static hls::stream<vec> int_stream_3_0("int_stream_3_0");
	#pragma HLS STREAM variable=int_stream_3_0 depth=13


	//RADIX_CONVERSIONS_TR #2
	static hls::stream<vec> int_stream_0_1("int_stream_0_1");
	#pragma HLS STREAM variable=int_stream_0_1 depth=13

	static hls::stream<vec> int_stream_1_1("int_stream_1_1");
	#pragma HLS STREAM variable=int_stream_1_1 depth=13

	static hls::stream<vec> int_stream_2_1("int_stream_2_1");
	#pragma HLS STREAM variable=int_stream_2_1 depth=13

	static hls::stream<vec> int_stream_3_1("int_stream_3_1");
	#pragma HLS STREAM variable=int_stream_3_1 depth=13

	static hls::stream<vec> out_stream_5_0_1("out_stream_5_0_1");
	#pragma HLS STREAM variable=out_stream_5_0_0 depth=13

	static hls::stream<vec> out_stream_5_1_1("out_stream_5_1_1");
	#pragma HLS STREAM variable=out_stream_5_1_1 depth=13

	static hls::stream<vec> out_stream_4_0_1("out_stream_4_0_1");
	#pragma HLS STREAM variable=out_stream_4_0_1 depth=13

	static hls::stream<vec> out_stream_4_1_1("out_stream_4_1_1");
	#pragma HLS STREAM variable=out_stream_4_1_1 depth=13

	static hls::stream<vec> out_stream_3_0_1("out_stream_3_0_1");
	#pragma HLS STREAM variable=out_stream_3_0_1 depth=13

	static hls::stream<vec> out_stream_3_1_1("out_stream_3_1_1");
	#pragma HLS STREAM variable=out_stream_3_1_1 depth=13

	static hls::stream<vec> out_stream_2_0_1("out_stream_2_0_1");
	#pragma HLS STREAM variable=out_stream_2_0_1 depth=13

	static hls::stream<vec> out_stream_2_1_1("out_stream_2_1_1");
	#pragma HLS STREAM variable=out_stream_2_1_1 depth=13

	static hls::stream<vec> out_stream_1_0_1("out_stream_1_0_1");
	#pragma HLS STREAM variable=out_stream_2_0_1 depth=13

	static hls::stream<vec> out_stream_1_1_1("out_stream_1_1_1");
	#pragma HLS STREAM variable=out_stream_1_1_1 depth=13

	static hls::stream<vec> s_stream_4_0_1("s_stream_4_0_1");
	#pragma HLS STREAM variable=s_stream_4_0_1 depth=144

	static hls::stream<vec> s_stream_4_1_1("s_stream_4_1_1");
	#pragma HLS STREAM variable=s_stream_4_1_1 depth=144

	static hls::stream<vec> s_stream_3_0_1("s_stream_3_0_1");
	#pragma HLS STREAM variable=s_stream_3_0_1 depth=144

	static hls::stream<vec> s_stream_3_1_1("s_stream_3_1_1");
	#pragma HLS STREAM variable=s_stream_3_1_1 depth=144

	static hls::stream<vec> s_stream_2_0_1("s_stream_2_0_1");
	#pragma HLS STREAM variable=s_stream_2_0_1 depth=144

	static hls::stream<vec> s_stream_2_1_1("s_stream_2_1_1");
	#pragma HLS STREAM variable=s_stream_2_1_1 depth=144

	static hls::stream<vec> s_stream_1_0_1("s_stream_1_0_1");
	#pragma HLS STREAM variable=s_stream_1_0_1 depth=144

	static hls::stream<vec> s_stream_1_1_1("s_stream_1_1_1");
	#pragma HLS STREAM variable=s_stream_1_1_1 depth=144

	static hls::stream<vec> s_stream_0_0_1("s_stream_0_0_1");
	#pragma HLS STREAM variable=s_stream_0_0_1 depth=144

	static hls::stream<vec> s_stream_0_1_1("s_stream_0_1_1");
	#pragma HLS STREAM variable=s_stream_0_1_1 depth=144

	static hls::stream<vec> mask_stream_5_0_1("mask_stream_5_0_1");
	#pragma HLS STREAM variable=mask_stream_5_0_1 depth=13

	static hls::stream<vec> mask_stream_5_1_1("mask_stream_5_1_1");
	#pragma HLS STREAM variable=mask_stream_5_1_1 depth=13

	static hls::stream<vec> mask_stream_4_0_1("mask_stream_4_0_1");
	#pragma HLS STREAM variable=mask_stream_4_0_1 depth=36

	static hls::stream<vec> mask_stream_4_1_1("mask_stream_4_1_1");
	#pragma HLS STREAM variable=mask_stream_4_1_1 depth=36

	static hls::stream<vec> mask_stream_3_0_1("mask_stream_3_0_1");
	#pragma HLS STREAM variable=mask_stream_3_0_1 depth=60

	static hls::stream<vec> mask_stream_3_1_1("mask_stream_3_1_1");
	#pragma HLS STREAM variable=mask_stream_3_1_1 depth=60

	static hls::stream<vec> mask_stream_2_0_1("mask_stream_2_0_1");
	#pragma HLS STREAM variable=mask_stream_2_0_1 depth=84

	static hls::stream<vec> mask_stream_2_1_1("mask_stream_2_1_1");
	#pragma HLS STREAM variable=mask_stream_2_1_1 depth=84

	static hls::stream<vec> mask_stream_1_0_1("mask_stream_1_0_1");
	#pragma HLS STREAM variable=mask_stream_1_0_1 depth=108

	static hls::stream<vec> mask_stream_1_1_1("mask_stream_1_1_1");
	#pragma HLS STREAM variable=mask_stream_1_1_1 depth=108

	static hls::stream<vec> mask_stream_0_0_1("mask_stream_0_0_1");
	#pragma HLS STREAM variable=mask_stream_0_0_1 depth=132

	static hls::stream<vec> mask_stream_0_1_1("mask_stream_0_1_1");
	#pragma HLS STREAM variable=mask_stream_0_1_1 depth=132

	//BUTTERFLIES
	//#0
	static hls::stream<mat_struct> struct_stream_1_0_0("struct_stream_1_0_0");
	#pragma HLS STREAM variable=struct_stream_1_0_0 depth=64

	static hls::stream<mat_struct> struct_stream_1_0_1("struct_stream_1_0_1");
	#pragma HLS STREAM variable=struct_stream_1_0_1 depth=64

	static hls::stream<mat_struct> struct_stream_1_1_0("struct_stream_1_1_0");
	#pragma HLS STREAM variable=struct_stream_1_1_0 depth=64

	static hls::stream<mat_struct> struct_stream_1_1_1("struct_stream_1_1_1");
	#pragma HLS STREAM variable=struct_stream_1_1_1 depth=64

	static hls::stream<mat_struct> struct_stream_1_2_0("struct_stream_1_2_0");
	#pragma HLS STREAM variable=struct_stream_1_2_0 depth=64

	static hls::stream<mat_struct> struct_stream_1_2_1("struct_stream_1_2_1");
	#pragma HLS STREAM variable=struct_stream_1_2_1 depth=64

	static hls::stream<mat_struct> struct_stream_1_3_0("struct_stream_1_3_0");
	#pragma HLS STREAM variable=struct_stream_1_3_0 depth=64

	static hls::stream<mat_struct> struct_stream_1_3_1("struct_stream_1_3_1");
	#pragma HLS STREAM variable=struct_stream_1_3_1 depth=64

	static hls::stream<mat_struct> struct_stream_1_4_0("struct_stream_1_4_0");
	#pragma HLS STREAM variable=struct_stream_1_4_0 depth=64

	static hls::stream<mat_struct> struct_stream_1_4_1("struct_stream_1_4_1");
	#pragma HLS STREAM variable=struct_stream_1_4_1 depth=64

	static hls::stream<mat_struct> struct_stream_1_5_0("struct_stream_1_5_0");
	#pragma HLS STREAM variable=struct_stream_1_5_0 depth=64

	static hls::stream<mat_struct> struct_stream_1_5_1("struct_stream_1_5_1");
	#pragma HLS STREAM variable=struct_stream_1_5_1 depth=64

	static hls::stream<mat_struct> struct_stream_1_6_0("struct_stream_1_6_0");
	#pragma HLS STREAM variable=struct_stream_1_6_0 depth=64

	static hls::stream<mat_struct> struct_stream_1_6_1("struct_stream_1_6_1");
	#pragma HLS STREAM variable=struct_stream_1_6_1 depth=64

	static hls::stream<vec> internal_stream_0_0("internal_stream_0_0");
	#pragma HLS STREAM variable=internal_stream_0_0 depth=64

	static hls::stream<vec> internal_stream_0_1("internal_stream_0_1");
	#pragma HLS STREAM variable=internal_stream_0_1 depth=64

	//#1
	static hls::stream<mat_struct> struct_stream_3_0_0("struct_stream_3_0_0");
	#pragma HLS STREAM variable=struct_stream_3_0_0 depth=64

	static hls::stream<mat_struct> struct_stream_3_0_1("struct_stream_3_0_1");
	#pragma HLS STREAM variable=struct_stream_3_0_1 depth=64

	static hls::stream<mat_struct> struct_stream_3_1_0("struct_stream_3_1_0");
	#pragma HLS STREAM variable=struct_stream_3_1_0 depth=64

	static hls::stream<mat_struct> struct_stream_3_1_1("struct_stream_3_1_1");
	#pragma HLS STREAM variable=struct_stream_3_1_1 depth=64

	static hls::stream<mat_struct> struct_stream_3_2_0("struct_stream_3_2_0");
	#pragma HLS STREAM variable=struct_stream_3_2_0 depth=64

	static hls::stream<mat_struct> struct_stream_3_2_1("struct_stream_3_2_1");
	#pragma HLS STREAM variable=struct_stream_3_2_1 depth=64

	static hls::stream<mat_struct> struct_stream_3_3_0("struct_stream_3_3_0");
	#pragma HLS STREAM variable=struct_stream_3_3_0 depth=64

	static hls::stream<mat_struct> struct_stream_3_3_1("struct_stream_3_3_1");
	#pragma HLS STREAM variable=struct_stream_3_3_1 depth=64

	static hls::stream<mat_struct> struct_stream_3_4_0("struct_stream_3_4_0");
	#pragma HLS STREAM variable=struct_stream_3_4_0 depth=64

	static hls::stream<mat_struct> struct_stream_3_4_1("struct_stream_3_4_1");
	#pragma HLS STREAM variable=struct_stream_3_4_1 depth=64

	static hls::stream<mat_struct> struct_stream_3_5_0("struct_stream_3_5_0");
	#pragma HLS STREAM variable=struct_stream_3_5_0 depth=64

	static hls::stream<mat_struct> struct_stream_3_5_1("struct_stream_3_5_1");
	#pragma HLS STREAM variable=struct_stream_3_5_1 depth=64

	static hls::stream<mat_struct> struct_stream_3_6_0("struct_stream_3_6_0");
	#pragma HLS STREAM variable=struct_stream_3_6_0 depth=64

	static hls::stream<mat_struct> struct_stream_3_6_1("struct_stream_3_6_1");
	#pragma HLS STREAM variable=struct_stream_3_6_1 depth=64

	static hls::stream<vec> internal_stream_1_0("internal_stream_1_0");
	#pragma HLS STREAM variable=internal_stream_1_0 depth=64

	static hls::stream<vec> internal_stream_1_1("internal_stream_1_1");
	#pragma HLS STREAM variable=internal_stream_1_1 depth=64

	static hls::stream<bit> check_synd_stream("check_synd_stream");
	#pragma HLS STREAM variable=check_synd_stream depth=1

	static hls::stream<bit> check_weight_stream("check_weight_stream");
	#pragma HLS STREAM variable=check_weight_stream depth=1

	//FFT
	//#0

	static hls::stream<vec> int_mat_0_stream("int_mat_0_stream");
	#pragma HLS STREAM variable=int_mat_0_stream depth=13

	static hls::stream<mat_struct> struct_stream_0_0_0("struct_stream_0_0_0");
	#pragma HLS STREAM variable=struct_stream_0_0_0 depth=64

	static hls::stream<mat_struct> struct_stream_0_0_1("struct_stream_0_0_1");
	#pragma HLS STREAM variable=struct_stream_0_0_1 depth=64

	static hls::stream<mat_struct> struct_stream_0_1_0("struct_stream_0_1_0");
	#pragma HLS STREAM variable=struct_stream_0_1_0 depth=64

	static hls::stream<mat_struct> struct_stream_0_1_1("struct_stream_0_1_1");
	#pragma HLS STREAM variable=struct_stream_0_1_1 depth=64

	static hls::stream<mat_struct> struct_stream_0_2_0("struct_stream_0_2_0");
	#pragma HLS STREAM variable=struct_stream_0_2_0 depth=64

	static hls::stream<mat_struct> struct_stream_0_2_1("struct_stream_0_2_1");
	#pragma HLS STREAM variable=struct_stream_0_2_1 depth=64

	static hls::stream<mat_struct> struct_stream_0_3_0("struct_stream_0_3_0");
	#pragma HLS STREAM variable=struct_stream_0_3_0 depth=64

	static hls::stream<mat_struct> struct_stream_0_3_1("struct_stream_0_3_1");
	#pragma HLS STREAM variable=struct_stream_0_3_1 depth=64

	static hls::stream<mat_struct> struct_stream_0_4_0("struct_stream_0_4_0");
	#pragma HLS STREAM variable=struct_stream_0_4_0 depth=64

	static hls::stream<mat_struct> struct_stream_0_4_1("struct_stream_0_4_1");
	#pragma HLS STREAM variable=struct_stream_0_4_1 depth=64

	static hls::stream<mat_struct> struct_stream_0_5_0("struct_stream_0_5_0");
	#pragma HLS STREAM variable=struct_stream_0_5_0 depth=64

	static hls::stream<mat_struct> struct_stream_0_5_1("struct_stream_0_5_1");
	#pragma HLS STREAM variable=struct_stream_0_5_1 depth=64

	static hls::stream<mat_struct> struct_stream_0_6_0("struct_stream_0_6_0");
	#pragma HLS STREAM variable=struct_stream_0_6_0 depth=64

	static hls::stream<mat_struct> struct_stream_0_6_1("struct_stream_0_6_1");
	#pragma HLS STREAM variable=struct_stream_0_6_1 depth=64

	static hls::stream<mat_struct> scaled_struct_stream_1("scaled_struct_stream_1");
	#pragma HLS STREAM variable=scaled_struct_stream_1 depth=64

	//#1

	static hls::stream<vec> int_mat_1_stream("int_mat_1_stream");
	#pragma HLS STREAM variable=int_mat_1_stream depth=13

	static hls::stream<mat_struct> struct_stream_2_0_0("struct_stream_2_0_0");
	#pragma HLS STREAM variable=struct_stream_2_0_0 depth=64

	static hls::stream<mat_struct> struct_stream_2_0_1("struct_stream_2_0_1");
	#pragma HLS STREAM variable=struct_stream_2_0_1 depth=64

	static hls::stream<mat_struct> struct_stream_2_1_0("struct_stream_2_1_0");
	#pragma HLS STREAM variable=struct_stream_2_1_0 depth=64

	static hls::stream<mat_struct> struct_stream_2_1_1("struct_stream_2_1_1");
	#pragma HLS STREAM variable=struct_stream_2_1_1 depth=64

	static hls::stream<mat_struct> struct_stream_2_2_0("struct_stream_2_2_0");
	#pragma HLS STREAM variable=struct_stream_2_2_0 depth=64

	static hls::stream<mat_struct> struct_stream_2_2_1("struct_stream_2_2_1");
	#pragma HLS STREAM variable=struct_stream_2_2_1 depth=64

	static hls::stream<mat_struct> struct_stream_2_3_0("struct_stream_2_3_0");
	#pragma HLS STREAM variable=struct_stream_2_3_0 depth=64

	static hls::stream<mat_struct> struct_stream_2_3_1("struct_stream_2_3_1");
	#pragma HLS STREAM variable=struct_stream_2_3_1 depth=64

	static hls::stream<mat_struct> struct_stream_2_4_0("struct_stream_2_4_0");
	#pragma HLS STREAM variable=struct_stream_2_4_0 depth=64

	static hls::stream<mat_struct> struct_stream_2_4_1("struct_stream_2_4_1");
	#pragma HLS STREAM variable=struct_stream_2_4_1 depth=64

	static hls::stream<mat_struct> struct_stream_2_5_0("struct_stream_2_5_0");
	#pragma HLS STREAM variable=struct_stream_2_5_0 depth=64

	static hls::stream<mat_struct> struct_stream_2_5_1("struct_stream_2_5_1");
	#pragma HLS STREAM variable=struct_stream_2_5_1 depth=64

	static hls::stream<mat_struct> struct_stream_2_6_0("struct_stream_2_6_0");
	#pragma HLS STREAM variable=struct_stream_2_6_0 depth=64

	static hls::stream<mat_struct> struct_stream_2_6_1("struct_stream_2_6_1");
	#pragma HLS STREAM variable=struct_stream_2_6_1 depth=64

	static hls::stream<mat_struct> scaled_struct_stream_2("scaled_struct_stream_2");
	#pragma HLS STREAM variable=scaled_struct_stream_2 depth=64


	#pragma HLS dataflow

	//OK
	load_s(s_mat_stream, s);

	//TODO maybe read 1k also for the sk rev and fwd to avoid delays in the AXI interface
	//OK
	load_sk_vec(sk_vec_rev_mat_stream, sk_rev);
	//OK
	load_sk_vec(sk_vec_fwd_mat_stream, sk_fwd);
	//OK
	load_irr_vec(sk_irr_vec_stream, sk_irr);
	//OK
	preprocess(recv_stream_preprocess, s_mat_stream);
	//OK
	benes(recv_stream_benes, recv_stream_preprocess, sk_vec_rev_mat_stream);

	//////////////////////////SCALING FUNCTION///////////////////////////////
	//OK
	irr_load_stream(irr_int_stream_0, irr_int_stream_1, sk_irr_vec_stream);


	//////////////////////////////FFT////////////////////////////
	//OK But look at the functions that feed it to make sure they are writig the 2d matrx in 1, 0 order. And for the consumers to consume the right thing
	radix_conversions_stream(int_mat_0_stream, irr_int_stream_0, irr_int_stream_1);

	//OK TODO later decouple the butterflies
	butterflies_stream(eval_in_stream, int_mat_0_stream);


//	 broadcast(scaled_struct_stream_0, int_mat_0_stream);
//	 split_struct_stream(struct_stream_0_0_0, struct_stream_0_0_1, scaled_struct_stream_0);
//	 fft_layer(struct_stream_0_1_0, struct_stream_0_1_1, struct_stream_0_0_0, struct_stream_0_0_1, 2, 32);
//	 fft_layer(struct_stream_0_2_0, struct_stream_0_2_1, struct_stream_0_1_0, struct_stream_0_1_1, 4, 16);
//	 fft_layer(struct_stream_0_3_0, struct_stream_0_3_1, struct_stream_0_2_0, struct_stream_0_2_1, 8, 8);
//	 fft_layer(struct_stream_0_4_0, struct_stream_0_4_1, struct_stream_0_3_0, struct_stream_0_3_1, 16, 4);
//	 fft_layer(struct_stream_0_5_0, struct_stream_0_5_1, struct_stream_0_4_0, struct_stream_0_4_1, 32, 2);
//	 fft_layer(struct_stream_0_6_0, struct_stream_0_6_1, struct_stream_0_5_0, struct_stream_0_5_1, 64, 1);
//	 write_output(eval_in_stream, struct_stream_0_6_0, struct_stream_0_6_1);


	////////////////////////////END FFT//////////////////////////
	//OK
    vec_sq_loop(eval_out_stream, eval_in_stream);
	//OK
    eval_stream_duplicate(eval_out_stream_fwd, eval_out_stream_rev, eval_out_stream);
	//OK
    mul_loop_1(inv_int_rev_stream, eval_out_stream_fwd);
	//OK
    mul_loop_2(inv_int_fwd_stream, inv_int_rev_stream, eval_out_stream_rev);
	//OK
    compute_outputs(scaled_struct_stream_1, inv_stream, inv_int_fwd_stream, recv_stream_benes);
	////////////////////////END SCALING FUNCTION/////////////////////////////



	//////////////////////////////FFT_TR////////////////////////////////////
		//OK
		butterflies_tr(int_stream_0_0, int_stream_1_0, int_stream_2_0, int_stream_3_0, scaled_struct_stream_1);
		// //////////////////////////////////BUTTERFLIES_TR//////////////////////////////////////
		// split_struct_stream(struct_stream_1_0_0, struct_stream_1_0_1, scaled_struct_stream_1);
		// fft_tr_layer(struct_stream_1_1_0, struct_stream_1_1_1, struct_stream_1_0_0, struct_stream_1_0_1, 31, 1);
		// fft_tr_layer(struct_stream_1_2_0, struct_stream_1_2_1, struct_stream_1_1_0, struct_stream_1_1_1, 15, 2);
		// fft_tr_layer(struct_stream_1_3_0, struct_stream_1_3_1, struct_stream_1_2_0, struct_stream_1_2_1, 7, 4);
		// fft_tr_layer(struct_stream_1_4_0, struct_stream_1_4_1, struct_stream_1_3_0, struct_stream_1_3_1, 3, 8);
		// fft_tr_layer(struct_stream_1_5_0, struct_stream_1_5_1, struct_stream_1_4_0, struct_stream_1_4_1, 1, 16);
		// fft_tr_layer(struct_stream_1_6_0, struct_stream_1_6_1, struct_stream_1_5_0, struct_stream_1_5_1, 0, 32);
		// merge_tr_struct_stream(internal_stream_0_0, internal_stream_0_1, struct_stream_1_6_0, struct_stream_1_6_1);

		// transpose_broadcast(int_stream_0_0, int_stream_1_0, internal_stream_0_0, internal_stream_0_1);
		// /////////////////////////////////END BUTTERFLIES_TR///////////////////////////////////

		//The last is always zero
		radix_conversions_tr(s_priv_stream_0, s_priv_stream_1, s_priv_stream_2, s_priv_stream_3, int_stream_0_0, int_stream_1_0, int_stream_2_0, int_stream_3_0);
		// //////////////////////////////RADIX_CONVERSIONS_TR////////////////////////////////////
		// mask_streaming(mask_stream_5_0_0, mask_stream_5_1_0, mask_stream_4_0_0, mask_stream_4_1_0, mask_stream_3_0_0, mask_stream_3_1_0, mask_stream_2_0_0, mask_stream_2_1_0, mask_stream_1_0_0, mask_stream_1_1_0, mask_stream_0_0_0,  mask_stream_0_1_0, mask_0);
		// s_streaming(s_stream_4_0_0, s_stream_4_1_0, s_stream_3_0_0, s_stream_3_1_0, s_stream_2_0_0, s_stream_2_1_0, s_stream_1_0_0, s_stream_1_1_0, s_stream_0_0_0, s_stream_0_1_0, s_radix_0);

		// loop_function_5(out_stream_5_0_0, out_stream_5_1_0, int_stream_0_0, int_stream_1_0, mask_stream_5_0_0, mask_stream_5_1_0);
		// loop_function_4(out_stream_4_0_0, out_stream_4_1_0, out_stream_5_0_0, out_stream_5_1_0, mask_stream_4_0_0, mask_stream_4_1_0, s_stream_4_0_0, s_stream_4_1_0);
		// loop_function_3(out_stream_3_0_0, out_stream_3_1_0, out_stream_4_0_0, out_stream_4_1_0, mask_stream_3_0_0, mask_stream_3_1_0, s_stream_3_0_0, s_stream_3_1_0);
		// loop_function_2(out_stream_2_0_0, out_stream_2_1_0, out_stream_3_0_0, out_stream_3_1_0, mask_stream_2_0_0, mask_stream_2_1_0, s_stream_2_0_0, s_stream_2_1_0);
		// loop_function_1(out_stream_1_0_0, out_stream_1_1_0, out_stream_2_0_0, out_stream_2_1_0, mask_stream_1_0_0, mask_stream_1_1_0, s_stream_1_0_0, s_stream_1_1_0);
		// loop_function_0(s_priv_stream_0, s_priv_stream_1, out_stream_1_0_0, out_stream_1_1_0, mask_stream_0_0_0, mask_stream_0_1_0, s_stream_0_0_0, s_stream_0_1_0);
		// ////////////////////////////END RADIX_CONVERSIONS_TR//////////////////////////////////

	////////////////////////////END FFT_TR//////////////////////////////////

	//OK
	priv_stream_duplicate(s_priv_stream_0_0, s_priv_stream_0_1, s_priv_stream_0_2, s_priv_stream_0_3, s_priv_stream_1_0, s_priv_stream_1_1, s_priv_stream_1_2, s_priv_stream_1_3, s_priv_stream_0, s_priv_stream_1, s_priv_stream_2, s_priv_stream_3);

	//OK
	bm(locator_stream_0, locator_stream_1, s_priv_stream_0_0, s_priv_stream_0_1, s_priv_stream_0_2, s_priv_stream_0_3);


	//////////////////////////////FFT////////////////////////////
	//OKK
	radix_conversions_stream(int_mat_1_stream, locator_stream_0, locator_stream_1);


	butterflies_stream(eval_stream, int_mat_1_stream);

//	 broadcast(scaled_struct_stream_2, int_mat_1_stream);
//	 split_struct_stream(struct_stream_2_0_0, struct_stream_2_0_1, scaled_struct_stream_2);
//	 fft_layer(struct_stream_2_1_0, struct_stream_2_1_1, struct_stream_2_0_0, struct_stream_2_0_1, 2, 32);
//	 fft_layer(struct_stream_2_2_0, struct_stream_2_2_1, struct_stream_2_1_0, struct_stream_2_1_1, 4, 16);
//	 fft_layer(struct_stream_2_3_0, struct_stream_2_3_1, struct_stream_2_2_0, struct_stream_2_2_1, 8, 8);
//	 fft_layer(struct_stream_2_4_0, struct_stream_2_4_1, struct_stream_2_3_0, struct_stream_2_3_1, 16, 4);
//	 fft_layer(struct_stream_2_5_0, struct_stream_2_5_1, struct_stream_2_4_0, struct_stream_2_4_1, 32, 2);
//	 fft_layer(struct_stream_2_6_0, struct_stream_2_6_1, struct_stream_2_5_0, struct_stream_2_5_1, 64, 1);
//	 write_output(eval_stream, struct_stream_2_6_0, struct_stream_2_6_1);

	////////////////////////////END FFT//////////////////////////

	// reencryption and weight check
	//OK
	error_compute(error_stream, eval_stream);
	//OK
	error_stream_duplicate(error_stream_0, error_stream_1, error_stream);
	//OK
	scaling_inv(scaled_struct_stream_3, inv_stream, error_stream_0);

	//////////////////////////////FFT_TR////////////////////////////////////
//	fft_tr(s_priv_cmp_stream_0, s_priv_cmp_stream_1, scaled_struct_stream_2);

		butterflies_tr(int_stream_0_1, int_stream_1_1, int_stream_2_1, int_stream_3_1, scaled_struct_stream_3);
		//////////////////////////////////BUTTERFLIES_TR//////////////////////////////////////
//		 split_struct_stream_tr(struct_stream_3_0_0, struct_stream_3_0_1, scaled_struct_stream_3);
//		 fft_tr_layer(struct_stream_3_1_0, struct_stream_3_1_1, struct_stream_3_0_0, struct_stream_3_0_1, 64, 1);
//		 fft_tr_layer(struct_stream_3_2_0, struct_stream_3_2_1, struct_stream_3_1_0, struct_stream_3_1_1, 32, 2);
//		 fft_tr_layer(struct_stream_3_3_0, struct_stream_3_3_1, struct_stream_3_2_0, struct_stream_3_2_1, 16, 4);
//		 fft_tr_layer(struct_stream_3_4_0, struct_stream_3_4_1, struct_stream_3_3_0, struct_stream_3_3_1, 8, 8);
//		 fft_tr_layer(struct_stream_3_5_0, struct_stream_3_5_1, struct_stream_3_4_0, struct_stream_3_4_1, 4, 16);
//		 fft_tr_layer(struct_stream_3_6_0, struct_stream_3_6_1, struct_stream_3_5_0, struct_stream_3_5_1, 2, 32);
//		 merge_tr_struct_stream(internal_stream_1_0, struct_stream_3_6_0, struct_stream_3_6_1);
//		 transpose_broadcast(int_stream_0_1, int_stream_1_1, int_stream_2_1, int_stream_3_1, internal_stream_1_0);
		/////////////////////////////////END BUTTERFLIES_TR///////////////////////////////////

		radix_conversions_tr(s_priv_cmp_stream_0, s_priv_cmp_stream_1, s_priv_cmp_stream_2, s_priv_cmp_stream_3, int_stream_0_1, int_stream_1_1, int_stream_2_1, int_stream_3_1);
		//////////////////////////////RADIX_CONVERSIONS_TR////////////////////////////////////
		// mask_streaming(mask_stream_5_0_1, mask_stream_5_1_1, mask_stream_4_0_1, mask_stream_4_1_1, mask_stream_3_0_1, mask_stream_3_1_1, mask_stream_2_0_1, mask_stream_2_1_1, mask_stream_1_0_1, mask_stream_1_1_1, mask_stream_0_0_1,  mask_stream_0_1_1, mask_1);
		// s_streaming(s_stream_4_0_1, s_stream_4_1_1, s_stream_3_0_1, s_stream_3_1_1, s_stream_2_0_1, s_stream_2_1_1, s_stream_1_0_1, s_stream_1_1_1, s_stream_0_0_1, s_stream_0_1_1, s_radix_1);

		// loop_function_5(out_stream_5_0_1, out_stream_5_1_1, int_stream_0_1, int_stream_1_1, mask_stream_5_0_1, mask_stream_5_1_1);
		// loop_function_4(out_stream_4_0_1, out_stream_4_1_1, out_stream_5_0_1, out_stream_5_1_1, mask_stream_4_0_1, mask_stream_4_1_1, s_stream_4_0_1, s_stream_4_1_1);
		// loop_function_3(out_stream_3_0_1, out_stream_3_1_1, out_stream_4_0_1, out_stream_4_1_1, mask_stream_3_0_1, mask_stream_3_1_1, s_stream_3_0_1, s_stream_3_1_1);
		// loop_function_2(out_stream_2_0_1, out_stream_2_1_1, out_stream_3_0_1, out_stream_3_1_1, mask_stream_2_0_1, mask_stream_2_1_1, s_stream_2_0_1, s_stream_2_1_1);
		// loop_function_1(out_stream_1_0_1, out_stream_1_1_1, out_stream_2_0_1, out_stream_2_1_1, mask_stream_1_0_1, mask_stream_1_1_1, s_stream_1_0_1, s_stream_1_1_1);
		// loop_function_0(s_priv_cmp_stream_0, s_priv_cmp_stream_1, out_stream_1_0_1, out_stream_1_1_1, mask_stream_0_0_1, mask_stream_0_1_1, s_stream_0_0_1, s_stream_0_1_1);
		////////////////////////////END RADIX_CONVERSIONS_TR//////////////////////////////////

	////////////////////////////END FFT_TR//////////////////////////////////
	//OK
	synd_cmp(check_synd_stream, s_priv_stream_1_0, s_priv_stream_1_1, s_priv_stream_1_2, s_priv_stream_1_3, s_priv_cmp_stream_0, s_priv_cmp_stream_1, s_priv_cmp_stream_2, s_priv_cmp_stream_3);
	//OK
	benes(error_stream_benes_out, error_stream_1, sk_vec_fwd_mat_stream);
	//OK
	error_stream_duplicate(error_stream_benes_out_0, error_stream_benes_out_1, error_stream_benes_out);
	//OK
	postprocess(e, e_stream_int, error_stream_benes_out_0);
	//OK
	weight_check(check_weight_stream, e_stream_int, error_stream_benes_out_1);
	//OK
	check_value_compute(check, check_synd_stream, check_weight_stream);



}
