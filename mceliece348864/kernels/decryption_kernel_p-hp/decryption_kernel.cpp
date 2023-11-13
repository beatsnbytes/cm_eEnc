#ifndef __SYNTHESIS__
#include <gmp.h>
#define __gmp_const const
#endif
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#define AP_INT_MAX_W 2048
#include"ap_int.h"
#include"hls_stream.h"
#include"params.h"


typedef ap_uint<12> gf;
typedef ap_uint<64> vec;
typedef ap_uint<12> gf_vec;
typedef ap_uint<2048> sk_2k_vec;
typedef ap_uint<1024> sk_1k_vec;
typedef ap_uint<1> bit;

typedef struct{
	vec mat[GFBITS];
}mat_struct;

typedef struct{
	vec mat[64];
}mat_64_struct;




//TODO maybe instead of 12 iterations play with mat[0] pointers
void stream_out(hls::stream<mat_struct> &out_stream, vec out[64][GFBITS]){

	int k, s, l;
	int stream_idx=0;

	mat_struct tmp_struct;

	for (k = 0; k < 32; k++)
	{
		for(l=0; l<12; l++){
			tmp_struct.mat[l] = out[k][l];
		}

		out_stream.write(tmp_struct);
		for(l=0; l<12; l++){
			tmp_struct.mat[l] = out[k+32][l];
		}
		out_stream.write(tmp_struct);
	}

}



//TODO read in 256B=2048b from the CPU
//Change accordingly the below
//let for later
void load_sk_vec(hls::stream<sk_2k_vec> &sk_wide_mat_stream, sk_2k_vec *sk_values){

	int i;

	//Reading in 256B=2048b
	LOOP_SK:
	for(i=0; i<23; i++){
		sk_wide_mat_stream.write(*(sk_values +i));
	}
}

/* bitsliced field squarings */
void vec_sq(vec *out, vec *in)
{
        int i;
        vec result[GFBITS];

        //

        result[0] = in[0]^in[6];
        result[1] = in[11];
        result[2] = in[1]^in[7];
        result[3] = in[6];
        result[4] = in[2] ^ in[11] ^ in[8];
        result[5] = in[7];
        result[6] = in[3]^in[9];
        result[7] = in[8];
        result[8] = in[4]^in[10];
        result[9] = in[9];
        result[10] = in[5] ^ in[11];
        result[11] = in[10];

        //

        for (i = 0; i < GFBITS; i++){
                out[i] = result[i];
        }
}


gf gf_mul(gf in0, gf in1)
{
	int i;

	uint32_t tmp;
	uint32_t t0;
	uint32_t t1;
	uint32_t t;

	t0 = in0;
	t1 = in1;

	tmp = t0 * (t1 & 1);

	//TODO maybe use vectors and bit retrieval accordingly?
	for (i = 1; i < GFBITS; i++)
		tmp ^= (t0 * (t1 & (1 << i)));

	t = tmp & 0x7FC000;
	tmp ^= t >> 9;
	tmp ^= t >> 12;

	t = tmp & 0x3000;
	tmp ^= t >> 9;
	tmp ^= t >> 12;

	return tmp & ((1 << GFBITS)-1);
}

/* input: field element in */
/* return: in^2 */
static inline gf gf_sq(gf in)
{
	const uint32_t B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};

	uint32_t x = in;
	uint32_t t;

	x = (x | (x << 8)) & B[3];
	x = (x | (x << 4)) & B[2];
	x = (x | (x << 2)) & B[1];
	x = (x | (x << 1)) & B[0];

	t = x & 0x7FC000;
	x ^= t >> 9;
	x ^= t >> 12;

	t = x & 0x3000;
	x ^= t >> 9;
	x ^= t >> 12;

	return x & ((1 << GFBITS)-1);
}



gf gf_inv(gf in)
{
	gf tmp_11;
	gf tmp_1111;

	gf out = in;

	out = gf_sq(out);
	tmp_11 = gf_mul(out, in); // 11

	out = gf_sq(tmp_11);
	out = gf_sq(out);
	tmp_1111 = gf_mul(out, tmp_11); // 1111

	out = gf_sq(tmp_1111);
	out = gf_sq(out);
	out = gf_sq(out);
	out = gf_sq(out);
	out = gf_mul(out, tmp_1111); // 11111111

	out = gf_sq(out);
	out = gf_sq(out);
	out = gf_mul(out, tmp_11); // 1111111111

	out = gf_sq(out);
	out = gf_mul(out, in); // 11111111111

	return gf_sq(out); // 111111111110
}

/* 2 field multiplications */

//TODO maybe perform two consecutively?
static inline vec gf_mul2(gf a, gf b0, gf b1)
{
	int i;

	vec tmp=0;
	vec t0;
	vec t1;
	vec t;
	vec mask = 0x0000000100000001;

	t0 = a;

	t1.range(31, 0) = b0;
	t1.range(63, 32) = b1;

	for (i = 0; i < GFBITS; i++)
	{
		tmp ^= t0 * (t1 & mask);
		mask = mask << 1;
	}

	//

	t = tmp & 0x007FC000007FC000;
	tmp ^= (t >> 9) ^ (t >> 12);

	t = tmp & 0x0000300000003000;
	tmp ^= (t >> 9) ^ (t >> 12);

	return tmp & 0x00000FFF00000FFF;
}


static inline vec vec_set1_16b(uint16_t v)
{
	vec ret;

	ret = v;
	ret |= ret << 16;
	ret |= ret << 32;

	return ret;
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
	vec int_mat[64];


	LOOP_TRANSPOSE:
	for(j=0; j<64; j++){
		#pragma HLS PIPELINE II=1
		for(i=0; i<64; i++){
			out_mat[j].bit(i) = in_mat[i].bit(j);
		}
	}
}


void transpose_64x64_stream_out(hls::stream<vec> &out_mat_stream, vec *in_mat){

	int i, j;
	vec out_mat;


	LOOP_TRANSPOSE:
	for(j=0; j<64; j++){
		#pragma HLS PIPELINE II=1
		for(i=0; i<64; i++){
			out_mat.bit(i) = in_mat[i].bit(j);
		}
		out_mat_stream.write(out_mat);
	}
}




static inline vec vec_setbits(vec b)
{
	vec ret = -b;

	return ret;
}

static inline vec vec_or_reduce(hls::stream<vec> &a)
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
static inline int vec_testz(vec a)
{
	a |= a >> 32;
	a |= a >> 16;
	a |= a >> 8;
	a |= a >> 4;
	a |= a >> 2;
	a |= a >> 1;

	return (a&1)^1;
}

static inline void vec_copy(vec * out, vec * in)
{
	int i;

	for (i = 0; i < GFBITS; i++)
		out[i] = in[i];
}


void vec_mul_bm_out(hls::stream<vec> &out_stream, const vec * f, const vec * g)
{
	int i, j;
	vec buf[ 2*GFBITS-1 ];
	#pragma HLS array partition variable=buf cyclic complete


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
		buf[i-GFBITS+3] ^= buf[i];
		buf[i-GFBITS+0] ^= buf[i];
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS pipeline
		out_stream.write(buf[i]);
	}
}

void vec_mul(vec * h, const vec * f, const vec * g)
{
	int i, j;
	vec buf[ 2*GFBITS-1 ];
	#pragma HLS array partition variable=buf cyclic complete


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
		buf[i-GFBITS+3] ^= buf[i];
		buf[i-GFBITS+0] ^= buf[i];
	}

	for (i = 0; i < GFBITS; i++){
	#pragma HLS unroll
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

static inline void radix_conversions_tr(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1)
{
	int i, j, k;

	vec in[2][GFBITS];

	const vec mask[6][2] =
	{
		{0x2222222222222222, 0x4444444444444444},
		{0x0C0C0C0C0C0C0C0C, 0x3030303030303030},
		{0x00F000F000F000F0, 0x0F000F000F000F00},
		{0x0000FF000000FF00, 0x00FF000000FF0000},
		{0x00000000FFFF0000, 0x0000FFFF00000000},
		{0xFFFFFFFF00000000, 0x00000000FFFFFFFF}
	};

	const vec s[5][2][GFBITS] =
	{
	{{
		0XF3CFC030FC30F003,
		0X3FCF0F003C00C00C,
		0X30033CC300C0C03C,
		0XCCFF0F3C0F30F0C0,
		0X0300C03FF303C3F0,
		0X3FFF3C0FF0CCCCC0,
		0XF3FFF0C00F3C3CC0,
		0X3003333FFFC3C000,
		0X0FF30FFFC3FFF300,
		0XFFC0F300F0F0CC00,
		0XC0CFF3FCCC3CFC00,
		0XFC3C03F0F330C000,
	},
	{
		0X000C03C0C3C0330C,
		0XF330CFFCC00F33C0,
		0XCCF330F00F3C0333,
		0XFF03FFF3FF0CF0C0,
		0X3CC3FCF00FCC303C,
		0X0F000C0FC30303F3,
		0XCF0FC3FF333CCF3C,
		0X003F3FC3C0FF333F,
		0X3CC3F0F3CF0FF00F,
		0XF3F33CC03FC30CC0,
		0X3CC330CFC333F33F,
		0X3CC0303FF3C3FFFC,
	}},
	{{
		0X000F00000000F00F,
		0X00000F00F00000F0,
		0X0F00000F00000F00,
		0XF00F00F00F000000,
		0X00F00000000000F0,
		0X0000000F00000000,
		0XF00000000F00F000,
		0X00F00F00000F0000,
		0X0000F00000F00F00,
		0X000F00F00F00F000,
		0X00F00F0000000000,
		0X0000000000F00000,
	},
	{
		0X0F00F00F00000000,
		0XF00000000000F000,
		0X00000F00000000F0,
		0X0F00F00000F00000,
		0X000F00000F00F00F,
		0X00F00F00F00F0000,
		0X0F00F00000000000,
		0X000000000F000000,
		0X00F00000000F00F0,
		0X0000F00F00000F00,
		0XF00000F00000F00F,
		0X00000F00F00F00F0,
	}},
	{{
		0X0000FF00FF0000FF,
		0X0000FF000000FF00,
		0XFF0000FF00FF0000,
		0XFFFF0000FF000000,
		0X00FF00FF00FF0000,
		0X0000FFFFFF000000,
		0X00FFFF00FF000000,
		0XFFFFFF0000FF0000,
		0XFFFF00FFFF00FF00,
		0X0000FF0000000000,
		0XFFFFFF00FF000000,
		0X00FF000000000000,
	},
	{
		0XFF00FFFFFF000000,
		0XFF0000FFFF000000,
		0XFFFF00FFFF000000,
		0XFF00FFFFFFFFFF00,
		0X00000000FF00FF00,
		0XFFFFFFFF00FF0000,
		0X00FFFFFF00FF0000,
		0XFFFF00FFFF00FFFF,
		0XFFFF0000FFFFFFFF,
		0XFF00000000FF0000,
		0X000000FF00FF00FF,
		0X00FF00FF00FFFF00,
	}},
	{{
		0X000000000000FFFF,
		0X00000000FFFF0000,
		0X0000000000000000,
		0XFFFF000000000000,
		0X00000000FFFF0000,
		0X0000FFFF00000000,
		0X0000000000000000,
		0X00000000FFFF0000,
		0X0000FFFF00000000,
		0X0000000000000000,
		0X0000000000000000,
		0X0000000000000000,
	},
	{
		0X0000000000000000,
		0XFFFF000000000000,
		0X0000000000000000,
		0X0000000000000000,
		0XFFFF00000000FFFF,
		0X0000000000000000,
		0X0000FFFF00000000,
		0XFFFF00000000FFFF,
		0X00000000FFFF0000,
		0X0000000000000000,
		0XFFFF00000000FFFF,
		0X00000000FFFF0000,
	}},
	{{
		0X00000000FFFFFFFF,
		0XFFFFFFFF00000000,
		0XFFFFFFFF00000000,
		0X0000000000000000,
		0X0000000000000000,
		0XFFFFFFFF00000000,
		0X0000000000000000,
		0X0000000000000000,
		0XFFFFFFFF00000000,
		0X0000000000000000,
		0X0000000000000000,
		0X0000000000000000,
	},
	{
		0X0000000000000000,
		0X0000000000000000,
		0X00000000FFFFFFFF,
		0XFFFFFFFF00000000,
		0XFFFFFFFF00000000,
		0X0000000000000000,
		0XFFFFFFFF00000000,
		0XFFFFFFFFFFFFFFFF,
		0XFFFFFFFF00000000,
		0X0000000000000000,
		0XFFFFFFFFFFFFFFFF,
		0XFFFFFFFF00000000,
	}}

	};

	//
	for(i=0; i<GFBITS; i++){
		in[0][i] = in_stream_0.read();
		in[1][i] = in_stream_1.read();
	}



	for (j = 5; j >= 0; j--)
	{
		if (j < 5)
		{
			vec_mul(in[0], in[0], s[j][0]); // scaling
			vec_mul(in[1], in[1], s[j][1]); // scaling
		}

		for (i = 0; i < GFBITS; i++){
			for (k = j; k <= 4; k++)
			#pragma HLS loop_tripcount min=1 max=5
			{
				in[0][i] ^= (in[0][i] & mask[k][0]) << (1 << k);
				in[0][i] ^= (in[0][i] & mask[k][1]) << (1 << k);

				in[1][i] ^= (in[1][i] & mask[k][0]) << (1 << k);
				in[1][i] ^= (in[1][i] & mask[k][1]) << (1 << k);
			}
		}

		for (i = 0; i < GFBITS; i++)
		{
			in[1][i] ^= (in[0][i] & mask[5][0]) >> 32;
			in[1][i] ^= (in[1][i] & mask[5][1]) << 32;
		}
	}

	for (i = 0; i < GFBITS; i++)
	{
		out_stream_0.write(in[0][i]);
		out_stream_1.write(in[1][i]);
//		printf("in[0][%d] = %lX\n", i, in[0][i].to_uint64());
//		printf("in[1][%d] = %lX\n", i, in[1][i].to_uint64());
	}

}

static inline void butterflies_tr(vec out[][ GFBITS ], hls::stream<mat_struct> &in_struct_stream)
{
	int i, j, k, s, b, l;

	vec tmp[ GFBITS ];
	vec pre[6][ GFBITS ];
	vec buf[64];
	vec buf_tr[64];

	#pragma HLS array partition variable=pre dim=2
	#pragma HLS array partition variable=tmp
	#pragma HLS array partition variable=buf_tr
	#pragma HLS array partition variable=buf

	const vec consts[ 63 ][ GFBITS ] =
	{
#include "consts.data"
	};

	vec consts_ptr = 63;

	const unsigned char reversal[64] =
	{
	  0, 32, 16, 48,  8, 40, 24, 56,
	  4, 36, 20, 52, 12, 44, 28, 60,
	  2, 34, 18, 50, 10, 42, 26, 58,
	  6, 38, 22, 54, 14, 46, 30, 62,
	  1, 33, 17, 49,  9, 41, 25, 57,
	  5, 37, 21, 53, 13, 45, 29, 61,
	  3, 35, 19, 51, 11, 43, 27, 59,
	  7, 39, 23, 55, 15, 47, 31, 63
	};

	#pragma HLS array partition variable=reversal

	const uint16_t beta[6] = {8, 1300, 3408, 1354, 2341, 1154};

	mat_struct in_struct;
	vec in_1[GFBITS];
	vec in_2[GFBITS];

	vec in[64][GFBITS];

	// butterflies

	for (i = 5; i >= 0; i--)
	{
		s = 1 << i;
		consts_ptr -= s;

		for (j = 0; j < 64; j += 2*s){
		#pragma HLS loop_tripcount min=1 max=32
			for (k = j; k < j+s; k++)
			#pragma HLS loop_tripcount min=1 max=32
			{
				if((i==5) && (j==0)){
					in_struct = in_struct_stream.read();
					for(l=0; l<12; l++){
						in[k][l] = in_struct.mat[l];
					}

					in_struct = in_struct_stream.read();
					for(l=0; l<12; l++){
						in[k+s][l] = in_struct.mat[l];
					}
				}

				vec_add(in[k], in[k], in[k+s]);
				vec_mul(tmp, in[k], consts[ consts_ptr + (k-j) ]);
				vec_add(in[k+s], in[k+s], tmp);
			}
		}
	}

	// transpose

	for (i = 0; i < GFBITS; i++)
	#pragma HLS pipeline
	{
		for (j = 0; j < 64; j++){
			buf[ reversal[j] ] = in[j][i];
		}

		transpose_64x64(buf_tr, buf);

		for (j = 0; j < 64; j++)
			in[j][i] = buf_tr[ j ];
	}

	// boradcast

	vec_copy(pre[0], in[32]); vec_add(in[33], in[33], in[32]);
	vec_copy(pre[1], in[33]); vec_add(in[35], in[35], in[33]);
	vec_add(pre[0], pre[0], in[35]); vec_add(in[34], in[34], in[35]);
	vec_copy(pre[2], in[34]); vec_add(in[38], in[38], in[34]);
	vec_add(pre[0], pre[0], in[38]); vec_add(in[39], in[39], in[38]);
	vec_add(pre[1], pre[1], in[39]); vec_add(in[37], in[37], in[39]);
	vec_add(pre[0], pre[0], in[37]); vec_add(in[36], in[36], in[37]);
	vec_copy(pre[3], in[36]); vec_add(in[44], in[44], in[36]);
	vec_add(pre[0], pre[0], in[44]); vec_add(in[45], in[45], in[44]);
	vec_add(pre[1], pre[1], in[45]); vec_add(in[47], in[47], in[45]);
	vec_add(pre[0], pre[0], in[47]); vec_add(in[46], in[46], in[47]);
	vec_add(pre[2], pre[2], in[46]); vec_add(in[42], in[42], in[46]);
	vec_add(pre[0], pre[0], in[42]); vec_add(in[43], in[43], in[42]);
	vec_add(pre[1], pre[1], in[43]); vec_add(in[41], in[41], in[43]);
	vec_add(pre[0], pre[0], in[41]); vec_add(in[40], in[40], in[41]);
	vec_copy(pre[4], in[40]); vec_add(in[56], in[56], in[40]);
	vec_add(pre[0], pre[0], in[56]); vec_add(in[57], in[57], in[56]);
	vec_add(pre[1], pre[1], in[57]); vec_add(in[59], in[59], in[57]);
	vec_add(pre[0], pre[0], in[59]); vec_add(in[58], in[58], in[59]);
	vec_add(pre[2], pre[2], in[58]); vec_add(in[62], in[62], in[58]);
	vec_add(pre[0], pre[0], in[62]); vec_add(in[63], in[63], in[62]);
	vec_add(pre[1], pre[1], in[63]); vec_add(in[61], in[61], in[63]);
	vec_add(pre[0], pre[0], in[61]); vec_add(in[60], in[60], in[61]);
	vec_add(pre[3], pre[3], in[60]); vec_add(in[52], in[52], in[60]);
	vec_add(pre[0], pre[0], in[52]); vec_add(in[53], in[53], in[52]);
	vec_add(pre[1], pre[1], in[53]); vec_add(in[55], in[55], in[53]);
	vec_add(pre[0], pre[0], in[55]); vec_add(in[54], in[54], in[55]);
	vec_add(pre[2], pre[2], in[54]); vec_add(in[50], in[50], in[54]);
	vec_add(pre[0], pre[0], in[50]); vec_add(in[51], in[51], in[50]);
	vec_add(pre[1], pre[1], in[51]); vec_add(in[49], in[49], in[51]);
	vec_add(pre[0], pre[0], in[49]); vec_add(in[48], in[48], in[49]);
	vec_copy(pre[5], in[48]); vec_add(in[16], in[16], in[48]);
	vec_add(pre[0], pre[0], in[16]); vec_add(in[17], in[17], in[16]);
	vec_add(pre[1], pre[1], in[17]); vec_add(in[19], in[19], in[17]);
	vec_add(pre[0], pre[0], in[19]); vec_add(in[18], in[18], in[19]);
	vec_add(pre[2], pre[2], in[18]); vec_add(in[22], in[22], in[18]);
	vec_add(pre[0], pre[0], in[22]); vec_add(in[23], in[23], in[22]);
	vec_add(pre[1], pre[1], in[23]); vec_add(in[21], in[21], in[23]);
	vec_add(pre[0], pre[0], in[21]); vec_add(in[20], in[20], in[21]);
	vec_add(pre[3], pre[3], in[20]); vec_add(in[28], in[28], in[20]);
	vec_add(pre[0], pre[0], in[28]); vec_add(in[29], in[29], in[28]);
	vec_add(pre[1], pre[1], in[29]); vec_add(in[31], in[31], in[29]);
	vec_add(pre[0], pre[0], in[31]); vec_add(in[30], in[30], in[31]);
	vec_add(pre[2], pre[2], in[30]); vec_add(in[26], in[26], in[30]);
	vec_add(pre[0], pre[0], in[26]); vec_add(in[27], in[27], in[26]);
	vec_add(pre[1], pre[1], in[27]); vec_add(in[25], in[25], in[27]);
	vec_add(pre[0], pre[0], in[25]); vec_add(in[24], in[24], in[25]);
	vec_add(pre[4], pre[4], in[24]); vec_add(in[8], in[8], in[24]);
	vec_add(pre[0], pre[0], in[8]); vec_add(in[9], in[9], in[8]);
	vec_add(pre[1], pre[1], in[9]); vec_add(in[11], in[11], in[9]);
	vec_add(pre[0], pre[0], in[11]); vec_add(in[10], in[10], in[11]);
	vec_add(pre[2], pre[2], in[10]); vec_add(in[14], in[14], in[10]);
	vec_add(pre[0], pre[0], in[14]); vec_add(in[15], in[15], in[14]);
	vec_add(pre[1], pre[1], in[15]); vec_add(in[13], in[13], in[15]);
	vec_add(pre[0], pre[0], in[13]); vec_add(in[12], in[12], in[13]);
	vec_add(pre[3], pre[3], in[12]); vec_add(in[4], in[4], in[12]);
	vec_add(pre[0], pre[0], in[4]); vec_add(in[5], in[5], in[4]);
	vec_add(pre[1], pre[1], in[5]); vec_add(in[7], in[7], in[5]);
	vec_add(pre[0], pre[0], in[7]); vec_add(in[6], in[6], in[7]);
	vec_add(pre[2], pre[2], in[6]); vec_add(in[2], in[2], in[6]);
	vec_add(pre[0], pre[0], in[2]); vec_add(in[3], in[3], in[2]);
	vec_add(pre[1], pre[1], in[3]); vec_add(in[1], in[1], in[3]);

	vec_add(pre[0], pre[0], in[1]); vec_add(out[0], in[0], in[1]);

	//

	for (j = 0; j < GFBITS; j++) {
		tmp[j] = (beta[0] >> j) & 1;
		tmp[j] = -tmp[j];
	}

	vec_mul(out[1], pre[0], tmp);

	for (i = 1; i < 6; i++)
	#pragma HLS pipeline
	{
		for (j = 0; j < GFBITS; j++) {
			tmp[j] = (beta[i] >> j) & 1;
			tmp[j] = -tmp[j];
		}

		vec_mul(tmp, pre[i], tmp);
		vec_add(out[1], out[1], tmp);
	}

}

//void fft_tr(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<mat_struct> &in_struct_stream)
//{
//		vec int_mat[2][GFBITS];
//
//		//TODO try to apply dataflow here
//		//Have to write in specific order since the radix reads the in like that
//		butterflies_tr(int_mat, in_struct_stream);
//		radix_conversions_tr_alt(out_stream_0, out_stream_1, int_mat);
//}



void split_struct_stream(hls::stream<mat_struct> &out_stream_0, hls::stream<mat_struct> &out_stream_1, hls::stream<mat_struct> &in_struct_stream){

	int i;
	mat_struct tmp;

	for(i=0; i<64; i++){
			tmp = in_struct_stream.read();
			if(i<32){
				out_stream_0.write(tmp);
			}else{
				out_stream_1.write(tmp);
			}
	}
}

void butterfly_tr_operation(mat_struct *op_a, mat_struct *op_b, vec consts[63][GFBITS], int index){

	vec tmp[GFBITS];
	//TODO apply dataflow
	//Maybe pass as streams? This will lower the FF too..current is 67 cycles per layer
	vec_add(op_a->mat, op_a->mat, op_b->mat);
	vec_mul(tmp, op_a->mat, consts[index]);
	vec_add(op_b->mat, op_b->mat, tmp);


}






static void fft_tr_layer(hls::stream<mat_struct> &sink_0, hls::stream<mat_struct> &sink_1, hls::stream<mat_struct> &src_0, hls::stream<mat_struct> &src_1, int start, int reps)
{

	vec consts[ 63 ][ GFBITS ] =
	{
#include "consts.data"
	};

	mat_struct wing1;
	mat_struct wing2;

	int iter=0;

	unsigned int index = start;

	for (unsigned int i = 0; i < 32; i++) {
		#pragma HLS pipeline
		wing1 = src_0.read();
		wing2 = src_1.read();

		//TODO apply dataflow inside the butterfly?
		butterfly_tr_operation(&wing1, &wing2, consts, index);
//		printf("idx=%d\n", index);

		if (i % 32 < 16) {
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


void merge_tr_struct_stream(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<mat_struct> &in_stream_0, hls::stream<mat_struct> &in_stream_1){

	int i, j;
	mat_struct tmp_0, tmp_1;

	int idx_0, idx_1;

	vec out_mat[64][GFBITS];

	#pragma HLS array_partition variable=out_mat dim=2 complete

	for(i=0; i<32; i++){
			tmp_0 = in_stream_0.read();
			tmp_1 = in_stream_1.read();
			for(j=0; j<GFBITS; j++){
				out_mat[i][j] = tmp_0.mat[j];
				out_mat[i+32][j] = tmp_1.mat[j];
			}
	}

	//TODO make more efficient
	for(i=0; i<GFBITS; i++){
	#pragma HLS pipeline
		for(j=0; j<32; j++){
			out_stream_0.write(out_mat[j][i]);
			out_stream_1.write(out_mat[j+32][i]);
		}
	}

}




void transpose_broadcast(hls::stream<vec> &stream_out_0, hls::stream<vec> &stream_out_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1){
	// transpose
	int i, j, k, s, b, l;

	vec tmp[ GFBITS ];
	vec pre[6][ GFBITS ];
	vec buf[64];
	vec buf_tr[64];
	vec in[64][GFBITS];
	vec out[2][GFBITS];

	const uint16_t beta[6] = {8, 1300, 3408, 1354, 2341, 1154};

	const unsigned char reversal[64] =
	{
	  0, 32, 16, 48,  8, 40, 24, 56,
	  4, 36, 20, 52, 12, 44, 28, 60,
	  2, 34, 18, 50, 10, 42, 26, 58,
	  6, 38, 22, 54, 14, 46, 30, 62,
	  1, 33, 17, 49,  9, 41, 25, 57,
	  5, 37, 21, 53, 13, 45, 29, 61,
	  3, 35, 19, 51, 11, 43, 27, 59,
	  7, 39, 23, 55, 15, 47, 31, 63
	};

	#pragma HLS array_partition variable=reversal
	#pragma HLS array_partition variable=pre complete dim=2
	#pragma HLS array_partition variable=tmp
	#pragma HLS array_partition variable=buf_tr
	#pragma HLS array_partition variable=buf
	#pragma HLS array_partition variable=in complete dim=2
	#pragma HLS array_partition variable=out

	for (i = 0; i < GFBITS; i++)
	#pragma HLS PIPELINE
	{
		for (j = 0; j < 32; j++){
			in[j][i] = in_stream_0.read();
			buf[ reversal[j] ] = in[j][i];
		}

		for (j = 32; j < 64; j++){
			in[j][i] = in_stream_1.read();
			buf[ reversal[j] ] = in[j][i];
		}

		transpose_64x64(buf_tr, buf);

		for (j = 0; j < 64; j++){
			in[j][i] = buf_tr[ j ];
		}
	}

	// broadcast

	vec_copy(pre[0], in[32]); vec_add(in[33], in[33], in[32]);
	vec_copy(pre[1], in[33]); vec_add(in[35], in[35], in[33]);
	vec_add(pre[0], pre[0], in[35]); vec_add(in[34], in[34], in[35]);
	vec_copy(pre[2], in[34]); vec_add(in[38], in[38], in[34]);
	vec_add(pre[0], pre[0], in[38]); vec_add(in[39], in[39], in[38]);
	vec_add(pre[1], pre[1], in[39]); vec_add(in[37], in[37], in[39]);
	vec_add(pre[0], pre[0], in[37]); vec_add(in[36], in[36], in[37]);
	vec_copy(pre[3], in[36]); vec_add(in[44], in[44], in[36]);
	vec_add(pre[0], pre[0], in[44]); vec_add(in[45], in[45], in[44]);
	vec_add(pre[1], pre[1], in[45]); vec_add(in[47], in[47], in[45]);
	vec_add(pre[0], pre[0], in[47]); vec_add(in[46], in[46], in[47]);
	vec_add(pre[2], pre[2], in[46]); vec_add(in[42], in[42], in[46]);
	vec_add(pre[0], pre[0], in[42]); vec_add(in[43], in[43], in[42]);
	vec_add(pre[1], pre[1], in[43]); vec_add(in[41], in[41], in[43]);
	vec_add(pre[0], pre[0], in[41]); vec_add(in[40], in[40], in[41]);
	vec_copy(pre[4], in[40]); vec_add(in[56], in[56], in[40]);
	vec_add(pre[0], pre[0], in[56]); vec_add(in[57], in[57], in[56]);
	vec_add(pre[1], pre[1], in[57]); vec_add(in[59], in[59], in[57]);
	vec_add(pre[0], pre[0], in[59]); vec_add(in[58], in[58], in[59]);
	vec_add(pre[2], pre[2], in[58]); vec_add(in[62], in[62], in[58]);
	vec_add(pre[0], pre[0], in[62]); vec_add(in[63], in[63], in[62]);
	vec_add(pre[1], pre[1], in[63]); vec_add(in[61], in[61], in[63]);
	vec_add(pre[0], pre[0], in[61]); vec_add(in[60], in[60], in[61]);
	vec_add(pre[3], pre[3], in[60]); vec_add(in[52], in[52], in[60]);
	vec_add(pre[0], pre[0], in[52]); vec_add(in[53], in[53], in[52]);
	vec_add(pre[1], pre[1], in[53]); vec_add(in[55], in[55], in[53]);
	vec_add(pre[0], pre[0], in[55]); vec_add(in[54], in[54], in[55]);
	vec_add(pre[2], pre[2], in[54]); vec_add(in[50], in[50], in[54]);
	vec_add(pre[0], pre[0], in[50]); vec_add(in[51], in[51], in[50]);
	vec_add(pre[1], pre[1], in[51]); vec_add(in[49], in[49], in[51]);
	vec_add(pre[0], pre[0], in[49]); vec_add(in[48], in[48], in[49]);
	vec_copy(pre[5], in[48]); vec_add(in[16], in[16], in[48]);
	vec_add(pre[0], pre[0], in[16]); vec_add(in[17], in[17], in[16]);
	vec_add(pre[1], pre[1], in[17]); vec_add(in[19], in[19], in[17]);
	vec_add(pre[0], pre[0], in[19]); vec_add(in[18], in[18], in[19]);
	vec_add(pre[2], pre[2], in[18]); vec_add(in[22], in[22], in[18]);
	vec_add(pre[0], pre[0], in[22]); vec_add(in[23], in[23], in[22]);
	vec_add(pre[1], pre[1], in[23]); vec_add(in[21], in[21], in[23]);
	vec_add(pre[0], pre[0], in[21]); vec_add(in[20], in[20], in[21]);
	vec_add(pre[3], pre[3], in[20]); vec_add(in[28], in[28], in[20]);
	vec_add(pre[0], pre[0], in[28]); vec_add(in[29], in[29], in[28]);
	vec_add(pre[1], pre[1], in[29]); vec_add(in[31], in[31], in[29]);
	vec_add(pre[0], pre[0], in[31]); vec_add(in[30], in[30], in[31]);
	vec_add(pre[2], pre[2], in[30]); vec_add(in[26], in[26], in[30]);
	vec_add(pre[0], pre[0], in[26]); vec_add(in[27], in[27], in[26]);
	vec_add(pre[1], pre[1], in[27]); vec_add(in[25], in[25], in[27]);
	vec_add(pre[0], pre[0], in[25]); vec_add(in[24], in[24], in[25]);
	vec_add(pre[4], pre[4], in[24]); vec_add(in[8], in[8], in[24]);
	vec_add(pre[0], pre[0], in[8]); vec_add(in[9], in[9], in[8]);
	vec_add(pre[1], pre[1], in[9]); vec_add(in[11], in[11], in[9]);
	vec_add(pre[0], pre[0], in[11]); vec_add(in[10], in[10], in[11]);
	vec_add(pre[2], pre[2], in[10]); vec_add(in[14], in[14], in[10]);
	vec_add(pre[0], pre[0], in[14]); vec_add(in[15], in[15], in[14]);
	vec_add(pre[1], pre[1], in[15]); vec_add(in[13], in[13], in[15]);
	vec_add(pre[0], pre[0], in[13]); vec_add(in[12], in[12], in[13]);
	vec_add(pre[3], pre[3], in[12]); vec_add(in[4], in[4], in[12]);
	vec_add(pre[0], pre[0], in[4]); vec_add(in[5], in[5], in[4]);
	vec_add(pre[1], pre[1], in[5]); vec_add(in[7], in[7], in[5]);
	vec_add(pre[0], pre[0], in[7]); vec_add(in[6], in[6], in[7]);
	vec_add(pre[2], pre[2], in[6]); vec_add(in[2], in[2], in[6]);
	vec_add(pre[0], pre[0], in[2]); vec_add(in[3], in[3], in[2]);
	vec_add(pre[1], pre[1], in[3]); vec_add(in[1], in[1], in[3]);
	vec_add(pre[0], pre[0], in[1]);
	//Here I can stream the 0 out. Does the next need it?
	vec_add(out[0], in[0], in[1]);



	for (j = 0; j < GFBITS; j++) {
		tmp[j] = (beta[0] >> j) & 1;
		tmp[j] = -tmp[j];
	}

	vec_mul(out[1], pre[0], tmp);

	for (i = 1; i < 6; i++)
	#pragma HLS pipeline
	{
		for (j = 0; j < GFBITS; j++) {
			tmp[j] = (beta[i] >> j) & 1;
			tmp[j] = -tmp[j];
		}

		vec_mul(tmp, pre[i], tmp);
		vec_add(out[1], out[1], tmp);
	}

	for (i = 0; i < GFBITS; i++) {
		stream_out_0.write(out[0][i]);
		stream_out_1.write(out[1][i]);
	}

//	for (i = 0; i < GFBITS; i++) {
//		printf("out[0][%d] = %lX\n", i, out[0][i].to_uint64());
//		printf("out[1][%d] = %lX\n", i, out[1][i].to_uint64());
//	}


}



void fft_tr_stream1(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<mat_struct> &in_struct_stream)
{

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

	static hls::stream<vec> int_stream_0_1("int_stream_0_1");
	#pragma HLS STREAM variable=int_stream_0_1 depth=12

	static hls::stream<vec> int_stream_1_1("int_stream_1_1");
	#pragma HLS STREAM variable=int_stream_1_1 depth=12

	static hls::stream<vec> int_stream_0_0("int_stream_0_0");
	#pragma HLS STREAM variable=int_stream_0_0 depth=12

	static hls::stream<vec> int_stream_1_0("int_stream_1_0");
	#pragma HLS STREAM variable=int_stream_1_0 depth=12

	static hls::stream<vec> internal_stream_1_0("internal_stream_1_0");
	#pragma HLS STREAM variable=internal_stream_1_0 depth=64

	static hls::stream<vec> internal_stream_1_1("internal_stream_1_1");
	#pragma HLS STREAM variable=internal_stream_1_1 depth=64

	static hls::stream<vec> s_priv_cmp_stream_0("s_priv_cmp_stream_0");
	#pragma HLS STREAM variable=s_priv_cmp_stream_0 depth=12

	static hls::stream<vec> s_priv_cmp_stream_1("s_priv_cmp_stream_1");
	#pragma HLS STREAM variable=s_priv_cmp_stream_1 depth=12

	#pragma HLS dataflow

	split_struct_stream(struct_stream_3_0_0, struct_stream_3_0_1, in_struct_stream);
	fft_tr_layer(struct_stream_3_1_0, struct_stream_3_1_1, struct_stream_3_0_0, struct_stream_3_0_1, 31, 1);
	fft_tr_layer(struct_stream_3_2_0, struct_stream_3_2_1, struct_stream_3_1_0, struct_stream_3_1_1, 15, 2);
	fft_tr_layer(struct_stream_3_3_0, struct_stream_3_3_1, struct_stream_3_2_0, struct_stream_3_2_1, 7, 4);
	fft_tr_layer(struct_stream_3_4_0, struct_stream_3_4_1, struct_stream_3_3_0, struct_stream_3_3_1, 3, 8);
	fft_tr_layer(struct_stream_3_5_0, struct_stream_3_5_1, struct_stream_3_4_0, struct_stream_3_4_1, 1, 16);
	fft_tr_layer(struct_stream_3_6_0, struct_stream_3_6_1, struct_stream_3_5_0, struct_stream_3_5_1, 0, 32);
	merge_tr_struct_stream(internal_stream_1_0, internal_stream_1_1, struct_stream_3_6_0, struct_stream_3_6_1);
	transpose_broadcast(int_stream_0_1, int_stream_1_1, internal_stream_1_0, internal_stream_1_1);
	radix_conversions_tr(out_stream_0, out_stream_1, int_stream_0_1, int_stream_1_1);
}




void fft_tr_stream2(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<mat_struct> &in_struct_stream)
{

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

	static hls::stream<vec> int_stream_0_1("int_stream_0_1");
	#pragma HLS STREAM variable=int_stream_0_1 depth=12

	static hls::stream<vec> int_stream_1_1("int_stream_1_1");
	#pragma HLS STREAM variable=int_stream_1_1 depth=12

	static hls::stream<vec> int_stream_0_0("int_stream_0_0");
	#pragma HLS STREAM variable=int_stream_0_0 depth=12

	static hls::stream<vec> int_stream_1_0("int_stream_1_0");
	#pragma HLS STREAM variable=int_stream_1_0 depth=12

	static hls::stream<vec> internal_stream_1_0("internal_stream_1_0");
	#pragma HLS STREAM variable=internal_stream_1_0 depth=64

	static hls::stream<vec> internal_stream_1_1("internal_stream_1_1");
	#pragma HLS STREAM variable=internal_stream_1_1 depth=64

	static hls::stream<vec> s_priv_cmp_stream_0("s_priv_cmp_stream_0");
	#pragma HLS STREAM variable=s_priv_cmp_stream_0 depth=12

	static hls::stream<vec> s_priv_cmp_stream_1("s_priv_cmp_stream_1");
	#pragma HLS STREAM variable=s_priv_cmp_stream_1 depth=12

	#pragma HLS dataflow

	split_struct_stream(struct_stream_3_0_0, struct_stream_3_0_1, in_struct_stream);
	fft_tr_layer(struct_stream_3_1_0, struct_stream_3_1_1, struct_stream_3_0_0, struct_stream_3_0_1, 31, 1);
	fft_tr_layer(struct_stream_3_2_0, struct_stream_3_2_1, struct_stream_3_1_0, struct_stream_3_1_1, 15, 2);
	fft_tr_layer(struct_stream_3_3_0, struct_stream_3_3_1, struct_stream_3_2_0, struct_stream_3_2_1, 7, 4);
	fft_tr_layer(struct_stream_3_4_0, struct_stream_3_4_1, struct_stream_3_3_0, struct_stream_3_3_1, 3, 8);
	fft_tr_layer(struct_stream_3_5_0, struct_stream_3_5_1, struct_stream_3_4_0, struct_stream_3_4_1, 1, 16);
	fft_tr_layer(struct_stream_3_6_0, struct_stream_3_6_1, struct_stream_3_5_0, struct_stream_3_5_1, 0, 32);
	merge_tr_struct_stream(internal_stream_1_0, internal_stream_1_1, struct_stream_3_6_0, struct_stream_3_6_1);
	transpose_broadcast(int_stream_0_1, int_stream_1_1, internal_stream_1_0, internal_stream_1_1);
	radix_conversions_tr(out_stream_0, out_stream_1, int_stream_0_1, int_stream_1_1);
}



static inline vec mask_nonzero(gf a)
{
	vec ret = a;

	ret -= 1;
	ret >>= 63;
	ret -= 1;

	return ret;
}

static inline vec mask_leq(uint16_t a, uint16_t b)
{
	vec a_tmp = a;
	vec b_tmp = b;
	vec ret = b_tmp - a_tmp;

	ret >>= 63;
	ret -= 1;

	return ret;
}

static inline void vec_cmov(vec * out, vec * in, uint16_t mask)
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

static inline void interleave(vec *in, int idx0, int idx1, vec *mask, int b)
{
	int s = 1 << b;

	vec x, y;

	x = (in[idx0] & mask[0]) | ((in[idx1] & mask[0]) << s);
	y = ((in[idx0] & mask[1]) >> s) | (in[idx1] & mask[1]);

	in[idx0] = x;
	in[idx1] = y;
}

/* input: in, field elements in bitsliced form */
/* output: out, field elements in non-bitsliced form */
static inline void get_coefs(gf *out, hls::stream<vec> &in_stream)
{
	int i, k;

	vec mask[4][2];
	vec buf[16];
	#pragma HLS array partition variable=buf

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
	interleave(buf,  2,  3, mask[0], 0);
	interleave(buf,  4,  5, mask[0], 0);
	interleave(buf,  6,  7, mask[0], 0);
	interleave(buf,  8,  9, mask[0], 0);
	interleave(buf, 10, 11, mask[0], 0);
	interleave(buf, 12, 13, mask[0], 0);
	interleave(buf, 14, 15, mask[0], 0);

	for (i = 0; i < 16; i++){
	#pragma HLS pipeline
		for (k = 0; k <  4; k++){
			out[ k*16 + i ] = buf[i].range((k+1)*16-5, k*16);
		}
	}
}


static inline gf vec_reduce(vec *in)
{
	int i;
	vec tmp;
	gf ret = 0;

//Give the d_vec immediately a value here?

	for (i = GFBITS-1; i >= 0; i--)
	#pragma HLS pipeline
	{
		ret.bit(i) = in[i].xor_reduce();
	}

	return ret;
}

static void update(vec *in, const gf e)
{
	int i;

	for (i = 0; i < GFBITS; i++)
	#pragma HLS pipeline
	{
		in[i].range(62, 0) = in[i].range(63, 1);
		in[i].bit(63) = e.bit(i);
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
void bm(hls::stream<vec> &out_stream, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1)
{
	uint16_t i;
	uint16_t N, L;

	vec prod[ GFBITS ];
	vec in_tmp[ GFBITS ];

	vec d_vec[ GFBITS ];
	vec b_vec[ GFBITS ];
	vec B[ GFBITS ], C[ GFBITS ];
	vec B_tmp[ GFBITS ], C_tmp[ GFBITS ];

	vec tmp[GFBITS];

	#pragma HLS array partition variable=B_tmp
	#pragma HLS array partition variable=C_tmp
	#pragma HLS array partition variable=B
	#pragma HLS array partition variable=C
	#pragma HLS array partition variable=in_tmp
	#pragma HLS array partition variable=prod

	#pragma HLS array partition variable=d_vec
	#pragma HLS array partition variable=b_vec


	vec mask, t;

	gf d, b, c0=1;

	gf coefs[SYS_T * 2];
	#pragma HLS array partition variable=coefs


	// init

	get_coefs(&coefs[  0], in_stream_0);
	get_coefs(&coefs[ 64], in_stream_1);

	C[0] = 0;
	B[0] = 1;
	B[0] <<= 63;

	for (i = 1; i < GFBITS; i++)
		B[i] = C[i] = 0;

	b = 1;
	L = 0;

	//

	for (i = 0; i < GFBITS; i++)
		in_tmp[i] = 0;

	//TODO this is the most cotlier part here
	for (N = 0; N < SYS_T * 2; N++)
	#pragma HLS pipeline
	{
		// computing d

		vec_mul(prod, in_tmp, C);

		update(in_tmp, coefs[N]);

		d = vec_reduce(prod);


		t = gf_mul2(c0, coefs[N], b);
		d ^= t & 0xFFFFFFFF;

		// 3 cases

		mask = mask_nonzero(d) & mask_leq(L*2, N);

		for (i = 0; i < GFBITS; i++)
		{
		#pragma HLS unroll
			//This will either be all 0's or all 1's
			d_vec[i] = (d.test(i)) ? (0xffffffffffffffff):(0);
			b_vec[i] = (b.test(i)) ? (0xffffffffffffffff):(0);
		}

		vec_mul(B_tmp, d_vec, B);
		vec_mul(C_tmp, b_vec, C);

		vec_cmov(B, C, mask);
		update(B, mask & c0);

		for (i = 0; i < GFBITS; i++){
		#pragma HLS unroll
			C[i] = B_tmp[i] ^ C_tmp[i];
		}

		c0 = t >> 32;
		b = (d & mask) | (b & ~mask);
		L = ((N+1-L) & mask) | (L & ~mask);

	}

	c0 = gf_inv(c0);

	for (i = 0; i < GFBITS; i++){
		tmp[i] = vec_setbits((c0 >> i) & 1);
	}

	vec_mul(tmp, tmp, C);

	store_to_stream(out_stream, tmp);


}



/* input: in, result of applying the radix conversions to the input polynomial */
/* output: out, evaluation results (by applying the FFT butterflies) */
static void butterflies(vec out[][ GFBITS ], vec *in)
{
	int i, j, k, s, b;

	vec tmp[ GFBITS ];
	vec consts[ 63 ][ GFBITS ] =
	{
#include "consts.data"
	};

	const vec powers[ 64 ][ GFBITS ] =
	{
#include "powers.data"
	};

	vec consts_ptr = 0;

	const unsigned char reversal[64] =
	{
	  0, 32, 16, 48,  8, 40, 24, 56,
	  4, 36, 20, 52, 12, 44, 28, 60,
	  2, 34, 18, 50, 10, 42, 26, 58,
	  6, 38, 22, 54, 14, 46, 30, 62,
	  1, 33, 17, 49,  9, 41, 25, 57,
	  5, 37, 21, 53, 13, 45, 29, 61,
	  3, 35, 19, 51, 11, 43, 27, 59,
	  7, 39, 23, 55, 15, 47, 31, 63
	};

	// boradcast

	for (j = 0; j < 64; j++){
	#pragma HLS pipeline
		for (i = 0; i < GFBITS; i++)
		{
				out[j][i] = (in[i] >> reversal[j]) & 1;
				out[j][i] = -out[j][i];
		}
	}

	// butterflies

	for (i = 0; i <= 5; i++)
	{
		s = 1 << i;

		for (j = 0; j < 64; j += 2*s)
		#pragma HLS loop_tripcount min=1 max=32
		{
			for (k = j; k < j+s; k++)
			#pragma HLS pipeline
			#pragma HLS loop_tripcount min=1 max=32
			{

				vec_mul(tmp, out[k+s], consts[ consts_ptr + (k-j) ]);

				for (b = 0; b < GFBITS; b++){
					out[k][b] ^= tmp[b];
				}
				for (b = 0; b < GFBITS; b++){
					out[k+s][b] ^= out[k][b];
				}
			}
		}

		consts_ptr += (1 << i);
	}

	//

	// adding the part contributed by x^64

	for (i = 0; i < 64; i++){
	#pragma HLS pipeline
		for (b = 0; b < GFBITS; b++){
			out[i][b] ^= powers[i][b];
		}

	}

}




/* input: in, polynomial in bitsliced form */
/* output: in, result of applying the radix conversions on in */
void radix_conversions_stream(hls::stream<vec> &out_stream, hls::stream<vec> &in_stream)
{
	int i, j, k;

	vec in[GFBITS];

	const vec mask[5][2] =
	{
		{0x8888888888888888, 0x4444444444444444},
		{0xC0C0C0C0C0C0C0C0, 0x3030303030303030},
		{0xF000F000F000F000, 0x0F000F000F000F00},
		{0xFF000000FF000000, 0x00FF000000FF0000},
		{0xFFFF000000000000, 0x0000FFFF00000000}
	};

	const vec s[5][GFBITS] =
	{
#include "scalars.data"
	};

	#pragma HLS array_partition variable=in complete
	#pragma HLS array_partition variable=mask complete

	//TODO maybe this could go back to the loop?
	for (i = 0; i < GFBITS; i++){
		in[i] = in_stream.read();
	}


	//TODO maybe decoupe in order to flatten loops?
	for (j = 0; j <= 4; j++)
	{
		for (i = 0; i < GFBITS; i++){
	#pragma HLS pipeline

			for (k = 4; k >= j; k--)
			//#pragma HLS loop_tripcount min=1 max=5 avg=3
			{
				in[i] ^= (in[i] & mask[k][0]) >> (1 << k);
				in[i] ^= (in[i] & mask[k][1]) >> (1 << k);
			}
		}

		vec_mul(in, in, s[j]); // scaling
	}


	for(i=0; i<64; i++){
		for(j=0; j<GFBITS; j++){
			out_stream.write(in[j]);
	//		printf("out[%d] = %lX\n", i, out[i].to_uint64());
		}
	}

}

/* input: in, result of applying the radix conversions to the input polynomial */
/* output: out, evaluation results (by applying the FFT butterflies) */
static void butterflies_stream(hls::stream<vec> &out_stream, vec *in)
{
	int i, j, k, s, b;

	vec out[64][GFBITS];
	vec tmp[ GFBITS ];
	vec consts[ 63 ][ GFBITS ] =
	{
#include "consts.data"
	};

	const vec powers[ 64 ][ GFBITS ] =
	{
#include "powers.data"
	};

	vec consts_ptr = 0;

	const unsigned char reversal[64] =
	{
	  0, 32, 16, 48,  8, 40, 24, 56,
	  4, 36, 20, 52, 12, 44, 28, 60,
	  2, 34, 18, 50, 10, 42, 26, 58,
	  6, 38, 22, 54, 14, 46, 30, 62,
	  1, 33, 17, 49,  9, 41, 25, 57,
	  5, 37, 21, 53, 13, 45, 29, 61,
	  3, 35, 19, 51, 11, 43, 27, 59,
	  7, 39, 23, 55, 15, 47, 31, 63
	};

	// boradcast

	for (j = 0; j < 64; j++){
	#pragma HLS pipeline
		for (i = 0; i < GFBITS; i++)
		{
				out[j][i] = (in[i] >> reversal[j]) & 1;
				out[j][i] = -out[j][i];
		}
	}

	// butterflies

	for (i = 0; i <= 5; i++)
	{
		s = 1 << i;

		for (j = 0; j < 64; j += 2*s)
		#pragma HLS loop_tripcount min=1 max=32
		{
			for (k = j; k < j+s; k++)
			#pragma HLS pipeline
			#pragma HLS loop_tripcount min=1 max=32
			{

				vec_mul(tmp, out[k+s], consts[ consts_ptr + (k-j) ]);

				for (b = 0; b < GFBITS; b++){
					out[k][b] ^= tmp[b];
				}
				for (b = 0; b < GFBITS; b++){
					out[k+s][b] ^= out[k][b];
				}
			}
		}

		consts_ptr += (1 << i);
	}

	//

	// adding the part contributed by x^64

	for (i = 0; i < 64; i++){
	#pragma HLS pipeline
		for (b = 0; b < GFBITS; b++){
			out[i][b] ^= powers[i][b];
			out_stream.write(out[i][b]);
//			printf("out[%d][%d] = %lX\n", i, b, out[i][b].to_uint64());
		}

	}

}


//void fft(vec out[][ GFBITS ], vec *in)
//{
//		radix_conversions(in);
//		butterflies(out, in);
//}


void broadcast(hls::stream<mat_struct> &out_stream, hls::stream<vec> &in_stream){

	int i, j;
	mat_struct tmp_struct;
	vec tmp_val;

	const unsigned char reversal[64] =
	{
	  0, 32, 16, 48,  8, 40, 24, 56,
	  4, 36, 20, 52, 12, 44, 28, 60,
	  2, 34, 18, 50, 10, 42, 26, 58,
	  6, 38, 22, 54, 14, 46, 30, 62,
	  1, 33, 17, 49,  9, 41, 25, 57,
	  5, 37, 21, 53, 13, 45, 29, 61,
	  3, 35, 19, 51, 11, 43, 27, 59,
	  7, 39, 23, 55, 15, 47, 31, 63
	};

#pragma HLS array_partition variable=reversal factor=4

	for (j = 0; j < 64; j++){
	#pragma HLS pipeline
		for (i = 0; i < GFBITS; i++)
		{
//				out[j][i] = (in_stream.read() >> reversal[j]) & 1;
//				out[j][i] = -out[j][i];
				tmp_val = (in_stream.read() >> reversal[j]) & 1;
//				tmp_val = -tmp_val;
				tmp_struct.mat[i] = -tmp_val;
		}
		out_stream.write(tmp_struct);
	}

}


void write_output(hls::stream<vec> &out_stream, hls::stream<mat_struct> &in_stream_0, hls::stream<mat_struct> &in_stream_1){

	int i, b;

	vec tmp_val_0, tmp_val_1;
	mat_struct tmp_struct;

	const vec powers[ 64 ][ GFBITS ] =
	{
#include "powers.data"
	};



	for (i = 0; i < 32; i++){
	#pragma HLS pipeline
		tmp_struct = in_stream_0.read();
		for (b = 0; b < GFBITS; b++){
			tmp_val_0 = tmp_struct.mat[b];
			tmp_val_0 ^= powers[i][b];
			out_stream.write(tmp_val_0);

		}

	}

	for (i = 32; i < 64; i++){
	#pragma HLS pipeline
		tmp_struct = in_stream_1.read();
		for (b = 0; b < GFBITS; b++){
			tmp_val_1 = tmp_struct.mat[b];
//			tmp_val_1 = in_stream_1.read();
			tmp_val_1 ^= powers[i][b];
			out_stream.write(tmp_val_1);

		}

	}
}


//TODO has to read the wing streams and write accordingly to the output
void butterfly_operation(mat_struct *op_a, mat_struct *op_b, vec consts[63][GFBITS], int index){

	vec tmp[GFBITS];
	vec_mul(tmp, op_b->mat, consts[index]);
	vec_add(op_a->mat, op_a->mat, tmp);
	vec_add(op_b->mat, op_b->mat, op_a->mat);

}



static void fft_layer(hls::stream<mat_struct> &sink_0, hls::stream<mat_struct> &sink_1, hls::stream<mat_struct> &src_0, hls::stream<mat_struct> &src_1, int start, int reps)
{

	vec consts[ 63 ][ GFBITS ] =
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

	for (unsigned int i = 0; i < 32; i++) {
		#pragma HLS pipeline

//		printf("Consts index=%d\n", index);

		if (i % 32 < 16) {
			wing1 = src_0.read();
			wing2 = src_0.read();

		}
		else {
			wing1 = src_1.read();
			wing2 = src_1.read();
		}

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


/* input: in, polynomial in bitsliced form */
/* output: in, result of applying the radix conversions on in */

//TODO decouple in and put an out variable to pipeline better?
static void radix_conversions(vec *in)
{
	int i, j, k;

	const vec mask[5][2] =
	{
		{0x8888888888888888, 0x4444444444444444},
		{0xC0C0C0C0C0C0C0C0, 0x3030303030303030},
		{0xF000F000F000F000, 0x0F000F000F000F00},
		{0xFF000000FF000000, 0x00FF000000FF0000},
		{0xFFFF000000000000, 0x0000FFFF00000000}
	};

	const vec s[5][GFBITS] =
	{
#include "scalars.data"
	};

	//

	//TODO maybe decoupe in order to flatten loops?
	for (j = 0; j <= 4; j++)
	{
		for (i = 0; i < GFBITS; i++){
	#pragma HLS pipeline
			for (k = 4; k >= j; k--)
			#pragma HLS loop_tripcount min=1 max=5
			{
				in[i] ^= (in[i] & mask[k][0]) >> (1 << k);
				in[i] ^= (in[i] & mask[k][1]) >> (1 << k);
			}
		}

		vec_mul(in, in, s[j]); // scaling
	}
}


void radix_broadcast(hls::stream<mat_struct> &out_stream, hls::stream<vec> &in_stream){

	int i, j, k;

	vec in[GFBITS];

	const vec mask[5][2] =
	{
		{0x8888888888888888, 0x4444444444444444},
		{0xC0C0C0C0C0C0C0C0, 0x3030303030303030},
		{0xF000F000F000F000, 0x0F000F000F000F00},
		{0xFF000000FF000000, 0x00FF000000FF0000},
		{0xFFFF000000000000, 0x0000FFFF00000000}
	};

	const vec s[5][GFBITS] =
	{
#include "scalars.data"
	};

	//Probably the reason it was failing!
//	#pragma HLS array_partition variable=in complete
//	#pragma HLS array_partition variable=mask complete

	//TODO maybe this could go back to the loop?
	for (i = 0; i < GFBITS; i++){
		in[i] = in_stream.read();
	}


	//TODO maybe decoupe in order to flatten loops?
	for (j = 0; j <= 4; j++)
	{
		for (i = 0; i < GFBITS; i++){
	#pragma HLS pipeline

			for (k = 4; k >= j; k--)
			//#pragma HLS loop_tripcount min=1 max=5 avg=3
			{
				in[i] ^= (in[i] & mask[k][0]) >> (1 << k);
				in[i] ^= (in[i] & mask[k][1]) >> (1 << k);
			}
		}

		vec_mul(in, in, s[j]); // scaling
	}


//	for(i=0; i<64; i++){
//		for(j=0; j<GFBITS; j++){
//			out_stream.write(in[j]);
//	//		printf("out[%d] = %lX\n", i, out[i].to_uint64());
//		}
//	}


	mat_struct tmp_struct;
	vec tmp_val;

	const unsigned char reversal[64] =
	{
	  0, 32, 16, 48,  8, 40, 24, 56,
	  4, 36, 20, 52, 12, 44, 28, 60,
	  2, 34, 18, 50, 10, 42, 26, 58,
	  6, 38, 22, 54, 14, 46, 30, 62,
	  1, 33, 17, 49,  9, 41, 25, 57,
	  5, 37, 21, 53, 13, 45, 29, 61,
	  3, 35, 19, 51, 11, 43, 27, 59,
	  7, 39, 23, 55, 15, 47, 31, 63
	};

	//Probably the reason it was failing!
//#pragma HLS array_partition variable=reversal factor=4

	for (j = 0; j < 64; j++){
	#pragma HLS pipeline
		for (i = 0; i < GFBITS; i++)
		{
//				out[j][i] = (in_stream.read() >> reversal[j]) & 1;
//				out[j][i] = -out[j][i];
				tmp_val = (in[i] >> reversal[j]) & 1;
//				tmp_val = -tmp_val;
				tmp_struct.mat[i] = -tmp_val;
		}
		out_stream.write(tmp_struct);
	}





}




void fft_stream_alt(hls::stream<vec> &out_stream, hls::stream<vec> &in_stream)
{
	static hls::stream<mat_struct> scaled_struct_stream_2("scaled_struct_stream_2");
	#pragma HLS STREAM variable=scaled_struct_stream_2 depth=64

	static hls::stream<vec> int_mat_1_stream("int_mat_1_stream");
	#pragma HLS STREAM variable=int_mat_1_stream depth=12

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

	#pragma HLS dataflow

//	radix_conversions_stream(int_mat_1_stream, in_stream);
//
//	broadcast(scaled_struct_stream_2, int_mat_1_stream);

	radix_broadcast(scaled_struct_stream_2, in_stream);

	split_struct_stream(struct_stream_2_0_0, struct_stream_2_0_1, scaled_struct_stream_2);
	fft_layer(struct_stream_2_1_0, struct_stream_2_1_1, struct_stream_2_0_0, struct_stream_2_0_1, 0, 32);
	fft_layer(struct_stream_2_2_0, struct_stream_2_2_1, struct_stream_2_1_0, struct_stream_2_1_1, 1, 16);
	fft_layer(struct_stream_2_3_0, struct_stream_2_3_1, struct_stream_2_2_0, struct_stream_2_2_1, 3, 8);
	fft_layer(struct_stream_2_4_0, struct_stream_2_4_1, struct_stream_2_3_0, struct_stream_2_3_1, 7, 4);
	fft_layer(struct_stream_2_5_0, struct_stream_2_5_1, struct_stream_2_4_0, struct_stream_2_4_1, 15, 2);
	fft_layer(struct_stream_2_6_0, struct_stream_2_6_1, struct_stream_2_5_0, struct_stream_2_5_1, 31, 1);
	write_output(out_stream, struct_stream_2_6_0, struct_stream_2_6_1);


}
//////////////////FFT//////////////////////


//void fft_stream(hls::stream<vec> &out_stream, hls::stream<vec> &in_stream)
//{
//	vec int_mat[GFBITS];
//	radix_conversions_stream(int_mat, in_stream);
//	butterflies_stream(out_stream, int_mat);
//
//}


static inline vec load8(const unsigned char * in)
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


/* one layer of the benes network */
/*parameters, r, cond , low*/
//TODO fix here to be able to pipeline and unroll...like 1st solution?
static void layer(vec * data, vec * bits, int lgs)
{
	int i, j, s;

	vec d;

	s = 1 << lgs;

	BENES_LAYER:
	for (i = 0; i < 64; i += s*2){
	#pragma HLS loop_tripcount min=1 max=32
		for (j = i; j < i+s; j++)
		#pragma HLS loop_tripcount min=1 max=32
		{

			d = (data[j+0] ^ data[j+s]);
			//first return bits value then increment
			d &= (*bits++);
			data[j+0] ^= d;
			data[j+s] ^= d;
		}
	}
}




static inline uint32_t load4(const unsigned char *src)
{
	uint32_t a;

	a  = src[3];
	a <<= 8;
	a |= src[2];
	a <<= 8;
	a |= src[1];
	a <<= 8;
	a |= src[0];

	return a;
}


/* bitsliced field inverses */
void vec_inv(vec *out, vec *in)
{
		vec tmp_11[GFBITS];
		vec tmp_1111[GFBITS];

#pragma HLS array partition variable=tmp_11 complete
#pragma HLS array partition variable=tmp_1111 complete

        vec_copy(out, in);

        vec_sq(out, out);
        vec_mul(tmp_11, out, in); // 11

        vec_sq(out, tmp_11);
        vec_sq(out, out);
        vec_mul(tmp_1111, out, tmp_11); // 1111

        vec_sq(out, tmp_1111);
        vec_sq(out, out);
        vec_sq(out, out);
        vec_sq(out, out);
        vec_mul(out, out, tmp_1111); // 11111111

        vec_sq(out, out);
        vec_sq(out, out);
        vec_mul(out, out, tmp_11); // 1111111111

        vec_sq(out, out);
        vec_mul(out, out, in); // 11111111111

        vec_sq(out, out); // 111111111110
}


void load_r(vec *r_mat, hls::stream<vec> &r_in_stream){
	int i, j;
	vec tmp_mat[64];

	for(i=0; i<64; i++){
		r_mat[i] = r_in_stream.read();
	}


}


//////////////////////////BENES LAYER///////////////////////////////


void benes_operation(vec *op_a, vec *op_b, vec mat, vec *d){

	*d = (*op_a ^ *op_b);
	*d &= mat;
	*op_a ^= *d;
	*op_b ^= *d;

}

static void benes_layer(hls::stream<vec> &sink_0, hls::stream<vec> &sink_1, hls::stream<vec> &src_0, hls::stream<vec> &src_1, hls::stream<mat_64_struct> &mat_stream, unsigned int layer_lvl)
{

	vec wing1;
	vec wing2;
	vec d=0;
	mat_64_struct tmp_struct;

	unsigned int index = 0;
	unsigned int start=0;
	unsigned int iter=0;

	unsigned int stride=1<<layer_lvl;

	tmp_struct = mat_stream.read();

	for (unsigned int i = 0; i < 32; i++) {
		#pragma HLS pipeline


		if (i % 32 < 16) {
			wing1 = src_0.read();
			wing2 = src_0.read();
		}
		else {
			wing1 = src_1.read();
			wing2 = src_1.read();

		}

		benes_operation(&wing1, &wing2, tmp_struct.mat[index], &d);

		sink_0.write(wing1);
		sink_1.write(wing2);

		index += stride;
		iter++;
		if(iter==(32>>layer_lvl)){
			start++;
			index=start;
			iter=0;
		}



	}

}

static void benes_tr_layer(hls::stream<vec> &sink_0, hls::stream<vec> &sink_1, hls::stream<vec> &src_0, hls::stream<vec> &src_1, hls::stream<mat_64_struct> &mat_stream, unsigned int layer_lvl)
{

	vec wing1;
	vec wing2;
	vec d=0;
	mat_64_struct tmp_struct;

	unsigned int index = 0;
	unsigned int start=0;
	unsigned int iter=0;

	unsigned int stride=1<<layer_lvl;

	tmp_struct = mat_stream.read();



	for (unsigned int i = 0; i < 32; i++) {
		#pragma HLS pipeline

		wing1 = src_0.read();
		wing2 = src_1.read();

		benes_operation(&wing1, &wing2, tmp_struct.mat[index], &d);

		if (i % 32 < 16) {
			sink_0.write(wing1);
			sink_0.write(wing2);

		}
		else {
			sink_1.write(wing1);
			sink_1.write(wing2);

		}

		index += stride;
		iter++;
		if(iter==(32>>layer_lvl)){
			start++;
			index=start;
			iter=0;
		}



	}

}


static void benes_short_layer(hls::stream<vec> &sink_0, hls::stream<vec> &sink_1, hls::stream<vec> &src_0, hls::stream<vec> &src_1, hls::stream<mat_64_struct> &mat_stream, unsigned int stride, unsigned int limit)
{

	vec wing1;
	vec wing2;
	vec d=0;
	mat_64_struct tmp_struct;

	unsigned int index = 0;
	unsigned int start=0;
	unsigned int iter=0;

	tmp_struct = mat_stream.read();


	for (unsigned int i = 0; i < 32; i++) {
		#pragma HLS pipeline


		wing1 = src_0.read();
		wing2 = src_1.read();


		benes_operation(&wing1, &wing2, tmp_struct.mat[index], &d);

		if ((i < 8)|(i>=16 && i<24)) {
			sink_0.write(wing1);
			sink_0.write(wing2);
		}
		else {
			sink_1.write(wing1);
			sink_1.write(wing2);
		}

		index += stride;
		iter++;

		if(iter%limit==0){
			start+=1;
			index=start;
		}
		if(iter==16){
			start=16;
			iter=0;
			index=start;
		}


	}

}


void split_benes_stream(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_struct_stream){

	int i;
	vec tmp;

	for(i=0; i<64; i++){
			tmp = in_struct_stream.read();
			if(i<32){
				out_stream_0.write(tmp);
			}else{
				out_stream_1.write(tmp);
			}
	}
}

void split_benes_short_stream(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1){

	int i;
	vec tmp;

	for(i=0; i<64; i++){
		if (i < 16) {
			tmp = in_stream_0.read();
			out_stream_0.write(tmp);
		}else if(i>=16 && i<32){
			tmp = in_stream_0.read();
			out_stream_1.write(tmp);
		}else if(i>=32 && i<48){
			tmp = in_stream_1.read();
			out_stream_0.write(tmp);
		}else{
			tmp = in_stream_1.read();
			out_stream_1.write(tmp);
		}
	}
}

void loop_cond_layer_1(hls::stream<mat_64_struct> &cond_out_stream_0, hls::stream<mat_64_struct> &cond_out_stream_1, hls::stream<mat_64_struct> &cond_out_stream_2, hls::stream<mat_64_struct> &cond_out_stream_3, hls::stream<mat_64_struct> &cond_out_stream_4, hls::stream<mat_64_struct> &cond_out_stream_5, sk_2k_vec *vec_mat){

	sk_2k_vec tmp;
	int low, i;
	vec cond[6][64];
	vec cond_tr[6][64];

	mat_64_struct tmp_mat_0_struct, tmp_mat_1_struct, tmp_mat_2_struct, tmp_mat_3_struct, tmp_mat_4_struct, tmp_mat_5_struct;

	BENES_LOOP_1:
	for (low = 0; low <= 5; low++)
	#pragma HLS pipeline
	{
		tmp = vec_mat[low];
//		tmp = vec_stream.read();
		for (i = 0; i < 64; i++){
			cond[low][i].range(31,0) = tmp.range(((i+1)*32)-1, (i*32));
			cond[low][i].range(63,32) = 0;
		}
		transpose_64x64(cond_tr[low], cond[low]);
	}

	for(i=0; i<64; i++){
		tmp_mat_0_struct.mat[i] = cond_tr[0][i];
		tmp_mat_1_struct.mat[i] = cond_tr[1][i];
		tmp_mat_2_struct.mat[i] = cond_tr[2][i];
		tmp_mat_3_struct.mat[i] = cond_tr[3][i];
		tmp_mat_4_struct.mat[i] = cond_tr[4][i];
		tmp_mat_5_struct.mat[i] = cond_tr[5][i];
	}

	cond_out_stream_0.write(tmp_mat_0_struct);
	cond_out_stream_1.write(tmp_mat_1_struct);
	cond_out_stream_2.write(tmp_mat_2_struct);
	cond_out_stream_3.write(tmp_mat_3_struct);
	cond_out_stream_4.write(tmp_mat_4_struct);
	cond_out_stream_5.write(tmp_mat_5_struct);

}

void loop_cond_layer_2(hls::stream<mat_64_struct> &cond_out_stream_0, hls::stream<mat_64_struct> &cond_out_stream_1, hls::stream<mat_64_struct> &cond_out_stream_2, hls::stream<mat_64_struct> &cond_out_stream_3, hls::stream<mat_64_struct> &cond_out_stream_4, hls::stream<mat_64_struct> &cond_out_stream_5, hls::stream<sk_2k_vec> &vec_stream){

	sk_2k_vec tmp;
	int low, i;
	vec cond[6][64];
	vec cond_tr[6][64];

	mat_64_struct tmp_mat_0_struct, tmp_mat_1_struct, tmp_mat_2_struct, tmp_mat_3_struct, tmp_mat_4_struct, tmp_mat_5_struct;

	BENES_LOOP_2:
	for (low = 0; low <= 5; low++)
	#pragma HLS pipeline
	{
//		tmp = vec_mat[low+6];
		tmp = vec_stream.read();
		for (i = 0; i < 32; i++){
			cond_tr[low][i] = tmp.range(((i+1)*64)-1, (i*64));
		}
	}

	for(i=0; i<64; i++){
		tmp_mat_0_struct.mat[i] = cond_tr[0][i];
		tmp_mat_1_struct.mat[i] = cond_tr[1][i];
		tmp_mat_2_struct.mat[i] = cond_tr[2][i];
		tmp_mat_3_struct.mat[i] = cond_tr[3][i];
		tmp_mat_4_struct.mat[i] = cond_tr[4][i];
		tmp_mat_5_struct.mat[i] = cond_tr[5][i];
	}

	cond_out_stream_0.write(tmp_mat_0_struct);
	cond_out_stream_1.write(tmp_mat_1_struct);
	cond_out_stream_2.write(tmp_mat_2_struct);
	cond_out_stream_3.write(tmp_mat_3_struct);
	cond_out_stream_4.write(tmp_mat_4_struct);
	cond_out_stream_5.write(tmp_mat_5_struct);

}


void preprocess_layer_3(hls::stream<mat_64_struct> &cond_out_stream_0, hls::stream<mat_64_struct> &cond_out_stream_1, hls::stream<mat_64_struct> &cond_out_stream_2, hls::stream<mat_64_struct> &cond_out_stream_3, hls::stream<mat_64_struct> &cond_out_stream_4, hls::stream<sk_2k_vec> &vec_stream){

	sk_2k_vec tmp;
	int low, i;
	vec cond[6][64];
	vec cond_tr[6][64];

	mat_64_struct tmp_mat_0_struct, tmp_mat_1_struct, tmp_mat_2_struct, tmp_mat_3_struct, tmp_mat_4_struct;


	BENES_LOOP_3:
	for (low = 4; low >= 0; low--)
	{
	#pragma HLS pipeline
//		tmp = vec_mat[16-low];
		tmp = vec_stream.read();
		for (i = 0; i < 32; i++){
			cond_tr[low][i] = tmp.range(((i+1)*64)-1, (i*64));

		}
	}

	for(i=0; i<64; i++){
		tmp_mat_0_struct.mat[i] = cond_tr[0][i];
		tmp_mat_1_struct.mat[i] = cond_tr[1][i];
		tmp_mat_2_struct.mat[i] = cond_tr[2][i];
		tmp_mat_3_struct.mat[i] = cond_tr[3][i];
		tmp_mat_4_struct.mat[i] = cond_tr[4][i];
	}

	cond_out_stream_0.write(tmp_mat_0_struct);
	cond_out_stream_1.write(tmp_mat_1_struct);
	cond_out_stream_2.write(tmp_mat_2_struct);
	cond_out_stream_3.write(tmp_mat_3_struct);
	cond_out_stream_4.write(tmp_mat_4_struct);

}


void loop_cond_layer_4(hls::stream<mat_64_struct> &cond_out_stream_0, hls::stream<mat_64_struct> &cond_out_stream_1, hls::stream<mat_64_struct> &cond_out_stream_2, hls::stream<mat_64_struct> &cond_out_stream_3, hls::stream<mat_64_struct> &cond_out_stream_4, hls::stream<mat_64_struct> &cond_out_stream_5, hls::stream<sk_2k_vec> &vec_stream){

	sk_2k_vec tmp;
	int low, i;
	vec cond[6][64];
	vec cond_tr[6][64];

	mat_64_struct tmp_mat_0_struct, tmp_mat_1_struct, tmp_mat_2_struct, tmp_mat_3_struct, tmp_mat_4_struct, tmp_mat_5_struct;


	BENES_LOOP_4:
	for (low = 5; low >= 0; low--)
	#pragma HLS pipeline
	{
//		tmp = vec_mat[-low+22];
		tmp = vec_stream.read();
		for (i = 0; i < 64; i++){
			cond[low][i].range(31,0) = tmp.range(((i+1)*32)-1, (i*32));
			cond[low][i].range(63,32) = 0;
		}
		transpose_64x64(cond_tr[low], cond[low]);
	}

	for(i=0; i<64; i++){
		tmp_mat_0_struct.mat[i] = cond_tr[0][i];
		tmp_mat_1_struct.mat[i] = cond_tr[1][i];
		tmp_mat_2_struct.mat[i] = cond_tr[2][i];
		tmp_mat_3_struct.mat[i] = cond_tr[3][i];
		tmp_mat_4_struct.mat[i] = cond_tr[4][i];
		tmp_mat_5_struct.mat[i] = cond_tr[5][i];
	}

	cond_out_stream_0.write(tmp_mat_0_struct);
	cond_out_stream_1.write(tmp_mat_1_struct);
	cond_out_stream_2.write(tmp_mat_2_struct);
	cond_out_stream_3.write(tmp_mat_3_struct);
	cond_out_stream_4.write(tmp_mat_4_struct);
	cond_out_stream_5.write(tmp_mat_5_struct);

}

void split_vec_stream(sk_2k_vec *vec_mat_out_0, hls::stream<sk_2k_vec> &vec_stream_out_1, hls::stream<sk_2k_vec> &vec_stream_out_2, hls::stream<sk_2k_vec> &vec_stream_out_3, hls::stream<sk_2k_vec> &vec_stream_in){

	int i;
	sk_2k_vec tmp_vec;

	for(i=0; i<6; i++){
		tmp_vec = vec_stream_in.read();
		vec_mat_out_0[i] = tmp_vec;
	}

	for(i=0; i<6; i++){
		tmp_vec = vec_stream_in.read();
		vec_stream_out_1.write(tmp_vec);
	}

	for(i=0; i<5; i++){
		tmp_vec = vec_stream_in.read();
		vec_stream_out_2.write(tmp_vec);
	}

	for(i=0; i<6; i++){
		tmp_vec = vec_stream_in.read();
		vec_stream_out_3.write(tmp_vec);
	}

}


void preprocess_layer_1(hls::stream<mat_64_struct> &cond_stream_1_0, hls::stream<mat_64_struct> &cond_stream_1_1, hls::stream<mat_64_struct> &cond_stream_1_2, hls::stream<mat_64_struct> &cond_stream_1_3, hls::stream<mat_64_struct> &cond_stream_1_4, hls::stream<mat_64_struct> &cond_stream_1_5, hls::stream<sk_2k_vec> &vec_mat_stream_2, hls::stream<sk_2k_vec> &vec_mat_stream_3, hls::stream<sk_2k_vec> &vec_mat_stream_4, hls::stream<vec> &benes_stream_1, hls::stream<vec> &r_in, hls::stream<sk_2k_vec> &vec_bits){

	vec r_mat[64];
	sk_2k_vec vec_mat_1[6];


	load_r(r_mat, r_in);
	transpose_64x64_stream_out(benes_stream_1, r_mat);
	split_vec_stream(vec_mat_1, vec_mat_stream_2, vec_mat_stream_3, vec_mat_stream_4, vec_bits);
	loop_cond_layer_1(cond_stream_1_0, cond_stream_1_1, cond_stream_1_2, cond_stream_1_3, cond_stream_1_4, cond_stream_1_5, vec_mat_1);

}


void preprocess_layer_2(hls::stream<mat_64_struct> &cond_stream_2_0, hls::stream<mat_64_struct> &cond_stream_2_1, hls::stream<mat_64_struct> &cond_stream_2_2, hls::stream<mat_64_struct> &cond_stream_2_3, hls::stream<mat_64_struct> &cond_stream_2_4, hls::stream<mat_64_struct> &cond_stream_2_5, hls::stream<vec> &benes_stream_2, hls::stream<sk_2k_vec> &vec_mat_stream_2, hls::stream<vec> &benes_stream_1_6_0, hls::stream<vec> &benes_stream_1_6_1){

	int i;
	vec r1[64];

	loop_cond_layer_2(cond_stream_2_0, cond_stream_2_1, cond_stream_2_2, cond_stream_2_3, cond_stream_2_4, cond_stream_2_5, vec_mat_stream_2);
	for(i=0; i<32; i++){
		r1[i] = benes_stream_1_6_0.read();
		r1[i+32] = benes_stream_1_6_1.read();
	}
	transpose_64x64_stream_out(benes_stream_2, r1);

}

void preprocess_layer_4(hls::stream<mat_64_struct> &cond_stream_4_0, hls::stream<mat_64_struct> &cond_stream_4_1, hls::stream<mat_64_struct> &cond_stream_4_2, hls::stream<mat_64_struct> &cond_stream_4_3, hls::stream<mat_64_struct> &cond_stream_4_4, hls::stream<mat_64_struct> &cond_stream_4_5, hls::stream<vec> &benes_stream_4, hls::stream<sk_2k_vec> &vec_mat_stream_4, hls::stream<vec> &benes_stream_3_5_0, hls::stream<vec> &benes_stream_3_5_1){

	int i;
	vec r2[64];

	for(i=0; i<32; i++){

		if(i<16){
			r2[i] = benes_stream_3_5_0.read();
			r2[i+16] = benes_stream_3_5_1.read();
		}else{
			r2[i+16] =  benes_stream_3_5_0.read();
			r2[i+32] =  benes_stream_3_5_1.read();
		}
	}
	transpose_64x64_stream_out(benes_stream_4, r2);
	loop_cond_layer_4(cond_stream_4_0, cond_stream_4_1, cond_stream_4_2, cond_stream_4_3, cond_stream_4_4, cond_stream_4_5, vec_mat_stream_4);

}

void postprocess_layer_4(hls::stream<vec> &r_out, hls::stream<vec> &benes_stream_4_6_0, hls::stream<vec> &benes_stream_4_6_1){

	int i;
	vec r3[64];

	for(i=0; i<32; i++){
		r3[i] = benes_stream_4_6_0.read();
		r3[i+32] = benes_stream_4_6_1.read();
	}
	transpose_64x64_stream_out(r_out, r3);

}

/* input: r, sequence of bits to be permuted */
/*        bits, condition bits of the Benes network */
/*        rev, 0 for normal application; !0 for inverse */
/* output: r, permuted bits */
void benes1(hls::stream<vec> &r_out, hls::stream<vec> &r_in, hls::stream<sk_2k_vec> &vec_bits)
{

	int i;


	const unsigned char *cond_ptr;
	int inc, low;

	vec r1[64], r2[64], r3[64], r4[64];

//	vec r_mat[64];
	sk_2k_vec tmp;

	sk_2k_vec vec_mat[23];

	vec cond[6][64];
	vec cond_tr[6][64];
	vec cond_tr_2[6][64];



static hls::stream<vec> cond_tr_stream("cond_tr_stream");
#pragma HLS STREAM variable=cond_tr_stream depth=2

static hls::stream<vec> benes_stream_1("benes_stream_1");
#pragma HLS STREAM variable=benes_stream_1 depth=64

static hls::stream<vec> benes_stream_1_0_0("benes_stream_1_0_0");
#pragma HLS STREAM variable=benes_stream_1_0_0 depth=64

static hls::stream<vec> benes_stream_1_0_1("benes_stream_1_0_1");
#pragma HLS STREAM variable=benes_stream_1_0_1 depth=64

static hls::stream<vec> benes_stream_1_1_0("benes_stream_1_1_0");
#pragma HLS STREAM variable=benes_stream_1_1_0 depth=64

static hls::stream<vec> benes_stream_1_1_1("benes_stream_1_1_1");
#pragma HLS STREAM variable=benes_stream_1_1_1 depth=64

static hls::stream<vec> benes_stream_1_2_0("benes_stream_1_2_0");
#pragma HLS STREAM variable=benes_stream_1_2_0 depth=64

static hls::stream<vec> benes_stream_1_2_1("benes_stream_1_2_1");
#pragma HLS STREAM variable=benes_stream_1_2_1 depth=64

static hls::stream<vec> benes_stream_1_3_0("benes_stream_1_3_0");
#pragma HLS STREAM variable=benes_stream_1_3_0 depth=64

static hls::stream<vec> benes_stream_1_3_1("benes_stream_1_3_1");
#pragma HLS STREAM variable=benes_stream_1_3_1 depth=64

static hls::stream<vec> benes_stream_1_4_0("benes_stream_1_4_0");
#pragma HLS STREAM variable=benes_stream_1_4_0 depth=64

static hls::stream<vec> benes_stream_1_4_1("benes_stream_1_4_1");
#pragma HLS STREAM variable=benes_stream_1_4_1 depth=64

static hls::stream<vec> benes_stream_1_5_0("benes_stream_1_5_0");
#pragma HLS STREAM variable=benes_stream_1_5_0 depth=64

static hls::stream<vec> benes_stream_1_5_1("benes_stream_1_5_1");
#pragma HLS STREAM variable=benes_stream_1_5_1 depth=64

static hls::stream<vec> benes_stream_1_6_0("benes_stream_1_6_0");
#pragma HLS STREAM variable=benes_stream_1_6_0 depth=64

static hls::stream<vec> benes_stream_1_6_1("benes_stream_1_6_1");
#pragma HLS STREAM variable=benes_stream_1_6_1 depth=64
//

static hls::stream<vec> benes_stream_2("benes_stream_2");
#pragma HLS STREAM variable=benes_stream_2 depth=64

static hls::stream<vec> benes_stream_2_0_0("benes_stream_2_0_0");
#pragma HLS STREAM variable=benes_stream_2_0_0 depth=64

static hls::stream<vec> benes_stream_2_0_1("benes_stream_2_0_1");
#pragma HLS STREAM variable=benes_stream_2_0_1 depth=64

static hls::stream<vec> benes_stream_2_1_0("benes_stream_2_1_0");
#pragma HLS STREAM variable=benes_stream_2_1_0 depth=64

static hls::stream<vec> benes_stream_2_1_1("benes_stream_2_1_1");
#pragma HLS STREAM variable=benes_stream_2_1_1 depth=64

static hls::stream<vec> benes_stream_2_2_0("benes_stream_2_2_0");
#pragma HLS STREAM variable=benes_stream_2_2_0 depth=64

static hls::stream<vec> benes_stream_2_2_1("benes_stream_2_2_1");
#pragma HLS STREAM variable=benes_stream_2_2_1 depth=64

static hls::stream<vec> benes_stream_2_3_0("benes_stream_2_3_0");
#pragma HLS STREAM variable=benes_stream_2_3_0 depth=64

static hls::stream<vec> benes_stream_2_3_1("benes_stream_2_3_1");
#pragma HLS STREAM variable=benes_stream_2_3_1 depth=64

static hls::stream<vec> benes_stream_2_4_0("benes_stream_2_4_0");
#pragma HLS STREAM variable=benes_stream_2_4_0 depth=64

static hls::stream<vec> benes_stream_2_4_1("benes_stream_2_4_1");
#pragma HLS STREAM variable=benes_stream_2_4_1 depth=64

static hls::stream<vec> benes_stream_2_5_0("benes_stream_2_5_0");
#pragma HLS STREAM variable=benes_stream_2_5_0 depth=64

static hls::stream<vec> benes_stream_2_5_1("benes_stream_2_5_1");
#pragma HLS STREAM variable=benes_stream_2_5_1 depth=64

static hls::stream<vec> benes_stream_2_6_0("benes_stream_2_6_0");
#pragma HLS STREAM variable=benes_stream_2_6_0 depth=64

static hls::stream<vec> benes_stream_2_6_1("benes_stream_2_6_1");
#pragma HLS STREAM variable=benes_stream_2_6_1 depth=64
//

static hls::stream<vec> benes_stream_3("benes_stream_3");
#pragma HLS STREAM variable=benes_stream_3 depth=64

static hls::stream<vec> benes_stream_3_0_0("benes_stream_3_0_0");
#pragma HLS STREAM variable=benes_stream_3_0_0 depth=64

static hls::stream<vec> benes_stream_3_0_1("benes_stream_3_0_1");
#pragma HLS STREAM variable=benes_stream_3_0_1 depth=64

static hls::stream<vec> benes_stream_3_1_0("benes_stream_3_1_0");
#pragma HLS STREAM variable=benes_stream_3_1_0 depth=64

static hls::stream<vec> benes_stream_3_1_1("benes_stream_3_1_1");
#pragma HLS STREAM variable=benes_stream_3_1_1 depth=64

static hls::stream<vec> benes_stream_3_2_0("benes_stream_3_2_0");
#pragma HLS STREAM variable=benes_stream_3_2_0 depth=64

static hls::stream<vec> benes_stream_3_2_1("benes_stream_3_2_1");
#pragma HLS STREAM variable=benes_stream_3_2_1 depth=64

static hls::stream<vec> benes_stream_3_3_0("benes_stream_3_3_0");
#pragma HLS STREAM variable=benes_stream_3_3_0 depth=64

static hls::stream<vec> benes_stream_3_3_1("benes_stream_3_3_1");
#pragma HLS STREAM variable=benes_stream_3_3_1 depth=64

static hls::stream<vec> benes_stream_3_4_0("benes_stream_3_4_0");
#pragma HLS STREAM variable=benes_stream_3_4_0 depth=64

static hls::stream<vec> benes_stream_3_4_1("benes_stream_3_4_1");
#pragma HLS STREAM variable=benes_stream_3_4_1 depth=64

static hls::stream<vec> benes_stream_3_5_0("benes_stream_3_5_0");
#pragma HLS STREAM variable=benes_stream_3_5_0 depth=64

static hls::stream<vec> benes_stream_3_5_1("benes_stream_3_5_1");
#pragma HLS STREAM variable=benes_stream_3_5_1 depth=64

static hls::stream<vec> benes_stream_3_6_0("benes_stream_3_6_0");
#pragma HLS STREAM variable=benes_stream_3_6_0 depth=64

static hls::stream<vec> benes_stream_3_6_1("benes_stream_3_6_1");
#pragma HLS STREAM variable=benes_stream_3_6_1 depth=64
//////

static hls::stream<vec> benes_stream_4("benes_stream_4");
#pragma HLS STREAM variable=benes_stream_4 depth=64

static hls::stream<vec> benes_stream_4_0_0("benes_stream_4_0_0");
#pragma HLS STREAM variable=benes_stream_4_0_0 depth=64

static hls::stream<vec> benes_stream_4_0_1("benes_stream_4_0_1");
#pragma HLS STREAM variable=benes_stream_4_0_1 depth=64

static hls::stream<vec> benes_stream_4_1_0("benes_stream_4_1_0");
#pragma HLS STREAM variable=benes_stream_4_1_0 depth=64

static hls::stream<vec> benes_stream_4_1_1("benes_stream_4_1_1");
#pragma HLS STREAM variable=benes_stream_4_1_1 depth=64

static hls::stream<vec> benes_stream_4_2_0("benes_stream_4_2_0");
#pragma HLS STREAM variable=benes_stream_4_2_0 depth=64

static hls::stream<vec> benes_stream_4_2_1("benes_stream_4_2_1");
#pragma HLS STREAM variable=benes_stream_4_2_1 depth=64

static hls::stream<vec> benes_stream_4_3_0("benes_stream_4_3_0");
#pragma HLS STREAM variable=benes_stream_4_3_0 depth=64

static hls::stream<vec> benes_stream_4_3_1("benes_stream_4_3_1");
#pragma HLS STREAM variable=benes_stream_4_3_1 depth=64

static hls::stream<vec> benes_stream_4_4_0("benes_stream_4_4_0");
#pragma HLS STREAM variable=benes_stream_4_4_0 depth=64

static hls::stream<vec> benes_stream_4_4_1("benes_stream_4_4_1");
#pragma HLS STREAM variable=benes_stream_4_4_1 depth=64

static hls::stream<vec> benes_stream_4_5_0("benes_stream_4_5_0");
#pragma HLS STREAM variable=benes_stream_4_5_0 depth=64

static hls::stream<vec> benes_stream_4_5_1("benes_stream_4_5_1");
#pragma HLS STREAM variable=benes_stream_4_5_1 depth=64

static hls::stream<vec> benes_stream_4_6_0("benes_stream_4_6_0");
#pragma HLS STREAM variable=benes_stream_4_6_0 depth=64

static hls::stream<vec> benes_stream_4_6_1("benes_stream_4_6_1");
#pragma HLS STREAM variable=benes_stream_4_6_1 depth=64

//
static hls::stream<mat_64_struct> cond_stream_1_0("cond_stream_1_0");
#pragma HLS STREAM variable=cond_stream_1_0 depth=2

static hls::stream<mat_64_struct> cond_stream_1_1("cond_stream_1_1");
#pragma HLS STREAM variable=cond_stream_1_1 depth=2

static hls::stream<mat_64_struct> cond_stream_1_2("cond_stream_1_2");
#pragma HLS STREAM variable=cond_stream_1_2 depth=2

static hls::stream<mat_64_struct> cond_stream_1_3("cond_stream_1_3");
#pragma HLS STREAM variable=cond_stream_1_3 depth=2

static hls::stream<mat_64_struct> cond_stream_1_4("cond_stream_1_4");
#pragma HLS STREAM variable=cond_stream_1_4 depth=2

static hls::stream<mat_64_struct> cond_stream_1_5("cond_stream_1_5");
#pragma HLS STREAM variable=cond_stream_1_5 depth=2

//
static hls::stream<mat_64_struct> cond_stream_2_0("cond_stream_2_0");
#pragma HLS STREAM variable=cond_stream_2_0 depth=2

static hls::stream<mat_64_struct> cond_stream_2_1("cond_stream_2_1");
#pragma HLS STREAM variable=cond_stream_2_1 depth=2

static hls::stream<mat_64_struct> cond_stream_2_2("cond_stream_2_2");
#pragma HLS STREAM variable=cond_stream_2_2 depth=2

static hls::stream<mat_64_struct> cond_stream_2_3("cond_stream_2_3");
#pragma HLS STREAM variable=cond_stream_2_3 depth=2

static hls::stream<mat_64_struct> cond_stream_2_4("cond_stream_2_4");
#pragma HLS STREAM variable=cond_stream_2_4 depth=2

static hls::stream<mat_64_struct> cond_stream_2_5("cond_stream_2_5");
#pragma HLS STREAM variable=cond_stream_2_5 depth=2

//
static hls::stream<mat_64_struct> cond_stream_3_0("cond_stream_3_0");
#pragma HLS STREAM variable=cond_stream_3_0 depth=2

static hls::stream<mat_64_struct> cond_stream_3_1("cond_stream_3_1");
#pragma HLS STREAM variable=cond_stream_3_1 depth=2

static hls::stream<mat_64_struct> cond_stream_3_2("cond_stream_3_2");
#pragma HLS STREAM variable=cond_stream_3_2 depth=2

static hls::stream<mat_64_struct> cond_stream_3_3("cond_stream_3_3");
#pragma HLS STREAM variable=cond_stream_3_3 depth=2

static hls::stream<mat_64_struct> cond_stream_3_4("cond_stream_3_4");
#pragma HLS STREAM variable=cond_stream_3_4 depth=2

//
static hls::stream<mat_64_struct> cond_stream_4_0("cond_stream_4_0");
#pragma HLS STREAM variable=cond_stream_4_0 depth=2

static hls::stream<mat_64_struct> cond_stream_4_1("cond_stream_4_1");
#pragma HLS STREAM variable=cond_stream_4_1 depth=2

static hls::stream<mat_64_struct> cond_stream_4_2("cond_stream_4_2");
#pragma HLS STREAM variable=cond_stream_4_2 depth=2

static hls::stream<mat_64_struct> cond_stream_4_3("cond_stream_4_3");
#pragma HLS STREAM variable=cond_stream_4_3 depth=2

static hls::stream<mat_64_struct> cond_stream_4_4("cond_stream_4_4");
#pragma HLS STREAM variable=cond_stream_4_4 depth=2

static hls::stream<mat_64_struct> cond_stream_4_5("cond_stream_4_5");
#pragma HLS STREAM variable=cond_stream_4_5 depth=2
//

static hls::stream<sk_2k_vec> vec_mat_stream_2("vec_mat_stream_2");
#pragma HLS STREAM variable=vec_mat_stream_2 depth=7

static hls::stream<sk_2k_vec> vec_mat_stream_3("vec_mat_stream_3");
#pragma HLS STREAM variable=vec_mat_stream_3 depth=6

static hls::stream<sk_2k_vec> vec_mat_stream_4("vec_mat_stream_4");
#pragma HLS STREAM variable=vec_mat_stream_4 depth=7

#pragma HLS dataflow



/////LAYER 1//////
	preprocess_layer_1(cond_stream_1_0, cond_stream_1_1, cond_stream_1_2, cond_stream_1_3, cond_stream_1_4, cond_stream_1_5, vec_mat_stream_2, vec_mat_stream_3, vec_mat_stream_4, benes_stream_1, r_in, vec_bits);

	split_benes_stream(benes_stream_1_0_0, benes_stream_1_0_1, benes_stream_1);
	benes_layer(benes_stream_1_1_0, benes_stream_1_1_1, benes_stream_1_0_0, benes_stream_1_0_1, cond_stream_1_0, 0);
	benes_layer(benes_stream_1_2_0, benes_stream_1_2_1, benes_stream_1_1_0, benes_stream_1_1_1, cond_stream_1_1, 1);
	benes_layer(benes_stream_1_3_0, benes_stream_1_3_1, benes_stream_1_2_0, benes_stream_1_2_1, cond_stream_1_2, 2);
	benes_layer(benes_stream_1_4_0, benes_stream_1_4_1, benes_stream_1_3_0, benes_stream_1_3_1, cond_stream_1_3, 3);
	benes_layer(benes_stream_1_5_0, benes_stream_1_5_1, benes_stream_1_4_0, benes_stream_1_4_1, cond_stream_1_4, 4);
	benes_layer(benes_stream_1_6_0, benes_stream_1_6_1, benes_stream_1_5_0, benes_stream_1_5_1, cond_stream_1_5, 5);



/////LAYER 2//////
	preprocess_layer_2(cond_stream_2_0, cond_stream_2_1, cond_stream_2_2, cond_stream_2_3, cond_stream_2_4, cond_stream_2_5, benes_stream_2, vec_mat_stream_2, benes_stream_1_6_0, benes_stream_1_6_1);

	split_benes_stream(benes_stream_2_0_0, benes_stream_2_0_1, benes_stream_2);
	benes_layer(benes_stream_2_1_0, benes_stream_2_1_1, benes_stream_2_0_0, benes_stream_2_0_1, cond_stream_2_0, 0);
	benes_layer(benes_stream_2_2_0, benes_stream_2_2_1, benes_stream_2_1_0, benes_stream_2_1_1, cond_stream_2_1, 1);
	benes_layer(benes_stream_2_3_0, benes_stream_2_3_1, benes_stream_2_2_0, benes_stream_2_2_1, cond_stream_2_2, 2);
	benes_layer(benes_stream_2_4_0, benes_stream_2_4_1, benes_stream_2_3_0, benes_stream_2_3_1, cond_stream_2_3, 3);
	benes_layer(benes_stream_2_5_0, benes_stream_2_5_1, benes_stream_2_4_0, benes_stream_2_4_1, cond_stream_2_4, 4);
	benes_layer(benes_stream_2_6_0, benes_stream_2_6_1, benes_stream_2_5_0, benes_stream_2_5_1, cond_stream_2_5, 5);


	/////LAYER 3//////
	preprocess_layer_3(cond_stream_3_0, cond_stream_3_1, cond_stream_3_2, cond_stream_3_3, cond_stream_3_4, vec_mat_stream_3);

	split_benes_short_stream(benes_stream_3_0_0, benes_stream_3_0_1, benes_stream_2_6_0, benes_stream_2_6_1);
	benes_short_layer(benes_stream_3_1_0, benes_stream_3_1_1, benes_stream_3_0_0, benes_stream_3_0_1, cond_stream_3_4, 1, 16);
	benes_short_layer(benes_stream_3_2_0, benes_stream_3_2_1, benes_stream_3_1_0, benes_stream_3_1_1, cond_stream_3_3, 8, 2);
	benes_short_layer(benes_stream_3_3_0, benes_stream_3_3_1, benes_stream_3_2_0, benes_stream_3_2_1, cond_stream_3_2, 4, 4);
	benes_short_layer(benes_stream_3_4_0, benes_stream_3_4_1, benes_stream_3_3_0, benes_stream_3_3_1, cond_stream_3_1, 2, 8);
	benes_short_layer(benes_stream_3_5_0, benes_stream_3_5_1, benes_stream_3_4_0, benes_stream_3_4_1, cond_stream_3_0, 1, 16);


/////LAYER 4//////
	preprocess_layer_4(cond_stream_4_0, cond_stream_4_1, cond_stream_4_2, cond_stream_4_3, cond_stream_4_4, cond_stream_4_5, benes_stream_4, vec_mat_stream_4, benes_stream_3_5_0, benes_stream_3_5_1);

	split_benes_stream(benes_stream_4_0_0, benes_stream_4_0_1, benes_stream_4);
	benes_tr_layer(benes_stream_4_1_0, benes_stream_4_1_1, benes_stream_4_0_0, benes_stream_4_0_1, cond_stream_4_5, 5);
	benes_tr_layer(benes_stream_4_2_0, benes_stream_4_2_1, benes_stream_4_1_0, benes_stream_4_1_1, cond_stream_4_4, 4);
	benes_tr_layer(benes_stream_4_3_0, benes_stream_4_3_1, benes_stream_4_2_0, benes_stream_4_2_1, cond_stream_4_3, 3);
	benes_tr_layer(benes_stream_4_4_0, benes_stream_4_4_1, benes_stream_4_3_0, benes_stream_4_3_1, cond_stream_4_2, 2);
	benes_tr_layer(benes_stream_4_5_0, benes_stream_4_5_1, benes_stream_4_4_0, benes_stream_4_4_1, cond_stream_4_1, 1);
	benes_tr_layer(benes_stream_4_6_0, benes_stream_4_6_1, benes_stream_4_5_0, benes_stream_4_5_1, cond_stream_4_0, 0);


	postprocess_layer_4(r_out, benes_stream_4_6_0, benes_stream_4_6_1);
}


void benes2(hls::stream<vec> &r_out, hls::stream<vec> &r_in, hls::stream<sk_2k_vec> &vec_bits)
{


	int i;


	const unsigned char *cond_ptr;
	int inc, low;

	vec r1[64], r2[64], r3[64], r4[64];

//	vec r_mat[64];
	sk_2k_vec tmp;

	sk_2k_vec vec_mat[23];

	vec cond[6][64];
	vec cond_tr[6][64];
	vec cond_tr_2[6][64];


static hls::stream<vec> benes_stream_1("benes_stream_1");
#pragma HLS STREAM variable=benes_stream_1 depth=64

static hls::stream<vec> benes_stream_1_0_0("benes_stream_1_0_0");
#pragma HLS STREAM variable=benes_stream_1_0_0 depth=64

static hls::stream<vec> benes_stream_1_0_1("benes_stream_1_0_1");
#pragma HLS STREAM variable=benes_stream_1_0_1 depth=64

static hls::stream<vec> benes_stream_1_1_0("benes_stream_1_1_0");
#pragma HLS STREAM variable=benes_stream_1_1_0 depth=64

static hls::stream<vec> benes_stream_1_1_1("benes_stream_1_1_1");
#pragma HLS STREAM variable=benes_stream_1_1_1 depth=64

static hls::stream<vec> benes_stream_1_2_0("benes_stream_1_2_0");
#pragma HLS STREAM variable=benes_stream_1_2_0 depth=64

static hls::stream<vec> benes_stream_1_2_1("benes_stream_1_2_1");
#pragma HLS STREAM variable=benes_stream_1_2_1 depth=64

static hls::stream<vec> benes_stream_1_3_0("benes_stream_1_3_0");
#pragma HLS STREAM variable=benes_stream_1_3_0 depth=64

static hls::stream<vec> benes_stream_1_3_1("benes_stream_1_3_1");
#pragma HLS STREAM variable=benes_stream_1_3_1 depth=64

static hls::stream<vec> benes_stream_1_4_0("benes_stream_1_4_0");
#pragma HLS STREAM variable=benes_stream_1_4_0 depth=64

static hls::stream<vec> benes_stream_1_4_1("benes_stream_1_4_1");
#pragma HLS STREAM variable=benes_stream_1_4_1 depth=64

static hls::stream<vec> benes_stream_1_5_0("benes_stream_1_5_0");
#pragma HLS STREAM variable=benes_stream_1_5_0 depth=64

static hls::stream<vec> benes_stream_1_5_1("benes_stream_1_5_1");
#pragma HLS STREAM variable=benes_stream_1_5_1 depth=64

static hls::stream<vec> benes_stream_1_6_0("benes_stream_1_6_0");
#pragma HLS STREAM variable=benes_stream_1_6_0 depth=64

static hls::stream<vec> benes_stream_1_6_1("benes_stream_1_6_1");
#pragma HLS STREAM variable=benes_stream_1_6_1 depth=64
//

static hls::stream<vec> benes_stream_2("benes_stream_2");
#pragma HLS STREAM variable=benes_stream_2 depth=64

static hls::stream<vec> benes_stream_2_0_0("benes_stream_2_0_0");
#pragma HLS STREAM variable=benes_stream_2_0_0 depth=64

static hls::stream<vec> benes_stream_2_0_1("benes_stream_2_0_1");
#pragma HLS STREAM variable=benes_stream_2_0_1 depth=64

static hls::stream<vec> benes_stream_2_1_0("benes_stream_2_1_0");
#pragma HLS STREAM variable=benes_stream_2_1_0 depth=64

static hls::stream<vec> benes_stream_2_1_1("benes_stream_2_1_1");
#pragma HLS STREAM variable=benes_stream_2_1_1 depth=64

static hls::stream<vec> benes_stream_2_2_0("benes_stream_2_2_0");
#pragma HLS STREAM variable=benes_stream_2_2_0 depth=64

static hls::stream<vec> benes_stream_2_2_1("benes_stream_2_2_1");
#pragma HLS STREAM variable=benes_stream_2_2_1 depth=64

static hls::stream<vec> benes_stream_2_3_0("benes_stream_2_3_0");
#pragma HLS STREAM variable=benes_stream_2_3_0 depth=64

static hls::stream<vec> benes_stream_2_3_1("benes_stream_2_3_1");
#pragma HLS STREAM variable=benes_stream_2_3_1 depth=64

static hls::stream<vec> benes_stream_2_4_0("benes_stream_2_4_0");
#pragma HLS STREAM variable=benes_stream_2_4_0 depth=64

static hls::stream<vec> benes_stream_2_4_1("benes_stream_2_4_1");
#pragma HLS STREAM variable=benes_stream_2_4_1 depth=64

static hls::stream<vec> benes_stream_2_5_0("benes_stream_2_5_0");
#pragma HLS STREAM variable=benes_stream_2_5_0 depth=64

static hls::stream<vec> benes_stream_2_5_1("benes_stream_2_5_1");
#pragma HLS STREAM variable=benes_stream_2_5_1 depth=64

static hls::stream<vec> benes_stream_2_6_0("benes_stream_2_6_0");
#pragma HLS STREAM variable=benes_stream_2_6_0 depth=64

static hls::stream<vec> benes_stream_2_6_1("benes_stream_2_6_1");
#pragma HLS STREAM variable=benes_stream_2_6_1 depth=64
//

static hls::stream<vec> benes_stream_3("benes_stream_3");
#pragma HLS STREAM variable=benes_stream_3 depth=64

static hls::stream<vec> benes_stream_3_0_0("benes_stream_3_0_0");
#pragma HLS STREAM variable=benes_stream_3_0_0 depth=64

static hls::stream<vec> benes_stream_3_0_1("benes_stream_3_0_1");
#pragma HLS STREAM variable=benes_stream_3_0_1 depth=64

static hls::stream<vec> benes_stream_3_1_0("benes_stream_3_1_0");
#pragma HLS STREAM variable=benes_stream_3_1_0 depth=64

static hls::stream<vec> benes_stream_3_1_1("benes_stream_3_1_1");
#pragma HLS STREAM variable=benes_stream_3_1_1 depth=64

static hls::stream<vec> benes_stream_3_2_0("benes_stream_3_2_0");
#pragma HLS STREAM variable=benes_stream_3_2_0 depth=64

static hls::stream<vec> benes_stream_3_2_1("benes_stream_3_2_1");
#pragma HLS STREAM variable=benes_stream_3_2_1 depth=64

static hls::stream<vec> benes_stream_3_3_0("benes_stream_3_3_0");
#pragma HLS STREAM variable=benes_stream_3_3_0 depth=64

static hls::stream<vec> benes_stream_3_3_1("benes_stream_3_3_1");
#pragma HLS STREAM variable=benes_stream_3_3_1 depth=64

static hls::stream<vec> benes_stream_3_4_0("benes_stream_3_4_0");
#pragma HLS STREAM variable=benes_stream_3_4_0 depth=64

static hls::stream<vec> benes_stream_3_4_1("benes_stream_3_4_1");
#pragma HLS STREAM variable=benes_stream_3_4_1 depth=64

static hls::stream<vec> benes_stream_3_5_0("benes_stream_3_5_0");
#pragma HLS STREAM variable=benes_stream_3_5_0 depth=64

static hls::stream<vec> benes_stream_3_5_1("benes_stream_3_5_1");
#pragma HLS STREAM variable=benes_stream_3_5_1 depth=64

static hls::stream<vec> benes_stream_3_6_0("benes_stream_3_6_0");
#pragma HLS STREAM variable=benes_stream_3_6_0 depth=64

static hls::stream<vec> benes_stream_3_6_1("benes_stream_3_6_1");
#pragma HLS STREAM variable=benes_stream_3_6_1 depth=64
//////

static hls::stream<vec> benes_stream_4("benes_stream_4");
#pragma HLS STREAM variable=benes_stream_4 depth=64

static hls::stream<vec> benes_stream_4_0_0("benes_stream_4_0_0");
#pragma HLS STREAM variable=benes_stream_4_0_0 depth=64

static hls::stream<vec> benes_stream_4_0_1("benes_stream_4_0_1");
#pragma HLS STREAM variable=benes_stream_4_0_1 depth=64

static hls::stream<vec> benes_stream_4_1_0("benes_stream_4_1_0");
#pragma HLS STREAM variable=benes_stream_4_1_0 depth=64

static hls::stream<vec> benes_stream_4_1_1("benes_stream_4_1_1");
#pragma HLS STREAM variable=benes_stream_4_1_1 depth=64

static hls::stream<vec> benes_stream_4_2_0("benes_stream_4_2_0");
#pragma HLS STREAM variable=benes_stream_4_2_0 depth=64

static hls::stream<vec> benes_stream_4_2_1("benes_stream_4_2_1");
#pragma HLS STREAM variable=benes_stream_4_2_1 depth=64

static hls::stream<vec> benes_stream_4_3_0("benes_stream_4_3_0");
#pragma HLS STREAM variable=benes_stream_4_3_0 depth=64

static hls::stream<vec> benes_stream_4_3_1("benes_stream_4_3_1");
#pragma HLS STREAM variable=benes_stream_4_3_1 depth=64

static hls::stream<vec> benes_stream_4_4_0("benes_stream_4_4_0");
#pragma HLS STREAM variable=benes_stream_4_4_0 depth=64

static hls::stream<vec> benes_stream_4_4_1("benes_stream_4_4_1");
#pragma HLS STREAM variable=benes_stream_4_4_1 depth=64

static hls::stream<vec> benes_stream_4_5_0("benes_stream_4_5_0");
#pragma HLS STREAM variable=benes_stream_4_5_0 depth=64

static hls::stream<vec> benes_stream_4_5_1("benes_stream_4_5_1");
#pragma HLS STREAM variable=benes_stream_4_5_1 depth=64

static hls::stream<vec> benes_stream_4_6_0("benes_stream_4_6_0");
#pragma HLS STREAM variable=benes_stream_4_6_0 depth=64

static hls::stream<vec> benes_stream_4_6_1("benes_stream_4_6_1");
#pragma HLS STREAM variable=benes_stream_4_6_1 depth=64

//
static hls::stream<mat_64_struct> cond_stream_1_0("cond_stream_1_0");
#pragma HLS STREAM variable=cond_stream_1_0 depth=2

static hls::stream<mat_64_struct> cond_stream_1_1("cond_stream_1_1");
#pragma HLS STREAM variable=cond_stream_1_1 depth=2

static hls::stream<mat_64_struct> cond_stream_1_2("cond_stream_1_2");
#pragma HLS STREAM variable=cond_stream_1_2 depth=2

static hls::stream<mat_64_struct> cond_stream_1_3("cond_stream_1_3");
#pragma HLS STREAM variable=cond_stream_1_3 depth=2

static hls::stream<mat_64_struct> cond_stream_1_4("cond_stream_1_4");
#pragma HLS STREAM variable=cond_stream_1_4 depth=2

static hls::stream<mat_64_struct> cond_stream_1_5("cond_stream_1_5");
#pragma HLS STREAM variable=cond_stream_1_5 depth=2

//
static hls::stream<mat_64_struct> cond_stream_2_0("cond_stream_2_0");
#pragma HLS STREAM variable=cond_stream_2_0 depth=2

static hls::stream<mat_64_struct> cond_stream_2_1("cond_stream_2_1");
#pragma HLS STREAM variable=cond_stream_2_1 depth=2

static hls::stream<mat_64_struct> cond_stream_2_2("cond_stream_2_2");
#pragma HLS STREAM variable=cond_stream_2_2 depth=2

static hls::stream<mat_64_struct> cond_stream_2_3("cond_stream_2_3");
#pragma HLS STREAM variable=cond_stream_2_3 depth=2

static hls::stream<mat_64_struct> cond_stream_2_4("cond_stream_2_4");
#pragma HLS STREAM variable=cond_stream_2_4 depth=2

static hls::stream<mat_64_struct> cond_stream_2_5("cond_stream_2_5");
#pragma HLS STREAM variable=cond_stream_2_5 depth=2

//
static hls::stream<mat_64_struct> cond_stream_3_0("cond_stream_3_0");
#pragma HLS STREAM variable=cond_stream_3_0 depth=2

static hls::stream<mat_64_struct> cond_stream_3_1("cond_stream_3_1");
#pragma HLS STREAM variable=cond_stream_3_1 depth=2

static hls::stream<mat_64_struct> cond_stream_3_2("cond_stream_3_2");
#pragma HLS STREAM variable=cond_stream_3_2 depth=2

static hls::stream<mat_64_struct> cond_stream_3_3("cond_stream_3_3");
#pragma HLS STREAM variable=cond_stream_3_3 depth=2

static hls::stream<mat_64_struct> cond_stream_3_4("cond_stream_3_4");
#pragma HLS STREAM variable=cond_stream_3_4 depth=2

//
static hls::stream<mat_64_struct> cond_stream_4_0("cond_stream_4_0");
#pragma HLS STREAM variable=cond_stream_4_0 depth=2

static hls::stream<mat_64_struct> cond_stream_4_1("cond_stream_4_1");
#pragma HLS STREAM variable=cond_stream_4_1 depth=2

static hls::stream<mat_64_struct> cond_stream_4_2("cond_stream_4_2");
#pragma HLS STREAM variable=cond_stream_4_2 depth=2

static hls::stream<mat_64_struct> cond_stream_4_3("cond_stream_4_3");
#pragma HLS STREAM variable=cond_stream_4_3 depth=2

static hls::stream<mat_64_struct> cond_stream_4_4("cond_stream_4_4");
#pragma HLS STREAM variable=cond_stream_4_4 depth=2

static hls::stream<mat_64_struct> cond_stream_4_5("cond_stream_4_5");
#pragma HLS STREAM variable=cond_stream_4_5 depth=2
//

static hls::stream<sk_2k_vec> vec_mat_stream_2("vec_mat_stream_2");
#pragma HLS STREAM variable=vec_mat_stream_2 depth=6

static hls::stream<sk_2k_vec> vec_mat_stream_3("vec_mat_stream_3");
#pragma HLS STREAM variable=vec_mat_stream_3 depth=5

static hls::stream<sk_2k_vec> vec_mat_stream_4("vec_mat_stream_4");
#pragma HLS STREAM variable=vec_mat_stream_4 depth=6

#pragma HLS dataflow



/////LAYER 1//////
	preprocess_layer_1(cond_stream_1_0, cond_stream_1_1, cond_stream_1_2, cond_stream_1_3, cond_stream_1_4, cond_stream_1_5, vec_mat_stream_2, vec_mat_stream_3, vec_mat_stream_4, benes_stream_1, r_in, vec_bits);

	split_benes_stream(benes_stream_1_0_0, benes_stream_1_0_1, benes_stream_1);
	benes_layer(benes_stream_1_1_0, benes_stream_1_1_1, benes_stream_1_0_0, benes_stream_1_0_1, cond_stream_1_0, 0);
	benes_layer(benes_stream_1_2_0, benes_stream_1_2_1, benes_stream_1_1_0, benes_stream_1_1_1, cond_stream_1_1, 1);
	benes_layer(benes_stream_1_3_0, benes_stream_1_3_1, benes_stream_1_2_0, benes_stream_1_2_1, cond_stream_1_2, 2);
	benes_layer(benes_stream_1_4_0, benes_stream_1_4_1, benes_stream_1_3_0, benes_stream_1_3_1, cond_stream_1_3, 3);
	benes_layer(benes_stream_1_5_0, benes_stream_1_5_1, benes_stream_1_4_0, benes_stream_1_4_1, cond_stream_1_4, 4);
	benes_layer(benes_stream_1_6_0, benes_stream_1_6_1, benes_stream_1_5_0, benes_stream_1_5_1, cond_stream_1_5, 5);



/////LAYER 2//////
	preprocess_layer_2(cond_stream_2_0, cond_stream_2_1, cond_stream_2_2, cond_stream_2_3, cond_stream_2_4, cond_stream_2_5, benes_stream_2, vec_mat_stream_2, benes_stream_1_6_0, benes_stream_1_6_1);

	split_benes_stream(benes_stream_2_0_0, benes_stream_2_0_1, benes_stream_2);
	benes_layer(benes_stream_2_1_0, benes_stream_2_1_1, benes_stream_2_0_0, benes_stream_2_0_1, cond_stream_2_0, 0);
	benes_layer(benes_stream_2_2_0, benes_stream_2_2_1, benes_stream_2_1_0, benes_stream_2_1_1, cond_stream_2_1, 1);
	benes_layer(benes_stream_2_3_0, benes_stream_2_3_1, benes_stream_2_2_0, benes_stream_2_2_1, cond_stream_2_2, 2);
	benes_layer(benes_stream_2_4_0, benes_stream_2_4_1, benes_stream_2_3_0, benes_stream_2_3_1, cond_stream_2_3, 3);
	benes_layer(benes_stream_2_5_0, benes_stream_2_5_1, benes_stream_2_4_0, benes_stream_2_4_1, cond_stream_2_4, 4);
	benes_layer(benes_stream_2_6_0, benes_stream_2_6_1, benes_stream_2_5_0, benes_stream_2_5_1, cond_stream_2_5, 5);


	/////LAYER 3//////
	preprocess_layer_3(cond_stream_3_0, cond_stream_3_1, cond_stream_3_2, cond_stream_3_3, cond_stream_3_4, vec_mat_stream_3);

	split_benes_short_stream(benes_stream_3_0_0, benes_stream_3_0_1, benes_stream_2_6_0, benes_stream_2_6_1);
	benes_short_layer(benes_stream_3_1_0, benes_stream_3_1_1, benes_stream_3_0_0, benes_stream_3_0_1, cond_stream_3_4, 1, 16);
	benes_short_layer(benes_stream_3_2_0, benes_stream_3_2_1, benes_stream_3_1_0, benes_stream_3_1_1, cond_stream_3_3, 8, 2);
	benes_short_layer(benes_stream_3_3_0, benes_stream_3_3_1, benes_stream_3_2_0, benes_stream_3_2_1, cond_stream_3_2, 4, 4);
	benes_short_layer(benes_stream_3_4_0, benes_stream_3_4_1, benes_stream_3_3_0, benes_stream_3_3_1, cond_stream_3_1, 2, 8);
	benes_short_layer(benes_stream_3_5_0, benes_stream_3_5_1, benes_stream_3_4_0, benes_stream_3_4_1, cond_stream_3_0, 1, 16);


/////LAYER 4//////
	preprocess_layer_4(cond_stream_4_0, cond_stream_4_1, cond_stream_4_2, cond_stream_4_3, cond_stream_4_4, cond_stream_4_5, benes_stream_4, vec_mat_stream_4, benes_stream_3_5_0, benes_stream_3_5_1);

	split_benes_stream(benes_stream_4_0_0, benes_stream_4_0_1, benes_stream_4);
	benes_tr_layer(benes_stream_4_1_0, benes_stream_4_1_1, benes_stream_4_0_0, benes_stream_4_0_1, cond_stream_4_5, 5);
	benes_tr_layer(benes_stream_4_2_0, benes_stream_4_2_1, benes_stream_4_1_0, benes_stream_4_1_1, cond_stream_4_4, 4);
	benes_tr_layer(benes_stream_4_3_0, benes_stream_4_3_1, benes_stream_4_2_0, benes_stream_4_2_1, cond_stream_4_3, 3);
	benes_tr_layer(benes_stream_4_4_0, benes_stream_4_4_1, benes_stream_4_3_0, benes_stream_4_3_1, cond_stream_4_2, 2);
	benes_tr_layer(benes_stream_4_5_0, benes_stream_4_5_1, benes_stream_4_4_0, benes_stream_4_4_1, cond_stream_4_1, 1);
	benes_tr_layer(benes_stream_4_6_0, benes_stream_4_6_1, benes_stream_4_5_0, benes_stream_4_5_1, cond_stream_4_0, 0);


	postprocess_layer_4(r_out, benes_stream_4_6_0, benes_stream_4_6_1);
}



void load_irr_vec(hls::stream<sk_1k_vec> &out, sk_1k_vec * sk){

	out.write(*sk);
}

static void irr_load(vec * out, hls::stream<sk_1k_vec> &in)
{
	int i, j;
	vec irr[ SYS_T + 1 ];
	sk_1k_vec tmp;

#pragma HLS array partition variable=irr complete

	tmp = in.read();
	for (i = 0; i < SYS_T; i++){
	#pragma HLS pipeline
		irr[i].range(63, 12) = 0;
		irr[i].range(11, 0) = tmp.range(((i+1)*16)-5, (i*16));
	}

	irr[ SYS_T ] = 1;

	for (i = 0; i < GFBITS; i++)
		out[i] = 0;

	//TODO fix here! The following loop with ap_uints was giving wrong results
	for (i = SYS_T; i >= 0; i--)
	#pragma HLS pipeline
	for (j = 0; j < GFBITS; j++)
	{
//		out[j].bit(i) = irr[i].bit(j);
		out[j] <<= 1;
		out[j] |= (irr[i] >> j) & 1;
	}
}







static void pre_fft_scaling(hls::stream<vec> &irr_out_stream, hls::stream<sk_1k_vec> &sk)
{
	vec irr_mat[GFBITS];
	int i;

	irr_load(irr_mat, sk);
	for (i=0; i<GFBITS; i++){
		irr_out_stream.write(irr_mat[i]);
	}

}


void post_fft_scaling(hls::stream<mat_struct> &out_struct,  hls::stream<vec> &inv_stream, hls::stream<vec> &eval_stream, hls::stream<vec> &recv_stream)
{
	int i, j;

	vec irr_int[ GFBITS ];
	vec tmp[ GFBITS ];
	vec eval[64][GFBITS];
	vec tmp_recv;

	vec inv[64][GFBITS];
	vec out[64][GFBITS];

	mat_struct tmp_struct;

#pragma HLS array partition variable=tmp complete

	for (i = 0; i < 64; i++){
	#pragma HLS pipeline
		for(j=0; j<GFBITS; j++){
			eval[i][j] = eval_stream.read();
		}
		vec_sq(eval[i], eval[i]);
	}

	vec_copy(inv[0], eval[0]);

	for (i = 1; i < 64; i++)
	#pragma HLS pipeline
		vec_mul(inv[i], inv[i-1], eval[i]);

	vec_inv(tmp, inv[63]);

	for (i = 62; i >= 0; i--)
	#pragma HLS pipeline
	{
		vec_mul(inv[i+1], tmp, inv[i]);
		vec_mul(tmp, tmp, eval[i+1]);
	}

	vec_copy(inv[0], tmp);

	//

	for (i = 0; i < 64; i++){
	#pragma HLS pipeline
		tmp_recv = recv_stream.read();
		for (j = 0; j < GFBITS; j++){
			out[i][j] = inv[i][j] & tmp_recv;
			inv_stream.write(inv[i][j]);

			tmp_struct.mat[j] = out[i][j];

//			printf("out[%d][%d] = %lX\n", i, j, out[i][j].to_uint64());
//			printf("inv[%d][%d] = %lX\n", i, j, inv[i][j].to_uint64());
		}

		out_struct.write(tmp_struct);
	}

//	stream_out(out_struct, out);

}



static void scaling_stream(hls::stream<mat_struct> &recv_stream_out, hls::stream<vec> &inv_stream, hls::stream<sk_1k_vec> &sk_irr_vec_stream, hls::stream<vec> &recv_stream)
{


//	static hls::stream<vec> int_mat_0_stream("int_mat_0_stream");
//	#pragma HLS STREAM variable=int_mat_0_stream depth=12

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

	static hls::stream<vec> irr_int_stream("irr_int_stream");
	#pragma HLS STREAM variable=irr_int_stream depth=12

	static hls::stream<mat_struct> scaled_struct_stream_0("scaled_struct_stream_0");
	#pragma HLS STREAM variable=scaled_struct_stream_0 depth=64

	static hls::stream<vec> eval_in_stream("eval_in_stream");
	#pragma HLS STREAM variable=eval_in_stream depth=768

	#pragma HLS dataflow

	pre_fft_scaling(irr_int_stream, sk_irr_vec_stream);

//	radix_conversions_stream(int_mat_0_stream, irr_int_stream);
//	broadcast(scaled_struct_stream_0, int_mat_0_stream);
	radix_broadcast(scaled_struct_stream_0, irr_int_stream);

	split_struct_stream(struct_stream_0_0_0, struct_stream_0_0_1, scaled_struct_stream_0);
	fft_layer(struct_stream_0_1_0, struct_stream_0_1_1, struct_stream_0_0_0, struct_stream_0_0_1, 0, 32);
	fft_layer(struct_stream_0_2_0, struct_stream_0_2_1, struct_stream_0_1_0, struct_stream_0_1_1, 1, 16);
	fft_layer(struct_stream_0_3_0, struct_stream_0_3_1, struct_stream_0_2_0, struct_stream_0_2_1, 3, 8);
	fft_layer(struct_stream_0_4_0, struct_stream_0_4_1, struct_stream_0_3_0, struct_stream_0_3_1, 7, 4);
	fft_layer(struct_stream_0_5_0, struct_stream_0_5_1, struct_stream_0_4_0, struct_stream_0_4_1, 15, 2);
	fft_layer(struct_stream_0_6_0, struct_stream_0_6_1, struct_stream_0_5_0, struct_stream_0_5_1, 31, 1);
	write_output(eval_in_stream, struct_stream_0_6_0, struct_stream_0_6_1);

	post_fft_scaling(recv_stream_out, inv_stream, eval_in_stream, recv_stream);
}

//static void scaling(hls::stream<mat_struct> &out_struct, hls::stream<vec> &inv_stream, hls::stream<sk_1k_vec> &sk, hls::stream<vec> &recv)
//{
//	int i, j;
//
//	vec irr_int[ GFBITS ];
//	vec eval[64][ GFBITS ];
//	vec tmp[ GFBITS ];
//	vec tmp_recv;
//
//	vec inv[64][GFBITS];
//	vec out[64][GFBITS];
//
//	mat_struct tmp_struct;
//
//#pragma HLS array partition variable=eval dim=2 complete
//#pragma HLS array partition variable=tmp complete
//#pragma HLS array partition variable=irr_int complete
//
//
//	irr_load(irr_int, sk);
//
//	fft(eval, irr_int);
//
//	for (i = 0; i < 64; i++){
//	#pragma HLS unroll
//		vec_sq(eval[i], eval[i]);
//	}
//
//	vec_copy(inv[0], eval[0]);
//
//	for (i = 1; i < 64; i++)
//	#pragma HLS pipeline
//		vec_mul(inv[i], inv[i-1], eval[i]);
//
//	vec_inv(tmp, inv[63]);
//
//	for (i = 62; i >= 0; i--)
//	#pragma HLS pipeline
//	{
//		vec_mul(inv[i+1], tmp, inv[i]);
//		vec_mul(tmp, tmp, eval[i+1]);
//	}
//
//	vec_copy(inv[0], tmp);
//
//	//
//
//	for (i = 0; i < 64; i++){
//	#pragma HLS pipeline
//		tmp_recv = recv.read();
//		for (j = 0; j < GFBITS; j++){
//			out[i][j] = inv[i][j] & tmp_recv;
//			inv_stream.write(inv[i][j]);
//
//			tmp_struct.mat[j] = out[i][j];
//
////			printf("out[%d][%d] = %lX\n", i, j, out[i][j].to_uint64());
////			printf("inv[%d][%d] = %lX\n", i, j, inv[i][j].to_uint64());
//		}
//
//		out_struct.write(tmp_struct);
//	}
//
////	stream_out(out_struct, out);
//
//}

static void preprocess(hls::stream<vec> &recv_stream, hls::stream<unsigned char> &s)
{
	int i;
	unsigned char r[ 512 ];
	vec tmp;

	for (i = 0; i < SYND_BYTES; i++){
		r[i] = s.read();
//		printf("r[%d]=%lX\n", i, r[i]);
	}

	for (i = SYND_BYTES; i < 512; i++)
		r[i] = 0;

	for (i = 0; i < 64; i++){
		tmp = load8(r + i*8);
		recv_stream.write(tmp);
//		printf("recv[%d]=%lX\n", i, tmp.to_uint64());
	}
}

//TODO Have extra unneeded iterations here (512-436)
//TODO maybe merge loops?
static void postprocess(unsigned char *e, hls::stream<unsigned char> &e_stream_out, hls::stream<vec> &err_stream_in)
{
	int i, j;
	unsigned char error8[ (1 << GFBITS)/8 ];
	vec tmp;

	for (i = 0; i < 64; i++){
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

//static void scaling_inv(hls::stream<mat_struct> &out_struct, hls::stream<vec> &inv, hls::stream<vec> &recv)
//{
//	int i, j;
//	vec recv_tmp;
//
//	vec out[64][GFBITS];
//
//	for (i = 0; i < 64; i++){
//		recv_tmp = recv.read();
//		for (j = 0; j < GFBITS; j++){
//			out[i][j] = inv.read() & recv_tmp;
////			printf("out[%d][%d] = %lX\n", i, j, out[i][j].to_uint64());
//		}
//	}
//
//	stream_out(out_struct, out);
//}

void scaling_inv(hls::stream<mat_struct> &out_struct, hls::stream<vec> &inv, hls::stream<vec> &recv_stream)
{
	int i, j;
	vec recv_tmp;

	vec out[64][GFBITS];

	for (i = 0; i < 64; i++){
		recv_tmp = recv_stream.read();
		for (j = 0; j < GFBITS; j++){
			out[i][j] = inv.read() & recv_tmp;
//			printf("out[%d][%d] = %lX\n", i, j, out[i][j].to_uint64());
		}
	}

//	stream_out(out_struct, out);
	int k, l;

	mat_struct tmp_struct;

	for (k = 0; k < 64; k++)
	{
		for(l=0; l<12; l++){
			tmp_struct.mat[l] = out[k][l];
		}
		out_struct.write(tmp_struct);
	}


}


static int weight_check(hls::stream<unsigned char> &e_stream, hls::stream<vec> &error_stream)
{
	int i, j;
	uint16_t w0 = 0;
	uint16_t w1 = 0;
	uint16_t check;

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
	check >>= 15;

	return check;
}

void synd_cmp(hls::stream<bit> &check_out, hls::stream<vec> &s_stream_0, hls::stream<vec> &s_stream_1, hls::stream<vec> &s_cmp_stream_0, hls::stream<vec> &s_cmp_stream_1)
{
	int i, j;
	vec diff = 0;
	vec diff_0 = 0;
	vec diff_1 = 0;

	vec out_bit;

	for (j = 0; j < GFBITS; j++){
	#pragma HLS pipeline
		diff_0 |= (s_stream_0.read() ^ s_cmp_stream_0.read());
		diff_1 |= (s_stream_1.read() ^ s_cmp_stream_1.read());
	}

	diff = diff_0 | diff_1;

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

void load_sk(unsigned char *sk_mat, unsigned char *sk){

	int i;
	LOOP_SK_READ:
	for (i = 0; i < 6492; i++){
		sk_mat[i] = *(sk+i);
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
	vec error_mat[64];
	int i, j;

	for (i = 0; i < 64; i++)
	{
		error_mat[i] = vec_or_reduce(in_stream);
		error_mat[i] ^= allone;
		out_stream.write(error_mat[i]);
//		printf("error[%d] = %lX\n", i, error_mat[i].to_uint64());
	}
}

void priv_stream_duplicate(hls::stream<vec> &out_stream_0_0, hls::stream<vec> &out_stream_0_1, hls::stream<vec> &out_stream_1_0, hls::stream<vec> &out_stream_1_1, hls::stream<vec> &in_stream_0, hls::stream<vec> &in_stream_1){

	int i, j;
	vec tmp_0, tmp_1;

	for (j = 0; j < GFBITS; j++){
		tmp_0 = in_stream_0.read();
		tmp_1 = in_stream_1.read();
		out_stream_0_0.write(tmp_0);
		out_stream_1_0.write(tmp_0);
		out_stream_0_1.write(tmp_1);
		out_stream_1_1.write(tmp_1);
	}

}

void error_stream_duplicate(hls::stream<vec> &out_stream_0, hls::stream<vec> &out_stream_1, hls::stream<vec> &in_stream){

	int i;
	vec tmp;

	for(i=0; i<64; i++){
		tmp = in_stream.read();
		out_stream_0.write(tmp);
		out_stream_1.write(tmp);
	}
}




void weight_check(hls::stream<bit> &check_weight, hls::stream<unsigned char> &e_stream, hls::stream<vec> &error_stream)
{
	int i, j;

	ap_uint<12> w0 = 0;
	ap_uint<12> w1 = 0;
	ap_uint<12> check;
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
	check_bit = check.bit(11);

	check_weight.write(check_bit);

}



void check_value_compute(int *check, hls::stream<bit> &check_synd, hls::stream<bit> &check_weight){

	bit check_value=0;
	bit one=1;

	check_value =  one - (check_synd.read() & check_weight.read());
	*check = (int)check_value;

}


void decryption_kernel_active(int *check, unsigned char *e, sk_1k_vec *sk_irr, sk_2k_vec *sk_rev, sk_2k_vec *sk_fwd, unsigned char *s){

#pragma HLS INTERFACE m_axi     port=e       offset=slave depth=436 bundle=gmem
#pragma HLS INTERFACE m_axi     port=sk_irr  offset=slave depth=1 bundle=gmem2
#pragma HLS INTERFACE m_axi     port=sk_rev  offset=slave depth=23 bundle=gmem3
#pragma HLS INTERFACE m_axi     port=sk_fwd  offset=slave depth=23 bundle=gmem4
#pragma HLS INTERFACE m_axi     port=s       offset=slave depth=96 bundle=gmem5
#pragma HLS INTERFACE m_axi     port=check   offset=slave depth=1 bundle=gmem6
#pragma HLS INTERFACE s_axilite port=e       bundle=control
#pragma HLS INTERFACE s_axilite port=sk_irr  bundle=control
#pragma HLS INTERFACE s_axilite port=sk_rev  bundle=control
#pragma HLS INTERFACE s_axilite port=sk_fwd  bundle=control
#pragma HLS INTERFACE s_axilite port=s       bundle=control
#pragma HLS INTERFACE s_axilite port=check   bundle=control
#pragma HLS INTERFACE s_axilite port=return  bundle=control



	static hls::stream<unsigned char> s_mat_stream("s_mat_stream");
	#pragma HLS STREAM variable=s_mat_stream depth=96

	static hls::stream<sk_2k_vec> sk_vec_rev_mat_stream("sk_vec_rev_mat_stream");
	#pragma HLS STREAM variable=sk_vec_rev_mat_stream depth=23

	static hls::stream<sk_2k_vec> sk_vec_fwd_mat_stream("sk_vec_fwd_mat_stream");
	#pragma HLS STREAM variable=sk_vec_fwd_mat_stream depth=23

	static hls::stream<sk_1k_vec> sk_irr_vec_stream("sk_irr_vec_stream");
	#pragma HLS STREAM variable=sk_irr_vec_stream depth=2

	static hls::stream<vec> recv_stream_in("recv_stream_in");
	#pragma HLS STREAM variable=recv_stream_in depth=64

	static hls::stream<vec> recv_stream_out("recv_stream_out");
	#pragma HLS STREAM variable=recv_stream_out depth=64

	static hls::stream<vec> inv_stream("inv_stream");
	#pragma HLS STREAM variable=inv_stream depth=768


	static hls::stream<mat_struct> scaled_struct_stream_2("scaled_struct_stream_2");
	#pragma HLS STREAM variable=scaled_struct_stream_2 depth=64

	static hls::stream<mat_struct> scaled_struct_stream_3("scaled_struct_stream_3");
	#pragma HLS STREAM variable=scaled_struct_stream_3 depth=64

	static hls::stream<vec> s_priv_stream_0("s_priv_stream_0");
	#pragma HLS STREAM variable=s_priv_stream_0 depth=12

	static hls::stream<vec> s_priv_stream_1("s_priv_stream_1");
	#pragma HLS STREAM variable=s_priv_stream_1 depth=12

	static hls::stream<vec> s_priv_stream_0_0("s_priv_stream_0_0");
	#pragma HLS STREAM variable=s_priv_stream_0_0 depth=12

	static hls::stream<vec> s_priv_stream_0_1("s_priv_stream_0_1");
	#pragma HLS STREAM variable=s_priv_stream_0_1 depth=12

	static hls::stream<vec> s_priv_stream_1_0("s_priv_stream_1_0");
	#pragma HLS STREAM variable=s_priv_stream_1_0 depth=12

	static hls::stream<vec> s_priv_stream_1_1("s_priv_stream_1_1");
	#pragma HLS STREAM variable=s_priv_stream_1_1 depth=12

	static hls::stream<vec> s_priv_cmp_stream_0("s_priv_cmp_stream_0");
	#pragma HLS STREAM variable=s_priv_cmp_stream_0 depth=12

	static hls::stream<vec> s_priv_cmp_stream_1("s_priv_cmp_stream_1");
	#pragma HLS STREAM variable=s_priv_cmp_stream_1 depth=12

	static hls::stream<vec> locator_stream("locator_stream");
	#pragma HLS STREAM variable=locator_stream depth=12

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

///////

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


		static hls::stream<vec> int_stream_comb("int_stream_comb");
		#pragma HLS STREAM variable=int_stream_comb depth=24

		static hls::stream<bit> check_synd_stream("check_synd_stream");
		#pragma HLS STREAM variable=check_synd_stream depth=1

		static hls::stream<bit> check_weight_stream("check_weight_stream");
		#pragma HLS STREAM variable=check_weight_stream depth=1

		static hls::stream<mat_struct> scaled_struct_stream_1("scaled_struct_stream_1");
		#pragma HLS STREAM variable=scaled_struct_stream_1 depth=64

		static hls::stream<vec> recv_stream_benes_in("recv_stream_benes_in");
		#pragma HLS STREAM variable=recv_stream_benes_in depth=64



	#pragma HLS dataflow

	load_s(s_mat_stream, s);

	load_sk_vec(sk_vec_rev_mat_stream, sk_rev);

	load_sk_vec(sk_vec_fwd_mat_stream, sk_fwd);

	load_irr_vec(sk_irr_vec_stream, sk_irr);

	preprocess(recv_stream_in, s_mat_stream);

	benes1(recv_stream_benes_in, recv_stream_in, sk_vec_rev_mat_stream);

///////////
//
//static hls::stream<vec> cond_tr_stream("cond_tr_stream");
//#pragma HLS STREAM variable=cond_tr_stream depth=2
//
//static hls::stream<vec> benes_stream_1("benes_stream_1");
//#pragma HLS STREAM variable=benes_stream_1 depth=64
//
//static hls::stream<vec> benes_stream_1_0_0("benes_stream_1_0_0");
//#pragma HLS STREAM variable=benes_stream_1_0_0 depth=64
//
//static hls::stream<vec> benes_stream_1_0_1("benes_stream_1_0_1");
//#pragma HLS STREAM variable=benes_stream_1_0_1 depth=64
//
//static hls::stream<vec> benes_stream_1_1_0("benes_stream_1_1_0");
//#pragma HLS STREAM variable=benes_stream_1_1_0 depth=64
//
//static hls::stream<vec> benes_stream_1_1_1("benes_stream_1_1_1");
//#pragma HLS STREAM variable=benes_stream_1_1_1 depth=64
//
//static hls::stream<vec> benes_stream_1_2_0("benes_stream_1_2_0");
//#pragma HLS STREAM variable=benes_stream_1_2_0 depth=64
//
//static hls::stream<vec> benes_stream_1_2_1("benes_stream_1_2_1");
//#pragma HLS STREAM variable=benes_stream_1_2_1 depth=64
//
//static hls::stream<vec> benes_stream_1_3_0("benes_stream_1_3_0");
//#pragma HLS STREAM variable=benes_stream_1_3_0 depth=64
//
//static hls::stream<vec> benes_stream_1_3_1("benes_stream_1_3_1");
//#pragma HLS STREAM variable=benes_stream_1_3_1 depth=64
//
//static hls::stream<vec> benes_stream_1_4_0("benes_stream_1_4_0");
//#pragma HLS STREAM variable=benes_stream_1_4_0 depth=64
//
//static hls::stream<vec> benes_stream_1_4_1("benes_stream_1_4_1");
//#pragma HLS STREAM variable=benes_stream_1_4_1 depth=64
//
//static hls::stream<vec> benes_stream_1_5_0("benes_stream_1_5_0");
//#pragma HLS STREAM variable=benes_stream_1_5_0 depth=64
//
//static hls::stream<vec> benes_stream_1_5_1("benes_stream_1_5_1");
//#pragma HLS STREAM variable=benes_stream_1_5_1 depth=64
//
//static hls::stream<vec> benes_stream_1_6_0("benes_stream_1_6_0");
//#pragma HLS STREAM variable=benes_stream_1_6_0 depth=64
//
//static hls::stream<vec> benes_stream_1_6_1("benes_stream_1_6_1");
//#pragma HLS STREAM variable=benes_stream_1_6_1 depth=64
////
//
//static hls::stream<vec> benes_stream_2("benes_stream_2");
//#pragma HLS STREAM variable=benes_stream_2 depth=64
//
//static hls::stream<vec> benes_stream_2_0_0("benes_stream_2_0_0");
//#pragma HLS STREAM variable=benes_stream_2_0_0 depth=64
//
//static hls::stream<vec> benes_stream_2_0_1("benes_stream_2_0_1");
//#pragma HLS STREAM variable=benes_stream_2_0_1 depth=64
//
//static hls::stream<vec> benes_stream_2_1_0("benes_stream_2_1_0");
//#pragma HLS STREAM variable=benes_stream_2_1_0 depth=64
//
//static hls::stream<vec> benes_stream_2_1_1("benes_stream_2_1_1");
//#pragma HLS STREAM variable=benes_stream_2_1_1 depth=64
//
//static hls::stream<vec> benes_stream_2_2_0("benes_stream_2_2_0");
//#pragma HLS STREAM variable=benes_stream_2_2_0 depth=64
//
//static hls::stream<vec> benes_stream_2_2_1("benes_stream_2_2_1");
//#pragma HLS STREAM variable=benes_stream_2_2_1 depth=64
//
//static hls::stream<vec> benes_stream_2_3_0("benes_stream_2_3_0");
//#pragma HLS STREAM variable=benes_stream_2_3_0 depth=64
//
//static hls::stream<vec> benes_stream_2_3_1("benes_stream_2_3_1");
//#pragma HLS STREAM variable=benes_stream_2_3_1 depth=64
//
//static hls::stream<vec> benes_stream_2_4_0("benes_stream_2_4_0");
//#pragma HLS STREAM variable=benes_stream_2_4_0 depth=64
//
//static hls::stream<vec> benes_stream_2_4_1("benes_stream_2_4_1");
//#pragma HLS STREAM variable=benes_stream_2_4_1 depth=64
//
//static hls::stream<vec> benes_stream_2_5_0("benes_stream_2_5_0");
//#pragma HLS STREAM variable=benes_stream_2_5_0 depth=64
//
//static hls::stream<vec> benes_stream_2_5_1("benes_stream_2_5_1");
//#pragma HLS STREAM variable=benes_stream_2_5_1 depth=64
//
//static hls::stream<vec> benes_stream_2_6_0("benes_stream_2_6_0");
//#pragma HLS STREAM variable=benes_stream_2_6_0 depth=64
//
//static hls::stream<vec> benes_stream_2_6_1("benes_stream_2_6_1");
//#pragma HLS STREAM variable=benes_stream_2_6_1 depth=64
////
//
//static hls::stream<vec> benes_stream_3("benes_stream_3");
//#pragma HLS STREAM variable=benes_stream_3 depth=64
//
//static hls::stream<vec> benes_stream_3_0_0("benes_stream_3_0_0");
//#pragma HLS STREAM variable=benes_stream_3_0_0 depth=64
//
//static hls::stream<vec> benes_stream_3_0_1("benes_stream_3_0_1");
//#pragma HLS STREAM variable=benes_stream_3_0_1 depth=64
//
//static hls::stream<vec> benes_stream_3_1_0("benes_stream_3_1_0");
//#pragma HLS STREAM variable=benes_stream_3_1_0 depth=64
//
//static hls::stream<vec> benes_stream_3_1_1("benes_stream_3_1_1");
//#pragma HLS STREAM variable=benes_stream_3_1_1 depth=64
//
//static hls::stream<vec> benes_stream_3_2_0("benes_stream_3_2_0");
//#pragma HLS STREAM variable=benes_stream_3_2_0 depth=64
//
//static hls::stream<vec> benes_stream_3_2_1("benes_stream_3_2_1");
//#pragma HLS STREAM variable=benes_stream_3_2_1 depth=64
//
//static hls::stream<vec> benes_stream_3_3_0("benes_stream_3_3_0");
//#pragma HLS STREAM variable=benes_stream_3_3_0 depth=64
//
//static hls::stream<vec> benes_stream_3_3_1("benes_stream_3_3_1");
//#pragma HLS STREAM variable=benes_stream_3_3_1 depth=64
//
//static hls::stream<vec> benes_stream_3_4_0("benes_stream_3_4_0");
//#pragma HLS STREAM variable=benes_stream_3_4_0 depth=64
//
//static hls::stream<vec> benes_stream_3_4_1("benes_stream_3_4_1");
//#pragma HLS STREAM variable=benes_stream_3_4_1 depth=64
//
//static hls::stream<vec> benes_stream_3_5_0("benes_stream_3_5_0");
//#pragma HLS STREAM variable=benes_stream_3_5_0 depth=64
//
//static hls::stream<vec> benes_stream_3_5_1("benes_stream_3_5_1");
//#pragma HLS STREAM variable=benes_stream_3_5_1 depth=64
//
//static hls::stream<vec> benes_stream_3_6_0("benes_stream_3_6_0");
//#pragma HLS STREAM variable=benes_stream_3_6_0 depth=64
//
//static hls::stream<vec> benes_stream_3_6_1("benes_stream_3_6_1");
//#pragma HLS STREAM variable=benes_stream_3_6_1 depth=64
////////
//
//static hls::stream<vec> benes_stream_4("benes_stream_4");
//#pragma HLS STREAM variable=benes_stream_4 depth=64
//
//static hls::stream<vec> benes_stream_4_0_0("benes_stream_4_0_0");
//#pragma HLS STREAM variable=benes_stream_4_0_0 depth=64
//
//static hls::stream<vec> benes_stream_4_0_1("benes_stream_4_0_1");
//#pragma HLS STREAM variable=benes_stream_4_0_1 depth=64
//
//static hls::stream<vec> benes_stream_4_1_0("benes_stream_4_1_0");
//#pragma HLS STREAM variable=benes_stream_4_1_0 depth=64
//
//static hls::stream<vec> benes_stream_4_1_1("benes_stream_4_1_1");
//#pragma HLS STREAM variable=benes_stream_4_1_1 depth=64
//
//static hls::stream<vec> benes_stream_4_2_0("benes_stream_4_2_0");
//#pragma HLS STREAM variable=benes_stream_4_2_0 depth=64
//
//static hls::stream<vec> benes_stream_4_2_1("benes_stream_4_2_1");
//#pragma HLS STREAM variable=benes_stream_4_2_1 depth=64
//
//static hls::stream<vec> benes_stream_4_3_0("benes_stream_4_3_0");
//#pragma HLS STREAM variable=benes_stream_4_3_0 depth=64
//
//static hls::stream<vec> benes_stream_4_3_1("benes_stream_4_3_1");
//#pragma HLS STREAM variable=benes_stream_4_3_1 depth=64
//
//static hls::stream<vec> benes_stream_4_4_0("benes_stream_4_4_0");
//#pragma HLS STREAM variable=benes_stream_4_4_0 depth=64
//
//static hls::stream<vec> benes_stream_4_4_1("benes_stream_4_4_1");
//#pragma HLS STREAM variable=benes_stream_4_4_1 depth=64
//
//static hls::stream<vec> benes_stream_4_5_0("benes_stream_4_5_0");
//#pragma HLS STREAM variable=benes_stream_4_5_0 depth=64
//
//static hls::stream<vec> benes_stream_4_5_1("benes_stream_4_5_1");
//#pragma HLS STREAM variable=benes_stream_4_5_1 depth=64
//
//static hls::stream<vec> benes_stream_4_6_0("benes_stream_4_6_0");
//#pragma HLS STREAM variable=benes_stream_4_6_0 depth=64
//
//static hls::stream<vec> benes_stream_4_6_1("benes_stream_4_6_1");
//#pragma HLS STREAM variable=benes_stream_4_6_1 depth=64
//
////
//static hls::stream<mat_64_struct> cond_stream_1_0("cond_stream_1_0");
//#pragma HLS STREAM variable=cond_stream_1_0 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_1_1("cond_stream_1_1");
//#pragma HLS STREAM variable=cond_stream_1_1 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_1_2("cond_stream_1_2");
//#pragma HLS STREAM variable=cond_stream_1_2 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_1_3("cond_stream_1_3");
//#pragma HLS STREAM variable=cond_stream_1_3 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_1_4("cond_stream_1_4");
//#pragma HLS STREAM variable=cond_stream_1_4 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_1_5("cond_stream_1_5");
//#pragma HLS STREAM variable=cond_stream_1_5 depth=2
//
////
//static hls::stream<mat_64_struct> cond_stream_2_0("cond_stream_2_0");
//#pragma HLS STREAM variable=cond_stream_2_0 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_2_1("cond_stream_2_1");
//#pragma HLS STREAM variable=cond_stream_2_1 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_2_2("cond_stream_2_2");
//#pragma HLS STREAM variable=cond_stream_2_2 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_2_3("cond_stream_2_3");
//#pragma HLS STREAM variable=cond_stream_2_3 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_2_4("cond_stream_2_4");
//#pragma HLS STREAM variable=cond_stream_2_4 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_2_5("cond_stream_2_5");
//#pragma HLS STREAM variable=cond_stream_2_5 depth=2
//
////
//static hls::stream<mat_64_struct> cond_stream_3_0("cond_stream_3_0");
//#pragma HLS STREAM variable=cond_stream_3_0 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_3_1("cond_stream_3_1");
//#pragma HLS STREAM variable=cond_stream_3_1 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_3_2("cond_stream_3_2");
//#pragma HLS STREAM variable=cond_stream_3_2 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_3_3("cond_stream_3_3");
//#pragma HLS STREAM variable=cond_stream_3_3 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_3_4("cond_stream_3_4");
//#pragma HLS STREAM variable=cond_stream_3_4 depth=2
//
////
//static hls::stream<mat_64_struct> cond_stream_4_0("cond_stream_4_0");
//#pragma HLS STREAM variable=cond_stream_4_0 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_4_1("cond_stream_4_1");
//#pragma HLS STREAM variable=cond_stream_4_1 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_4_2("cond_stream_4_2");
//#pragma HLS STREAM variable=cond_stream_4_2 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_4_3("cond_stream_4_3");
//#pragma HLS STREAM variable=cond_stream_4_3 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_4_4("cond_stream_4_4");
//#pragma HLS STREAM variable=cond_stream_4_4 depth=2
//
//static hls::stream<mat_64_struct> cond_stream_4_5("cond_stream_4_5");
//#pragma HLS STREAM variable=cond_stream_4_5 depth=2
////
//
//static hls::stream<sk_2k_vec> vec_mat_stream_2("vec_mat_stream_2");
//#pragma HLS STREAM variable=vec_mat_stream_2 depth=7
//
//static hls::stream<sk_2k_vec> vec_mat_stream_3("vec_mat_stream_3");
//#pragma HLS STREAM variable=vec_mat_stream_3 depth=6
//
//static hls::stream<sk_2k_vec> vec_mat_stream_4("vec_mat_stream_4");
//#pragma HLS STREAM variable=vec_mat_stream_4 depth=7
//
//
//
//
///////LAYER 1//////
//	preprocess_layer_1(cond_stream_1_0, cond_stream_1_1, cond_stream_1_2, cond_stream_1_3, cond_stream_1_4, cond_stream_1_5, vec_mat_stream_2, vec_mat_stream_3, vec_mat_stream_4, benes_stream_1, recv_stream_in, sk_vec_rev_mat_stream);
//
//	split_benes_stream(benes_stream_1_0_0, benes_stream_1_0_1, benes_stream_1);
//	benes_layer(benes_stream_1_1_0, benes_stream_1_1_1, benes_stream_1_0_0, benes_stream_1_0_1, cond_stream_1_0, 0);
//	benes_layer(benes_stream_1_2_0, benes_stream_1_2_1, benes_stream_1_1_0, benes_stream_1_1_1, cond_stream_1_1, 1);
//	benes_layer(benes_stream_1_3_0, benes_stream_1_3_1, benes_stream_1_2_0, benes_stream_1_2_1, cond_stream_1_2, 2);
//	benes_layer(benes_stream_1_4_0, benes_stream_1_4_1, benes_stream_1_3_0, benes_stream_1_3_1, cond_stream_1_3, 3);
//	benes_layer(benes_stream_1_5_0, benes_stream_1_5_1, benes_stream_1_4_0, benes_stream_1_4_1, cond_stream_1_4, 4);
//	benes_layer(benes_stream_1_6_0, benes_stream_1_6_1, benes_stream_1_5_0, benes_stream_1_5_1, cond_stream_1_5, 5);
//
//
//
///////LAYER 2//////
//	preprocess_layer_2(cond_stream_2_0, cond_stream_2_1, cond_stream_2_2, cond_stream_2_3, cond_stream_2_4, cond_stream_2_5, benes_stream_2, vec_mat_stream_2, benes_stream_1_6_0, benes_stream_1_6_1);
//
//	split_benes_stream(benes_stream_2_0_0, benes_stream_2_0_1, benes_stream_2);
//	benes_layer(benes_stream_2_1_0, benes_stream_2_1_1, benes_stream_2_0_0, benes_stream_2_0_1, cond_stream_2_0, 0);
//	benes_layer(benes_stream_2_2_0, benes_stream_2_2_1, benes_stream_2_1_0, benes_stream_2_1_1, cond_stream_2_1, 1);
//	benes_layer(benes_stream_2_3_0, benes_stream_2_3_1, benes_stream_2_2_0, benes_stream_2_2_1, cond_stream_2_2, 2);
//	benes_layer(benes_stream_2_4_0, benes_stream_2_4_1, benes_stream_2_3_0, benes_stream_2_3_1, cond_stream_2_3, 3);
//	benes_layer(benes_stream_2_5_0, benes_stream_2_5_1, benes_stream_2_4_0, benes_stream_2_4_1, cond_stream_2_4, 4);
//	benes_layer(benes_stream_2_6_0, benes_stream_2_6_1, benes_stream_2_5_0, benes_stream_2_5_1, cond_stream_2_5, 5);
//
//
//	/////LAYER 3//////
//	preprocess_layer_3(cond_stream_3_0, cond_stream_3_1, cond_stream_3_2, cond_stream_3_3, cond_stream_3_4, vec_mat_stream_3);
//
//	split_benes_short_stream(benes_stream_3_0_0, benes_stream_3_0_1, benes_stream_2_6_0, benes_stream_2_6_1);
//	benes_short_layer(benes_stream_3_1_0, benes_stream_3_1_1, benes_stream_3_0_0, benes_stream_3_0_1, cond_stream_3_4, 1, 16);
//	benes_short_layer(benes_stream_3_2_0, benes_stream_3_2_1, benes_stream_3_1_0, benes_stream_3_1_1, cond_stream_3_3, 8, 2);
//	benes_short_layer(benes_stream_3_3_0, benes_stream_3_3_1, benes_stream_3_2_0, benes_stream_3_2_1, cond_stream_3_2, 4, 4);
//	benes_short_layer(benes_stream_3_4_0, benes_stream_3_4_1, benes_stream_3_3_0, benes_stream_3_3_1, cond_stream_3_1, 2, 8);
//	benes_short_layer(benes_stream_3_5_0, benes_stream_3_5_1, benes_stream_3_4_0, benes_stream_3_4_1, cond_stream_3_0, 1, 16);
//
//
///////LAYER 4//////
//	preprocess_layer_4(cond_stream_4_0, cond_stream_4_1, cond_stream_4_2, cond_stream_4_3, cond_stream_4_4, cond_stream_4_5, benes_stream_4, vec_mat_stream_4, benes_stream_3_5_0, benes_stream_3_5_1);
//
//	split_benes_stream(benes_stream_4_0_0, benes_stream_4_0_1, benes_stream_4);
//	benes_tr_layer(benes_stream_4_1_0, benes_stream_4_1_1, benes_stream_4_0_0, benes_stream_4_0_1, cond_stream_4_5, 5);
//	benes_tr_layer(benes_stream_4_2_0, benes_stream_4_2_1, benes_stream_4_1_0, benes_stream_4_1_1, cond_stream_4_4, 4);
//	benes_tr_layer(benes_stream_4_3_0, benes_stream_4_3_1, benes_stream_4_2_0, benes_stream_4_2_1, cond_stream_4_3, 3);
//	benes_tr_layer(benes_stream_4_4_0, benes_stream_4_4_1, benes_stream_4_3_0, benes_stream_4_3_1, cond_stream_4_2, 2);
//	benes_tr_layer(benes_stream_4_5_0, benes_stream_4_5_1, benes_stream_4_4_0, benes_stream_4_4_1, cond_stream_4_1, 1);
//	benes_tr_layer(benes_stream_4_6_0, benes_stream_4_6_1, benes_stream_4_5_0, benes_stream_4_5_1, cond_stream_4_0, 0);
//
//
//	postprocess_layer_4(recv_stream_benes_in, benes_stream_4_6_0, benes_stream_4_6_1);
//
//
//
//
//
//
//
///////////



//	scaling_stream(scaled_struct_stream_1, inv_stream, sk_irr_vec_stream, recv_stream_benes_in);

	///////////////////////////////////////////////
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

	static hls::stream<vec> irr_int_stream("irr_int_stream");
	#pragma HLS STREAM variable=irr_int_stream depth=12

	static hls::stream<mat_struct> scaled_struct_stream_0("scaled_struct_stream_0");
	#pragma HLS STREAM variable=scaled_struct_stream_0 depth=64

	static hls::stream<vec> eval_in_stream("eval_in_stream");
	#pragma HLS STREAM variable=eval_in_stream depth=768


	pre_fft_scaling(irr_int_stream, sk_irr_vec_stream);

	radix_broadcast(scaled_struct_stream_0, irr_int_stream);

	split_struct_stream(struct_stream_0_0_0, struct_stream_0_0_1, scaled_struct_stream_0);
	fft_layer(struct_stream_0_1_0, struct_stream_0_1_1, struct_stream_0_0_0, struct_stream_0_0_1, 0, 32);
	fft_layer(struct_stream_0_2_0, struct_stream_0_2_1, struct_stream_0_1_0, struct_stream_0_1_1, 1, 16);
	fft_layer(struct_stream_0_3_0, struct_stream_0_3_1, struct_stream_0_2_0, struct_stream_0_2_1, 3, 8);
	fft_layer(struct_stream_0_4_0, struct_stream_0_4_1, struct_stream_0_3_0, struct_stream_0_3_1, 7, 4);
	fft_layer(struct_stream_0_5_0, struct_stream_0_5_1, struct_stream_0_4_0, struct_stream_0_4_1, 15, 2);
	fft_layer(struct_stream_0_6_0, struct_stream_0_6_1, struct_stream_0_5_0, struct_stream_0_5_1, 31, 1);
	write_output(eval_in_stream, struct_stream_0_6_0, struct_stream_0_6_1);

	post_fft_scaling(scaled_struct_stream_1, inv_stream, eval_in_stream, recv_stream_benes_in);
	///////////////////////////////////////////////

	fft_tr_stream1(s_priv_stream_0, s_priv_stream_1, scaled_struct_stream_1);

	priv_stream_duplicate(s_priv_stream_0_0, s_priv_stream_0_1, s_priv_stream_1_0, s_priv_stream_1_1, s_priv_stream_0, s_priv_stream_1);

	bm(locator_stream, s_priv_stream_0_0, s_priv_stream_0_1);

	fft_stream_alt(eval_stream, locator_stream);

	error_compute(error_stream, eval_stream);

	error_stream_duplicate(error_stream_0, error_stream_1, error_stream);

	scaling_inv(scaled_struct_stream_3, inv_stream, error_stream_0);

	fft_tr_stream2(s_priv_cmp_stream_0, s_priv_cmp_stream_1, scaled_struct_stream_3);

	synd_cmp(check_synd_stream, s_priv_stream_1_0, s_priv_stream_1_1, s_priv_cmp_stream_0, s_priv_cmp_stream_1);

//	benes2(error_stream_benes_out, error_stream_1, sk_vec_fwd_mat_stream);

	/////////

	static hls::stream<vec> benes_stream_5("benes_stream_5");
	#pragma HLS STREAM variable=benes_stream_5 depth=64

	static hls::stream<vec> benes_stream_5_0_0("benes_stream_5_0_0");
	#pragma HLS STREAM variable=benes_stream_5_0_0 depth=64

	static hls::stream<vec> benes_stream_5_0_1("benes_stream_5_0_1");
	#pragma HLS STREAM variable=benes_stream_5_0_1 depth=64

	static hls::stream<vec> benes_stream_5_1_0("benes_stream_5_1_0");
	#pragma HLS STREAM variable=benes_stream_5_1_0 depth=64

	static hls::stream<vec> benes_stream_5_1_1("benes_stream_5_1_1");
	#pragma HLS STREAM variable=benes_stream_5_1_1 depth=64

	static hls::stream<vec> benes_stream_5_2_0("benes_stream_5_2_0");
	#pragma HLS STREAM variable=benes_stream_5_2_0 depth=64

	static hls::stream<vec> benes_stream_5_2_1("benes_stream_5_2_1");
	#pragma HLS STREAM variable=benes_stream_5_2_1 depth=64

	static hls::stream<vec> benes_stream_5_3_0("benes_stream_5_3_0");
	#pragma HLS STREAM variable=benes_stream_5_3_0 depth=64

	static hls::stream<vec> benes_stream_5_3_1("benes_stream_5_3_1");
	#pragma HLS STREAM variable=benes_stream_5_3_1 depth=64

	static hls::stream<vec> benes_stream_5_4_0("benes_stream_5_4_0");
	#pragma HLS STREAM variable=benes_stream_5_4_0 depth=64

	static hls::stream<vec> benes_stream_5_4_1("benes_stream_5_4_1");
	#pragma HLS STREAM variable=benes_stream_5_4_1 depth=64

	static hls::stream<vec> benes_stream_5_5_0("benes_stream_5_5_0");
	#pragma HLS STREAM variable=benes_stream_5_5_0 depth=64

	static hls::stream<vec> benes_stream_5_5_1("benes_stream_5_5_1");
	#pragma HLS STREAM variable=benes_stream_5_5_1 depth=64

	static hls::stream<vec> benes_stream_5_6_0("benes_stream_5_6_0");
	#pragma HLS STREAM variable=benes_stream_5_6_0 depth=64

	static hls::stream<vec> benes_stream_5_6_1("benes_stream_5_6_1");
	#pragma HLS STREAM variable=benes_stream_5_6_1 depth=64
	//

	static hls::stream<vec> benes_stream_6("benes_stream_6");
	#pragma HLS STREAM variable=benes_stream_6 depth=64

	static hls::stream<vec> benes_stream_6_0_0("benes_stream_6_0_0");
	#pragma HLS STREAM variable=benes_stream_6_0_0 depth=64

	static hls::stream<vec> benes_stream_6_0_1("benes_stream_6_0_1");
	#pragma HLS STREAM variable=benes_stream_6_0_1 depth=64

	static hls::stream<vec> benes_stream_6_1_0("benes_stream_6_1_0");
	#pragma HLS STREAM variable=benes_stream_6_1_0 depth=64

	static hls::stream<vec> benes_stream_6_1_1("benes_stream_6_1_1");
	#pragma HLS STREAM variable=benes_stream_6_1_1 depth=64

	static hls::stream<vec> benes_stream_6_2_0("benes_stream_6_2_0");
	#pragma HLS STREAM variable=benes_stream_6_2_0 depth=64

	static hls::stream<vec> benes_stream_6_2_1("benes_stream_6_2_1");
	#pragma HLS STREAM variable=benes_stream_6_2_1 depth=64

	static hls::stream<vec> benes_stream_6_3_0("benes_stream_6_3_0");
	#pragma HLS STREAM variable=benes_stream_6_3_0 depth=64

	static hls::stream<vec> benes_stream_6_3_1("benes_stream_6_3_1");
	#pragma HLS STREAM variable=benes_stream_6_3_1 depth=64

	static hls::stream<vec> benes_stream_6_4_0("benes_stream_6_4_0");
	#pragma HLS STREAM variable=benes_stream_6_4_0 depth=64

	static hls::stream<vec> benes_stream_6_4_1("benes_stream_6_4_1");
	#pragma HLS STREAM variable=benes_stream_6_4_1 depth=64

	static hls::stream<vec> benes_stream_6_5_0("benes_stream_6_5_0");
	#pragma HLS STREAM variable=benes_stream_6_5_0 depth=64

	static hls::stream<vec> benes_stream_6_5_1("benes_stream_6_5_1");
	#pragma HLS STREAM variable=benes_stream_6_5_1 depth=64

	static hls::stream<vec> benes_stream_6_6_0("benes_stream_6_6_0");
	#pragma HLS STREAM variable=benes_stream_6_6_0 depth=64

	static hls::stream<vec> benes_stream_6_6_1("benes_stream_6_6_1");
	#pragma HLS STREAM variable=benes_stream_6_6_1 depth=64
	//

	static hls::stream<vec> benes_stream_7("benes_stream_7");
	#pragma HLS STREAM variable=benes_stream_7 depth=64

	static hls::stream<vec> benes_stream_7_0_0("benes_stream_7_0_0");
	#pragma HLS STREAM variable=benes_stream_7_0_0 depth=64

	static hls::stream<vec> benes_stream_7_0_1("benes_stream_7_0_1");
	#pragma HLS STREAM variable=benes_stream_7_0_1 depth=64

	static hls::stream<vec> benes_stream_7_1_0("benes_stream_7_1_0");
	#pragma HLS STREAM variable=benes_stream_7_1_0 depth=64

	static hls::stream<vec> benes_stream_7_1_1("benes_stream_7_1_1");
	#pragma HLS STREAM variable=benes_stream_7_1_1 depth=64

	static hls::stream<vec> benes_stream_7_2_0("benes_stream_7_2_0");
	#pragma HLS STREAM variable=benes_stream_7_2_0 depth=64

	static hls::stream<vec> benes_stream_7_2_1("benes_stream_7_2_1");
	#pragma HLS STREAM variable=benes_stream_7_2_1 depth=64

	static hls::stream<vec> benes_stream_7_3_0("benes_stream_7_3_0");
	#pragma HLS STREAM variable=benes_stream_7_3_0 depth=64

	static hls::stream<vec> benes_stream_7_3_1("benes_stream_7_3_1");
	#pragma HLS STREAM variable=benes_stream_7_3_1 depth=64

	static hls::stream<vec> benes_stream_7_4_0("benes_stream_7_4_0");
	#pragma HLS STREAM variable=benes_stream_7_4_0 depth=64

	static hls::stream<vec> benes_stream_7_4_1("benes_stream_7_4_1");
	#pragma HLS STREAM variable=benes_stream_7_4_1 depth=64

	static hls::stream<vec> benes_stream_7_5_0("benes_stream_7_5_0");
	#pragma HLS STREAM variable=benes_stream_7_5_0 depth=64

	static hls::stream<vec> benes_stream_7_5_1("benes_stream_7_5_1");
	#pragma HLS STREAM variable=benes_stream_7_5_1 depth=64

	static hls::stream<vec> benes_stream_7_6_0("benes_stream_7_6_0");
	#pragma HLS STREAM variable=benes_stream_7_6_0 depth=64

	static hls::stream<vec> benes_stream_7_6_1("benes_stream_7_6_1");
	#pragma HLS STREAM variable=benes_stream_7_6_1 depth=64
	//////

	static hls::stream<vec> benes_stream_8("benes_stream_8");
	#pragma HLS STREAM variable=benes_stream_8 depth=64

	static hls::stream<vec> benes_stream_8_0_0("benes_stream_8_0_0");
	#pragma HLS STREAM variable=benes_stream_8_0_0 depth=64

	static hls::stream<vec> benes_stream_8_0_1("benes_stream_8_0_1");
	#pragma HLS STREAM variable=benes_stream_8_0_1 depth=64

	static hls::stream<vec> benes_stream_8_1_0("benes_stream_8_1_0");
	#pragma HLS STREAM variable=benes_stream_8_1_0 depth=64

	static hls::stream<vec> benes_stream_8_1_1("benes_stream_8_1_1");
	#pragma HLS STREAM variable=benes_stream_8_1_1 depth=64

	static hls::stream<vec> benes_stream_8_2_0("benes_stream_8_2_0");
	#pragma HLS STREAM variable=benes_stream_8_2_0 depth=64

	static hls::stream<vec> benes_stream_8_2_1("benes_stream_8_2_1");
	#pragma HLS STREAM variable=benes_stream_8_2_1 depth=64

	static hls::stream<vec> benes_stream_8_3_0("benes_stream_8_3_0");
	#pragma HLS STREAM variable=benes_stream_8_3_0 depth=64

	static hls::stream<vec> benes_stream_8_3_1("benes_stream_8_3_1");
	#pragma HLS STREAM variable=benes_stream_8_3_1 depth=64

	static hls::stream<vec> benes_stream_8_4_0("benes_stream_8_4_0");
	#pragma HLS STREAM variable=benes_stream_8_4_0 depth=64

	static hls::stream<vec> benes_stream_8_4_1("benes_stream_8_4_1");
	#pragma HLS STREAM variable=benes_stream_8_4_1 depth=64

	static hls::stream<vec> benes_stream_8_5_0("benes_stream_8_5_0");
	#pragma HLS STREAM variable=benes_stream_8_5_0 depth=64

	static hls::stream<vec> benes_stream_8_5_1("benes_stream_8_5_1");
	#pragma HLS STREAM variable=benes_stream_8_5_1 depth=64

	static hls::stream<vec> benes_stream_8_6_0("benes_stream_8_6_0");
	#pragma HLS STREAM variable=benes_stream_8_6_0 depth=64

	static hls::stream<vec> benes_stream_8_6_1("benes_stream_8_6_1");
	#pragma HLS STREAM variable=benes_stream_8_6_1 depth=64



	//
	static hls::stream<mat_64_struct> cond_stream_5_0("cond_stream_5_0");
	#pragma HLS STREAM variable=cond_stream_5_0 depth=2

	static hls::stream<mat_64_struct> cond_stream_5_1("cond_stream_5_1");
	#pragma HLS STREAM variable=cond_stream_5_1 depth=2

	static hls::stream<mat_64_struct> cond_stream_5_2("cond_stream_5_2");
	#pragma HLS STREAM variable=cond_stream_5_2 depth=2

	static hls::stream<mat_64_struct> cond_stream_5_3("cond_stream_5_3");
	#pragma HLS STREAM variable=cond_stream_5_3 depth=2

	static hls::stream<mat_64_struct> cond_stream_5_4("cond_stream_5_4");
	#pragma HLS STREAM variable=cond_stream_5_4 depth=2

	static hls::stream<mat_64_struct> cond_stream_5_5("cond_stream_5_5");
	#pragma HLS STREAM variable=cond_stream_5_5 depth=2

	//
	static hls::stream<mat_64_struct> cond_stream_6_0("cond_stream_6_0");
	#pragma HLS STREAM variable=cond_stream_6_0 depth=2

	static hls::stream<mat_64_struct> cond_stream_6_1("cond_stream_6_1");
	#pragma HLS STREAM variable=cond_stream_6_1 depth=2

	static hls::stream<mat_64_struct> cond_stream_6_2("cond_stream_6_2");
	#pragma HLS STREAM variable=cond_stream_6_2 depth=2

	static hls::stream<mat_64_struct> cond_stream_6_3("cond_stream_6_3");
	#pragma HLS STREAM variable=cond_stream_6_3 depth=2

	static hls::stream<mat_64_struct> cond_stream_6_4("cond_stream_6_4");
	#pragma HLS STREAM variable=cond_stream_6_4 depth=2

	static hls::stream<mat_64_struct> cond_stream_6_5("cond_stream_6_5");
	#pragma HLS STREAM variable=cond_stream_6_5 depth=2

	//
	static hls::stream<mat_64_struct> cond_stream_7_0("cond_stream_7_0");
	#pragma HLS STREAM variable=cond_stream_7_0 depth=2

	static hls::stream<mat_64_struct> cond_stream_7_1("cond_stream_7_1");
	#pragma HLS STREAM variable=cond_stream_7_1 depth=2

	static hls::stream<mat_64_struct> cond_stream_7_2("cond_stream_7_2");
	#pragma HLS STREAM variable=cond_stream_7_2 depth=2

	static hls::stream<mat_64_struct> cond_stream_7_3("cond_stream_7_3");
	#pragma HLS STREAM variable=cond_stream_7_3 depth=2

	static hls::stream<mat_64_struct> cond_stream_7_4("cond_stream_7_4");
	#pragma HLS STREAM variable=cond_stream_7_4 depth=2

	static hls::stream<mat_64_struct> cond_stream_7_5("cond_stream_7_5");
	#pragma HLS STREAM variable=cond_stream_7_5 depth=2

	//
	static hls::stream<mat_64_struct> cond_stream_8_0("cond_stream_8_0");
	#pragma HLS STREAM variable=cond_stream_8_0 depth=2

	static hls::stream<mat_64_struct> cond_stream_8_1("cond_stream_8_1");
	#pragma HLS STREAM variable=cond_stream_8_1 depth=2

	static hls::stream<mat_64_struct> cond_stream_8_2("cond_stream_8_2");
	#pragma HLS STREAM variable=cond_stream_8_2 depth=2

	static hls::stream<mat_64_struct> cond_stream_8_3("cond_stream_8_3");
	#pragma HLS STREAM variable=cond_stream_8_3 depth=2

	static hls::stream<mat_64_struct> cond_stream_8_4("cond_stream_8_4");
	#pragma HLS STREAM variable=cond_stream_8_4 depth=2

	static hls::stream<mat_64_struct> cond_stream_8_5("cond_stream_8_5");
	#pragma HLS STREAM variable=cond_stream_8_5 depth=2

	//

	static hls::stream<sk_2k_vec> vec_mat_stream_5("vec_mat_stream_5");
	#pragma HLS STREAM variable=vec_mat_stream_5 depth=7

	static hls::stream<sk_2k_vec> vec_mat_stream_6("vec_mat_stream_6");
	#pragma HLS STREAM variable=vec_mat_stream_6 depth=6

	static hls::stream<sk_2k_vec> vec_mat_stream_7("vec_mat_stream_7");
	#pragma HLS STREAM variable=vec_mat_stream_7 depth=7



	/////LAYER 1//////
		preprocess_layer_1(cond_stream_5_0, cond_stream_5_1, cond_stream_5_2, cond_stream_5_3, cond_stream_5_4, cond_stream_5_5, vec_mat_stream_5, vec_mat_stream_6, vec_mat_stream_7, benes_stream_5, error_stream_1, sk_vec_fwd_mat_stream);

		split_benes_stream(benes_stream_5_0_0, benes_stream_5_0_1, benes_stream_5);
		benes_layer(benes_stream_5_1_0, benes_stream_5_1_1, benes_stream_5_0_0, benes_stream_5_0_1, cond_stream_5_0, 0);
		benes_layer(benes_stream_5_2_0, benes_stream_5_2_1, benes_stream_5_1_0, benes_stream_5_1_1, cond_stream_5_1, 1);
		benes_layer(benes_stream_5_3_0, benes_stream_5_3_1, benes_stream_5_2_0, benes_stream_5_2_1, cond_stream_5_2, 2);
		benes_layer(benes_stream_5_4_0, benes_stream_5_4_1, benes_stream_5_3_0, benes_stream_5_3_1, cond_stream_5_3, 3);
		benes_layer(benes_stream_5_5_0, benes_stream_5_5_1, benes_stream_5_4_0, benes_stream_5_4_1, cond_stream_5_4, 4);
		benes_layer(benes_stream_5_6_0, benes_stream_5_6_1, benes_stream_5_5_0, benes_stream_5_5_1, cond_stream_5_5, 5);



	/////LAYER 2//////
		preprocess_layer_2(cond_stream_6_0, cond_stream_6_1, cond_stream_6_2, cond_stream_6_3, cond_stream_6_4, cond_stream_6_5, benes_stream_6, vec_mat_stream_5, benes_stream_5_6_0, benes_stream_5_6_1);

		split_benes_stream(benes_stream_6_0_0, benes_stream_6_0_1, benes_stream_6);
		benes_layer(benes_stream_6_1_0, benes_stream_6_1_1, benes_stream_6_0_0, benes_stream_6_0_1, cond_stream_6_0, 0);
		benes_layer(benes_stream_6_2_0, benes_stream_6_2_1, benes_stream_6_1_0, benes_stream_6_1_1, cond_stream_6_1, 1);
		benes_layer(benes_stream_6_3_0, benes_stream_6_3_1, benes_stream_6_2_0, benes_stream_6_2_1, cond_stream_6_2, 2);
		benes_layer(benes_stream_6_4_0, benes_stream_6_4_1, benes_stream_6_3_0, benes_stream_6_3_1, cond_stream_6_3, 3);
		benes_layer(benes_stream_6_5_0, benes_stream_6_5_1, benes_stream_6_4_0, benes_stream_6_4_1, cond_stream_6_4, 4);
		benes_layer(benes_stream_6_6_0, benes_stream_6_6_1, benes_stream_6_5_0, benes_stream_6_5_1, cond_stream_6_5, 5);


		/////LAYER 3//////
		preprocess_layer_3(cond_stream_7_0, cond_stream_7_1, cond_stream_7_2, cond_stream_7_3, cond_stream_7_4, vec_mat_stream_6);

		split_benes_short_stream(benes_stream_7_0_0, benes_stream_7_0_1, benes_stream_6_6_0, benes_stream_6_6_1);
		benes_short_layer(benes_stream_7_1_0, benes_stream_7_1_1, benes_stream_7_0_0, benes_stream_7_0_1, cond_stream_7_4, 1, 16);
		benes_short_layer(benes_stream_7_2_0, benes_stream_7_2_1, benes_stream_7_1_0, benes_stream_7_1_1, cond_stream_7_3, 8, 2);
		benes_short_layer(benes_stream_7_3_0, benes_stream_7_3_1, benes_stream_7_2_0, benes_stream_7_2_1, cond_stream_7_2, 4, 4);
		benes_short_layer(benes_stream_7_4_0, benes_stream_7_4_1, benes_stream_7_3_0, benes_stream_7_3_1, cond_stream_7_1, 2, 8);
		benes_short_layer(benes_stream_7_5_0, benes_stream_7_5_1, benes_stream_7_4_0, benes_stream_7_4_1, cond_stream_7_0, 1, 16);


	/////LAYER 4//////
		preprocess_layer_4(cond_stream_8_0, cond_stream_8_1, cond_stream_8_2, cond_stream_8_3, cond_stream_8_4, cond_stream_8_5, benes_stream_8, vec_mat_stream_7, benes_stream_7_5_0, benes_stream_7_5_1);

		split_benes_stream(benes_stream_8_0_0, benes_stream_8_0_1, benes_stream_8);
		benes_tr_layer(benes_stream_8_1_0, benes_stream_8_1_1, benes_stream_8_0_0, benes_stream_8_0_1, cond_stream_8_5, 5);
		benes_tr_layer(benes_stream_8_2_0, benes_stream_8_2_1, benes_stream_8_1_0, benes_stream_8_1_1, cond_stream_8_4, 4);
		benes_tr_layer(benes_stream_8_3_0, benes_stream_8_3_1, benes_stream_8_2_0, benes_stream_8_2_1, cond_stream_8_3, 3);
		benes_tr_layer(benes_stream_8_4_0, benes_stream_8_4_1, benes_stream_8_3_0, benes_stream_8_3_1, cond_stream_8_2, 2);
		benes_tr_layer(benes_stream_8_5_0, benes_stream_8_5_1, benes_stream_8_4_0, benes_stream_8_4_1, cond_stream_8_1, 1);
		benes_tr_layer(benes_stream_8_6_0, benes_stream_8_6_1, benes_stream_8_5_0, benes_stream_8_5_1, cond_stream_8_0, 0);


		postprocess_layer_4(error_stream_benes_out, benes_stream_8_6_0, benes_stream_8_6_1);


	/////////

	error_stream_duplicate(error_stream_benes_out_0, error_stream_benes_out_1, error_stream_benes_out);

	postprocess(e, e_stream_int, error_stream_benes_out_0);

	weight_check(check_weight_stream, e_stream_int, error_stream_benes_out_1);

	check_value_compute(check, check_synd_stream, check_weight_stream);

}
