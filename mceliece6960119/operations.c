#include "operations.h"

#include "controlbits.h"
#include "randombytes.h"
#include "crypto_hash.h"
#include "encrypt.h"
#include "decrypt.h"
#include "params.h"
#include "sk_gen.h"
#include "pk_gen.h"
#include "util.h"
#include "kat_kem.h"



#include <stdint.h>
#include <string.h>

//Added libraries
#include <sys/time.h>
#include "../common/custom_util.h"

#include <CL/opencl.h>
#include <CL/cl_ext.h>
#include "params.h"

int times_event_enq=0;
double sum_list_event_enq[1];

//Global variables to report timing of KEM parts
double sum_encrypt=0.0;
int times_encrypt=0;
double sum_decrypt=0.0;
int times_decrypt=0;
double sum_keyop=0.0;
int times_keyop=0;


/* check if the padding bits of pk are all zero */
static int check_pk_padding(const unsigned char * pk)
{
	unsigned char b;
	int i, ret;

	b = 0;
	for (i = 0; i < PK_NROWS; i++)
		b |= pk[i*PK_ROW_BYTES + PK_ROW_BYTES-1];

	b >>= (PK_NCOLS % 8);
	b -= 1;
	b >>= 7;
	ret = b;

	return ret-1;
}

int crypto_kem_enc(
       unsigned char *c,
       unsigned char *key,
       const unsigned char *pk
)
{
	unsigned char e[ SYS_N/8 ];
	unsigned char one_ec[ 1 + SYS_N/8 + SYND_BYTES ] = {1};
	unsigned char mask;
	int i, j, padding_ok;

	//

	padding_ok = check_pk_padding(pk);


	cl_event event_enq, event_migr_tokern, event_migr_tohost;

	#ifdef ENCRYPTION_KERNEL
	{

		//The kernel inputs
		memcpy(ptr_v_in, DRBG_ctx.V, sizeof(unsigned char)*16);
		memcpy(ptr_key_in, DRBG_ctx.Key, sizeof(unsigned char)*32);

		// for(i=0; i<PARALLELIZATION_FACTOR; i++){
		// 	memcpy(ptr_pk_in[i], pk+(i*(PK_NROWS*PK_ROW_BYTES)/PARALLELIZATION_FACTOR), sizeof(unsigned char)*(PK_NROWS*PK_ROW_BYTES)/PARALLELIZATION_FACTOR);
		// }

		for(i=0; i<PARALLELIZATION_FACTOR; i++){
			for(j=0; j<(PK_NROWS); j++){
				memcpy( ptr_pk_in[0]+j*(PK_ROW_BYTES+91), pk + (j*PK_ROW_BYTES), sizeof(unsigned char)*(PK_ROW_BYTES));
			}
		}

		clEnqueueMigrateMemObjects(commands, PARALLELIZATION_FACTOR+2, buffer_list_enc_in, 0, 0, NULL, &event_migr_tokern);

		clEnqueueTask(commands, encryption_kernel, 1, &event_migr_tokern, &event_enq);

		clEnqueueMigrateMemObjects(commands, PARALLELIZATION_FACTOR+1, buffer_list_enc_out, CL_MIGRATE_MEM_OBJECT_HOST, 1, &event_enq, &event_migr_tohost);

		clWaitForEvents(1, &event_migr_tohost);

		cl_profile_print(&event_enq, 1, sum_list_event_enq, &times_event_enq);

		for(i=0; i<PARALLELIZATION_FACTOR; i++){
			memcpy(c+(i*(SYND_BYTES)/PARALLELIZATION_FACTOR), ptr_s_out[i], sizeof(unsigned char)*(SYND_BYTES)/PARALLELIZATION_FACTOR);
		}
		memcpy(e, ptr_e_out, sizeof(unsigned char)*(SYS_N/8));

		clReleaseEvent(event_enq);
		clReleaseEvent(event_migr_tokern);
		clReleaseEvent(event_migr_tohost);

	}
	#endif

	#ifndef  ENCRYPTION_KERNEL
	{
	struct timeval start_encrypt, end_encrypt;
	gettimeofday(&start_encrypt, NULL);
	
	encrypt(c, pk, e);

	gettimeofday(&end_encrypt, NULL);
	get_event_time(&start_encrypt, &end_encrypt, &sum_encrypt, &times_encrypt);
	}
	#endif

	#if KAT==1
	{
	int k;
	printf("encrypt e: positions");
	for (k = 0;k < SYS_N;++k)
		if (e[k/8] & (1 << (k&7)))
		printf(" %d",k);
	printf("\n");
	}
	#endif

	///////////////////////////

	memcpy(one_ec + 1, e, SYS_N/8);
	memcpy(one_ec + 1 + SYS_N/8, c, SYND_BYTES);

	crypto_hash_32b(key, one_ec, sizeof(one_ec));

	// clear outputs (set to all 0's) if padding bits are not all zero

	mask = padding_ok;
	mask ^= 0xFF;

	for (i = 0; i < SYND_BYTES; i++)
		c[i] &= mask;

	for (i = 0; i < 32; i++)
		key[i] &= mask;

	return padding_ok;
}

/* check if the padding bits of c are all zero */
static int check_c_padding(const unsigned char * c)
{
	unsigned char b;
	int ret;

	b = c[ SYND_BYTES-1 ] >> (PK_NROWS % 8);
	b -= 1;
	b >>= 7;
	ret = b;

	return ret-1;
}

int crypto_kem_dec(
       unsigned char *key,
       const unsigned char *c,
       const unsigned char *sk
)
{
	int i, padding_ok;

	unsigned char mask;
	unsigned char ret_decrypt = 0;

	uint16_t m;

	unsigned char e[ SYS_N/8 ];
	unsigned char preimage[ 1 + SYS_N/8 + SYND_BYTES ];
	unsigned char *x = preimage;
	const unsigned char *s = sk + 40 + IRR_BYTES + COND_BYTES;

	//

	padding_ok = check_c_padding(c);

	struct timeval start_decrypt, end_decrypt;
	gettimeofday(&start_decrypt, NULL);

	ret_decrypt = decrypt(e, sk + 40, c);

	gettimeofday(&end_decrypt, NULL);
	get_event_time(&start_decrypt, &end_decrypt, &sum_decrypt, &times_decrypt);


	m = ret_decrypt;
	m -= 1;
	m >>= 8;

	*x++ = m & 1;
	for (i = 0; i < SYS_N/8; i++) 
		*x++ = (~m & s[i]) | (m & e[i]);

	for (i = 0; i < SYND_BYTES; i++) 
		*x++ = c[i];

	crypto_hash_32b(key, preimage, sizeof(preimage)); 

	// clear outputs (set to all 1's) if padding bits are not all zero

	mask = padding_ok;

	for (i = 0; i < 32; i++)
		key[i] |= mask;

	return padding_ok;
}

int crypto_kem_keypair
(
       unsigned char *pk,
       unsigned char *sk 
)
{
	int i;
	unsigned char seed[ 33 ] = {64};
	unsigned char r[ SYS_N/8 + (1 << GFBITS)*sizeof(uint32_t) + SYS_T*2 + 32 ];
	unsigned char *rp, *skp;

	gf f[ SYS_T ]; // element in GF(2^mt)
	gf irr[ SYS_T ]; // Goppa polynomial
	uint32_t perm[ 1 << GFBITS ]; // random permutation as 32-bit integers
	int16_t pi[ 1 << GFBITS ]; // random permutation

	randombytes(seed+1, 32);

	while (1)
	{
		rp = &r[ sizeof(r)-32 ];
		skp = sk;

		// expanding and updating the seed

		shake(r, sizeof(r), seed, 33);
		memcpy(skp, seed+1, 32);
		skp += 32 + 8;
		memcpy(seed+1, &r[ sizeof(r)-32 ], 32);

		// generating irreducible polynomial

		rp -= sizeof(f); 

		for (i = 0; i < SYS_T; i++) 
			f[i] = load_gf(rp + i*2); 

		if (genpoly_gen(irr, f)) 
			continue;

		for (i = 0; i < SYS_T; i++)
			store_gf(skp + i*2, irr[i]);

		skp += IRR_BYTES;

		// generating permutation

		rp -= sizeof(perm);

		for (i = 0; i < (1 << GFBITS); i++) 
			perm[i] = load4(rp + i*4); 

		if (pk_gen(pk, skp - IRR_BYTES, perm, pi))
			continue;

		controlbitsfrompermutation(skp, pi, GFBITS, 1 << GFBITS);
		skp += COND_BYTES;

		// storing the random string s

		rp -= SYS_N/8;
		memcpy(skp, rp, SYS_N/8);

		// storing positions of the 32 pivots

		store8(sk + 32, 0xFFFFFFFF);

		break;
	}

	return 0;
}

