#include <gmp.h>
#define __gmp_const const

#include <stdint.h>
#include <stdlib.h>
#define AP_INT_MAX_W 4096
#include"ap_int.h"

typedef ap_uint<64> vec;
typedef ap_uint<4096> sk_4k_vec;
typedef ap_uint<2048> sk_2k_vec;
typedef ap_uint<1024> sk_1k_vec;
typedef ap_uint<512> sk_05k_vec;

#define GFBITS 13
#define SYS_N 4608
#define SYS_T 96

#define COND_BYTES ((1 << (GFBITS-4))*(2*GFBITS - 1))
#define IRR_BYTES (SYS_T * 2)

#define PK_NROWS (SYS_T*GFBITS)
#define PK_NCOLS (SYS_N - PK_NROWS)
#define PK_ROW_BYTES ((PK_NCOLS + 7)/8)

#define SYND_BYTES ((PK_NROWS + 7)/8)

#define GFMASK ((1 << GFBITS) - 1)

void decryption_kernel(int *, unsigned char *, sk_05k_vec *, sk_4k_vec *, sk_4k_vec *, unsigned char *);


