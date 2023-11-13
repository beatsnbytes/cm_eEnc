#include <gmp.h>
#define __gmp_const const

#include <stdint.h>
#include <stdlib.h>
#define AP_INT_MAX_W 4096
#include"ap_int.h"


#define BIT_WIDTH 840
#define ALIGNED_WIDTH 1024
#define ALIGNED_BYTES 128
#define PACK_FACTOR BIT_WIDTH/8
typedef ap_uint<BIT_WIDTH> custom_width;
typedef ap_uint<ALIGNED_WIDTH> aligned_width;

typedef ap_uint<64> vec;
typedef ap_uint<2048> sk;

#define GFBITS 12
#define SYS_N 3488
#define SYS_T 64

#define COND_BYTES ((1 << (GFBITS-4))*(2*GFBITS - 1))
#define IRR_BYTES (SYS_T * 2)

#define PK_NROWS (SYS_T*GFBITS)
#define PK_NCOLS (SYS_N - PK_NROWS)
#define PK_ROW_BYTES ((PK_NCOLS + 7)/8)

#define SYND_BYTES ((PK_NROWS + 7)/8)

#define GFMASK ((1 << GFBITS) - 1)

void benes_kernel(vec *r, unsigned char *sk_values, unsigned char *s);


