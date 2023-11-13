/*
  This file is for Niederreiter encryption
*/

#ifndef ENCRYPT_H
#define ENCRYPT_H
#define encrypt CRYPTO_NAMESPACE(encrypt)

void encrypt(unsigned char *, const unsigned char *, unsigned char *);

extern double sum_encryption_func;
extern int times_encryption_func;

#endif

