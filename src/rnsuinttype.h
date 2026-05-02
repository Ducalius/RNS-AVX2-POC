#ifndef RNSUINTTYPE_H_INCLUDED
#define RNSUINTTYPE_H_INCLUDED

#include <immintrin.h>
#include <stdint.h>
#include "consts.h"
#include "bn.h"
#include "rnsblocktype.h"




typedef struct _rns4096_uint {
    int blocks;
    _rns256_block b_st[RNS4096_BLOCKSIZE];
    _rns256_block b_mg[RNS4096_BLOCKSIZE];
} _rns4096_uint;

void rns4096_init(_rns4096_uint* n);

void rns4096_from_hex(_rns4096_uint* val, char* str, int n);

void rns4096_from_dec(_rns4096_uint* val, char* str, int n);

void rns4096_copy(_rns4096_uint* src, _rns4096_uint* dest);

void rns4096_copy1(_rns4096_uint* src, _rns4096_uint* dest, __RNS_BASE base);

void rns4096_to_dec(_rns4096_uint* val, mpz_t dest, __RNS_BASE base);

void rns4096_to_dec_mod(_rns4096_uint* val, char* pstr, mpz_t dest, __RNS_BASE base);

void rns4096_from_arr(_rns4096_uint* val, uint32_t* arr, __RNS_BASE base);

void rns4096_base_copy(_rns4096_uint* val, __RNS_BASE base);

void rns4096_add(_rns4096_uint* a, _rns4096_uint* b, _rns4096_uint* dest, __RNS_BASE base);

void rns4096_mul(_rns4096_uint* a, _rns4096_uint* b, _rns4096_uint* dest, __RNS_BASE base);

void rns4096_mul_1mod(_rns4096_uint* a, _rns4096_uint* b, _rns4096_uint* dest, int mod_id, __RNS_BASE base);
#endif // RNSUINTTYPE_H_INCLUDED
