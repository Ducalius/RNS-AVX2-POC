#ifndef RNSBLOCKTYPE_H_INCLUDED
#define RNSBLOCKTYPE_H_INCLUDED

#include <immintrin.h>
#include <stdint.h>
#include "consts.h"
#include "bn.h"
#include <stdio.h>

#define RNS4096_BLOCKSIZE 16

typedef enum{
    BASE_ST,
    BASE_MG
} __RNS_BASE;

typedef struct _rns256_block {
    int id;
    __m256i data;
} _rns256_block;

void _rns256_block_init(_rns256_block* b, int id);

void _rns256_block_add(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int id, __RNS_BASE base);

void _rns256_block_mul(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int id, __RNS_BASE base);

void _rns256_block_trunc(_rns256_block* a, _rns256_block* dest, uint32_t roffset, uint32_t loffset);

void _rns256_block_mul_1mod(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int mod_id, __RNS_BASE base);

void _rns256_block_swap(_rns256_block* a, _rns256_block* b, int apos, int bpos);

#endif // RNSBLOCKTYPE_H_INCLUDED
