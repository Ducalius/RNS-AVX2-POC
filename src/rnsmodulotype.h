#ifndef RNSMODULOTYPE_H_INCLUDED
#define RNSMODULOTYPE_H_INCLUDED

#include "rnsblocktype.h"

typedef struct _rns4096_modulo {
    int blocks;
    _rns256_block p_inv[RNS4096_BLOCKSIZE];
    _rns256_block p_mg[RNS4096_BLOCKSIZE];
    _rns256_block p_M_st[RNS4096_BLOCKSIZE];
    _rns256_block p_M_mg[RNS4096_BLOCKSIZE];
    _rns256_block p_M2_st[RNS4096_BLOCKSIZE];
    _rns256_block p_M2_mg[RNS4096_BLOCKSIZE];
} _rns4096_modulo;

void rns4096_modulo_init(_rns4096_modulo* n);

void rns4096_modulo_from_dec(_rns4096_modulo* dest, char* str, int n);



#endif // RNSMODULOTYPE_H_INCLUDED
