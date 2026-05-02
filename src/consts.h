#ifndef CONSTS_H_INCLUDED
#define CONSTS_H_INCLUDED
#include <stdint.h>
#include <gmp.h>
#include "bn.h"

#define RNSBASE_SIZE 128
#define RNS4096_BLOCKSIZE 16


#define CORRECTION_ROFFSET 16
#define CORRECTION_LOFFSET 8

// Base
extern const uint32_t __RNSBASE_ci[];
extern const uint32_t __RNSBASE_mi[];
// M
extern const char __RNSBASE_M_plain[];
extern mpz_t __RNSBASE_M;

// Alt base
extern const uint32_t __RNSBASE_MG_ci[];
extern const uint32_t __RNSBASE_MG_mi[];
// M'
extern const char __RNSBASE_MG_M_plain[];
extern mpz_t __RNSBASE_MG_M;


// Step 2
extern const uint32_t __RNSBASE_Mi_inv[];

// Step 3
extern const uint32_t __RNSBASE_to_MG_matrix[RNSBASE_SIZE][RNSBASE_SIZE];
extern const uint32_t __RNSBASE_Mneg_mod_MG_mi[];

// Step 5
extern const uint32_t __RNSBASE_MG_M_inv[];

// Step 6
extern const uint32_t __RNSBASE_MG_Mi_inv_mg[];

// Step 7
extern const uint32_t __RNSBASE_from_MG_matrix[RNSBASE_SIZE][RNSBASE_SIZE];
extern const uint32_t __RNSBASE_MG_Mneg_mod_mi[];

extern const uint32_t testing[];

void __RNSBASE_init_consts();

#endif // CONSTS_H_INCLUDED
