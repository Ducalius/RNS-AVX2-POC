#ifndef MONTGOMERY_H_INCLUDED
#define MONTGOMERY_H_INCLUDED

#include "rnsmodulotype.h"
#include "rnsblocktype.h"
#include "rnsuinttype.h"
#include "consts.h"


void rns4096_modulo_mul(_rns4096_uint* a, _rns4096_modulo* p, _rns4096_uint* dest, __RNS_BASE base);

void rns4096_base_conversion(_rns4096_uint* a, _rns4096_uint* dest, __RNS_BASE base);

void rns4096_mred(_rns4096_uint* a, _rns4096_modulo *p, _rns4096_uint* dest);

uint32_t rns4096_get_corr_factor_alpha(_rns4096_uint* a, __RNS_BASE base, uint32_t roffset, uint32_t loffset, float alpha);

uint32_t rns4096_get_corr_factor(_rns4096_uint* a, __RNS_BASE base, uint32_t roffset, uint32_t loffset);

uint32_t get_corr_factor_long(_rns4096_uint* a, __RNS_BASE base);

void rns4096_mred_init2(_rns4096_uint* a, _rns4096_uint* dest, _rns4096_modulo* p);

void rns4096_mred_normalize(_rns4096_uint* a, _rns4096_modulo* p, _rns4096_uint* dest, __RNS_BASE base);

void rns4096_mred_finalize(_rns4096_uint* a, _rns4096_uint* dest, _rns4096_modulo* p);

void rns4096_MM(_rns4096_uint* a, _rns4096_uint* b, _rns4096_modulo* p, _rns4096_uint* dest);

void rns4096_MExp(_rns4096_uint* a, uint32_t e, _rns4096_modulo* p, _rns4096_uint* dest);

void rns4096_MExp_str(_rns4096_uint* a, char* ebuf, _rns4096_modulo* p, _rns4096_uint* dest);

void rns4096_mat_transpose(_rns4096_uint** mat, __RNS_BASE base);

#endif