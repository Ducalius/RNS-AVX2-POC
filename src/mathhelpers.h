#ifndef MATHHELPERS_H_INCLUDED
#define MATHHELPERS_H_INCLUDED

#include <stdint.h>

uint32_t mod_invert(uint32_t num, uint32_t mod);

uint32_t neg_mod(uint32_t num, uint32_t mod);

uint32_t fast_mod(char *str, uint32_t mod, int len);

void hex_to_dec(char* hex, char* dec);

void dec_to_hex(char* dec, char* hex);

#endif // MATHHELPERS_H_INCLUDED
