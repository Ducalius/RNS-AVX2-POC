#include "mathhelpers.h"
#include <gmp.h>

uint32_t mod_invert(uint32_t num, uint32_t mod) {
    uint64_t A = (uint64_t)num;
    uint64_t M = (uint64_t)mod;
    uint64_t m0 = M;
    long y = 0, x = 1;
    if (M == 1)
        return 0;
    while (A > 1) {
        // q is quotient
        uint64_t q = A / M;
        uint64_t t = M;
        // m is remainder now, process same as
        // Euclid's algo
        M = A % M, A = t;
        t = y;

        // Update y and x
        y = x - q * y;
        x = t;
    }
    // Make x positive
    if (x < 0)
        x += m0;
    return (uint32_t)x;
}

uint32_t neg_mod(uint32_t num, uint32_t mod) {
    long long A = -(uint64_t)num;
    long long M = (uint64_t)mod;
    return (uint32_t)(((A % M) + M) % M);
}

uint32_t fast_mod(char* str, uint32_t mod, int len) {
    long long res = 0;
    for (int i = 0; i < len; i++) {
        if (str[i] == 0)
            break;
        res = (res * 10 + str[i] - '0') % mod;
    }

    return (uint32_t)res;
}

void dec_to_hex(char* dec, char* hex) {
    mpz_t num;
    mpz_init(num);
    mpz_set_str(num, dec, 10);
    mpz_get_str(hex, 16, num);
}

void hex_to_dec(char* hex, char* dec) {
    mpz_t num;
    mpz_init(num);
    mpz_set_str(num, hex, 16);
    mpz_get_str(dec, 10, num);
}