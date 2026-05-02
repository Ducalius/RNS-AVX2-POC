#include "rnsmodulotype.h"
#include "mathhelpers.h"
#include <gmp.h>

void rns4096_modulo_init(_rns4096_modulo* n){
    n->blocks = RNS4096_BLOCKSIZE;
    for (int d = 0; d < n->blocks; d++){
        _rns256_block_init(&(n->p_inv[d]), d);
        _rns256_block_init(&(n->p_mg[d]), d);
        _rns256_block_init(&(n->p_M_st[d]), d);
        _rns256_block_init(&(n->p_M_mg[d]), d);
        _rns256_block_init(&(n->p_M2_st[d]), d);
        _rns256_block_init(&(n->p_M2_mg[d]), d);
    }
}


void rns4096_modulo_from_dec(_rns4096_modulo* dest, char* str, int n) {
    mpz_t Msq, Pdec, mpztemp, Msq2;
    mpz_init(Pdec);
    mpz_init(Msq2);
    mpz_init(mpztemp);
    mpz_set_str(Pdec, str, 10);

    mpz_init_set(Msq, __RNSBASE_M);

    mpz_mod(Msq, Msq, Pdec);

    mpz_mul(Msq2, Msq, Msq);
    mpz_mod(Msq2, Msq2, Pdec);

    for (int b = 0; b < dest->blocks; b++) {
        uint32_t vec_inv[8];
        uint32_t vec_mg[8];
        uint32_t vec_M_st[8];
        uint32_t vec_M_mg[8];
        uint32_t vec_M2_st[8];
        uint32_t vec_M2_mg[8];
        for (int i = 0; i < 8; i++) {
            int res = 0;
            uint32_t mi = __RNSBASE_mi[b * 8 + i];
            uint32_t temp = fast_mod(str, mi, n);
            temp = mod_invert(temp,mi);
            vec_inv[i] = neg_mod(temp,mi);
            //vec_inv[i] = temp;

            uint32_t mi_mg = __RNSBASE_MG_mi[b * 8 + i];
            vec_mg[i] = fast_mod(str, mi_mg, n);

            mpz_mod_ui(mpztemp, Msq, mi);
            
            vec_M_st[i] = (uint32_t)mpz_get_ui(mpztemp);

            mpz_mod_ui(mpztemp, Msq, mi_mg);
            vec_M_mg[i] = (uint32_t)mpz_get_ui(mpztemp);

            mpz_mod_ui(mpztemp, Msq2, mi);
            vec_M2_st[i] = (uint32_t)mpz_get_ui(mpztemp);
            mpz_mod_ui(mpztemp, Msq2, mi_mg);
            vec_M2_mg[i] = (uint32_t)mpz_get_ui(mpztemp);

        };
        dest->p_inv[b].data = _mm256_loadu_si256(vec_inv);
        dest->p_mg[b].data = _mm256_loadu_si256(vec_mg);
        dest->p_M_st[b].data = _mm256_loadu_si256(vec_M_st);
        dest->p_M_mg[b].data = _mm256_loadu_si256(vec_M_mg);
        dest->p_M2_st[b].data = _mm256_loadu_si256(vec_M2_st);
        dest->p_M2_mg[b].data = _mm256_loadu_si256(vec_M2_mg);
    };
};
