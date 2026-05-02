#include "rnsuinttype.h"
#include "mathhelpers.h"
#include <stdlib.h>

void rns4096_init(_rns4096_uint* n){
    n->blocks = RNS4096_BLOCKSIZE;
    for (int d = 0; d < n->blocks; d++){
        _rns256_block_init(&(n->b_st[d]), d);
        _rns256_block_init(&(n->b_mg[d]), d);
    }
}

void rns4096_from_hex(_rns4096_uint* val, char* hexstr, int n) {
    char* str = (char*)malloc(n*2*sizeof(char));
    hex_to_dec(hexstr, str);
    for (int b = 0; b < val->blocks; b++) {
        uint32_t vec_st[8];
        uint32_t vec_mg[8];
        for (int i = 0; i < 8; i++) {
            int res = 0;
            uint32_t mi = __RNSBASE_mi[b * 8 + i];
            vec_st[i] = fast_mod(str, mi, n);

            uint32_t mi_mg = __RNSBASE_MG_mi[b * 8 + i];
            vec_mg[i] = fast_mod(str, mi_mg, n);
        };
        val->b_st[b].data = _mm256_loadu_si256(vec_st);
        val->b_mg[b].data = _mm256_loadu_si256(vec_mg);
    };
};

void rns4096_from_dec(_rns4096_uint* val, char* str, int n) {
    for (int b = 0; b < val->blocks; b++) {
        uint32_t vec_st[8];
        uint32_t vec_mg[8];
        for (int i = 0; i < 8; i++) {
            int res = 0;
            uint32_t mi = __RNSBASE_mi[b * 8 + i];
            vec_st[i] = fast_mod(str, mi, n);

            uint32_t mi_mg = __RNSBASE_MG_mi[b * 8 + i];
            vec_mg[i] = fast_mod(str, mi_mg, n);
        };
        val->b_st[b].data = _mm256_loadu_si256(vec_st);
        val->b_mg[b].data = _mm256_loadu_si256(vec_mg);
    };
};

void rns4096_from_arr(_rns4096_uint* val, uint32_t* arr, __RNS_BASE base) {
    switch (base) {
    case BASE_ST:
        //printf("arr load 1st element %u\n", arr[0]);
        for (int b = 0; b < val->blocks; b++) {
            val->b_st[b].data = _mm256_loadu_si256(&(arr[b * 8]));
        }
        break;
    case BASE_MG:
        for (int b = 0; b < val->blocks; b++)
            val->b_mg[b].data = _mm256_loadu_si256(&(arr[b * 8]));
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };
}


void rns4096_base_copy(_rns4096_uint* val, __RNS_BASE base) {
    switch (base) {
    case BASE_ST:
        for (int b = 0; b < val->blocks; b++)
            val->b_mg[b].data = val->b_st[b].data;
        break;
    case BASE_MG:
        for (int b = 0; b < val->blocks; b++)
            val->b_st[b].data = val->b_mg[b].data;
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };
}


void rns4096_add(_rns4096_uint* a, _rns4096_uint* b, _rns4096_uint* dest, __RNS_BASE base) {
    for (int i = 0; i < a->blocks; i++) {
        switch (base) {
        case BASE_ST:
            _rns256_block_add(&(a->b_st[i]), &(b->b_st[i]), &(dest->b_st[i]), i, base);
            break;
        case BASE_MG:
            _rns256_block_add(&(a->b_mg[i]), &(b->b_mg[i]), &(dest->b_mg[i]), i, base);
            break;
        };
    }
}


void rns4096_mul(_rns4096_uint* a, _rns4096_uint* b, _rns4096_uint* dest, __RNS_BASE base) {
    for (int i = 0; i < a->blocks; i++) {
        switch (base) {
        case BASE_ST:
            _rns256_block_mul(&(a->b_st[i]), &(b->b_st[i]), &(dest->b_st[i]), i, base);
            break;
        case BASE_MG:
            _rns256_block_mul(&(a->b_mg[i]), &(b->b_mg[i]), &(dest->b_mg[i]), i, base);
            break;
        };
    }

}


void rns4096_to_dec(_rns4096_uint* val, mpz_t dest, __RNS_BASE base) {
    _rns4096_uint w, res;
    rns4096_init(&w);
    rns4096_init(&res);
    mpz_t Mi, x, temp;
    mpz_set_ui(dest, 0);
    switch (base) {
    case BASE_ST:
        rns4096_from_arr(&w, __RNSBASE_Mi_inv, BASE_ST);
        rns4096_mul(val, &w, &res, BASE_ST);

        for (int b = 0; b < val->blocks; b++) {
            uint32_t res_vec[8];
            _mm256_storeu_epi32(res_vec, res.b_st[b].data);
            for (int i = 0; i < 8; i++) {
                mpz_init_set_ui(x, res_vec[i]);

                mpz_init_set(Mi, __RNSBASE_M);

                mpz_init(temp);
                mpz_divexact_ui(Mi, Mi, __RNSBASE_mi[b * 8 + i]);

                mpz_mul(temp, x, Mi);
                mpz_add(dest, dest, temp);
            }
        }
        mpz_mod(dest, dest, __RNSBASE_M);
        break;
    case BASE_MG:
        rns4096_from_arr(&w, __RNSBASE_MG_Mi_inv_mg, BASE_MG);
        rns4096_mul(val, &w, &res, BASE_MG);


        for (int b = 0; b < val->blocks; b++) {
            uint32_t res_vec[8];
            _mm256_storeu_epi32(res_vec, res.b_mg[b].data);
            for (int i = 0; i < 8; i++) {
                mpz_init_set_ui(x, res_vec[i]);

                mpz_init_set(Mi, __RNSBASE_MG_M);

                mpz_init(temp);
                mpz_divexact_ui(Mi, Mi, __RNSBASE_MG_mi[b * 8 + i]);

                mpz_mul(temp, x, Mi);
                mpz_add(dest, dest, temp);
            }
        }
        mpz_mod(dest, dest, __RNSBASE_MG_M);
        break;
    };
        mpz_clear(Mi);
        mpz_clear(x);
        mpz_clear(temp);
    
}

void rns4096_to_dec_mod(_rns4096_uint* val, char* pstr, mpz_t dest, __RNS_BASE base) {
    _rns4096_uint w, res;
    rns4096_init(&w);
    rns4096_init(&res);
    mpz_t Mi, x, temp, p;
    mpz_init(p);
    mpz_set_str(p, pstr, 10);
    mpz_set_ui(dest, 0);
    switch (base) {
    case BASE_ST:
        rns4096_from_arr(&w, __RNSBASE_Mi_inv, BASE_ST);
        rns4096_mul(val, &w, &res, BASE_ST);

        for (int b = 0; b < val->blocks; b++) {
            uint32_t res_vec[8];
            _mm256_storeu_epi32(res_vec, res.b_st[b].data);
            for (int i = 0; i < 8; i++) {
                mpz_init_set_ui(x, res_vec[i]);

                mpz_init_set(Mi, __RNSBASE_M);

                mpz_init(temp);
                mpz_divexact_ui(Mi, Mi, __RNSBASE_mi[b * 8 + i]);

                mpz_mul(temp, x, Mi);
                mpz_add(dest, dest, temp);
            }
        }
        mpz_mod(dest, dest, __RNSBASE_M);
        break;
    case BASE_MG:
        rns4096_from_arr(&w, __RNSBASE_MG_Mi_inv_mg, BASE_MG);
        rns4096_mul(val, &w, &res, BASE_MG);


        for (int b = 0; b < val->blocks; b++) {
            uint32_t res_vec[8];
            _mm256_storeu_epi32(res_vec, res.b_mg[b].data);
            for (int i = 0; i < 8; i++) {
                mpz_init_set_ui(x, res_vec[i]);

                mpz_init_set(Mi, __RNSBASE_MG_M);

                mpz_init(temp);
                mpz_divexact_ui(Mi, Mi, __RNSBASE_MG_mi[b * 8 + i]);

                mpz_mul(temp, x, Mi);
                mpz_add(dest, dest, temp);
            }
        }
        mpz_mod(dest, dest, __RNSBASE_MG_M);
        break;
    };
    if (mpz_cmp(p, dest) < 0) {
        mpz_sub(dest, dest, p);
    }
    mpz_clear(p);
    mpz_clear(Mi);
    mpz_clear(x);
    mpz_clear(temp);

}

void rns4096_copy(_rns4096_uint* src, _rns4096_uint* dest) {
    for (int b = 0; b < RNS4096_BLOCKSIZE; b++) {
        //printf("block %d\n", b);
        dest->b_st[b].data = src->b_st[b].data;
        dest->b_mg[b].data = src->b_mg[b].data;
    }
}

void rns4096_copy1(_rns4096_uint* src, _rns4096_uint* dest, __RNS_BASE base) {

    for (int b = 0; b < src->blocks; b++) {
        switch (base) {
        case BASE_ST:
            dest->b_st[b].data = src->b_st[b].data;
            break;
        case BASE_MG:
            dest->b_mg[b].data = src->b_mg[b].data;
            break;
        };
    }
}

void rns4096_mul_1mod(_rns4096_uint* a, _rns4096_uint* b, _rns4096_uint* dest, int mod_id,  __RNS_BASE base) {
    for (int i = 0; i < a->blocks; i++) {
        switch (base) {
        case BASE_ST:
            _rns256_block_mul_1mod(&(a->b_st[i]), &(b->b_st[i]), &(dest->b_st[i]), mod_id, base);
            break;
        case BASE_MG:
            _rns256_block_mul_1mod(&(a->b_mg[i]), &(b->b_mg[i]), &(dest->b_mg[i]), mod_id, base);
            break;
        };
    }

}

void rns4096_mul_1mod_debug(_rns4096_uint* a, _rns4096_uint* b, _rns4096_uint* dest, int mod_id, __RNS_BASE base, int iter) {
    for (int i = 0; i < a->blocks; i++) {
        switch (base) {
            _rns256_block_mul_1mod(&(a->b_st[i]), &(b->b_st[i]), &(dest->b_st[i]), mod_id, base);
            break;
        case BASE_MG:
            _rns256_block_mul_1mod(&(a->b_mg[i]), &(b->b_mg[i]), &(dest->b_mg[i]), mod_id, base);
            break;
        };
    }

}