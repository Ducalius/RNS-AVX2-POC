#include "montgomery.h"
#include <math.h>
#include <time.h>

void rns4096_modulo_mul(_rns4096_uint* a, _rns4096_modulo* p, _rns4096_uint* dest, __RNS_BASE base) {

    for (int i = 0; i < a->blocks; i++) {
        switch (base) {
        case BASE_ST:
            _rns256_block_mul(&(a->b_st[i]), &(p->p_inv[i]), &(dest->b_st[i]), i, base);
            break;
        case BASE_MG:
            _rns256_block_mul(&(a->b_mg[i]), &(p->p_mg[i]), &(dest->b_mg[i]), i, base);
            break;
        };

    }
}

void rns4096_mred_normalize(_rns4096_uint* a, _rns4096_modulo* p, _rns4096_uint* dest, __RNS_BASE base) {

    for (int i = 0; i < a->blocks; i++) {
        switch (base) {
        case BASE_ST:
            _rns256_block_mul(&(a->b_st[i]), &(p->p_M_st[i]), &(dest->b_st[i]), i, base);
            break;
        case BASE_MG:
            _rns256_block_mul(&(a->b_mg[i]), &(p->p_M_mg[i]), &(dest->b_mg[i]), i, base);
            
            break;
        };

    }
}

void rns4096_mred_init2(_rns4096_uint* a, _rns4096_uint* dest, _rns4096_modulo* p) {
    for (int i = 0; i < a->blocks; i++) {
        _rns256_block_mul(&(a->b_st[i]), &(p->p_M2_st[i]), &(dest->b_st[i]), i, BASE_ST);
        _rns256_block_mul(&(a->b_mg[i]), &(p->p_M2_mg[i]), &(dest->b_mg[i]), i, BASE_MG);
    }
}

void rns4096_mred_finalize(_rns4096_uint* a, _rns4096_uint* dest, _rns4096_modulo* p) {
    _rns4096_uint unityvec;
    rns4096_init(&unityvec);
    rns4096_from_dec(&unityvec, "1", 1);
    rns4096_MM(&a, &unityvec, p, &a);
}

uint32_t get_corr_factor_long(_rns4096_uint* a, __RNS_BASE base) {
    double res = 0.0;
    switch (base) {
    case BASE_ST:
        for (int b = 0; b < 16; b++) {
            uint32_t theta_vec[8];
            _mm256_storeu_si256(theta_vec, a->b_st[b].data);
            for (int i = 0; i < 8; i++) {
                //printf("%f / %f=", (double)(uint64_t)theta_vec[i], (double)(uint64_t)__RNSBASE_mi[b * 8 + i]);
                res += (double)(uint64_t)theta_vec[i] / (double)(uint64_t)__RNSBASE_mi[b * 8 + i];
                //printf("%f res=%f\n", (double)(uint64_t)theta_vec[i] / (double)(uint64_t)__RNSBASE_mi[b * 8 + i],res);
            }
        }
        //printf("--------%.9f\n", res);
        break;
    case BASE_MG:
        for (int b = 0; b < 16; b++) {
            uint32_t theta_vec[8];
            _mm256_storeu_si256(theta_vec, a->b_mg[b].data);
            for (int i = 0; i < 8; i++) {
                //printf("%f / %f=", (double)(uint64_t)theta_vec[i], (double)(uint64_t)__RNSBASE_mi[b * 8 + i]);
                res += (double)(uint64_t)theta_vec[i] / (double)(uint64_t)__RNSBASE_MG_mi[b * 8 + i];
                //printf("%f res=%f\n", (double)(uint64_t)theta_vec[i] / (double)(uint64_t)__RNSBASE_mi[b * 8 + i],res);
            }
        }
        break;
    };
    return (uint32_t)(int)res;

}

uint32_t rns4096_get_corr_factor(_rns4096_uint* a, __RNS_BASE base, uint32_t roffset, uint32_t loffset) {
    return rns4096_get_corr_factor_alpha(a, base, roffset, loffset, 0);
}

uint32_t rns4096_get_corr_factor_alpha(_rns4096_uint* a, __RNS_BASE base, uint32_t roffset, uint32_t loffset, float alpha) {
    _rns4096_uint temp;
    rns4096_init(&temp);
    for (int i = 0; i < a->blocks; i++) {
        switch (base) {
        case BASE_ST:
            //printf("trunc:\n");
            //print_block(a->b_st[i]);
            _rns256_block_trunc(&(a->b_st[i]), &(temp.b_st[i]), roffset, loffset);
            break;
        case BASE_MG:
            _rns256_block_trunc(&(a->b_mg[i]), &(temp.b_st[i]), roffset, loffset);
            break;
        };
    }
    
    __m256i resvec = _mm256_setzero_si256();
    for (int i = 0; i < a->blocks; i++) 
        resvec = _mm256_add_epi32(resvec, temp.b_st[i].data);
    uint32_t resvec_buf[8];
    _mm256_storeu_si256(resvec_buf, resvec);
    uint32_t res = resvec_buf[0];
    for (int i = 1; i < 8; i++)
        res += resvec_buf[i];
    res = (uint32_t)(((double)res / pow(2, 24)) + alpha);
    return res;
}

void rns4096_mat_transpose(_rns4096_uint** mat, __RNS_BASE base) {

}


void rns4096_base_conversion(_rns4096_uint* a, _rns4096_uint* dest, __RNS_BASE base, uint32_t fact) {
    _rns4096_uint mat, temp, res, corr_vec, temp2, _a;
    rns4096_init(&mat);
    rns4096_init(&res);
    rns4096_init(&corr_vec);

    uint32_t res_raw[128][128];

    _rns4096_uint res_mat[128];

    _rns256_block fact_block, lookup_block;


    switch (base) {
    case BASE_ST:
        rns4096_base_copy(a, BASE_ST);

        // Multiply
        for (int i = 0; i < RNSBASE_SIZE; i++) {
            rns4096_init(&temp);
            rns4096_from_arr(&mat, __RNSBASE_to_MG_matrix[i], BASE_MG);
            rns4096_mul_1mod(a, &mat, &temp, i, BASE_MG);
            for (int b = 0; b < RNS4096_BLOCKSIZE; b++) {
                _mm256_storeu_si256((__m256i const*)(&(res_raw[i][b * 8])), temp.b_mg[b].data);
            }
        };
        //Add vectors to get new terms

        for (int i = 0; i < RNSBASE_SIZE; i++) {
            rns4096_init(&temp);
            rns4096_init(&temp2);
            for (int b = 0; b < RNS4096_BLOCKSIZE; b++) {
                uint32_t vec[8];
                for (int j = 0; j < 8; j++) 
                    vec[j] = res_raw[b * 8 + j][i];
                //print_vec_test_dec(vec);
                temp.b_mg[b].data = _mm256_loadu_si256(vec);
            }

            rns4096_add(&res, &temp, &temp2, BASE_MG);
            rns4096_copy(&temp2, &res);


        };



        //Get correction vector 
        rns4096_init(&temp);
        
        _rns256_block_init(&fact_block, 0);
        fact_block.data = _mm256_set1_epi32(fact);
        
        for (int b = 0; b < RNS4096_BLOCKSIZE; b++) {

            _rns256_block_init(&lookup_block, b);
            lookup_block.data = _mm256_loadu_si256(&__RNSBASE_Mneg_mod_MG_mi[b * 8]);
            _rns256_block_mul(&lookup_block, &fact_block, &(corr_vec.b_mg[b]), b, BASE_MG);
        }

        rns4096_init(&temp);
        //Correction step
        rns4096_add(&res, &corr_vec, &temp, BASE_MG);

        //Load result
        for (int b = 0; b < RNS4096_BLOCKSIZE; b++)
            dest->b_mg[b].data = temp.b_mg[b].data;
        
        break;

    case BASE_MG:
        
        rns4096_init(&_a);
        rns4096_copy(a, &_a);
        rns4096_base_copy(&_a, BASE_MG);
        
        // Multiply
        for (int i = 0; i < RNSBASE_SIZE; i++) {
            rns4096_init(&temp);
            rns4096_from_arr(&mat, __RNSBASE_from_MG_matrix[i], BASE_ST);
            rns4096_mul_1mod(&_a, &mat, &temp, i, BASE_ST, i);

            for (int b = 0; b < RNS4096_BLOCKSIZE; b++) {
                _mm256_storeu_si256((__m256i const*)(&(res_raw[i][b * 8])), temp.b_st[b].data);
            }
        };
        //Add vectors to get new terms
        
        for (int i = 0; i < RNSBASE_SIZE; i++) {
            rns4096_init(&temp);
            rns4096_init(&temp2);

            for (int b = 0; b < RNS4096_BLOCKSIZE; b++) {
                uint32_t vec[8];
                for (int j = 0; j < 8; j++)
                    vec[j] = res_raw[b * 8 + j][i];

                temp.b_st[b].data = _mm256_loadu_si256(vec);
            }

      
            rns4096_add(&res, &temp, &res, BASE_ST);

        };


        //Get correction vector 
        rns4096_init(&temp);

        _rns256_block_init(&fact_block, 0);
        fact_block.data = _mm256_set1_epi32(fact);

        for (int b = 0; b < RNS4096_BLOCKSIZE; b++) {
            
            _rns256_block_init(&lookup_block, b);
            lookup_block.data = _mm256_loadu_si256(&__RNSBASE_MG_Mneg_mod_mi[b * 8]);
            _rns256_block_mul(&lookup_block, &fact_block, &(corr_vec.b_st[b]), b, BASE_ST);
        }

        rns4096_init(&temp);
        //Correction step
        rns4096_add(&res, &corr_vec, &temp, BASE_ST);

        //Load result
        for (int b = 0; b < RNS4096_BLOCKSIZE; b++)
            dest->b_st[b].data = temp.b_st[b].data;

        break;
    };
}


void rns4096_mred(_rns4096_uint* a, _rns4096_modulo* p, _rns4096_uint* dest) {
    _rns4096_uint _q, _theta, _r, _s, _sigma, _x;
    rns4096_init(&_q);
    rns4096_init(&_theta);
    rns4096_init(&_r);
    rns4096_init(&_s);
    rns4096_init(&_sigma);
    rns4096_init(&_x);
    uint32_t corr_fact;

    //Step 0
    rns4096_mred_normalize(a, p, &_x, BASE_ST);
    rns4096_mred_normalize(a, p, &_x, BASE_MG);

    //Step 1
    rns4096_modulo_mul(&_x, p, &_q, BASE_ST);

    //Step 2
    rns4096_from_arr(&_theta, __RNSBASE_Mi_inv, BASE_ST);
    rns4096_mul(&_q, &_theta, &_theta, BASE_ST);

    //Step 3
    corr_fact = rns4096_get_corr_factor(&_theta, BASE_ST, CORRECTION_ROFFSET, CORRECTION_LOFFSET);
    rns4096_base_conversion(&_theta, &_q, BASE_ST, corr_fact);

    //Step 4 
    rns4096_modulo_mul(&_q, p, &_r, BASE_MG);
    rns4096_add(&_x, &_r, &_r, BASE_MG);
    

    //Step 5
    rns4096_from_arr(&_s, __RNSBASE_MG_M_inv, BASE_MG);
    rns4096_mul(&_s, &_r, &_s, BASE_MG);
    rns4096_copy1(&_s, dest, BASE_MG);

    //Step 6
    rns4096_from_arr(&_sigma, __RNSBASE_MG_Mi_inv_mg, BASE_MG);
    rns4096_mul(&_s, &_sigma, &_sigma, BASE_MG);


    //Step 7
    corr_fact = rns4096_get_corr_factor_alpha(&_sigma, BASE_MG, CORRECTION_ROFFSET, CORRECTION_LOFFSET, 0.5);
    rns4096_base_conversion(&_sigma, dest, BASE_MG, corr_fact);
}

void rns4096_MM(_rns4096_uint* a, _rns4096_uint* b, _rns4096_modulo* p, _rns4096_uint* dest) {
    _rns4096_uint _q, _theta, _r, _s, _sigma,_x, _xt;
    rns4096_init(&_q);
    rns4096_init(&_theta);
    rns4096_init(&_r);
    rns4096_init(&_s);
    rns4096_init(&_x);
    rns4096_init(&_xt);
    rns4096_init(&_sigma);
    uint32_t corr_fact;

    //Multiply
    rns4096_mul(a, b, &_xt, BASE_ST);
    rns4096_mul(a, b, &_xt, BASE_MG);

    //Step 0
    rns4096_mred_normalize(&_xt, p, &_x, BASE_ST);
    rns4096_mred_normalize(&_xt, p, &_x, BASE_MG);


    //Step 1
    rns4096_modulo_mul(&_x, p, &_q, BASE_ST);


    //Step 2
    rns4096_from_arr(&_theta, __RNSBASE_Mi_inv, BASE_ST);
    rns4096_mul(&_q, &_theta, &_theta, BASE_ST);


    //Step 3
    corr_fact = rns4096_get_corr_factor(&_theta, BASE_ST, CORRECTION_ROFFSET, CORRECTION_LOFFSET);
    rns4096_base_conversion(&_theta, &_q, BASE_ST, corr_fact);


    //Step 4 
    rns4096_modulo_mul(&_q, p, &_r, BASE_MG);
    rns4096_add(&_x, &_r, &_r, BASE_MG);


    //Step 5
    rns4096_from_arr(&_s, __RNSBASE_MG_M_inv, BASE_MG);
    rns4096_mul(&_s, &_r, &_s, BASE_MG);
    rns4096_copy1(&_s, dest, BASE_MG);

    //Step 6
    rns4096_from_arr(&_sigma, __RNSBASE_MG_Mi_inv_mg, BASE_MG);
    rns4096_mul(&_s, &_sigma, &_sigma, BASE_MG);


    //Step 7
    corr_fact = rns4096_get_corr_factor_alpha(&_sigma, BASE_MG, CORRECTION_ROFFSET, CORRECTION_LOFFSET, 0.5);
    rns4096_base_conversion(&_sigma, dest, BASE_MG, corr_fact);

    
}



void rns4096_MExp(_rns4096_uint* a, uint32_t e, _rns4096_modulo* p, _rns4096_uint* dest) {
    _rns4096_uint _r;
    rns4096_init(&_r);
    rns4096_from_dec(&_r, "1", 1);
    for (int i = floor(log2(e)); i >= 0; i--) {
        rns4096_MM(&_r, &_r, p, &_r);
        if ((e >> i) & 1) {
            rns4096_MM(&_r, a, p, &_r);
        }
    }
    rns4096_copy(&_r, dest);
}

void rns4096_MExp_str(_rns4096_uint* a, char* ebuf, _rns4096_modulo* p, _rns4096_uint* dest) {
    mpz_t e;
    mpz_init(e);
    mpz_set_str(e, ebuf, 10);
    size_t len = mpz_sizeinbase(e,2);
    char* bin_e = (char*)malloc(sizeof(char) * len);
    
    mpz_get_str(bin_e, 2, e);

    clock_t begin = clock();
    _rns4096_uint _r;
    rns4096_init(&_r);
    rns4096_from_dec(&_r, "1", 1);
    for (int i = 0; i < len; i++) {

        rns4096_MM(&_r, &_r, p, &_r);

        if (bin_e[i] == '1') {
            rns4096_MM(&_r, a, p, &_r);

        }
        
    }
    rns4096_copy(&_r, dest);
    clock_t end = clock();
    printf("rns exp time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
}