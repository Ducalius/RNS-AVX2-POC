#include "rnsblocktype.h"
#include <string.h>

void print_block(_rns256_block b) {
    uint32_t output[8];
    _mm256_storeu_epi32(output, b.data);
    for(int i = 0; i < 8; i++)
        printf("%8x ",output[i]);
    printf("\n");
}

void print_block_item(_rns256_block b, int i) {
    uint32_t output[8];
    _mm256_storeu_epi32(output, b.data);
    printf("%8u ", output[i]);
    printf("\n");
}

void print_block_diff(_rns256_block b, _rns256_block c) {
    uint32_t output1[8];
    uint32_t output2[8];
    _mm256_storeu_epi32(output1, b.data);
    _mm256_storeu_epi32(output2, c.data);
    for (int i = 0; i < 8; i++)
        printf("%10d ", output1[i] - output2[i]);
    printf("\n");
}

void print_block_dec(_rns256_block b) {
    uint32_t output[8];
    _mm256_storeu_epi32(output, b.data);
    for (int i = 0; i < 8; i++)
        printf("%10u ", output[i]);
    printf("\n");
}

void print_block_dec_format(_rns256_block b) {
    uint32_t output[8];
    _mm256_storeu_epi32(output, b.data);
    for (int i = 0; i < 8; i++)
        printf("%10u, ", output[i]);
    printf("\n");
}

void print_vec(__m256i* v) {
    uint32_t output[8];
    
    _mm256_storeu_epi32(output, *v);
    for (int i = 0; i < 8; i++)
        printf("%8x ", output[i]);
    printf("\n");
}

void print_vec_dec(__m256i* v) {
    uint32_t output[8];

    _mm256_storeu_epi32(output, *v);
    for (int i = 0; i < 8; i++)
        printf("%10u ", output[i]);
    printf("\n");
}

void print_vec_test(uint32_t* v) {
    for(int i = 0; i < 8; i++)
        printf("%8x ",v[i]);
    printf("\n");
}

void print_vec_test_dec(uint32_t* v) {
    for (int i = 0; i < 8; i++)
        printf("%10u ", v[i]);
    printf("\n");
}

void _rns256_block_init(_rns256_block* b, int id) {
    b->id = id;
    b->data = _mm256_setzero_si256();
}

void _rns256_block_add(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int id, __RNS_BASE base){
    uint32_t* bp;
    if (a->id != b->id) {
        fprintf(stderr, "Operands are not the same base. %d %d \n", a->id, b->id);
        return;
    }
    switch(base) {
    case BASE_ST: 
        bp = &__RNSBASE_mi[id * 8];
        break;
    case BASE_MG:
        bp = &__RNSBASE_MG_mi[id * 8];
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };

    /*
    __m256i m = _mm256_loadu_si256(bp);

    __m256i c = _mm256_add_epi32(a->data, b->data);
    __m256i cb = _mm256_cmpgt_epi32(c,a->data);
    __m256i d = _mm256_and_si256(cb,m);
    c = _mm256_sub_epi32(c,d);

    dest->data = c;
    */

    __m256i m = _mm256_loadu_si256(bp);

    __m256i c = _mm256_add_epi32(a->data, b->data);
    __m256i opm = _mm256_max_epu32(c, a->data);
    __m256i cb = _mm256_cmpeq_epi32(a->data, opm);

    __m256i d = _mm256_and_si256(cb, m);

    c = _mm256_sub_epi32(c, d);


    dest->data = c;
}

void __rns256_block_add_debug(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int id, __RNS_BASE base) {
    uint32_t* bp;
    if (a->id != b->id) {
        fprintf(stderr, "Operands are not the same base. %d %d \n", a->id, b->id);
        return;
    }
    switch (base) {
    case BASE_ST:
        bp = &__RNSBASE_mi[id * 8];
        break;
    case BASE_MG:
        bp = &__RNSBASE_MG_mi[id * 8];
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };

    __m256i m = _mm256_loadu_si256(bp);

    __m256i c = _mm256_add_epi32(a->data, b->data);
    print_vec_dec(&(a->data));
    print_vec_dec(&(b->data));
    print_vec_dec(&c);
    //__m256i cba = _mm256_cmpgt_epi32(c, a->data);
    //__m256i cbb = _mm256_cmpgt_epi32(c, b->data);
    __m256i opm = _mm256_max_epu32(c, a->data);
    print_vec_dec(&(a->data));
    __m256i cb = _mm256_cmpeq_epi32(a->data, opm);
    //__m256i cb = _mm256_or_si256(cbb, cba);
    print_vec_dec(&cb);
    __m256i d = _mm256_and_si256(cb, m);
    print_vec_dec(&d);
    c = _mm256_sub_epi32(c, d);
    print_vec_dec(&c);

    dest->data = c;
}

void _rns256_block_mul(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int id, __RNS_BASE base) {
    uint32_t* mup;
    switch (base) {
    case BASE_ST:
        mup = &__RNSBASE_ci[id * 8];
        break;
    case BASE_MG:
        mup = &__RNSBASE_MG_ci[id * 8];
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };
    __m256i mu = _mm256_loadu_si256(mup);
    //print_vec(&mu);
    __m256i mus = _mm256_shuffle_epi32(mu, 0xb1);

    __m256i d0 = _mm256_mul_epu32(a->data, b->data);

    __m256i as = _mm256_shuffle_epi32(a->data, 0xb1);
    __m256i bs = _mm256_shuffle_epi32(b->data, 0xb1);
    __m256i d1 = _mm256_mul_epu32(as, bs);

    for (int i = 0; i < 3; i++) {

        d0 = _mm256_shuffle_epi32(d0, 0xb1);
        d1 = _mm256_shuffle_epi32(d1, 0xb1);

        __m256i e0 = _mm256_mul_epu32(d0, mu);
        __m256i e1 = _mm256_mul_epu32(d1, mus);

        __m256i f0 = _mm256_srli_epi64(d0, 32);
        __m256i f1 = _mm256_srli_epi64(d1, 32);

        d0 = _mm256_add_epi64(e0, f0);
        d1 = _mm256_add_epi64(e1, f1);
    }
    d1 = _mm256_slli_epi64(d1, 32);
    dest->data = _mm256_xor_epi32(d0, d1);
}

void __rns256_block_mul_debug(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int id, __RNS_BASE base) {
    uint32_t* mup;
    switch (base) {
    case BASE_ST:
        mup = &__RNSBASE_ci[id * 8];
        break;
    case BASE_MG:
        mup = &__RNSBASE_MG_ci[id * 8];
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };
    __m256i mu = _mm256_loadu_si256(mup);
    
    __m256i mus = _mm256_shuffle_epi32(mu, 0xb1);
    
    __m256i d0 = _mm256_mul_epu32(a->data, b->data);
    __m256i as = _mm256_shuffle_epi32(a->data, 0xb1);
    __m256i bs = _mm256_shuffle_epi32(b->data, 0xb1);
    __m256i d1 = _mm256_mul_epu32(as, bs);
    printf("=======================d vecs:\n");
    print_vec_dec(&d0);
    print_vec_dec(&d1);
    
    for (int i = 0; i < 3; i++) {
        printf("=======================iteration   %d  =======================\n", i);
        printf("=======================d vecs bef:\n");
        
        d0 = _mm256_shuffle_epi32(d0, 0xb1);
        d1 = _mm256_shuffle_epi32(d1, 0xb1);
        print_vec(&d0);
        print_vec(&d1);
        __m256i e0 = _mm256_mul_epu32(d0, mu);
        __m256i e1 = _mm256_mul_epu32(d1, mus);
        printf("=======================e vecs:\n");
        print_vec(&e0);
        print_vec(&e1);
        __m256i f0 = _mm256_srli_epi64(d0, 32);
        __m256i f1 = _mm256_srli_epi64(d1, 32);
        printf("=======================f vecs:\n");
        print_vec(&f0);
        print_vec(&f1);
        d0 = _mm256_add_epi64(e0, f0);
        d1 = _mm256_add_epi64(e1, f1);
        printf("=======================d vecs:\n");
        print_vec(&d0);
        print_vec(&d1);
    }
    d1 = _mm256_slli_epi64(d1, 32);
    
    dest->data = _mm256_xor_epi32(d0, d1);
}

void _rns256_block_mul_1mod(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int mod_id, __RNS_BASE base) {
    uint32_t* mup;
    switch (base) {
    case BASE_ST:
        mup = __RNSBASE_ci[mod_id];
        break;
    case BASE_MG:
        mup = __RNSBASE_MG_ci[mod_id];
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };
    __m256i mu = _mm256_set1_epi32(mup);
    //if(mod_id == 81) print_vec(&mu);
    __m256i mus = _mm256_shuffle_epi32(mu, 0xb1);
    __m256i as = a->data;
    __m256i bs = b->data;

    __m256i d0 = _mm256_mul_epu32(as, bs);

    as = _mm256_shuffle_epi32(as, 0xb1);
    bs = _mm256_shuffle_epi32(bs, 0xb1);
    __m256i d1 = _mm256_mul_epu32(as, bs);

    for (int i = 0; i < 3; i++) {

        d0 = _mm256_shuffle_epi32(d0, 0xb1);
        d1 = _mm256_shuffle_epi32(d1, 0xb1);

        __m256i e0 = _mm256_mul_epu32(d0, mu);
        __m256i e1 = _mm256_mul_epu32(d1, mus);

        __m256i f0 = _mm256_srli_epi64(d0, 32);
        __m256i f1 = _mm256_srli_epi64(d1, 32);

        d0 = _mm256_add_epi64(e0, f0);
        d1 = _mm256_add_epi64(e1, f1);
    }
    d1 = _mm256_slli_epi64(d1, 32);
    dest->data = _mm256_xor_epi32(d0, d1);
}


void __rns256_block_mul_1mod_debug(_rns256_block* a, _rns256_block* b, _rns256_block* dest, int mod_id, __RNS_BASE base) {
    uint32_t* mup;
    switch (base) {
    case BASE_ST:
        mup = __RNSBASE_ci[mod_id];
        break;
    case BASE_MG:
        mup = __RNSBASE_MG_ci[mod_id];
        break;
    default:
        fprintf(stderr, "Base selection error \n");
        return;
    };
    __m256i mu = _mm256_set1_epi32(mup);
    printf("block 10, debug\n");
    __m256i mus = _mm256_shuffle_epi32(mu, 0xb1);
    print_block_dec(*a);
    print_block_dec(*b);
    print_vec(&mu);
    printf("mul\n");
    __m256i d0 = _mm256_mul_epu32(a->data, b->data);

    __m256i as = _mm256_shuffle_epi32(a->data, 0xb1);
    __m256i bs = _mm256_shuffle_epi32(b->data, 0xb1);
    __m256i d1 = _mm256_mul_epu32(as, bs);

    for (int i = 0; i < 3; i++) {

        d0 = _mm256_shuffle_epi32(d0, 0xb1);
        d1 = _mm256_shuffle_epi32(d1, 0xb1);

        __m256i e0 = _mm256_mul_epu32(d0, mu);
        __m256i e1 = _mm256_mul_epu32(d1, mus);

        __m256i f0 = _mm256_srli_epi64(d0, 32);
        __m256i f1 = _mm256_srli_epi64(d1, 32);

        d0 = _mm256_add_epi64(e0, f0);
        d1 = _mm256_add_epi64(e1, f1);
    }
    d1 = _mm256_slli_epi64(d1, 32);
    printf("res\n");
    dest->data = _mm256_xor_epi32(d0, d1);
    print_block_dec(*dest);
}

void _rns256_block_trunc(_rns256_block* a, _rns256_block* dest, uint32_t roffset, uint32_t loffset) {
    
    //print_block(*a);
    __m256i res = a->data;
    //print_vec(&res);
    res = _mm256_srli_epi32(res, roffset);
    //print_vec(&res);
    res = _mm256_slli_epi32(res, loffset);
    //print_vec(&res);
    //printf("\n");
    dest->data = res;
}

void _rns256_block_swap(_rns256_block* a, _rns256_block* b, int apos, int bpos) {
    return;
}
