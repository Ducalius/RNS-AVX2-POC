/* Construct a 256-bit vector from 4 64-bit doubles. Add it to itself
 * and print the result.
 */

#include <stdio.h>
#include <immintrin.h>
#include <stdint.h>
#include <gmp.h>
#include <math.h>
#include "consts.h"
#include "rnsuinttype.h"
#include "rnsmodulotype.h"
#include <time.h>
#include <stdlib.h>
#include "rsa.h"

int main() {

    __RNSBASE_init_consts();
    char as[] = "11af39d123ad42f3123c41b9";
    char bs[] = "866af481d95791ec0ff8123b";
    char as_dec[] = "7b8b071a2e4b88b";
    char bs_dec[] = "41600292105652342275874230843";
    char t_dec[] = "15921986908262798289";
    char p_dec[] = "1443492886803489649";
    //printf("%d_\n", mod_invert(32, 11));

    _rns4096_uint nl;
    _rns4096_uint nl2;
    _rns4096_uint ttest;
    _rns4096_uint res;
    _rns4096_modulo p;

    rns4096_modulo_init(&p);
    rns4096_modulo_from_dec(&p, p_dec, RSA_BIT_SIZE);
    rns4096_init(&nl);
    rns4096_from_hex(&nl, as_dec, 28);
    rns4096_init(&nl2);
    rns4096_from_dec(&nl2, bs_dec, 29);
    rns4096_init(&res);
    rns4096_init(&ttest);
    rns4096_from_dec(&ttest, t_dec, 20);
    
  
    //print_block(nl.b_st[0]);
    //print_block(nl2.b_st[0]);
    //rns4096_modulo_from_hex(&p,as,24);

    /*
    __m256i t = _mm256_setr_epi32(8,7,6,5,4,3,2,1);
    t = _mm256_shuffle_epi32(t, 0xb1);
    print_vec(&t);
    */
    //print_vec(&r);



    
    //_rns256_block d;
    //_rns256_block_init(&d,0);
    //_rns256_block_mul(&(nl.b_st[0]),&(nl2.b_st[0]),&d,0,BASE_ST);
    _rns4096_uint _q, _theta, _r, _s, _sigma;
    rns4096_init(&_q);
    rns4096_init(&_theta);
    rns4096_init(&_r);
    rns4096_init(&_s);
    rns4096_init(&_sigma);
    uint32_t corr_fact;


    
    /* MG -> ST */
    /*
    rns4096_from_arr(&_theta, __RNSBASE_MG_Mi_inv_mg, BASE_MG);
    rns4096_mul(&nl, &_theta, &_theta, BASE_MG);
    print_block(_theta.b_mg[0]);

    //Step 3
    corr_fact = rns4096_get_corr_factor(&_theta, BASE_MG, CORRECTION_ROFFSET, CORRECTION_LOFFSET);
    printf("correction factor: %u\n", corr_fact);
    printf("correction factor test: %u\n", get_corr_factor_long(&_theta, BASE_MG));
    
    rns4096_base_conversion(&_theta, &res, BASE_MG, corr_fact);
    
    printf("test:\n");
    print_block(nl.b_st[0]);
    //print_block(nl.b_mg[0]);
    //rns4096_mred(&nl, &p, &res);
    printf("res:\n");
    print_block(res.b_st[0]);
    
    printf("compare blocks:\n");
    for (int i = 0; i < 16; i++) {
        printf("%d:\n", i);
        print_block(nl.b_st[i]);
        print_block(res.b_st[i]);
        print_block_diff(nl.b_st[i], res.b_st[i]);
    }
    */
    


    rns4096_MExp(&nl, "56", &p, &res);
    //rns4096_MExp(&ttest, 2, &p, &res);
    /*
    rns4096_mul(&ttest, &ttest, &ttest, BASE_ST);
    rns4096_mul(&ttest, &ttest, &ttest, BASE_MG);
    rns4096_mred(&ttest, &p, &res);
    */
    //printf("cool results: %%%%%%%%%%%%%%%%%%%%%%%%%5\n");
    mpz_t resdec;
    mpz_init(resdec);
    rns4096_to_dec(&res, resdec, BASE_ST);
    char buf[8192];
    //mpz_out_str(stdout,10,resdec);
    //gmp_printf("result: %Zd\n", resdec);
    //printf("%s\n", buf);
    /*
    for (int i = 0; i < 16; i++) {
        print_block_dec_format(&(res.b_st[i]));
    }
    
    _rns4096_uint vecset1;
    rns4096_init(&vecset1);
    rns4096_from_dec(&vecset1, "1", 1);
    rns4096_mul(&res, &vecset1, &res, BASE_ST);
    rns4096_mred(&res, &p, &res);
    //rns4096_mul(&nl, &nl2, &res, BASE_ST);
    for (int i = 0; i < 16; i++) {
        print_block_dec_format(&(res.b_st[i]));
    }
    _rns4096_uint num, temp, _theta2;
    rns4096_init(&num);
    rns4096_init(&temp);
    rns4096_init(&_theta2);
    mpz_t resdec2, nummpz;
    mpz_init(resdec2);
    mpz_init(nummpz);
    int count = 0;
    for (int i = 0; i < 20; i++) {
        char str[256];
        long long num_dec = rand() % (long long)pow(2,48);
        num_dec = num_dec * 54232;
        long long num_dec2 = rand() % (long long)pow(2, 48);
        num_dec2 = num_dec2 * 54232;
        sprintf(str, "%d%d%d%d", num_dec, num_dec2,num_dec, num_dec2);
        rns4096_from_dec(&num, str, 256);
        
        rns4096_from_arr(&_theta2, __RNSBASE_Mi_inv, BASE_ST);
        rns4096_mul(&num, &_theta2, &_theta2, BASE_ST);
        corr_fact = rns4096_get_corr_factor(&_theta2, BASE_ST, CORRECTION_ROFFSET, CORRECTION_LOFFSET);
        rns4096_base_conversion(&_theta2, &temp, BASE_ST, corr_fact);
        
        rns4096_from_arr(&_theta2, __RNSBASE_MG_Mi_inv_mg, BASE_MG);
        rns4096_mul(&temp, &_theta2, &_theta2, BASE_MG);
        corr_fact = rns4096_get_corr_factor(&_theta2, BASE_MG, CORRECTION_ROFFSET, CORRECTION_LOFFSET);
        rns4096_base_conversion(&_theta2, &temp, BASE_MG, corr_fact);

        rns4096_to_dec(&temp, resdec2, BASE_ST);
        //mpz_set_si(nummpz, num_dec);
        rns4096_to_dec(&num, nummpz, BASE_ST);
        printf("testval: %s\n", str);
        gmp_printf("convresult: %Zd\n", nummpz);
        gmp_printf("result: %Zd\n", resdec2);
        if (mpz_cmp(resdec2, nummpz) == 0) {
            count += 1;
        }
    }
    printf("\n\n$$$$$$$$$$$$$$$$$$$$ test result = %d\n\n", count);
    */
    /* ST -> MG 
    rns4096_from_arr(&_theta, __RNSBASE_Mi_inv, BASE_ST);
    rns4096_mul(&nl, &_theta, &_theta, BASE_ST);
    print_block(_theta.b_st[0]);

    //Step 3
    corr_fact = rns4096_get_corr_factor(&_theta, BASE_ST, CORRECTION_ROFFSET, CORRECTION_LOFFSET);
    printf("correction factor: %u\n", corr_fact);
    printf("correction factor test: %u\n", get_corr_factor_long(&_theta, BASE_ST));

    rns4096_base_conversion(&_theta, &res, BASE_ST, corr_fact);

    printf("testing base ext 1%%%%%%%:\n");
    print_block(nl.b_mg[0]);
    //print_block(nl.b_mg[0]);
    //rns4096_mred(&nl, &p, &res);
    printf("res:\n");
    print_block(res.b_mg[0]);

    printf("compare blocks:\n");
    for (int i = 0; i < 16; i++) {
        //printf("%d:\n", i);
        //print_block(nl.b_mg[i]);
        print_block_dec(res.b_mg[i]);
        //print_block_diff(nl.b_mg[i], res.b_mg[i]);
    }
    printf("bef =\n");
    for (int i = 0; i < 16; i++)
        print_block_dec(nl.b_mg[i]);

    
    
    _rns256_block at, bt, rt, tt;
    _rns256_block_init(&at, 0);
    
    _rns256_block_init(&bt, 0);
    _rns256_block_init(&tt, 0);
    */
    //61
    /* 
    at.data = _mm256_set_epi32(2578992039, 1541926255, 3867230482, 2138179790, 3251436267, 3591926808, 1993078732, 3656345472);
    bt.data = _mm256_set_epi32(2475138070, 1224111614, 1753827002, 704190883, 653331144, 697569530, 404848244, 3593196155);
    */

    //60
    
    /*
    at.data = _mm256_set_epi32(3325573944, 2131637670, 4216536803, 3208995677, 1036704945, 1000360676, 3513344523, 1202744040);
    bt.data = _mm256_set_epi32(3548385716, 3705255311, 3945660840, 3224151607, 2214731214, 2591566200, 2774701469, 2453601432);
    tt.data = _mm256_set_epi32(2578992429, 1541925742, 3867230392, 2138180021, 3251436159, 3591926876, 1993078705, 3656345472);
    */

    //
    /*
    uint32_t tempvec[8] = { 83912483, 3512335608, 1779452436, 2464379211, 2660805789, 1877526636, 3739729191, 229108740 };
    uint32_t xsumvec[8] = { 325103101, 1473361039, 226312823, 940774333, 258764089, 106818172, 2197242611, 1357010493 };
    uint32_t testvec[8] = { 409015584, 690729360, 2005765259, 3405153544, 2919569878, 1984344808, 1642004563, 1586119233 };

    at.data = _mm256_loadu_epi32(tempvec);
    bt.data = _mm256_loadu_epi32(xsumvec);
    tt.data = _mm256_loadu_epi32(testvec);
    */
    /*
    _rns256_block_init(&rt, 0);
    printf("\n");
    __rns256_block_add_debug(&at, &bt, &rt, 0, BASE_MG);
    print_block_dec(rt);
    print_block_dec(tt);
    
    
    //MUL test
    _rns256_block at, bt, rt, tt;
    _rns256_block_init(&at, 0);

    _rns256_block_init(&bt, 0);
    _rns256_block_init(&tt, 0);
    _rns256_block_init(&rt, 0);
    at.data = _mm256_setr_epi32(3335310125, 4070850807, 2238218534, 1494409479, 2078093172, 892053154, 1531806292, 3832503055);
    bt.data = _mm256_setr_epi32(3447428592, 3974847155, 855327748, 3894446523, 3246423782, 1950380532, 3125832992, 2437619545);
    tt.data = _mm256_setr_epi32(1834075761, 668318, 1153224235, 2185915206, 2312605691, 72381581, 1652684862, 3456134352);
    print_block_dec(at);
    print_block_dec(bt);
    printf("\n");
    __rns256_block_mul_debug(&at, &bt, &rt, 10, BASE_MG);
    print_block_dec(rt);
    print_block_dec(tt);
    print_block_diff(rt,tt);
    //printf("compare blocks:\n");
    */
    char pbuf[8192];
    char qbuf[8192];
    char Nbuf[8192];
    char dbuf[8192];

    char plain[8192] = "8008135322800813532280081353228008135322800813532280081353228008135322800813532280081353228008135322";
    char cipher[8192] = "";
    char dec[8192] = "";
    RSA_keygen(pbuf, qbuf, Nbuf, dbuf, RSA_BIT_SIZE);
    RSA_enc_mpz(Nbuf, plain,100,cipher,16,100);
    printf("cipher: %s\n", cipher);
    RSA_dec_mpz(Nbuf, dbuf, cipher, 12, dec, 16);
    printf("dec: %s\n", dec);

    RSA_enc_rns4096(Nbuf, plain, cipher, 16);
    printf("cipher: %s\n", cipher);
    RSA_dec_rns4096(Nbuf, dbuf, cipher, dec, 16);
    printf("dec: %s\n", dec);

    //print_vec_test(__RNSBASE_ci);
    mpz_clear(resdec);

  return 0;
}
