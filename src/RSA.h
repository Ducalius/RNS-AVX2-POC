#pragma once
#include "montgomery.h"
#include <gmp.h>
#include "consts.h"
#include "rnsuinttype.h"
#include "rnsmodulotype.h"
#include <openssl/rsa.h>

#define RSA_PUB_EXP 65537
#define RSA_BIT_SIZE 2047

void RSA_keygen(char* p_buf, char* q_buf, char* N_buf, char* d_buf, int bitc);

void RSA_enc_mpz(char* Nbuf, char* msg, size_t msgsize, char* out, int base);

void RSA_dec_mpz(char* Nbuf, char* dbuf, char* msg, size_t outsize, char* out, int base);

void RSA_enc_rns4096(char* Nbuf, char* msg, char* out, int base);

void RSA_dec_rns4096(char* Nbuf, char* dbuf, char* msg, char* out, int base);
