#include "RSA.h"
#include <strsafe.h>
#include <time.h>

void RSA_keygen(char* p_buf, char* q_buf, char* N_buf, char* d_buf, int bitc) {
	mpz_t p, q, N, d, phi, e;

	gmp_randstate_t state;
	gmp_randinit_mt(state);

	mpz_init(p);
	mpz_init(q);
	mpz_init(N);
	mpz_init(d);
	mpz_init(phi);
	mpz_init_set_ui(e, RSA_PUB_EXP);

	mpz_urandomb(p, state, bitc / 2);
	mpz_nextprime(p, p);
	gmp_snprintf(p_buf, RSA_BIT_SIZE / 2 + 1,"%Zd",p);

	mpz_urandomb(q, state, bitc / 2);
	mpz_nextprime(q, q);
	gmp_snprintf(q_buf, RSA_BIT_SIZE / 2 + 1, "%Zd", q);

	mpz_mul(N, p, q);
	gmp_snprintf(N_buf, RSA_BIT_SIZE + 1, "%Zd", N);
	gmp_printf("N: %Zd\n", N);

	mpz_sub_ui(p, p, 1);
	mpz_sub_ui(q, q, 1);

	mpz_mul(phi, p, q);
	mpz_invert(d, e, phi);
	gmp_snprintf(d_buf, RSA_BIT_SIZE + 1, "%Zd", d);
	
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(N);
	mpz_clear(d);
	mpz_clear(phi);
	mpz_clear(e);
}

void RSA_enc_mpz(char* Nbuf, char* msg, size_t msgsize, char* out, int base) {
	clock_t begin = clock();
	mpz_t N, msgi, e;
	mpz_init(N);
	mpz_init(msgi);
	mpz_init_set_ui(e, RSA_PUB_EXP);
	mpz_set_str(N, Nbuf, 10);
	char* msghex;
	if (base == 16)
		msghex = msg;
	else if (base == 10){
		msghex = (char*)malloc(2 * msgsize * sizeof(char));
		dec_to_hex(msg, msghex);
	}

	mpz_set_str(msgi, msghex, 16);
	mpz_powm(msgi, msgi, e, N);
	clock_t end = clock();
	printf("mpz enc time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
	gmp_snprintf(out, (RSA_BIT_SIZE) + 10, "%Zx", msgi);

}

void RSA_dec_mpz(char* Nbuf, char* dbuf, char* msg, size_t outsize, char* out, int base) {
	clock_t begin = clock();
	mpz_t N, msgi, d;
	mpz_init(N);
	mpz_init(msgi);
	mpz_init(d);
	mpz_set_str(d, dbuf, 10);
	mpz_set_str(N, Nbuf, 10);

	char* msghex;
	if (base == 16)
		msghex = msg;
	else
		return;
	mpz_set_str(msgi, msghex, 16);
	mpz_powm(msgi, msgi, d, N);
	clock_t end = clock();
	printf("mpz dec time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
	gmp_snprintf(out, (RSA_BIT_SIZE / 8 + 2), "%Zx", msgi);
}

void RSA_enc_rns4096(char* Nbuf, char* msg, char* out, int base) {
	clock_t begin = clock();
	_rns4096_uint m, c;
	_rns4096_modulo N;
	mpz_t res;
	mpz_init(res);
	rns4096_modulo_init(&N);
	rns4096_modulo_from_dec(&N, Nbuf, RSA_BIT_SIZE);

	rns4096_init(&c);
	rns4096_init(&m);
	if (base == 16)
		rns4096_from_hex(&m, msg, RSA_BIT_SIZE);
	else if (base == 10)
		rns4096_from_dec(&m, msg, RSA_BIT_SIZE);
	else return;

	rns4096_MExp(&m, RSA_PUB_EXP, &N, &c);
	
	rns4096_to_dec_mod(&c,Nbuf, res, BASE_ST);
	clock_t end = clock();
	printf("rns enc time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
	gmp_snprintf(out, RSA_BIT_SIZE, "%Zx", res);
}

void RSA_dec_rns4096(char* Nbuf, char* dbuf, char* msg, char* out, int base) {
	clock_t begin = clock();
	_rns4096_uint m, c;
	_rns4096_modulo N;
	mpz_t res;
	mpz_init(res);
	rns4096_modulo_init(&N);
	rns4096_modulo_from_dec(&N, Nbuf, RSA_BIT_SIZE);

	rns4096_init(&c);
	rns4096_init(&m);
	if (base == 16)
		rns4096_from_hex(&m, msg, RSA_BIT_SIZE);
	else if (base == 10)
		rns4096_from_dec(&m, msg, RSA_BIT_SIZE);
	else return;

	rns4096_MExp_str(&m, dbuf, &N, &c);
	rns4096_to_dec_mod(&c,Nbuf, res, BASE_ST);
	clock_t end = clock();
	printf("rns dec time: %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
	gmp_snprintf(out, RSA_BIT_SIZE, "%Zx", res);
}