#include "mpfr.h" 
#include "context_variables.h"
#include "hor_recursion.h"
#include <stdlib.h>
#include <stdio.h>

#ifndef CONFORMAL_BLOCK_KTYPE
#define CONFORMAL_BLOCK_KTYPE
void set_k_coeffs(mpfr_t a[4][3], mpfr_t h, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd);
//void initialize_k_coeffs(mpfr_t a[4][3], mpfr_prec_t prec);
//void deallocate_k_coeffs(mpfr_t a[4][3], mpfr_prec_t prec);
void k_evaluate_at_n(mpfr_t a[4], mpfr_t[4][3],long n, mpfr_prec_t prec, mp_rnd_t rnd);
mpfr_t* recursion_k(unsigned long nMax, mpfr_t h, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd);

mpfr_t* k_table_c(mpfr_t h,mpfr_t S, mpfr_t P, cb_context context);
mpfr_t* chiral_h_times_rho_to_n_c(unsigned long n, mpfr_t h,mpfr_t S, mpfr_t P, cb_context context);
mpfr_t* chiral_h_asymptotic_c(mpfr_t S, cb_context context);

#endif
