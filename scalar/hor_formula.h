#include "stdlib.h"
#include "mpfr.h" 

#ifndef __RHOREC_MPFR_2_
#define __RHOREC_MPFR_2_
/*
 * Hogervorst-Osborn-Rycykov recursion relation takes the form 
 * 
 *  \sum_{i=0}^{7} b[n-i] *  a[i](n,parameters (epsilon, spin ,Delta, S, P)) =0
 *
 * and each coefficients a[i] is like
 *
 *  a[i][n,params] = a[i,4](params)*n^4 + a[i,3](params)*n*^3 + ... + a[i,0](params)
 *
 * These parameters are multi-precision objects and we want to minimize the number of operations among them. Our trick is to pre-compute a[i,j] (i=0 to 7)(j=0 to 4). Since they do not include n, we can re-use them at every stage of recursion process (that is, computation of next b[n]).
 *  
 * The process for obtaining a[i,j] is implemented in 
 * set_zero_spin_rec_coeffs (for ell >0) 
 * and
 * set_zero_spin_rec_coeffs 
 *
 * */

void set_nonzero_spin_rec_coeffs(mpfr_t result[8][5], mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd); 
void set_zero_spin_rec_coeffs(mpfr_t result[6][4], mpfr_t epsilon, mpfr_t Delta, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd);

void set_nonzero_spin_rec_coeffs_deriv(mpfr_t result[8][5], mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd); 
void set_zero_spin_rec_coeffs_deriv(mpfr_t result[6][4], mpfr_t epsilon, mpfr_t Delta, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd);

void initialize_spin_nonzero_coeffs_folder(mpfr_t a[8][5], mpfr_prec_t prec);
void initialize_spin_zero_coeffs_folder(mpfr_t a[6][4], mpfr_prec_t prec);

void deallocate_spin_nonzero_coeffs_folder(mpfr_t a[8][5]);
void deallocate_spin_zero_coeffs_folder(mpfr_t a[6][4]);

void spin_nonzero_evaluate_at_n(mpfr_t a[8], mpfr_t rho[8][5], long n, mpfr_prec_t prec, mp_rnd_t rnd);
void spin_zero_evaluate_at_n(mpfr_t a[6], mpfr_t rho[6][4], long n, mpfr_prec_t prec, mp_rnd_t rnd);
#endif

