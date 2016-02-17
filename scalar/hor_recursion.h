#include "context_variables.h"
#include "hor_formula.h"
#include "stdio.h"

#ifndef CONFORMAL_BLOCK_COMPUTE
#define CONFORMAL_BLOCK_COMPUTE
mpfr_t* real_axis_result(mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P, cb_context context);
mpfr_t* gBlock_full(mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P, cb_context context);
mpfr_t* hBlock_times_rho_n(unsigned long n, mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P,cb_context context);
mpfr_t* h_asymptotic(mpfr_t epsilon, mpfr_t S, cb_context context);
#endif
