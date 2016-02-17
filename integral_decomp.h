#include "stdlib.h"
#include "stdio.h"
#include "mpfr.h"

mpfr_t* simple_pole_case_c(long pole_order_max, mpfr_t base, mpfr_t pole_position, mpfr_t incomplete_gamma_factor, mpfr_prec_t prec);
mpfr_t* double_pole_case_c(long pole_order_max, mpfr_t base, mpfr_t pole_position, mpfr_t incomplete_gamma_factor, mpfr_prec_t prec);
