#include "stdlib.h"
#include "stdio.h"
#include "mpfr.h"
mpfr_t* fast_partial_fraction_c(mpfr_t* pole_locations, int* double_or_single, int expected_result_length, mpfr_prec_t prec);
