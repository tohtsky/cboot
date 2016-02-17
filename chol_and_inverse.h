#include "mpfr.h"

mpfr_t* mpfr_triangular_inverse(mpfr_t* A, int dim,mpfr_prec_t prec);

mpfr_t* mpfr_cholesky(mpfr_t* A, int dim,mpfr_prec_t prec);

mpfr_t* form_anti_band(mpfr_t* ab_vector, int dim, mpfr_prec_t prec);
