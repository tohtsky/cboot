#include "stdlib.h"
#include "mpfr.h"

#ifndef CONFORMAL_BOOTSTRAP_CONTEXT
#define CONFORMAL_BOOTSTRAP_CONTEXT
mpfr_t* compute_rho_to_z_matrix(unsigned long Lambda_arg, long prec);
typedef struct _context {
	long n_Max;
	mpfr_prec_t prec;
	mp_rnd_t rnd;
	int lambda;
	mpfr_t* rho_to_z_matrix;	
	mpfr_t rho;
} cb_context; 

cb_context context_construct(long n_Max, mpfr_prec_t prec, int lambda); 
void clear_cb_context(cb_context context);
#endif
