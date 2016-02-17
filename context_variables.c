#include "context_variables.h" 


/* basic constructor for cb_context */
cb_context context_construct(long n_Max, mpfr_prec_t prec, int lambda){
    cb_context result;
    result.n_Max = n_Max;
    result.prec = prec;
    result.rnd = MPFR_RNDN;
    result.rho_to_z_matrix = compute_rho_to_z_matrix(lambda,prec);
    
    result.lambda = lambda;
    mpfr_init2(result.rho, prec);
    mpfr_set_ui(result.rho,8, MPFR_RNDN);
    mpfr_sqrt(result.rho,result.rho,MPFR_RNDN);
    mpfr_ui_sub(result.rho,3,result.rho,MPFR_RNDN); 
    return result;
}

void clear_cb_context(cb_context context){
	for(int i=0;i<=context.lambda;i++){
		for(int j=0;j<=context.lambda;j++){
			mpfr_clear(context.rho_to_z_matrix[i*(context.lambda+1)+j]);
		}
	}
	free(context.rho_to_z_matrix);
	mpfr_clear(context.rho);
}

mpfr_t* compute_rho_to_z_matrix(unsigned long Lambda_arg, long prec){
	/* To avoid writing lambda + 1 so many times...*/
	unsigned long Lambda=Lambda_arg+1;
	mpfr_t* temps=malloc(sizeof(mpfr_t)*(Lambda));
	mpfr_init2(temps[0],prec);
	mpfr_set_ui(temps[0],8,MPFR_RNDN);
	mpfr_sqrt(temps[0],temps[0],MPFR_RNDN);
	mpfr_neg(temps[0],temps[0],MPFR_RNDN);
	for(unsigned long j=1;j<Lambda;j++){
		mpfr_init2(temps[j],prec);
		mpfr_mul_si(temps[j],temps[j-1],2*j-3,MPFR_RNDN);
		mpfr_div_ui(temps[j],temps[j],j,MPFR_RNDN); 
	}
	mpfr_sub_ui(temps[1],temps[1],2,MPFR_RNDN);
	mpfr_add_ui(temps[0],temps[0],3,MPFR_RNDN);
	mpfr_t temp;
	mpfr_init2(temp,prec);
	mpfr_t temp2;
	mpfr_init2(temp2,prec);

	mpfr_t* result=malloc(sizeof(mpfr_t)*(Lambda)*(Lambda));
	mpfr_init2(result[0],prec);
	mpfr_set_ui(result[0],1,MPFR_RNDN);
	for(unsigned long j=1; j<(Lambda*Lambda); j++){
		mpfr_init2(result[j],prec);
		mpfr_set_zero(result[j],1);
	}
	for(unsigned long j=1;j<Lambda;j++){
		mpfr_set_ui(temp,1,MPFR_RNDN);
		for(unsigned long k=0;k<=j;k++){
			mpfr_mul(temp2,temps[j-k],temp,MPFR_RNDN);
			mpfr_add(result[j+Lambda],result[j+Lambda],temp2,MPFR_RNDN);
			mpfr_mul_si(temp,temp,-2,MPFR_RNDN); 
		} 
	}
	for(unsigned long i=2;i<Lambda;i++){
		for(unsigned long j=1;j<Lambda;j++){
			for(unsigned long k=i-1;k<Lambda-j;k++){
				mpfr_mul(temp,result[Lambda*(i-1)+k],result[j+Lambda],MPFR_RNDN);
				mpfr_add(result[Lambda*i+k+j],result[Lambda*i+k+j],temp,MPFR_RNDN);
			}

		} 
	} 

	/* transposition */
	for(unsigned long i=0;i<Lambda;i++){
		for(unsigned long j=0;j<i;j++){
		mpfr_swap(result[i+Lambda*j],result[j+Lambda*i]);
		}
	}
	for(unsigned long j=0;j<Lambda;j++){
		mpfr_clear(temps[j]);
	}
	free(temps);
	mpfr_clear(temp);
	mpfr_clear(temp2);
	return result; 
}
