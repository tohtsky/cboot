#include "chol_and_inverse.h"
mpfr_t* mpfr_triangular_inverse(mpfr_t* A, int dim,mpfr_prec_t prec){

	mpfr_t *res=malloc(sizeof(mpfr_t)*dim*dim);

	mpfr_t s;
	mpfr_t mul_temp;

	mpfr_init2(s,prec);
	mpfr_init2(mul_temp,prec); 
/*
 *
 *  \sum_{j=0}^{dim-1} A[i,k] res[k,j] = delta _ {i,j}
 *  e.g. i=0 -> only A[0+0] can contribute, res[0,j]*A[0,0] = delta_{0,j}
 *  
 *  i=1 -> A[1,0]*res[0,k] + A[1,1]*res[1,k]=delta_{1,k}
 *  res[1,k]= 0 for k>1
 *  res[1,1] = A[1,1]^{-1}
 *  res[1,0] = (-A[1,0]*res[0,0] )A[1,1]^{-1} 
 *  
 *  i=2 -> A[2,0]*res[0,k] + A[2,1]*res[1,k] + A[2,2]*res[2,k] = delta(2,k)
 *         k=2 -> ok
 *         k=1 -> A[2,i]*res[0,k] 
 *         
 *         
 *         A[i,0]*res[0,1] + A[i,k] * res[k,j] +...+ A[i,i]res[i,k] = delta(i,k)
 *	   
 *	   A[i,i]*res[i,j]  =  -( A[i,0]*res[0,j] +...+ A[i,k]*res[k,j]
 *	                      + A[i,i-1] * res[i-1,j])
 *	                    = sum_{k= j}^{i-1} A[i,k]*res[k,j]
 *
 * */
	for(int i=0;i<dim;i++){
		for(int j=i+1;j<dim;j++){
			mpfr_init2(res[i*dim+j],prec);
			//mpfr_set_ui(res[i*dim+j],0,MPFR_RNDN);
			mpfr_set_zero(res[i*dim+j],1);
		}
		mpfr_init2(res[i*dim+i],prec);
		mpfr_ui_div(res[i*dim+i],1,A[i*dim+i],MPFR_RNDN);
		for(int j=0;j<i;j++){
			mpfr_init2(res[i*dim+j],prec);
			//mpfr_set_ui(s,0,MPFR_RNDN);
			mpfr_set_zero(s,1);
			for(int k=j;k<i;k++){
				mpfr_mul(mul_temp,A[i*dim+k],res[k*dim+j],MPFR_RNDN);
				mpfr_add(s,s,mul_temp,MPFR_RNDN);
				}
			mpfr_neg(s,s,MPFR_RNDN);
			mpfr_div(res[i*dim+j],s,A[i*dim+i],MPFR_RNDN);
			} 
		} 
	mpfr_clear(mul_temp);
	mpfr_clear(s);
	return res;

}

mpfr_t* mpfr_cholesky(mpfr_t* A, int dim,mpfr_prec_t prec){

	mpfr_t *res=malloc(sizeof(mpfr_t)*dim*dim);

	mpfr_t s;
	mpfr_t mul_temp;

	mpfr_init2(s,prec);
	mpfr_init2(mul_temp,prec);

	for (int i = 0; i < dim; i++){
		for (int j=i+1; j<dim;j++){
			mpfr_init2(res[i*dim+j],prec);
			//mpfr_set_ui(res[i*dim+j],0,MPFR_RNDN); 
			mpfr_set_zero(res[i*dim+j],1);
		}

		for (int j = 0; j < (i+1); j++) {
			//mpfr_set_ui(s,0,MPFR_RNDN);
			mpfr_set_zero(s,1);
			mpfr_init2(res[i*dim+j],prec);
			for (int k = 0; k < j; k++){
				mpfr_mul(mul_temp,res[i*dim+k],res[j*dim+k],MPFR_RNDN);	
				mpfr_add(s,s,mul_temp,MPFR_RNDN);
			}
				mpfr_sub(res[i*dim+j],A[i*dim+j],s,MPFR_RNDN);
			if(i==j){
				mpfr_sqrt(res[i*dim+j],res[i*dim+j],MPFR_RNDN);
			}
			else{ 
				mpfr_div(res[i*dim+j],res[i*dim+j],res[j*dim+j],MPFR_RNDN);
			}

		} 
	}
	return res;
}

mpfr_t* form_anti_band(mpfr_t* ab_vector, int dim, mpfr_prec_t prec){
	int dim_res=dim*dim;
	mpfr_t *res=malloc(sizeof(mpfr_t)*dim_res);
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			mpfr_init2(res[i*dim+j],prec);
			mpfr_set(res[i*dim+j],ab_vector[i+j],MPFR_RNDN);
		}
	}
	return res;
}
