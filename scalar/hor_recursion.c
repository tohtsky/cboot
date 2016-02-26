#include "hor_recursion.h"
#include "stdio.h"
#define deb(x)\
	printf("debug at position %f\n",x)


mpfr_t* recursionNonZeroVector(unsigned long nMax, mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd){ 
	if(nMax<=7){
		printf("error: too small order of expansion in \"recursionSpinZeroVector\" function.");
		return NULL;
	}
	unsigned long order;
	order = nMax+1;
	mpfr_t *result=malloc(sizeof(mpfr_t)*order);
	mpfr_t *result_deriv=malloc(sizeof(mpfr_t)*order);

	mpfr_init2(result[0],prec);
	mpfr_set_ui(result[0],1,rnd);
	mpfr_t smallNumber;
	mpfr_init2(smallNumber,prec);
	mpfr_set_si_2exp(smallNumber,1,-(prec*3)/8+15,rnd);
	mpfr_t temp;
	mpfr_t temp2;
	
	mpfr_init2(temp,prec);
	mpfr_init2(temp2,prec);
	
	mpfr_t recCoeffs[8][5];
	mpfr_t recCoeffs_deriv[8][5];

	int has_a_zero_been_zero=0;

	initialize_spin_nonzero_coeffs_folder(recCoeffs, prec);
	set_nonzero_spin_rec_coeffs(recCoeffs, epsilon,ell,Delta,S,P,prec,rnd);

	mpfr_t as[8]; 

	for(int i=0;i<=7;i++){
		mpfr_init2(as[i],prec);
	}
	for(unsigned long i=1;i<=7;i++){
        	mpfr_set_si(temp,0,rnd);
		spin_nonzero_evaluate_at_n(as, recCoeffs, i, prec, rnd);
        	if(mpfr_cmp_ui(as[0],0)==0){
			for(unsigned long k=0;k<i;k++){
				mpfr_clear(result[k]);
			}
			free(result);
			mpfr_clear(temp);
			deallocate_spin_nonzero_coeffs_folder(recCoeffs); 

			mpfr_t shiftedDelta;
			mpfr_init2(shiftedDelta,prec);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_set(shiftedDelta,smallNumber,rnd);	
			}
			else{
				mpfr_add_ui(temp2,smallNumber,1,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd);
			}
			result=recursionNonZeroVector(nMax, epsilon, ell, shiftedDelta, S, P, prec, rnd);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_neg(smallNumber,smallNumber,rnd);
				mpfr_set(shiftedDelta,smallNumber,rnd);	
			}
			else{
				mpfr_ui_sub(temp2,1,smallNumber,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd);
			} 

			mpfr_t* result2=recursionNonZeroVector(nMax, epsilon, ell, shiftedDelta, S, P, prec, rnd);

			for(unsigned long k=0;k<=nMax;k++){
				mpfr_add(result[k],result[k],result2[k],rnd);
				mpfr_div_ui(result[k],result[k],2,rnd);
				mpfr_clear(result2[k]);
			}

			for(int i=0;i<=7;++i){
				mpfr_clear(as[i]);
			} 
			mpfr_clear(temp2);
			mpfr_clear(shiftedDelta);
			mpfr_clear(smallNumber);
			free(result2);
			return result;
		}
		for(unsigned long j=1;j<=i;j++){
		    mpfr_mul(temp2,as[j],result[i-j],rnd);
		    mpfr_add(temp,temp,temp2,rnd);
		} 
		mpfr_neg(temp,temp,rnd);
		mpfr_init2(result[i],prec);
		mpfr_div(result[i],temp,as[0],rnd);
	}

	for(unsigned long i=8;i<=nMax;i++){ 
		mpfr_set_si(temp,0,rnd);
		spin_nonzero_evaluate_at_n(as, recCoeffs, i, prec, rnd);
        	if(mpfr_cmp_ui(as[0],0)==0){
			for(unsigned long k=0;k<i;k++){
				mpfr_clear(result[k]);
			}
			free(result);
			mpfr_clear(temp);
			deallocate_spin_nonzero_coeffs_folder(recCoeffs); 

			mpfr_t shiftedDelta;
			mpfr_init2(shiftedDelta,prec);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_set(shiftedDelta,smallNumber,rnd);	
			}
			else{
				mpfr_add_ui(temp2,smallNumber,1,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd);
			}
			result=recursionNonZeroVector(nMax, epsilon, ell, shiftedDelta, S, P, prec, rnd);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_ui_sub(shiftedDelta,0,smallNumber,rnd);	
			}
			else{
				mpfr_ui_sub(temp2,1,smallNumber,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd);
			}

			mpfr_t* result2=recursionNonZeroVector(nMax, epsilon, ell, shiftedDelta, S, P, prec, rnd);

			for(unsigned long k=0;k<=nMax;k++){
				mpfr_add(result[k],result[k],result2[k],rnd);
				mpfr_div_ui(result[k],result[k],2,rnd);
				mpfr_clear(result2[k]);
			}

			for(int i=0;i<=7;++i){
				mpfr_clear(as[i]);
			} 
			mpfr_clear(temp2);
			mpfr_clear(shiftedDelta);
			mpfr_clear(smallNumber);
			free(result2);
			return result;
		}

		for(unsigned long j=1;j<=7;j++){
			mpfr_mul(temp2,as[j],result[i-j],rnd);
			mpfr_add(temp,temp,temp2,rnd);
		} 
		mpfr_neg(temp,temp,rnd);
		mpfr_init2(result[i],prec);
		mpfr_div(result[i],temp,as[0],rnd);
	}

	for(int i=0;i<=7;++i){
		mpfr_clear(as[i]);
		for(int j=0;j<=4;j++){
			mpfr_clear(recCoeffs[i][j]);
		}
	} 
	mpfr_clear(temp);
	mpfr_clear(temp2);
	mpfr_clear(smallNumber);
	return result;
}

mpfr_t* recursionSpinZeroVector(unsigned long nMax, mpfr_t epsilon, mpfr_t Delta, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd){ 
	if(nMax<=5){
		printf("error: too small order of expansion in \"recursionSpinZeroVector\" function.");
		return NULL;
	}
	unsigned long order;
	order = nMax+1;
	mpfr_t *result=malloc(sizeof(mpfr_t)*order);
	mpfr_init2(result[0],prec);
	mpfr_set_ui(result[0],1,rnd);
	mpfr_t smallNumber;
	mpfr_init2(smallNumber,prec);
	mpfr_set_si_2exp(smallNumber,1,-(3*prec)/8+15,rnd);
	
	mpfr_t temp;
	mpfr_t temp2;
	
	mpfr_init2(temp,prec);
	mpfr_init2(temp2,prec);
	
	mpfr_t recCoeffs[6][4];
	
	initialize_spin_zero_coeffs_folder(recCoeffs, prec);

	set_zero_spin_rec_coeffs(recCoeffs, epsilon,Delta,S,P,prec,rnd);
	mpfr_t as[6];
	for(int i=0;i<=5;i++){
		mpfr_init2(as[i],prec);
	}
	for(unsigned long i=1;i<=5;i++){
        	mpfr_set_si(temp,0,rnd);
		spin_zero_evaluate_at_n(as, recCoeffs, i, prec, rnd); 
        	if(mpfr_cmp_ui(as[0],0)==0){
			for(unsigned long k=0;k<i;k++){
				mpfr_clear(result[k]);
			}
			free(result);

			mpfr_clear(temp);
			deallocate_spin_zero_coeffs_folder(recCoeffs); 

			mpfr_t shiftedDelta;
			mpfr_init2(shiftedDelta,prec);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_set(Delta,smallNumber,rnd);	
			}
			else{
				mpfr_add_ui(temp2,smallNumber,1,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd);
			}
			result=recursionSpinZeroVector(nMax, epsilon, shiftedDelta, S, P, prec, rnd);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_ui_sub(Delta,0,smallNumber,rnd);	
			}
			else{
				mpfr_ui_sub(temp2,1,smallNumber,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd); 
			}

			mpfr_t* result2=recursionSpinZeroVector(nMax, epsilon,shiftedDelta, S, P, prec, rnd);

			for(unsigned long k=0;k<=nMax;k++){
				mpfr_add(result[k],result[k],result2[k],rnd);
				mpfr_div_ui(result[k],result[k],2,rnd);
				mpfr_clear(result2[k]);
			}

			for(int i=0;i<=5;++i){
				mpfr_clear(as[i]);
			} 
			mpfr_clear(temp2);
			mpfr_clear(shiftedDelta);
			mpfr_clear(smallNumber);
			free(result2);
			return result;
		}
        	
		for(unsigned long j=1;j<=i;j++){
		    mpfr_mul(temp2,as[j],result[i-j],rnd);
		    mpfr_add(temp,temp,temp2,rnd);
		} 
		mpfr_neg(temp,temp,rnd);
		mpfr_init2(result[i],prec);
		mpfr_div(result[i],temp,as[0],rnd);
	}
	for(unsigned long i=6;i<=nMax;i++){
		mpfr_set_si(temp,0,rnd);
		spin_zero_evaluate_at_n(as, recCoeffs, i, prec, rnd);
        	if(mpfr_cmp_ui(as[0],0)==0){
			for(unsigned long k=0;k<i;k++){
				mpfr_clear(result[k]);
			}
			free(result);

			mpfr_clear(temp);
			deallocate_spin_zero_coeffs_folder(recCoeffs); 

			mpfr_t shiftedDelta;
			mpfr_init2(shiftedDelta,prec);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_set(Delta,smallNumber,rnd);	
			}
			else{
				mpfr_add_ui(temp2,smallNumber,1,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd);
			}
			result=recursionSpinZeroVector(nMax, epsilon, shiftedDelta, S, P, prec, rnd);
			if(mpfr_cmp_ui(Delta,0)==0){
				mpfr_ui_sub(Delta,0,smallNumber,rnd);	
			}
			else{
				mpfr_ui_sub(temp2,1,smallNumber,rnd);
				mpfr_mul(shiftedDelta,Delta,temp2,rnd); 
			}

			mpfr_t* result2=recursionSpinZeroVector(nMax, epsilon,shiftedDelta, S, P, prec, rnd);

			for(unsigned long k=0;k<=nMax;k++){
				mpfr_add(result[k],result[k],result2[k],rnd);
				mpfr_div_ui(result[k],result[k],2,rnd);
				mpfr_clear(result2[k]);
			}

			for(int i=0;i<=5;++i){
				mpfr_clear(as[i]);
			} 
			mpfr_clear(temp2);
			mpfr_clear(shiftedDelta);
			mpfr_clear(smallNumber);
			free(result2);
			return result;
		}

		for(unsigned long j=1;j<=5;j++){
			mpfr_mul(temp2,as[j],result[i-j],rnd);
			mpfr_add(temp,temp,temp2,rnd);
		} 
		mpfr_neg(temp,temp,rnd);
		mpfr_init2(result[i],prec);
		mpfr_div(result[i],temp,as[0],rnd);
	}
	for(int i=0;i<=5;i++){
		mpfr_clear(as[i]);
		for(int j=0;j<=3;j++){
			mpfr_clear(recCoeffs[i][j]);
		}
	}
	mpfr_clear(temp);
	mpfr_clear(temp2);
	mpfr_clear(smallNumber);
	return result;
}

mpfr_t* real_axis_result(mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P,cb_context context){
	mpfr_t* hBlock;
	if(mpfr_cmp_ui(ell,0)==0){
		hBlock = recursionSpinZeroVector(context.n_Max,epsilon,Delta,S,P,context.prec, context.rnd);
	}
	else{
		hBlock = recursionNonZeroVector(context.n_Max,epsilon,ell,Delta,S,P,context.prec, context.rnd); 
	}
	mpfr_t* result_in_rho=malloc(sizeof(mpfr_t)*(context.lambda+1));
	mpfr_t temp1;
	mpfr_t temp2;
	mpfr_t temp3;
	mpfr_init2(temp1,context.prec);
	mpfr_init2(temp2,context.prec);
	mpfr_init2(temp3,context.prec); 

	mpfr_init2(result_in_rho[0],context.prec);
	mpfr_set_si(result_in_rho[0],0,context.rnd); 
	mpfr_mul_ui(temp1,context.rho,4,context.rnd);
	mpfr_pow(temp1,temp1,Delta, context.rnd);

	for(unsigned long j=0;j<=context.n_Max;j++){ 
		mpfr_mul(hBlock[j],hBlock[j],temp1,context.rnd);
		mpfr_add(result_in_rho[0],result_in_rho[0],hBlock[j],context.rnd);
		if(j<context.n_Max){
			mpfr_mul(temp1,temp1,context.rho,context.rnd);
		}
	}
	for(int i=1; i<=context.lambda;i++){
		mpfr_init2(result_in_rho[i],context.prec);
		mpfr_set_si(result_in_rho[i],0,context.rnd); 
		for(unsigned long j=0;j<=context.n_Max;j++){
			mpfr_add_si(temp1,Delta,(j-i+1),context.rnd);
			mpfr_mul(hBlock[j],hBlock[j],temp1,context.rnd);
			mpfr_add(result_in_rho[i],result_in_rho[i],hBlock[j],context.rnd);
		}
		mpfr_pow_ui(temp1,context.rho,i,MPFR_RNDN); 
		mpfr_fac_ui(temp2,i,MPFR_RNDN);
		mpfr_mul(temp2,temp2,temp1,MPFR_RNDN);
		mpfr_div(result_in_rho[i],result_in_rho[i],temp2,MPFR_RNDN);
	}

	for(long i=0;i<=context.n_Max;i++){
		mpfr_clear(hBlock[i]);
	} 

	free(hBlock);
	mpfr_t* result=malloc(sizeof(mpfr_t)*(context.lambda+1));
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_init2(result[j],context.prec);
		mpfr_set_zero(result[j],1); 
		for(unsigned long k=0;k<=context.lambda;k++){
			mpfr_mul(temp1,result_in_rho[k],context.rho_to_z_matrix[k+(context.lambda+1)*j],MPFR_RNDN);
			mpfr_add(result[j],result[j],temp1,MPFR_RNDN); 
		}
	} 
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_clear(result_in_rho[j]);
	}

	free(result_in_rho); 
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	mpfr_clear(temp3);

	return result; 
}
long indexOfConformalBlock(cb_context context, int n, int m){
	if(m>=0&&n>=0){
		return (context.lambda+2-n)*n+m;
	}
	else{
		return -1;
	}
}


void element_helper(cb_context context, mpfr_t* array,mpfr_t r,int m, int n){
	long j=indexOfConformalBlock(context,m,n);
	if(j>=0){
		mpfr_set(r,*(array+j),MPFR_RNDN);
	}
	else{
		mpfr_set_zero(r,1);
	}
} 
 

mpfr_t* casimirExpander(mpfr_t* realAxisResult, mpfr_t epsilon, mpfr_t ell,  mpfr_t Delta, mpfr_t S, mpfr_t P, cb_context context){
	mpfr_t* result;
	if(context.lambda%2){
		result=malloc(sizeof(mpfr_t)*((context.lambda+5)*(context.lambda-1)/4 +2));
	}
	else{
		result=malloc(sizeof(mpfr_t)*((context.lambda+4)*(context.lambda)/4 +1)); 
	}
	for(int i=0;i<=context.lambda;++i){
		mpfr_init2(result[i],context.prec);
		mpfr_set(result[i],realAxisResult[i],MPFR_RNDN);
	}

	mpfr_t Casimir;
	mpfr_t temp1, temp2, temp3, r;
	mpfr_init2(temp1,context.prec);
	mpfr_init2(temp2,context.prec);
	mpfr_init2(temp3,context.prec);
	mpfr_init2(temp3,context.prec);
	mpfr_init2(r,context.prec);
	mpfr_init2(Casimir,context.prec);

	/*computing quadratic Casimir times 2 */
	mpfr_mul_ui(temp3,epsilon,2,MPFR_RNDN);
	mpfr_sub(temp2,Delta,temp3,MPFR_RNDN);
	mpfr_sub_ui(temp2,temp2,2,MPFR_RNDN);
	mpfr_mul(temp2,temp2,Delta,MPFR_RNDN); 

	mpfr_add(temp3,temp3,ell,MPFR_RNDN);
	mpfr_mul(temp3,temp3,ell,MPFR_RNDN); 
	mpfr_add(Casimir,temp3,temp2,MPFR_RNDN); 

	for(int i=1;i<=(context.lambda/2);i++){
		for(int j=0;j<=context.lambda-2*i;++j){ 
			mpfr_init2(result[indexOfConformalBlock(context,i,j)],context.prec);
			mpfr_set_zero(temp1,1);
			/* The first line in arXiv 1203.6064 (C.1)
			 * note: n(there) = i (here), and m(there)=j (here)
			 * (a/2)=x
			 * (b/2)=sqrt(y)
			 * d/da ^m d/db ^n (g) = h_{n,m} = (1/2)^(m+n)d/dx ^m d/ dy ^n g
			 * and h_{m,n}=(1/2)^(i+j) * i!j!h_{i,j}
			 * */
			element_helper(context,result,r,i,j-2);
			mpfr_mul_ui(temp2,r,4,MPFR_RNDN);
			element_helper(context,result,r,i,j-3);
			mpfr_mul_ui(temp3,r,8,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN);

			element_helper(context,result,r,i,j-1);
			mpfr_mul_ui(temp3,r,2,MPFR_RNDN);
			mpfr_sub(temp2,temp2,temp3,MPFR_RNDN);
			
			mpfr_mul_ui(temp3,epsilon,2,MPFR_RNDN);
			mpfr_add_si(temp3,temp3,2*i-1,MPFR_RNDN);
			mpfr_mul_ui(temp3,temp3,2,MPFR_RNDN);
			mpfr_mul(temp2,temp2,temp3,MPFR_RNDN);
			mpfr_add(temp1,temp1,temp2,MPFR_RNDN);

			/* The second line */
			element_helper(context,result,r,i-1,j+2);
			mpfr_mul_si(temp2,r,-(j+2)*(j+1),MPFR_RNDN);
			mpfr_div_ui(temp2,temp2,i,MPFR_RNDN);
			mpfr_add(temp1,temp1,temp2,MPFR_RNDN);

			mpfr_mul_ui(temp2,epsilon,2,MPFR_RNDN);
			mpfr_add_si(temp2,temp2, -j-4*i +6,MPFR_RNDN);


			mpfr_mul_ui(temp3,S,2,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN);
			/*
			 *
			 *add 2*alpha + 2*beta = 2*S to the above!
			 *  -done
			 * */
			element_helper(context,result,r,i-1,j+1);
			mpfr_mul(temp2,temp2,r,MPFR_RNDN);
			mpfr_mul_si(temp2,temp2,2*(j+1),MPFR_RNDN);
			mpfr_div_ui(temp2,temp2,i,MPFR_RNDN);
			mpfr_add(temp1,temp1,temp2,MPFR_RNDN);

			/* The third line */
			mpfr_mul_si(temp2,epsilon,4*(j+i-1),MPFR_RNDN);
			mpfr_add_si(temp2,temp2,j*j+8*j*i-5*j+4*i*i-2*i-2,MPFR_RNDN); 

			mpfr_mul_ui(temp3,P,2,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN);
			mpfr_mul_si(temp3,S,8*i+4*j-8,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN);

			/*
			 *
			 *add 4*alpha*beta + alpha*(8*i +4* j-8) + beta*(8*i+4*j - 8) to the above!
			 * 2*P + S*(8*i+4*j - 8)
			 * -done
			 *
			 * */

			mpfr_mul_ui(temp3,Casimir,2,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN);
			element_helper(context,result,r,i-1,j);
			mpfr_mul(temp2,temp2,r,MPFR_RNDN);
			mpfr_mul_ui(temp2,temp2,4,MPFR_RNDN);
			mpfr_div_ui(temp2,temp2,i,MPFR_RNDN);
			mpfr_add(temp1,temp1,temp2,MPFR_RNDN);

			/* The fourth line */
			mpfr_mul_si(temp2,epsilon,2*(j-2*i+1),MPFR_RNDN);
			mpfr_add_si(temp2,temp2,j*j+12*j*i-13*j+12*i*i-34*i+22,MPFR_RNDN);

			mpfr_mul_ui(temp3,P,2,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN);
			mpfr_mul_si(temp3,S,8*i+2*j-10,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN);

			/*
			 *
			 *add 4*alpha*beta + (alpha+beta)*(8*i+2*j-10) to the above!!
			 * =2*P+S*(8*i+2*j -10) 
			 *
			 * */

			element_helper(context,result,r,i-1,j-1);
			mpfr_mul(temp2,temp2,r,MPFR_RNDN);
			mpfr_mul_ui(temp2,temp2,8,MPFR_RNDN);
			mpfr_div_ui(temp2,temp2,i,MPFR_RNDN);
			mpfr_add(temp1,temp1,temp2,MPFR_RNDN);

			/* The last line */
			mpfr_mul_ui(temp2,epsilon,2,MPFR_RNDN);
			mpfr_add_si(temp2,temp2,-3*j-4*i+6,MPFR_RNDN);
			
			mpfr_mul_si(temp3, S, -2,MPFR_RNDN);
			mpfr_add(temp2,temp2,temp3,MPFR_RNDN); 
			/*
			 *
			 * add -2*(alpha + beta ) to this..
			 * = -2*S
			 * -done
			 * */
			element_helper(context,result,r,i-2,j+1);
			mpfr_mul(temp2,temp2,r,MPFR_RNDN);
			mpfr_mul_ui(temp2,temp2,8*(j+1),MPFR_RNDN);

			element_helper(context,result,r,i-2,j+2);
			mpfr_mul_ui(temp3,r,4*(j+1)*(j+2),MPFR_RNDN);
			mpfr_sub(temp2,temp3,temp2,MPFR_RNDN); 
			mpfr_div_ui(temp2,temp2,i,MPFR_RNDN);
			mpfr_add(temp1,temp1,temp2,MPFR_RNDN);

			/* finally division by 2(D+2n-3) */
			mpfr_mul_ui(temp2,epsilon,4,MPFR_RNDN);
			mpfr_add_si(temp2,temp2,4*i-2,MPFR_RNDN); 
			mpfr_div(result[indexOfConformalBlock(context,i,j)],temp1,temp2,MPFR_RNDN);
		}
	}
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	mpfr_clear(temp3);
	mpfr_clear(Casimir);
	return result;

}


mpfr_t* gBlock_full(mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P,cb_context context){

	mpfr_t* realAxisResult=real_axis_result(epsilon,ell,Delta,S,P,context); 
	mpfr_t* result=casimirExpander(realAxisResult,epsilon,ell,Delta,S,P,context);
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_clear(realAxisResult[j]);
	}
	free(realAxisResult);
	return result;
}

mpfr_t* hBlock_times_rho_n(unsigned long n, mpfr_t epsilon, mpfr_t ell, mpfr_t Delta, mpfr_t S, mpfr_t P, cb_context context){
	/* *
	 * gives (4*rho)^{n} * h(\Delta, l,...) = (4*rho)^{n-\Delta}g*(\Delta, l, ...)
	 * , evaluated in x-y coordinate
	 * */
	mpfr_t* hBlock;
	if(mpfr_cmp_ui(ell,0)==0){
		hBlock = recursionSpinZeroVector(context.n_Max,epsilon,Delta,S,P,context.prec, context.rnd);
	}
	else{
		hBlock = recursionNonZeroVector(context.n_Max,epsilon,ell,Delta,S,P,context.prec, context.rnd); 
	}
	mpfr_t* result_in_rho=malloc(sizeof(mpfr_t)*(context.lambda+1));
	mpfr_t temp1;
	mpfr_t temp2;
	mpfr_t temp3;
	mpfr_init2(temp1,context.prec);
	mpfr_init2(temp2,context.prec);
	mpfr_init2(temp3,context.prec); 

	mpfr_init2(result_in_rho[0],context.prec);
	mpfr_set_si(result_in_rho[0],0,context.rnd); 
	mpfr_mul_ui(temp1,context.rho,4,context.rnd);
	mpfr_pow_ui(temp1,temp1,n, context.rnd);
	for(unsigned long j=0;j<=context.n_Max;j++){ 
		mpfr_mul(hBlock[j],hBlock[j],temp1,context.rnd);
		mpfr_add(result_in_rho[0],result_in_rho[0],hBlock[j],context.rnd);
		if(j<context.n_Max){
			mpfr_mul(temp1,temp1,context.rho,context.rnd);
		}
	}
	for(int i=1; i<=context.lambda;i++){
		mpfr_init2(result_in_rho[i],context.prec);
		mpfr_set_si(result_in_rho[i],0,context.rnd); 
		for(unsigned long j=0;j<=context.n_Max;j++){
			mpfr_mul_si(hBlock[j],hBlock[j],n+j-i+1,context.rnd);
			mpfr_add(result_in_rho[i],result_in_rho[i],hBlock[j],context.rnd);
		}
		mpfr_pow_ui(temp1,context.rho,i,MPFR_RNDN); 
		mpfr_fac_ui(temp2,i,MPFR_RNDN);
		mpfr_mul(temp2,temp2,temp1,MPFR_RNDN);
		mpfr_div(result_in_rho[i],result_in_rho[i],temp2,MPFR_RNDN);
	}

	for(long i=0;i<=context.n_Max;i++){
		mpfr_clear(hBlock[i]);
	} 

	free(hBlock);
	mpfr_t* result=malloc(sizeof(mpfr_t)*(context.lambda+1));
	for(unsigned long i=0;i<=context.lambda;i++){
		mpfr_init2(result[i],context.prec);
		mpfr_set_zero(result[i],1); 
		for(unsigned long j=0;j<=context.lambda;j++){
			mpfr_mul(temp1,result_in_rho[j],context.rho_to_z_matrix[j+(context.lambda+1)*i],MPFR_RNDN);
			mpfr_add(result[i],result[i],temp1,MPFR_RNDN); 
		}
	} 
	for(unsigned long i=0;i<=context.lambda;i++){
		mpfr_clear(result_in_rho[i]);
	}

	free(result_in_rho); 
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	mpfr_clear(temp3);

	return result; 
}


mpfr_t* h_asymptotic(mpfr_t epsilon, mpfr_t S, cb_context context){
	mpfr_t temp1, temp2, temp3;
	mpfr_init2(temp1,context.prec);
	mpfr_init2(temp2,context.prec);
	mpfr_init2(temp3,context.prec);
	/* first factor */
	mpfr_mul_ui(temp1,S,2,context.rnd);
	mpfr_add_si(temp1,temp1,-1,context.rnd);
	mpfr_sub(temp1,temp1,epsilon,context.rnd);
    /* temp1 = 2S -1 - epsilon*/
	mpfr_t* firstFactor=malloc(sizeof(mpfr_t)*(context.lambda+1)); 
	mpfr_init2(firstFactor[0],context.prec); 
	mpfr_add_ui(temp2,context.rho,1,context.rnd);
    /* temp2 = 1 + rho */
	mpfr_pow(firstFactor[0],temp2,temp1,context.rnd);
    /* firstFactor[0] = (1+rho)**(-1-epsilon + 2S)*/
	mpfr_ui_div(temp2,1,temp2,context.rnd);
	for(unsigned long j=1;j<=context.lambda;j++){
		mpfr_init2(firstFactor[j],context.prec); 
		mpfr_add_si(temp3,temp1,-j+1,context.rnd);
		mpfr_mul(firstFactor[j],firstFactor[j-1],temp3,context.rnd);
		mpfr_mul(firstFactor[j],firstFactor[j],temp2,context.rnd);
		mpfr_div_ui(firstFactor[j],firstFactor[j],j,context.rnd); 
	}

	/*second factor*/
	mpfr_mul_si(temp1,S,-2,context.rnd);
	mpfr_add_si(temp1,temp1,-1,context.rnd);
	mpfr_sub(temp1,temp1,epsilon,context.rnd);
	mpfr_t* secondFactor=malloc(sizeof(mpfr_t)*(context.lambda+1)); 
	mpfr_init2(secondFactor[0],context.prec); 
	mpfr_ui_sub(temp2,1,context.rho,context.rnd);
	mpfr_pow(secondFactor[0],temp2,temp1,context.rnd);
	mpfr_ui_div(temp2,1,temp2,context.rnd);
	mpfr_neg(temp2,temp2,context.rnd);
	for(unsigned long j=1;j<=context.lambda;j++){
		mpfr_init2(secondFactor[j],context.prec); 
		mpfr_add_si(temp3,temp1,-j+1,context.rnd);
		mpfr_mul(secondFactor[j],secondFactor[j-1],temp3,context.rnd);
		mpfr_mul(secondFactor[j],secondFactor[j],temp2,context.rnd);
		mpfr_div_ui(secondFactor[j],secondFactor[j],j,context.rnd); 
	}

    /* Convolve firstFactor and secondFactor */
	mpfr_t* result_in_rho=malloc(sizeof(mpfr_t)*(context.lambda+1)); 
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_init2(result_in_rho[j],context.prec);
		mpfr_set_zero(result_in_rho[j],1);
		for(unsigned long k=0;k<=j;k++){
			mpfr_mul(temp1,firstFactor[k],secondFactor[j-k],context.rnd);
			mpfr_add(result_in_rho[j],result_in_rho[j],temp1,context.rnd);	
		}
	}
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_clear(firstFactor[j]);	
	};
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_clear(secondFactor[j]);	
	};
	free(firstFactor);
	free(secondFactor);

	mpfr_t* result=(mpfr_t *) malloc(sizeof(mpfr_t)*(context.lambda+1));
	for(int i=0;i<=context.lambda;i++){
		mpfr_init2(result[i],context.prec);
		mpfr_set_zero(result[i],1); 
		for(unsigned long j=0;j<=context.lambda;j++){
			mpfr_mul(temp1,result_in_rho[j],context.rho_to_z_matrix[j+(context.lambda+1)*i],context.rnd);
			mpfr_add(result[i],result[i],temp1,context.rnd); 
		}
	} 

	for(unsigned long i=0;i<=context.lambda;i++){
		mpfr_clear(result_in_rho[i]);
	}

	free(result_in_rho); 
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	mpfr_clear(temp3); 
	return result; 
}
