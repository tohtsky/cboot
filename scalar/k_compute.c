#include "k_compute.h"

#define deb(x)\
	printf("deb %f\n",x)

void set_k_coeffs(mpfr_t a[4][3], mpfr_t h, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd){
    mpfr_t _t_0;
    mpfr_init2(_t_0, prec);
    mpfr_t _t_1;
    mpfr_init2(_t_1, prec);
    mpfr_set_si(a[0][0],0,rnd);
    mpfr_mul_si(_t_0,h,-8,rnd);
    mpfr_add_si(a[0][1],_t_0,4,rnd);
    mpfr_set_si(a[0][2],-4,rnd);
    mpfr_mul_si(_t_0,S,16,rnd);
    mpfr_mul(_t_0,_t_0,h,rnd);
    mpfr_mul_si(_t_1,S,16,rnd);
    mpfr_sub(_t_0,_t_0,_t_1,rnd);
    mpfr_mul_si(_t_1,P,8,rnd);
    mpfr_add(_t_0,_t_0,_t_1,rnd);
    mpfr_mul_si(_t_1,h,8,rnd);
    mpfr_add(_t_0,_t_0,_t_1,rnd);
    mpfr_sub_si(a[1][0],_t_0,8,rnd);
    mpfr_mul_si(_t_0,S,16,rnd);
    mpfr_mul_si(_t_1,h,8,rnd);
    mpfr_sub(_t_0,_t_0,_t_1,rnd);
    mpfr_add_si(a[1][1],_t_0,12,rnd);
    mpfr_set_si(a[1][2],-4,rnd);
    mpfr_mul_si(_t_0,S,16,rnd);
    mpfr_mul(_t_0,_t_0,h,rnd);
    mpfr_mul_si(_t_1,S,32,rnd);
    mpfr_sub(_t_0,_t_0,_t_1,rnd);
    mpfr_mul_si(_t_1,P,8,rnd);
    mpfr_sub(_t_0,_t_0,_t_1,rnd);
    mpfr_mul_si(_t_1,h,8,rnd);
    mpfr_sub(_t_0,_t_0,_t_1,rnd);
    mpfr_add_si(a[2][0],_t_0,8,rnd);
    mpfr_mul_si(_t_0,S,16,rnd);
    mpfr_mul_si(_t_1,h,8,rnd);
    mpfr_add(_t_0,_t_0,_t_1,rnd);
    mpfr_sub_si(a[2][1],_t_0,12,rnd);
    mpfr_set_si(a[2][2],4,rnd);
    mpfr_mul_si(_t_0,h,-16,rnd);
    mpfr_add_si(a[3][0],_t_0,24,rnd);
    mpfr_mul_si(_t_0,h,8,rnd);
    mpfr_sub_si(a[3][1],_t_0,20,rnd);
    mpfr_set_si(a[3][2],4,rnd);
    mpfr_clear(_t_0);
    mpfr_clear(_t_1); 
}
void initialize_k_coeffs(mpfr_t a[4][3], mpfr_prec_t prec){
	for(int i=0;i<=3;i++){
		for(int j=0;j<=2;j++){
			mpfr_init2(a[i][j],prec);
		}
	}
}
void deallocate_k_coeffs(mpfr_t a[4][3]){
	for(int i=0;i<=3;i++){
		for(int j=0;j<=2;j++){
			mpfr_clear(a[i][j]);
		}
	}
}

void k_evaluate_at_n(mpfr_t a[4], mpfr_t rho[4][3],long n, mpfr_prec_t prec, mp_rnd_t rnd){
	for(int i=0;i<=3;i++){
		mpfr_set(a[i],rho[i][2],rnd);
		for(int j=0;j<=1;j++){
			mpfr_mul_si(a[i],a[i],n,rnd);
			mpfr_add(a[i],a[i],rho[i][1-j],rnd);
		}
	}
}
//
//mpfr_t* rhoInZ(unsigned long Lambda, long prec);
//
//static int lambda_at_Module = -1;
//static mp_rnd_t rnd = MPFR_RNDN;
//static mpfr_prec_t precision_at_Module = -1;
//static long nMax_at_Module=250; 
//static mpfr_t zeroAtPrecision;
//static mpfr_t rho;
//static mpfr_t rhoPrime;
//
mpfr_t* rhoInZ(unsigned long Lambda_arg, long prec){
	/* This time I was bothered to write lambda + 1 so many times...*/
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
		//mpfr_set_ui(result[j],0,MPFR_RNDN); 
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

mpfr_t* recursion_k(unsigned long nMax, mpfr_t h, mpfr_t S, mpfr_t P, mpfr_prec_t prec, mp_rnd_t rnd){
	if(nMax<=3){
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
	mpfr_set_si_2exp(smallNumber,1,-(prec*3)/8+15,rnd);
	mpfr_t temp;
	mpfr_t temp2;
	mpfr_init2(temp,prec);
	mpfr_init2(temp2,prec);
	
	mpfr_t recCoeffs[4][3];
	initialize_k_coeffs(recCoeffs,prec);

	set_k_coeffs(recCoeffs,h,S,P,prec,rnd);
	mpfr_t as[4];
	for(int i=0;i<=3;i++){
		mpfr_init2(as[i],prec);
	}
	for(unsigned long i=1;i<=2;i++){
        	mpfr_set_si(temp,0,rnd);
		k_evaluate_at_n(as, recCoeffs, i, prec, rnd);
        	if(mpfr_cmp_ui(as[0],0)==0){
			for(unsigned long k=0;k<i;k++){
				mpfr_clear(result[k]);
			}
			free(result);
			mpfr_clear(temp);
			deallocate_k_coeffs(recCoeffs); 
			
			mpfr_t shiftedh;
			mpfr_init2(shiftedh,prec);
			if(mpfr_cmp_ui(h,0)==0){
				mpfr_set(shiftedh,smallNumber,rnd);	
			}
			else{
				mpfr_add_ui(temp2,smallNumber,1,rnd);
				mpfr_mul(shiftedh,h,temp2,rnd);
			}
			result=recursion_k(nMax, shiftedh, S, P, prec, rnd);
			if(mpfr_cmp_ui(h,0)==0){
				mpfr_neg(smallNumber,smallNumber,rnd);
				mpfr_set(shiftedh,smallNumber,rnd);	
				//mpfr_ui_sub(shiftedh,0,smallNumber,rnd);	
			}
			else{
				mpfr_ui_sub(temp2,1,smallNumber,rnd);
				mpfr_mul(shiftedh,h,temp2,rnd);
			} 

			mpfr_t* result2= recursion_k(nMax,shiftedh, S, P, prec, rnd);

			for(unsigned long k=0;k<=nMax;k++){
				mpfr_add(result[k],result[k],result2[k],rnd);
				mpfr_div_ui(result[k],result[k],2,rnd);
				mpfr_clear(result2[k]);
			}

			for(int i=0;i<=3;++i){
				mpfr_clear(as[i]);
			} 
			mpfr_clear(temp2);
			mpfr_clear(shiftedh);
			mpfr_clear(smallNumber);
			free(result2);
			return result;
		}
		//mpfr_printf("result[%d] = %.32RNf\n",i-1,result[i-1]);
		for(unsigned long j=1;j<=i;j++){
		    mpfr_mul(temp2,as[j],result[i-j],rnd);
		    mpfr_add(temp,temp,temp2,rnd);
		} 
		mpfr_neg(temp,temp,rnd);
		mpfr_init2(result[i],prec);
		mpfr_div(result[i],temp,as[0],rnd);
	}
	for(unsigned long i=3;i<=nMax;i++){
        	mpfr_set_si(temp,0,rnd);
		k_evaluate_at_n(as, recCoeffs, i, prec, rnd);
        	if(mpfr_cmp_ui(as[0],0)==0){
			printf("possible 0-division at loop %ld \n",i);
			for(unsigned long k=0;k<i;k++){
				mpfr_clear(result[k]);
			}
			free(result);
			mpfr_clear(temp);
			deallocate_k_coeffs(recCoeffs); 

			mpfr_t shiftedh;
			mpfr_init2(shiftedh,prec);
			if(mpfr_cmp_ui(h,0)==0){
				mpfr_set(shiftedh,smallNumber,rnd);	
			}
			else{
				mpfr_add_ui(temp2,smallNumber,1,rnd);
				mpfr_mul(shiftedh,h,temp2,rnd);
			}
			result=recursion_k(nMax, shiftedh, S, P, prec, rnd);
			if(mpfr_cmp_ui(h,0)==0){
				mpfr_neg(smallNumber,smallNumber,rnd);
				mpfr_set(shiftedh,smallNumber,rnd);	
				//mpfr_ui_sub(shiftedh,0,smallNumber,rnd);	
			}
			else{
				mpfr_ui_sub(temp2,1,smallNumber,rnd);
				mpfr_mul(shiftedh,h,temp2,rnd);
			} 

			mpfr_t* result2= recursion_k(nMax,shiftedh, S, P, prec, rnd);

			for(unsigned long k=0;k<=nMax;k++){
				mpfr_add(result[k],result[k],result2[k],rnd);
				mpfr_div_ui(result[k],result[k],2,rnd);
				mpfr_clear(result2[k]);
			}

			for(int i=0;i<=3;++i){
				mpfr_clear(as[i]);
			} 
			mpfr_clear(temp2);
			mpfr_clear(shiftedh);
			mpfr_clear(smallNumber);
			free(result2);
			return result;
		}
		for(unsigned long j=1;j<=3;j++){
		    mpfr_mul(temp2,as[j],result[i-j],rnd);
		    mpfr_add(temp,temp,temp2,rnd);
		}

		mpfr_neg(temp,temp,rnd);
		mpfr_init2(result[i],prec);
		mpfr_div(result[i],temp,as[0],rnd);
	}

	for(int i=0;i<=3;++i){
		mpfr_clear(as[i]);
		for(int j=0;j<=2;j++){
			mpfr_clear(recCoeffs[i][j]);
		}
	} 
	mpfr_clear(temp);
	mpfr_clear(temp2);
	mpfr_clear(smallNumber);
	return result;
}

mpfr_t* k_table_c(mpfr_t h,mpfr_t S, mpfr_t P, cb_context context){
	mpfr_t* hBlock = recursion_k(context.n_Max, h, S, P, context.prec, context.rnd);
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
	mpfr_pow(temp1,temp1,h,context.rnd);

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
			mpfr_add_si(temp1,h,(j-i+1),context.rnd);
			mpfr_mul(hBlock[j],hBlock[j],temp1,context.rnd);
			mpfr_add(result_in_rho[i],result_in_rho[i],hBlock[j],context.rnd);
		}
		mpfr_pow_ui(temp1,context.rho,i,context.rnd); 
		mpfr_fac_ui(temp2,i,context.rnd);
		mpfr_mul(temp2,temp2,temp1,context.rnd);
		mpfr_div(result_in_rho[i],result_in_rho[i],temp2,context.rnd);
	}

	for(long i=0;i<=context.n_Max;i++){
		mpfr_clear(hBlock[i]);
	} 
	free(hBlock);

	mpfr_t* result=malloc(sizeof(mpfr_t)*(context.lambda+1));
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_init2(result[j],context.prec);
		//mpfr_set_ui(result[j],0,context.rnd);
		mpfr_set_zero(result[j],1);
		for(unsigned long k=0;k<=context.lambda;k++){
			mpfr_mul(temp1,result_in_rho[k],context.rho_to_z_matrix[k+(context.lambda+1)*j],context.rnd);
			mpfr_add(result[j],result[j],temp1,context.rnd); 
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

mpfr_t* chiral_h_times_rho_to_n_c(unsigned long n, mpfr_t h,mpfr_t S, mpfr_t P, cb_context context){
	mpfr_t* hBlock = recursion_k(context.n_Max, h, S, P, context.prec, context.rnd);
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
	mpfr_pow_ui(temp1,temp1,n,context.rnd);

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
		mpfr_pow_ui(temp1,context.rho,i,context.rnd); 
		mpfr_fac_ui(temp2,i,context.rnd);
		mpfr_mul(temp2,temp2,temp1,context.rnd);
		mpfr_div(result_in_rho[i],result_in_rho[i],temp2,context.rnd);
	}

	for(long i=0;i<=context.n_Max;i++){
		mpfr_clear(hBlock[i]);
	}
	free(hBlock);

	mpfr_t* result=malloc(sizeof(mpfr_t)*(context.lambda+1));
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_init2(result[j],context.prec);
		//mpfr_set_ui(result[j],0,context.rnd);
		mpfr_set_zero(result[j],1);
		for(unsigned long k=0;k<=context.lambda;k++){
			mpfr_mul(temp1,result_in_rho[k],context.rho_to_z_matrix[k+(context.lambda+1)*j],context.rnd);
			mpfr_add(result[j],result[j],temp1,context.rnd); 
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

mpfr_t* chiral_h_asymptotic_c(mpfr_t S, cb_context context){
	mpfr_t temp1, temp2, temp3;
	mpfr_init2(temp1,context.prec);
	mpfr_init2(temp2,context.prec);
	mpfr_init2(temp3,context.prec);
	/* first factor */
	mpfr_mul_ui(temp1,S,2,context.rnd);
	mpfr_add_si(temp1,temp1,-1,context.rnd);
	mpfr_div_ui(temp1,temp1,2,context.rnd);
	mpfr_t* firstFactor=malloc(sizeof(mpfr_t)*(context.lambda+1)); 
	mpfr_init2(firstFactor[0],context.prec); 
	mpfr_add_ui(temp2,context.rho,1,context.rnd);
	mpfr_pow(firstFactor[0],temp2,temp1,context.rnd);
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
	mpfr_div_ui(temp1,temp1,2,context.rnd);
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
	mpfr_t* result_in_rho=malloc(sizeof(mpfr_t)*(context.lambda+1)); 
	for(unsigned long j=0;j<=context.lambda;j++){
		mpfr_init2(result_in_rho[j],context.prec);
		//mpfr_set_ui(result_in_rho[j],0,context.rnd);
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

	mpfr_t* result=malloc(sizeof(mpfr_t)*(context.lambda+1));
	for(unsigned long i=0;i<=context.lambda;i++){
		mpfr_init2(result[i],context.prec);
		//mpfr_set_ui(result[i],0,context.rnd);
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
