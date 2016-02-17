#include "integral_decomp.h"
#define deb(d)\
	printf("deb %f\n",d)
#undef debug_mode
//#define debug_mode 1

mpfr_t* simple_pole_case_c(long pole_order_max, mpfr_t base, mpfr_t pole_position, mpfr_t incomplete_gamma_factor, mpfr_prec_t prec){ 

	mpfr_t* result=malloc(sizeof(mpfr_t)*(pole_order_max+1)); 

	mpfr_t temp1;
	mpfr_init2(temp1,prec);
	mpfr_t temp2;
	mpfr_init2(temp2,prec); 
	
	mpfr_t temp3;
	mpfr_init2(temp3,prec); 

	mpfr_t minus_pole_position;
	mpfr_init2(minus_pole_position,prec); 
	mpfr_neg(minus_pole_position,pole_position,MPFR_RNDN);

	mpfr_t factorial;
	mpfr_init2(factorial,prec); 
	mpfr_set_ui(factorial,1,MPFR_RNDN);

	mpfr_t minus_log_base;
	mpfr_init2(minus_log_base,prec); 
	mpfr_log(minus_log_base,base,prec); 
	mpfr_neg(minus_log_base,minus_log_base,MPFR_RNDN);
	mpfr_ui_div(minus_log_base,1,minus_log_base,MPFR_RNDN);

	mpfr_t log_base_power;
	mpfr_init2(log_base_power,prec); 
	mpfr_set(log_base_power,minus_log_base,MPFR_RNDN);

	mpfr_set_ui(temp1,0,MPFR_RNDN);
	mpfr_mul(temp2,pole_position,incomplete_gamma_factor,MPFR_RNDN); 

	mpfr_init2(result[0],prec);
	mpfr_set(result[0],incomplete_gamma_factor,MPFR_RNDN);
	for(long j=1;j<=pole_order_max;j++){
		mpfr_init2(result[j],prec);
		mpfr_mul(temp1,temp1,pole_position,MPFR_RNDN);	
		mpfr_mul(temp3,factorial,log_base_power,MPFR_RNDN);
		mpfr_add(temp1,temp3,temp1,MPFR_RNDN);
		mpfr_add(result[j],temp1,temp2,MPFR_RNDN);

		if(j<pole_order_max){ 
			mpfr_mul(temp2,temp2,pole_position,MPFR_RNDN);
			mpfr_mul(log_base_power,log_base_power,minus_log_base,MPFR_RNDN);
			mpfr_mul_si(factorial,factorial,j,MPFR_RNDN);
		}
	}

	mpfr_clear(temp1);
	mpfr_clear(temp2);
	mpfr_clear(temp3);
	mpfr_clear(minus_pole_position);
	mpfr_clear(factorial);
	mpfr_clear(minus_log_base);
	mpfr_clear(log_base_power);
	return result;
}

mpfr_t* double_pole_case_c(long pole_order_max, mpfr_t base, mpfr_t pole_position, mpfr_t incomplete_gamma_factor, mpfr_prec_t prec){ 

	mpfr_t* result=malloc(sizeof(mpfr_t)*(pole_order_max+1)); 

	mpfr_t temp1;
	mpfr_init2(temp1,prec);

	mpfr_t double_pole_integral;
	mpfr_init2(double_pole_integral,prec); 

	mpfr_t temp2;
	mpfr_init2(temp2,prec); 

	mpfr_t minus_pole_position;
	mpfr_init2(minus_pole_position,prec); 
	mpfr_neg(minus_pole_position,pole_position,MPFR_RNDN);

//	mpfr_t factorial;
//	mpfr_init2(factorial,prec); 
//	mpfr_set_ui(factorial,1,MPFR_RNDN);

	mpfr_t minus_log_base;
	mpfr_init2(minus_log_base,prec); 
	mpfr_log(minus_log_base,base,prec); 
	mpfr_neg(minus_log_base,minus_log_base,MPFR_RNDN);
	mpfr_ui_div(minus_log_base,1,minus_log_base,MPFR_RNDN);

	mpfr_t log_base_power;
	mpfr_init2(log_base_power,prec); 
	mpfr_set(log_base_power,minus_log_base,MPFR_RNDN);

	//mpfr_set_ui(temp1,0,MPFR_RNDN);
	
	//mpfr_mul(double_pole_integral,temp1,incomplete_gamma_factor,MPFR_RNDN); 
	mpfr_log(temp1,base,MPFR_RNDN); 
	mpfr_mul(double_pole_integral,incomplete_gamma_factor,temp1,MPFR_RNDN);
	mpfr_ui_div(temp1,1,minus_pole_position,MPFR_RNDN);
	mpfr_add(double_pole_integral,double_pole_integral,temp1,MPFR_RNDN); 
	//mpfr_printf("before loop: double_pole_integral = %.32RNf\n",double_pole_integral);


	for(int i=0;i<=pole_order_max;i++){
		mpfr_init2(result[i],prec);
		mpfr_set(result[i],double_pole_integral,MPFR_RNDN); 
		if(i<pole_order_max){
			mpfr_mul(double_pole_integral,double_pole_integral,pole_position,MPFR_RNDN); 
		} 
	}

#ifdef debug_mode
	for(int i=0;i<=pole_order_max;i++){
		mpfr_printf("result[%d] = %.16RNf\n",i,result[i],MPFR_RNDN);
	} 
#endif

	mpfr_t * factorial_times_power_lnb;
	mpfr_t * single_pole_coeffs;
	/*
	 * x**1/(x+a)**2 case
	 * */
	
	if(pole_order_max >= 1){ 
		single_pole_coeffs=malloc(sizeof(mpfr_t)*(pole_order_max));
		/*  x^(j+1)=( 
		 *  single_pole_coeffs[0] x^(j-1) + 
		 *  single_pole_coeffs[1] x^(j-2) + ...  +
		 *  single_pole_coeffs[j-1] x^0
		 *  )*(x-a)^2  +
		 *
		 *  single_pole_coeffs[j](x-a) * ((x-a) + a)
		 *
		 *  + a^(j+1)
		 *
		 *  => single_pole_coeffs[j+1]
		 *
		 *  single_pole_coeffs[j+1] = single_pole_coeffs[j]*a + a^j+1
		 * single_pole_coeffs[0] =  
		 * */
		if(pole_order_max>=2){
			factorial_times_power_lnb=malloc(sizeof(mpfr_t)*(pole_order_max-1));
			mpfr_init2(factorial_times_power_lnb[0],prec);
			mpfr_set(factorial_times_power_lnb[0],minus_log_base,MPFR_RNDN);
		}

		mpfr_set(temp1,pole_position,MPFR_RNDN);		
		/* below temp1 is used as pole_position^j*/

		mpfr_init2(single_pole_coeffs[0],prec);
		mpfr_set_ui(single_pole_coeffs[0], 1 ,MPFR_RNDN);
		mpfr_add(result[1],result[1],incomplete_gamma_factor,MPFR_RNDN);
		
		for(int j=1;j<=pole_order_max-1;j++){
			mpfr_init2(single_pole_coeffs[j],prec);
			mpfr_mul(single_pole_coeffs[j],single_pole_coeffs[j-1],pole_position,MPFR_RNDN);
			mpfr_add(single_pole_coeffs[j],single_pole_coeffs[j],temp1,MPFR_RNDN);

			mpfr_mul(temp2,single_pole_coeffs[j],incomplete_gamma_factor,MPFR_RNDN);
			mpfr_add(result[j+1],result[j+1],temp2,MPFR_RNDN);
			if(j<=pole_order_max-2){
				mpfr_init2(factorial_times_power_lnb[j],prec);
				mpfr_mul(temp1,temp1,pole_position,MPFR_RNDN);

				mpfr_mul(factorial_times_power_lnb[j],factorial_times_power_lnb[j-1],minus_log_base,MPFR_RNDN); 
				mpfr_mul_ui(factorial_times_power_lnb[j],factorial_times_power_lnb[j],j,MPFR_RNDN);
			}
		}
	}
#ifdef debug_mode
	for(int j=0;j<=pole_order_max-1;j++){
		mpfr_printf("single_pole_coeffs[%d] = %.16RNf\n",j,single_pole_coeffs[j]);
	}
#endif
#ifdef debug_mode
	for(int j=0;j<=pole_order_max-2;j++){
		mpfr_printf("factorial_times_power_lnb[%d] = %.16RNf\n",j,factorial_times_power_lnb[j]);
	}
#endif
	for(int j=0;j<=pole_order_max-2;j++){
		for(int k=0;k<=j;k++){
			mpfr_mul(temp2,factorial_times_power_lnb[k],single_pole_coeffs[j-k],MPFR_RNDN);
			mpfr_add(result[j+2],result[j+2],temp2,MPFR_RNDN);
		}
	}

	for(int i=0;i<=pole_order_max-1;i++){
		mpfr_clear(single_pole_coeffs[i]);
	}
	for(int i=0;i<=pole_order_max-2;i++){
		mpfr_clear(factorial_times_power_lnb[i]);
	}
	if(pole_order_max > 0){
		free(single_pole_coeffs);
	}
	if(pole_order_max > 1){
		free(factorial_times_power_lnb);
	}
	mpfr_clear(temp1);
	mpfr_clear(double_pole_integral);
	mpfr_clear(temp2);
	mpfr_clear(minus_pole_position);
	//mpfr_clear(factorial);
	mpfr_clear(minus_log_base);
	mpfr_clear(log_base_power);
	return result;
} 
