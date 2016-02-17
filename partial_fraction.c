#include "partial_fraction.h"
#define deb(num)\
	printf("deb %f\n",num)

mpfr_t* fast_partial_fraction_c(mpfr_t* pole_locations, int* double_or_single, int n_poles,  mpfr_prec_t prec){
	int expected_result_length=n_poles;
	for(int i=0;i<n_poles;i++){
		//mpfr_printf("given pole_locations[%d]=%.16RNf\n",i,pole_locations[i]);
		//printf("given double_or_single[%d]=%d\n",i,double_or_single[i]);
		if(double_or_single[i]){
			++expected_result_length;	
		}	
	}
	mpfr_t* result=malloc(sizeof(mpfr_t)*expected_result_length);
	mpfr_t temp1;
	mpfr_init2(temp1,prec);
	int count_result_location=0;
	for(int i=0;i<n_poles;i++){
		mpfr_init2(result[count_result_location],prec);
		mpfr_set_ui(result[count_result_location],1,MPFR_RNDN);
		for(int j=0;j<n_poles;j++){
			if(i!=j){
				mpfr_sub(temp1,pole_locations[i],pole_locations[j],MPFR_RNDN);	
				mpfr_mul(result[count_result_location],result[count_result_location],temp1,MPFR_RNDN);
				if(double_or_single[j]){
					mpfr_mul(result[count_result_location],result[count_result_location],temp1,MPFR_RNDN);	
				}
			}
		}
		mpfr_ui_div(result[count_result_location],1,result[count_result_location],MPFR_RNDN); 
		if(double_or_single[i]){
			++count_result_location; 
			mpfr_init2(result[count_result_location],prec);
			//mpfr_set_ui(result[count_result_location],0,MPFR_RNDN);
			mpfr_set_zero(result[count_result_location],1);
			for(int j=0;j<n_poles;j++){
				if(i!=j){ 
					mpfr_sub(temp1,pole_locations[i],pole_locations[j],MPFR_RNDN);	
					if(double_or_single[j]){ 
						mpfr_ui_div(temp1,2,temp1,MPFR_RNDN);
					}
					else{
						mpfr_ui_div(temp1,1,temp1,MPFR_RNDN);
					}
					mpfr_add(result[count_result_location], result[count_result_location], temp1,MPFR_RNDN);
				}
			
			}
			mpfr_mul(result[count_result_location],result[count_result_location],result[count_result_location-1],MPFR_RNDN);
			mpfr_neg(result[count_result_location],result[count_result_location],MPFR_RNDN);
		}
		++count_result_location;
	}
	mpfr_clear(temp1);
	return result;
}
