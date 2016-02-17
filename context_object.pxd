from sage.libs.mpfr cimport *
from sage.rings.real_mpfr cimport *
cimport numpy as np

cdef extern from "stdlib.h":
    void* malloc(size_t size)
    void free (void* ptr)

cdef extern from "sage/cboot/integral_decomp.h":
    mpfr_t* simple_pole_case_c(long pole_order_max, mpfr_t base, mpfr_t pole_position, mpfr_t incomplete_gamma_factor, mp_prec_t prec)
    mpfr_t* double_pole_case_c(long pole_order_max, mpfr_t base, mpfr_t pole_position, mpfr_t incomplete_gamma_factor, mp_prec_t prec)

cdef extern from "sage/cboot/partial_fraction.h":
    mpfr_t* fast_partial_fraction_c(mpfr_t* pole_locations, int* double_or_single, int expected_result_length, mp_prec_t prec)
   
cdef extern from "sage/cboot/chol_and_inverse.h":
    mpfr_t* mpfr_triangular_inverse(mpfr_t* A, int dim,mp_prec_t prec)
    mpfr_t* mpfr_cholesky(mpfr_t* A, int dim,mp_prec_t prec) 
    mpfr_t* form_anti_band(mpfr_t* ab_vector, int dim, mp_prec_t prec)

cdef extern from "sage/cboot/context_variables.h":
    ctypedef struct cb_context:
        mpfr_t* rho_to_z_matrix
    cb_context context_construct(long nMax, mp_prec_t prec, int Lambda)
    void clear_cb_context(cb_context context)
#
#cdef class clevel_context:
#    cdef cb_context context
#    cdef public np.ndarray zrmat
#
cdef class cb_universal_context:
    cdef cb_context c_context
    cdef public mp_prec_t precision
    cdef public RealField_class field 
    cdef public object Delta_Field
    cdef public object Delta
    cdef public int Lambda
    cdef public RealNumber rho
    cdef public long maxExpansionOrder
    cdef public object polynomial_vector_shift
    cdef public object polynomial_vector_evaluate
    cdef public object convert_to_polynomial_vector
    cdef public object convert_to_real_vector
    cdef public object rho_to_z_matrix 
    cdef public object zzbar_to_xy_marix
    cdef public object index_list
    cdef public object rho_to_delta
    cdef public object null_ftype
    cdef public object null_htype 

cdef class damped_rational:
    cdef public object poles
    cdef public RealNumber base
    cdef public RealNumber pref_constant
    cdef public cb_universal_context context 

cdef class positive_matrix_with_prefactor:
    cdef public damped_rational prefactor 
    cdef public cb_universal_context context 
    cdef public object matrix

cdef class prefactor_numerator(positive_matrix_with_prefactor):
    pass
