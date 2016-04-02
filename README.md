This tiny module for [SageMath](http://www.sagemath.org) contains helper functions to generate a bootstrap problem to be solved by [SDPB](https://github.com/davidsd/sdpb), written originally for Ph.D. thesis of the author (who was almost driven mad by a buggy proprietary software). Almost the same functionality is offered by [PyCFTBoot](https://github.com/cbehan/pycftboot), but the implementation detail and notation is somewhat different. This makes extensive use of "rings.RealNumber" and "rings.RealDensePolynomials" classes contained in Sage, to handle arbitrary precision number/polynomials.

###Install

1. Install sage from source. This module has been originally written in the environment of sage-6.8 on OSX Yosemite, but confirmed to work up to sage-7.0.

2. Place cboot directory under 
	`your/sage/src/sage`.
The actual location of `your/sage/src` can be checked by running sage and enter 
	`SAGE_ENV['SAGE_SRC']` 
In my environment, for example, this is '/Users/tomoki/sources/sage-6.8/src'.

3. Edit `your/sage/src/module_list.py` to add the following items in the variable `ext_modules`.  
    ```
    Extension('sage.cboot.context_object',sources=['sage/cboot/context_object.pyx','sage/cboot/partial_fraction.c','sage/cboot/integral_decomp.c','sage/cboot/chol_and_inverse.c','sage/cboot/context_variables.c'], extra_compile_args=['-std=c99'],libraries=['gmp','mpfr'],language = 'c'),

    Extension('sage.cboot.scalar_context',sources=['sage/cboot/scalar/scalar_context.pyx','sage/cboot/scalar/hor_formula.c','sage/cboot/scalar/hor_recursion.c','sage/cboot/scalar/k_compute.c'],include_dirs=['sage/cboot','sage/cboot/scalar'],depends=['sage/cboot/context_variables.h'], extra_compile_args=['-std=c99'],libraries=['gmp','mpfr'],language = 'c'),
    ```
4. Then run from your terminal `sage -b`.


###Examples 
Examples are provided in ipython notebooks. These include
1. The 3d Ising model (single correlator)
2. The 3d XY model (single correlator)
3. The 3d Ising model (mixed correlator)

Supersymmetric examples are under construction.
