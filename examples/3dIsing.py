from sage.all import 
import re
# The purpose of this tutorial is to let you see around basic functions in the module.
from subprocess import Popen, PIPE # To execute sdpb
from sage.cboot.scalar import context_for_scalar

# Define the context for the bootstrap problems.
# Lambda controls the maximal number of derivatives acting on conformal blocks.
# epsilon is (d-2)/2, so 0 for 2d, 0.5 for 3d, and 1 for 4d.
# nMax represents the order of expansion in the Hogervorst-Osborn-Rychkov series
context=context_for_scalar(Lambda=15,epsilon=0.5,nMax=250)

def fixed_spin_contributuion(delta, spin,shift,nu_max=15):
    # rational_approx_data specifies a data to approximate a conformal block of given spin.
    # nu_max is the cutoff for the included pole.
    # The last two argumensts are to be used for 4pt functions with non-identical dimensions
    q=context.rational_approx_data(nu_max,spin,0,0)
    
    # Compute the approximation.
    # Returned instance is "prefactor_numerator" instance, which stores
    # the positive prefactor ("damped_rational") and the numerator polynomial (numpy ndarray)
    approx=q.approx()
    
    # Shift by the unitarity bound / hypothetical gap
    # The method "shift" of the prefactor_numerator shift both the prefactor and the body.
    approx=approx.shift(shift)
    
    # Apply the matrix which represents v^{\Delta _ \phi}-multplication and anti-symmetrization
    # Attributes "matrix" is the numerator polynomial
    F_body=context.make_F_minus_matrix(delta).dot(approx.matrix)
    
    # Reconstruct with convolved conformal block numerator.
    return context.prefactor_numerator(approx.prefactor,F_body.reshape((1,1,len(F_body))), )

def spin_vs_shift(lmax,modify_list):
    # Returns a dictionary giving the amount of shift for the conformal block.
    # Default amount of shift should be the unitarity bound.
    ubs=dict([(x,x+2*context.epsilon) for x in range(0,lmax+1,2)])    
    # Modify unitarity bounds to hypothetical gaps.
    updates=dict([(x[0],context(x[1])) for x in modify_list.items()])
    ubs.update(updates)
    return ubs

def ising_singlet_bound_SDP(delta,modify_list,lmax=30,nu_max=15):
    delta=context(delta) # convert delta into a multi-precision number.
    shifts=spin_vs_shift(lmax,modify_list)
    polynomial_Vector_Matrices=\
    [fixed_spin_contributuion(delta,x,shifts[x],nu_max=nu_max) for x in shifts]
    
    # context.gBlock returns the conformal block derivative table 
    # as the 1d numpy array.
    normalization=context.make_F_minus_matrix(delta).dot(context.gBlock(0,0,0,0))
    obj=normalization*0
    return context.SDP(normalization,obj,polynomial_Vector_Matrices)

m=ising_singlet_bound_SDP(0.518,{0:1.5}) 
sdpbres=Popen(["sdpb","-s","test.xml","--findPrimalFeasible",\
       "--findDualFeasible","--noFinalCheckpoint"],\
                 stderr=PIPE,stdout=PIPE).communicate()
