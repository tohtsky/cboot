from sage.all import *
import sage.cboot.scalar as cb
import numpy as np
import re


#Store the conformal blocks
def store_conformal_blocks(lmax,nu_max,context):
    res={}
    for spin in range(0,lmax+1):
        g=context.rational_approx_data(nu_max,spin).approx()
        res.update({spin:g})
    return res


def g_to_F(delta,sector,g,NSO,shift=None):
    # Setup a function which returns a vector included in O(N) sum-rule.

    # We should distinguish between constant vector
    #(used e.g. for c-minimization) and rational-approximated conformal block
    if shift==None:
        g_body=g
    else:
        g_body=g.matrix
        
    # Multiply with a matrix representing a convolution with v^{\delta}
    # accompanied by (anti)-symmetrization    
    F=context.make_F_minus_matrix(delta).dot(g_body)
    H=context.make_F_plus_matrix(delta).dot(g_body)

    if sector=="S":
        # context.null_ftype is a null vector with dimension that of F^{(-)}
        body=np.concatenate((context.null_ftype,F,H))
        
    elif sector=="T":
        body=np.concatenate((F,(1-2/context(NSO))*F,-(1+2/context(NSO))*H))

    elif sector=="A":
        body=np.concatenate((-F,F,-H))
        
    if shift==None:
        return body
    else:
        return context.prefactor_numerator(g.prefactor,body).shift(shift)

def O_n_problem(delta,modify_list,block_store,NSO=2,norm_point=("S",0,0),
        obj_point=None):
    """
    modify_list argument should be a dictionary to specify where we insert hypothetical gaps.
    For example, if we assume a gap 1.8 in O(N)-singlet scalar operator, we set
        modify_list = {("S",0):1.8}        
    """
    delta=context(delta)
    pvms=[]
    for sector in ("S","T","A"):
        if sector is not "A":
            spins=[spin for spin in block_store.keys() if not spin%2]
        else:
            spins=[spin for spin in block_store.keys() if spin%2] 
        for spin in spins:
            try:
                shift=context(modify_list[(sector,spin)])
            except KeyError:
                if spin==0:
                    shift=context.epsilon
                else:
                    shift=spin+context.epsilon*2
            pvms.append(g_to_F(delta,sector,block_store[spin],NSO,shift=shift))

    norm=g_to_F(delta,norm_point[0],
            context.gBlock(norm_point[1],norm_point[2],0,0),NSO)

    if not obj_point:
        obj=norm*0
    else:
        obj=g_to_F(delta,obj_point[0],
                context.gBlock(obj_point[1],obj_point[2],0,0),NSO) 

    return context.SDP(norm,obj,pvms) 

if __name__=="__main__":
    context=cb.context_for_scalar(epsilon=0.5,)
    g_stored=store_conformal_blocks(25,15,context) 
    prob=O_n_problem(0.52,{("S",0):2.0},g_stored)
    prob.write("on.xml") 

