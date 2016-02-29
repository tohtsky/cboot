import sage.cboot.scalar as cb
import numpy as np
import re

context=cb.context_for_scalar(epsilon=0.5,)

#Store the conformal blocks
def store_conformal_blocks(lmax,nu_max,context):
    res={}
    for spin in range(0,lmax+1):
        g=context.rational_approx_data(nu_max,spin).approx()
        res.update({spin:g})
    return res

g_stored=store_conformal_blocks(25,15,context)


def g_to_F(delta,sector,g,NSO,shift=None):
    # We should distinguish between constant vector (used e.g. for c-minimization)
    # and rational-approximated conformal block
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

def single_continuous_sector(delta,sector,spin,modify_list,block_store,NSO):
    """
    modify_list argument should be a dictionary to specify where we insert hypothetical gaps.
    For example, if we assume a gap 1.8 in O(N)-singlet scalar operator, we set
        modify_list = {("S",0):1.8}        
    """
    try:
        shift=context(modify_list[(sector,spin)])
    except KeyError:
        if spin==0:
            shift=context.epsilon
        else:
            shift=spin+context.epsilon*2
    return g_to_F(delta,sector,block_store[spin],NSO,shift=shift)


def O_n_problem(delta,modify_list,block_store,NSO=2):
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
            spins=[spin for spin in block_store.keys() if not spin%2]
        
        for spin in spins:
            try:
                shift=context(modify_list[(sector,spin)])
            except KeyError:
                if spin==0:
                    shift=context.epsilon
                else:
                    shift=spin+context.epsilon*2
            pvms.append(g_to_F(delta,sector,block_store[spin],NSO,shift=shift))
    norm=g_to_F(delta,"S",context.gBlock(0,0,0,0),NSO)
    return context.SDP(norm,norm*0,pvms)



