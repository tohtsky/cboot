import sage.cboot as cb
import numpy as np
from subprocess import Popen, PIPE
import re

context=cb.context_for_scalar(epsilon=0.5,Lambda=13) 
lmax=25
nu_max=12
cbs={}
for spin in range(0,lmax):
    g=context.approx_cb(nu_max,spin)
    cbs.update({spin:g}) 

def make_F(delta,sector,spin,gap_dict,NSO):
    delta=context(delta)
    try:
        gap=context(gap_dict[(sector,spin)])
    except KeyError:
        if spin==0:
            gap=context.epsilon
        else:
            gap=2*context.epsilon+spin
    g_shift=cbs[spin].shift(gap)

    g_num=g_shift.matrix
    g_pref=g_shift.prefactor
    F=context.make_F_minus_matrix(delta).dot(g_num)
    H=context.make_F_plus_matrix(delta).dot(g_num) 

    if sector=="S":
        num=np.concatenate((context.null_ftype,F,H))
        
    elif sector=="T":
        num=np.concatenate((F,(1-2/context(NSO))*F,-(1+2/context(NSO))*H))

    elif sector=="A":
        num=np.concatenate((-F,F,-H))
        
    return context.prefactor_numerator(g_pref,num)

def make_SDP(delta,gap_dict,NSO=2):
    delta=context(delta)
    pvms=[]
    for sector in ("S","T","A"):
        if sector is not "A":
            spins=[spin for spin in cbs.keys() if not spin%2]
        else:
            spins=[spin for spin in cbs.keys() if spin%2] 
        for spin in spins:
            pvms.append(make_F(delta,sector,spin,gap_dict,NSO))

    norm_F=context.make_F_minus_matrix(delta).dot(context.gBlock(0,0,0,0))
    norm_H=context.make_F_plus_matrix(delta).dot(context.gBlock(0,0,0,0))
    norm=np.concatenate((context.null_ftype,norm_F,norm_H))

    obj=norm*0
    return context.SDP(norm,obj,pvms) 

sdpb="sdpb"
sdpbparams=["--findPrimalFeasible","--findDualFeasible","--noFinalCheckpoint"]

def bs(delta,upper=3,lower=1,sector="S",sdp_method=make_SDP,NSO=2):
    upper=context(upper)
    lower=context(lower)
    while upper - lower > 0.001:
        D_try=(upper+lower)/2
        prob=sdp_method(delta,{(sector,0):D_try},NSO=NSO)
        prob.write("3d_Ising_binary.xml")
        sdpbargs=[sdpb,"-s","3d_Ising_binary.xml"]+sdpbparams
        out, err=Popen(sdpbargs,stdout=PIPE,stderr=PIPE).communicate()
        sol=re.compile(r'found ([^ ]+) feasible').search(out).groups()[0] 
        if sol=="dual":
            print("(Delta_phi, Delta_{1})={0} is excluded."\
            .format((float(delta),float(D_try)),sector)) 
            upper=D_try
        elif sol=="primal":
            print("(Delta_phi, Delta_{1})={0} is not excluded."\
            .format((float(delta),float(D_try)),sector)) 
            lower=D_try
        else:
            raise RuntimeError("Unexpected return from sdpb") 
    return upper


if __name__=="__main__":
    # The default example
    print bs(0.52)

    # ======================================
    # if you want to derive the bound on Delta_T 
    #print bs(0.52,sector="T")


