from sage.all import *
import sage.cboot as cb
import numpy as np
from subprocess import Popen, PIPE
import re


context=cb.context_for_scalar(epsilon=0.5,Lambda=13) 
lmax=20
nu_max=12

sdpb="sdpb"
sdpbparams=["--findPrimalFeasible","--findDualFeasible","--noFinalCheckpoint"]

def cSCblock(nu_max,spin,Delta=None):
    d=context.epsilon*2+2
    epsilon=context.epsilon 
    pole_add_data=[((1,1),[-spin-1],[-spin],(spin+2*epsilon)/4/(spin+epsilon)),\
            ((0,2),[epsilon,epsilon-1,context(-1-spin),spin+d-3],\
            [0,d-3,context(-spin),context(spin+2*epsilon)],1/context(16))]
    if spin > 0:
        pole_add_data.append(((-1,1),[spin+d-3],[spin+2*epsilon]\
                ,context(spin)/context(4*(spin+epsilon))))
    if Delta ==None:
        res_g=context.approx_cb(nu_max,spin,0,0)
        res_g_tilde=res_g
        bosonic = res_g 
    else:
        Delta=context(Delta)
        res_g=context.gBlock(spin,Delta,0,0)
        res_g_tilde=np.copy(res_g)
        bosonic=np.copy(res_g)
        if Delta==0:
            return (res_g,res_g_tilde,bosonic)
    for shift,poles,factors,C in pole_add_data:
        if Delta==None:
            g=context.approx_cb(nu_max,spin+shift[0]).shift(shift[1])\
                    .multiply_factored_rational(poles,factors,C)
        else:
            g=context.gBlock(spin+shift[0],Delta+shift[1],0,0)
            for x in poles:
                g=g/(Delta-context(x))
            for x in factors:
                g*=(Delta-context(x))
            g=g*C
        res_g+=g
        if shift[0]%2:
            res_g_tilde=res_g_tilde-g
        else:
            res_g_tilde=res_g_tilde+g
    return (res_g, res_g_tilde, bosonic)

cbs={}
for spin in range(0,lmax):
    cbs.update({spin:cSCblock(nu_max,spin)})


def make_F(delta,sector,spin,gap_dict,Delta=None):
    delta=context(delta)
    if Delta==None:
        try:
            gap=context(gap_dict[(sector,spin)])
        except KeyError:
            if sector=="0+" or sector=="0-":
                gap=spin+context.epsilon*2
            elif sector=="2":
                gap=abs(2*delta-2*context.epsilon-1)+spin+(2*context.epsilon+1) 
        gs=[x.shift(gap) for x in cbs[spin]]
    else:
        Delta=context(Delta)
        gs=cSCblock(0,spin,Delta=Delta) 

    F=context.dot(context.F_minus_matrix(delta),gs[0])
    H=context.dot(context.F_plus_matrix(delta),gs[0])
    Ftilde=context.dot(context.F_minus_matrix(delta),gs[1])
    Htilde=context.dot(context.F_plus_matrix(delta),gs[1])
    Fbos=context.dot(context.F_minus_matrix(delta),gs[2])
    Hbos=context.dot(context.F_plus_matrix(delta),gs[2])

    if sector=="0+":
        return [F,Ftilde,Htilde]
    elif sector=="0-":
        return [F,-Ftilde,-Htilde]
    elif sector=="2":
        return [0,Fbos,-Hbos]

def make_SDP(delta,gap_dict,norm_point=("0+",0,0),obj_point=None):
    delta=context(delta)
    pvms=[]
    for spin in cbs:
        if not spin%2:
            pvms.append(make_F(delta,"0+",spin,gap_dict))
            pvms.append(make_F(delta,"2",spin,gap_dict))
            pvms.append(make_F(delta,"2",spin,gap_dict,Delta=2*delta+spin))

        else:
            pvms.append(make_F(delta,"0-",spin,gap_dict))
    
    if delta < (2*context.epsilon + 2)/4:
        pvms.append(make_F(delta,"2",0,gap_dict,\
                Delta=2*context.epsilon+2-2*delta))

    norm=make_F(delta,norm_point[0],norm_point[1],{},Delta=norm_point[2])
    if not obj_point:
        obj=0 
    else:
        obj=make_F(delta,obj_point[0],obj_point[1],{},Delta=obj_point[2])

    return context.sumrule_to_SDP(norm,obj,pvms)

def bs(delta,upper=3,lower=1,sector="0+",sdp_method=make_SDP):
    upper=context(upper)
    lower=context(lower)
    while upper - lower > 0.001:
        D_try=(upper+lower)/2
        prob=sdp_method(delta,{(sector,0):D_try})
        prob.write("3d_sc_binary.xml")
        sdpbargs=[sdpb,"-s","3d_sc_binary.xml"]+sdpbparams
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

def cc(delta):
    prob=make_SDP(delta,{},norm_point=("0-",1,1+2*context.epsilon),\
            obj_point=("0+",0,0)) 
    prob.write("3dsc_cc.xml")
    sdpbargs=[sdpb,"-s","3dsc_cc.xml","--noFinalCheckpoint"]

    out,err=Popen(sdpbargs,stdout=PIPE,stderr=PIPE).communicate()
    sol=re.compile(r'primalObjective *= *([^ ]+) *$',re.MULTILINE)\
            .search(out).groups()[0]
    return -delta**2/float(sol)


if __name__=='__main__': 

    # binary search for the N=2 super Ising model
    print bs(0.6666)

    # ===============================================
    # If you want to derive the central charge lower bound,
    # uncomment the following.
    #print cc(0.6666)
