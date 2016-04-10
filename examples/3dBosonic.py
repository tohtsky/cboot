import sage.cboot as cb
context=cb.context_for_scalar(epsilon=0.5,Lambda=13)

lmax=25
nu_max=12
cbs={}
for spin in range(0,lmax,2):
    g=context.approx_cb(nu_max,spin)
    cbs.update({spin:g}) 

def make_F(delta,spin,gap_dict):
    mat_F=context.F_minus_matrix(delta)
    try:
        gap=context(gap_dict[spin])
    except KeyError:
        if spin==0:
            gap=context.epsilon
        else:
            gap=2*context.epsilon+spin
    g_shift=cbs[spin].shift(gap)
    F=context.dot(mat_F,g_shift)
    return F 

def make_SDP(delta,gap_dict): 
    delta=context(delta)
    Fs=[make_F(delta,spin,gap_dict) for spin in cbs.keys()]
    mat_F=context.make_F_minus_matrix(delta)
    norm=context.dot(mat_F,context.gBlock(0,0,0,0))
    obj=norm*0
    return context.SDP(norm,obj,Fs) 

from subprocess import Popen, PIPE
import re

sdpb="sdpb"
sdpbparams=["--findPrimalFeasible","--findDualFeasible","--noFinalCheckpoint"] 
def bs(delta,upper=3,lower=1,sdp_method=make_SDP):
    upper=context(upper)
    lower=context(lower)
    while upper - lower > 0.001:
        D_try=(upper+lower)/2
        prob=sdp_method(delta,{0:D_try})
        prob.write("3d_Ising_binary.xml")
        sdpbargs=[sdpb,"-s","3d_Ising_binary.xml"]+sdpbparams
        out, err=Popen(sdpbargs,stdout=PIPE,stderr=PIPE).communicate()
        sol=re.compile(r'found ([^ ]+) feasible').search(out).groups()[0] 
        if sol=="dual":
            print("(Delta_phi, Delta_epsilon)={0} is excluded."\
            .format((float(delta),float(D_try)))) 
            upper=D_try
        elif sol=="primal":
            print("(Delta_phi, Delta_epsilon)={0} is not excluded."\
            .format((float(delta),float(D_try)))) 
            lower=D_try
        else:
            raise RuntimeError("Unexpected return from sdpb") 
    return upper

def make_SDP_for_cc(delta,gap_dict={0:1}):
    delta=context(delta)
    Fs=[make_F(delta,spin,gap_dict) for spin in cbs.keys()]
    mat_F=context.make_F_minus_matrix(delta)
    norm=context.dot(mat_F,context.gBlock(2,3,0,0))
    obj=context.dot(mat_F,context.gBlock(0,0,0,0))
    return context.SDP(norm,obj,Fs)

def cc(delta):
    prob=make_SDP_for_cc(delta)
    prob.write("3d_Ising_cc.xml")
    sdpbargs=[sdpb,"-s","3d_Ising_cc.xml","--noFinalCheckpoint"]
    out, err=Popen(sdpbargs,stdout=PIPE,stderr=PIPE).communicate()
    sol=re.compile(r'primalObjective *= *([^ ]+) *$',re.MULTILINE)\
            .search(out).groups()[0]
    return -delta**2/float(sol)

def make_SDP_epsilon_prime(delta,gap_dict): 
    delta=context(delta)
    Fs=[make_F(delta,spin,gap_dict) for spin in cbs.keys()]
    mat_F=context.make_F_minus_matrix(delta)
    Fs+=[context.dot(mat_F,context.gBlock(0,delta,0,0))]
    norm=context.dot(mat_F,context.gBlock(0,0,0,0))
    obj=norm*0
    return context.SDP(norm,obj,Fs) 

if __name__=='__main__':
    # The default example
    delta=0.518
    print(bs(delta))

    # ===============================================
    # if you want to derive the central charge lower bound,
    # uncomment the follwing lines.
    #delta=0.518
    #print("central charge lower bound at delta={0} is {1}"\
    #        .format(delta,cc(delta)))

    # ===============================================
    # The upper bound on epsilon' dimension.
    #Delta_epsilon=0.8
    #print(bs(Delta_epsilon,sdp_method=make_SDP_epsilon_prime))
