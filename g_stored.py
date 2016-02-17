from sage.all import *
import sage.cboot.scalar as cb
import numpy as np
from sympy.abc import F,H
from sympy import expand
from copy import copy

class cached_cb:
    def __init__(self,context):
        self.context=context
        self.is_initialized = 0
        self.g=[]
        self.spins=[]
    def compute_table(self,l_max_arg,nu_max):
        l_max=2*(l_max_arg//2)+1
        spins=(range(0,l_max,2)+range(l_max+3,2*l_max,4))+ range(1,l_max,2)+range(l_max+2,2*l_max,4)
        spins=sorted(spins)
        def func(ell):
            g=self.context.rational_approx_data(nu_max-(ell//6),ell,is_correlator_multiple=False)
            return [ell,g.prefactor(),g.approx_g()]
        self.g=map(func,spins) 
        print "initialization done"
        self.spins=spins
        self.is_initialized=1
    def index_for_spin(self,spin):
        if self.is_initialized==0:
            print "g not initialized"
            return -1
        else:
            temp_list=[[abs(self.spins[i]-spin),i] for i in range(0,len(self.spins))]
            temp_list=sorted(temp_list)
            if not temp_list[0][0]==0:
                print "conformal block of such a pin is not stored"
                return -1
            else:
                return temp_list[0][1] 
    def unitarity_bounds(self,option,modify_list=[],scalar_default=None):
        if scalar_default==None:
            scalar_default=self.context.epsilon
        special_spins=[x[0] for x in modify_list]
        default_ub=lambda x: (x+2*self.context.epsilon) if not x==0 else scalar_default
        return sorted([[x, default_ub(x)] for x in self.spins if (not (x+option)%2 and not x in special_spins)] + [x for x in modify_list if (not (x[0]+option)%2)]) 
    def __call__(self,spin):
        if self.is_initialized==0:
            raise RuntimeError("g has not been initialized!")
        else:
            index=self.index_for_spin(spin)
            return (self.g)[index]

class sumrule:
    def __init__(self,coeffs,types,sector_names,context,sector_wise_option=None):
        self.context=context
        self.types=types
        self.coeffs=[[context.field(x) for x in y] for y in coeffs]
        self.coeffs=dict(zip(sector_names,self.coeffs))
        if not (len(sector_names)==len(set(sector_names))):
            raise RuntimeError("duplicate in sector_names") 
        self.num_sectors=len(sector_names)
        if not (len(coeffs)==self.num_sectors):
            raise RuntimeError("number of sector_names in appropriate?")
        self.num_equalities=len(coeffs[0])
        self.sector_names=sector_names
        self.nullMinus=context.null_ftype
        self.nullPlus=context.null_htype
        self.dim=sum(map(lambda x: len(context.null_ftype) if (x=="f_minus" or x=="F") else len(context.null_htype), self.types))
        if sector_wise_option==None:
            self.sector_wise_option=dict(zip(self.sector_names,[None for x in range(0,self.num_sectors)]))
        else:
            if not len(sector_wise_option)==self.num_sectors:
                raise RuntimeError("option ")
            else:
                self.sector_wise_option=dict(zip(self.sector_names,sector_wise_option))

    def __call__(self,sector,pref,g,delta):
        fplus=self.context.v_to_d_and_symmetrizing_matrix(delta).dot(g)
        fminus=self.context.v_to_d_and_anti_symmetrizing_matrix(delta).dot(g)
        coeffs=self.coeffs[sector] 
        res=[]
        for x in zip(coeffs,self.types):
            if x[0]==0:
                if x[1]=="f_plus" or x[1]=="H":
                    res.append(self.context.null_htype)
                 #   res=np.concatenate((res,self.context.null_ftype))
                elif x[1]=="f_minus" or x[1]=="F":
                    res.append(self.context.null_ftype)
            else:
                if x[1]=="f_plus" or x[1]=="H":
                    res.append(fplus*x[0])
                elif x[1]=="f_minus" or x[1]=="F":
                    res.append(fminus*x[0])

        return self.context.positive_matrix_with_prefactor(pref, np.concatenate(res).reshape([1,1,self.dim])) 

    def sector_wise_problem(self,sector,delta,gs,modify_list=[],**kwargs): 
        """
        elements of modify_list takes the form: 
        [spin, modified_lowerbound ]
        """
        #print ("sector="+sector)
        #print ("modify_list="+repr(modify_list))
        ubs=gs.unitarity_bounds(self.sector_wise_option[sector],modify_list,**kwargs)
        res=[]
        for ub in ubs:
            #print (sector+repr(ub))
            gdata=gs(ub[0])
            pref=gdata[1].shift(ub[1])
            g=gs.context.polynomial_vector_shift(gdata[2],ub[1])
            res.append(self(sector,pref,g,delta))
        return res 
    def vector_at_sector(self,sector,delta,g):
        fplus=self.context.v_to_d_and_symmetrizing_matrix(delta).dot(g)
        fminus=self.context.v_to_d_and_anti_symmetrizing_matrix(delta).dot(g)
        coeffs=self.coeffs[sector] 
        res=[]
        for x in zip(coeffs,self.types):
            if x[0]==0:
                if x[1]=="f_plus" or x[1]=="H":
                    res.append(self.context.null_htype)
                 #   res=np.concatenate((res,self.context.null_ftype))
                elif x[1]=="f_minus" or x[1]=="F":
                    res.append(self.context.null_ftype)
                    #res=np.concatenate((res,self.context.null_htype))
            else:
                if x[1]=="f_plus" or x[1]=="H":
                    res.append(fplus*x[0])
                    #res=np.concatenate((res,fplus*x[0]))
                elif x[1]=="f_minus" or x[1]=="F":
                    res.append(fminus*x[0])
                    #res=np.concatenate((res,fminus*x[0]))
        #res=np.concatenate(res)
        return np.concatenate(res)

    def entire_continuous(self,delta,gs,modify_list_all,**kwargs):
        """
        elements of modify_list takes the form: 
        [sector_name, spin, modified_lowerbound ]
        """
        #print repr(modify_list_all)
        res=[]
        for y in self.sector_wise_option.items():
            modify_list=[x[1:3] for x in modify_list_all if x[0]==y[0]]
            #print ("y[0]="+repr(y[0]))
            #print ("modify_list at \"entire_continuous\"="+repr(modify_list))
            res=res+self.sector_wise_problem(y[0],delta,gs,modify_list,**kwargs)
        return res 


def ON_sumrule(N,context):
    return sumrule([[0,1,1],[1,(1-2/context.field(N)),-(1+2/context.field(N)),],[-1,1,-1]],["f_minus","f_minus","f_plus"],["S","T","A"],context,sector_wise_option=[0,0,1])

def Ising_sumrule(context):
    return sumrule([[1]],["f_minus"],["S"],context,sector_wise_option=[0])

def ONOM_sumrule(N,M,context):
    NSO1=context.field(N)
    NSO2=context.field(M)
    coeffs=\
    [[1, 0, 0, 0, 0, 1, 0, 0, 0],\
     [-2/NSO2, 1, 1, 0, 0, -2/NSO2, 1, 1, 0],\
      [0, -1, 1, 0, 0, 0, -1, 1, 0],\
       [-2/NSO1, 0, 1, 1, 0, -2/NSO1, 0, -1, 1],\
        [1 + 4/(NSO1*NSO2),\
                  1 - 2/NSO1,\
                    -2/NSO2 - 2/NSO1,\
                      1 - 2/NSO2,\
                        1,\
                          -1 + 4/(NSO1*NSO2),\
                            -1 - 2/NSO1,\
                              2/NSO2 - 2/NSO1,\
                                -1 - 2/NSO2],\
         [1, -1 + 2/NSO1, -2/NSO1, 1, -1, -1, 1 + 2/NSO1, -2/NSO1, -1],\
          [0, 0, 1, -1, 0, 0, 0, -1, -1],\
           [1, 1, -2/NSO2, -1 + 2/NSO2, -1, -1, -1, 2/NSO2, 1 + 2/NSO2],\
            [1, -1, 0, -1, 1, -1, 1, 0, 1]]
    types=["f_minus","f_minus","f_minus","f_minus","f_minus","f_plus","f_plus","f_plus","f_plus"]
    sector_names=[x+y for x in ["S","T","A"] for y in ["S","T","A"]]
    return sumrule(coeffs,types,sector_names,context,sector_wise_option=[0,0,1,0,0,1,1,1,0]) 

def determine_coeffs(x_arg):
    try:
        x=x_arg.expand()
        coeff_f=x.coeff(F)
        type_defined=False
        if not coeff_f ==0:
            parity_type = "F" 
            type_defined=True
            coeff=coeff_f
        coeff_h=x.coeff(H)
        if not coeff_h==0:
            if type_defined: 
                raise RuntimeError("F^{(+)} and F^{(-)} exist simultaneously!")
            parity_type="H"
            coeff=coeff_h
            type_defined=True
        if type_defined:
            return (coeff,parity_type)
        else:
            return (0,None)
    except AttributeError:
        return (0,None)

def fh_check_and_collect(fhtypes,known_types=None):
    if not known_types:
        confirmed=[None for x in fhtypes]
    else:
        confirmed=copy(known_types)
    for (x,y) in enumerate(zip(fhtypes,confirmed)): 
        if y[1]==None:
            if not y[0]==None:
                confirmed[x]=y[0]
        else:
            if (not y[0]==None) and (not y[0]==y[1]):
                raise RuntimeError("The sum rule contains an inconsistent (anti)-symmetrization.")
    return confirmed
#    print(fhtypes)

def sumrule_vector_parse(v,fhtypes):
    coeff_and_fhtype=[determine_coeffs(x) for x in v]
    coeffs,check_fhtypes=list(zip(*coeff_and_fhtype))
    new_fh=fh_check_and_collect(check_fhtypes,known_types=fhtypes)
    return list(coeffs),new_fh
    #print("coeffs={0}".format(coeffs))
    #print("type={0}".format(check_fhtypes))
   
def sumrule_parse(sv):
    #num_constraints=len(sv[0])
    fhtypes=None#[None for x in range(num_constraints)]
    coeffs=[]
    for sector in sv:
        new_coeff, fhtypes=sumrule_vector_parse(sv[sector],fhtypes) 
        coeffs.append(new_coeff) 
    if len({len(x) for x in coeffs})>1:
        raise RuntimeError("The sum rule vector has an non-identical dimension among some sector.") 

    sector_names=sv.keys()
    return (sector_names,coeffs,fhtypes)
    
def single_correlator_sumrule_parse(sv,spins,context):
    (sector_names,coeffs,fhtypes)=sumrule_parse(sv)
    if not sv.keys()==spins.keys():
        raise RuntimeError("The spins and the sum rule have incompatible sector name:\nthe one from spins       :{0}\nthe one from sector names:{1}".format(spins.keys(),sv.keys()))
    spins_bin=[1 if spins[x]=="odd" else 0 for x in spins]
    return sumrule(coeffs,fhtypes,sector_names,context,sector_wise_option=spins_bin) 

if __name__=='__main__':
    p=cb.context_for_scalar(Lambda=15,epsilon=0.5)
    gs=cached_cb(p)
    gs.compute_table(3,13)
    rule=single_correlator_sumrule_parse({"S":[F,H,0],"T":[-F,H,F],"A":[-F,-H,F]},{"S":"even","T":"even","A":"odd"},p)

    file_path="test.xml"
    spin=0
    sector="S"
    delta=0.518
    dtry=2.5
    pms=rule.entire_continuous(p.field(delta),gs,[[sector,spin,dtry]],scalar_default=p.field(0.5))+[\
            rule.vector_at_sector("S",p.field(delta),p.gBlock(0,1.5,0,0))]
    norm=rule.vector_at_sector("S",p.field(delta),p.identity_vector())
    obj=norm*0
    problem=p.SDP(norm,obj,pms)
    problem.write(file_path)
