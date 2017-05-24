#!/usr/bin/python
################################################################################
# Simplified Test Case for TP functions
################################################################################
# Units     P       T      rho
#===============================================================================
# Input     bar     K      g/cc   
# Output    atm     K      g/cc    #kg/m**3 
##############################################################################80
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
from scipy import optimize as op
import utils as ut
###############################################################################


def SRK(R,P,T,x,Pc,Tc,omega,kij):
    
    # Psuedoreduced values
    Tr = T/Tc
    Pr = P/Pc

    # SRK Parameters
    u  = 1.
    w  = 0.
    b  = 0.08664 * R * Tc / Pc
    fw = 0.48 + 1.574*omega - 0.176*omega**2
    alpha = (1.0 + fw * (1 - np.sqrt(Tr)))**2
    a  = (0.42748 *  (R * Tc)**2  / Pc ) * alpha


    # Incorporate Mixture 
    if np.size(Pc) > 1:
        am = np.sum( (x.reshape([np.size(x),1]) * x)  *          \
                      np.sqrt(a.reshape([np.size(a),1]) * a) *   \
                      (1 - kij) )
        bm = np.sum(x * b)
#        print "am, ",am
#        print "bm, ",bm
#        print "P, ",P
#        print "R, ",R
#        print "T, ",T

    # Set Coefficients
    
    if np.size(Pc) > 1:
        AC = am * P / R**2 / T**2
        BC = bm * P / R / T
    elif np.size(Pc) == 1:
        AC = a * P / R**2 / T**2
        BC = b * P / R / T  
#    print "AC, ",AC
#    print "BC, ",BC
    coeff = np.zeros(4)
        # build cubic polynomial
    def g(z):
        """
        Cubic polynomial in z from EOS. This should be zero.
        :param z: float compressibility factor
        """
        return z**3 - (1 + BC - u*BC)  * z**2 +           \
                      (AC + w*BC**2 - u*BC - u*BC**2) * z    -           \
                      AC*BC   - w*BC**2 - w*BC**3
                
                
    coeff[0] =  1.0
    coeff[1] = -1.0         - BC      + u*BC
    coeff[2] =  AC          + w*BC**2 - u*BC - u*BC**2
    coeff[3] = -1.0*AC*BC   - w*BC**2 - w*BC**3
    z = np.real( np.roots(coeff)[np.flatnonzero(np.iscomplex(np.roots(coeff))==False)])

#    alfa = (3*coeff[2] - coeff[1]**2)/3
#    beta = (2*(coeff[1]**3) - 9*coeff[1]*coeff[2] + 27*coeff[3])/27
#    rtest = ( (beta**2) /4 ) + ( (alfa**3) /27 )
#    
#    z = np.zeros(3)
#    if rtest < 0:
#        costheta = (-1.*beta) / (2 * np.sqrt( (-1.*alfa**3)/27))
#        theta = np.arccos(costheta)
#        
#        z[0] = 2 * np.sqrt( -1.0*alfa / 3) * np.cos(theta/3)
#        z[1] = 2 * np.sqrt( -1.0*alfa / 3) * np.cos(theta/3 + 2*np.pi / 3)
#        z[2] = 2 * np.sqrt( -1.0*alfa / 3) * np.cos(theta/3 + 4*np.pi / 3)

#        z = z - coeff[1]/3
#        
#    else:
#        
#        A_an = ( -1.*(beta/2) + np.sqrt(beta**2/4 + alfa**3/27) )**(1/3)        
#        B_an = ( -1.*(beta/2) - np.sqrt(beta**2/4 + alfa**3/27) )**(1/3)        
#        
#        z[0] = A_an + B_an
#        z[1] = -1.0*(A_an + B_an)/2 + ((A_an-B_an)/2)*np.sqrt(-3)
#        z[2] = -1.0*(A_an + B_an)/2 - ((A_an-B_an)/2)*np.sqrt(-3)        
#        z = z - coeff[1]/3
#    z = z[~np.isnan(z)]
#    z = np.unique(z)        
    #cubic = lambda Z:   Z**3 -                                                    \
    #                    Z**2 * (1.0 + BC - u*BC) +                                \
    #                    Z    * (AC + w*BC**2 - u*BC - u*BC**2) -                  \
    #                            AC*BC - w*BC**2 - w*BC**3
    #z = np.roots(coeff).real.clip(min=0) #op.newton(g, 1.0)  # compressibility factor
    #print z
    # V = maximum root,
    # l = minimum root
    # EOS is run per phase, thus not specific here
    if np.size(Pc) > 1:     
        lnfug_v = np.nan_to_num(lnfugcoef(np.max(z),x,Tc,Pc,u,w,a,b,am,bm,AC,BC,kij))
        lnfug_l = np.nan_to_num(lnfugcoef(np.min(z),x,Tc,Pc,u,w,a,b,am,bm,AC,BC,kij))
        dG = np.sum(x*lnfug_l) - np.sum(x*lnfug_v)
        if dG > 0:
            Zn = np.max(z)
            lnfug = lnfug_v
        elif dG < 0:
            Zn = np.min(z)
            lnfug = lnfug_l
        elif dG == 0:
            Zn = z
            lnfug = lnfug_v    
    elif np.size(Pc) == 1:
        lnfug = np.exp(z - 1 - np.log(z - BC) - AC / np.sqrt(8) / BC * np.log(
                (z + (1 + np.sqrt(2)) * BC) / (z + (1 - np.sqrt(2)) * BC)))
        z.clip(min=0)
        Zn = np.zeros(2)
        # Liquid
        Zn[0] = z.min()
        lnfug = lnfug_l
        # Vapor
        Zn[1] = z.max()
    rho = P / (R * T * Zn)  # density

    return Zn,rho,lnfug

def lnfugcoef(Zc,x,Tc,Pc,u,w,a,b,am,bm,AC,BC,kij):
    # Fugacity Coefficient
    
    CC      = np.sqrt(u**2 - 4.*w)

    gamma_i = ( (2 * np.sqrt(a)) / am) * np.sum(x * np.sqrt(a) * (1-kij))
    able    = (b/bm) * (Zc - 1.0)
    baker   = -1.0*np.log(Zc - BC)
    charlie = (AC/(BC*CC)) * (b/bm - gamma_i)
    
    lnfug = able + baker + charlie
    return lnfug        


def dlfdp(R,P,T,x,Pc,Tc,omega,kij,eptol):
    Zc_p,rho_p,lf_p = SRK(R,P+eptol,T,x,Pc,Tc,omega,kij)
    Zc_m,rho_m,lf_m = SRK(R,P-eptol,T,x,Pc,Tc,omega,kij)
    dlfdp = (lf_p - lf_m) / (2*eptol)
    return dlfdp
    
def WilsonEqn(Pc,Tc,omega,T,P):
    K = (Pc/P) * np.exp(5.377 * (1. + omega)*(1. - Tc/T))
    return K

def P_Bubble_Ideal(Pc,Tc,omega,T,z):
    A = 5.377*(1 + omega) * (1 - Tc / T)
    B = Pc * np.exp(A)
    P = np.sum(B*z)
    K = WilsonEqn(Pc,Tc,omega,T,P)    
    y = K*z    
    return P,y


def P_Bubble_Calc(R,T,z,Pc,Tc,omega,kij,maxiter):
    P,y = P_Bubble_Ideal(Pc,Tc,omega,T,z)
 #   print P
    f = 10
    df = 10
    eptol = 1e-6
    for i in np.arange(0,maxiter):
        Zc_y,rho_y,lf_y = SRK(R,P,T,y,Pc,Tc,omega,kij)
        Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
        K = np.exp(lf_z - lf_y)
        f = np.sum(z*K) - 1
        
        # Bubble Point / Psat
        dlf_z = dlfdp(R,P,T,z,Pc,Tc,omega,kij,eptol)
        dlf_y = dlfdp(R,P,T,y,Pc,Tc,omega,kij,eptol)
        df = np.sum(np.sum(z*K) * (dlf_z - dlf_y))
        P = P + f/df
        y = K*z
#        print P, f, df
        if f**2 < 1e-12:
            print("solution converged")
            print P, f, df            
            break
	
    return P,y    

def RR_Vapor(z,VF,K):
    f_beta = np.sum((z * (K - 1) ) / (VF*(K-1) + 1 ))
    return f_beta
    
def RRPrime_Vapor(z,VF,K):
    fprime_beta = np.sum((z * (K - 1)**2 ) / (VF*(K-1) + 1 )**2)
    return -1.0*fprime_beta

def RR_Vapor_Calc(z,VF,K,eptol):
    eps = 1
    while True: 
        f = RR_Vapor(z,VF,K) 
        df = RRPrime_Vapor(z,VF,K)
        VF = VF - f/df
        # Stabilize
        print VF,f      
        if f**2 < eptol:
            break            
        if VF > 1:
            break
            VF = (0.75+1)/2 - f/df
        elif VF < 0:
            VF = 0.25/2 - f/df                    

    x = z / (1 + VF*(K - 1))
    y = K*x        
    print VF,f            
    return VF,x,y

def StabilityVapor(R,P,T,z,Pc,Tc,omega,kij):
    # 1) Calculate fugacity coefficients for overall composition z and calculate di
    Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
    d = np.log(z) + lf_z
    fug_z = np.exp(lf_z)*z*P
    # 2) Generate initial estimates for vapor like and liquid like trial phases
    K = WilsonEqn(Pc,Tc,omega,T,P)

    # 3) Run Stability Test
    while True:

        W = z*K
        y = W/np.sum(W)

        #print "K, ", K
        #print "y, ", y
        # First test: Vapor like trial phase, components: y
        Zc,rho,lf = SRK(R,P,T,y,Pc,Tc,omega,kij)    
        fug = np.exp(lf)*y*P
        Res = (fug_z/fug) * (1/np.sum(W))
        K = K*Res

        # Update W
#        W = np.exp(d - lf)
 #       y = W / np.sum(W)
        # Reduced TPD
        tpd = np.sum(y*(np.log(y) + lf - d))
        # TPD 
        TPD = tpd * R * T
        # Modified Tangent Plane Function
        tm = 1 + np.sum(W*(np.log(W) + lf - d - 1) )
        # Check for convergence
  #      r = np.sum(np.log(W) + lf - d)**2
        r = np.sum(Res-1)**2
        trivial = np.sum(np.log(K))**2
        #print r, trivial
        if r < 1e-10:
            break
        if trivial < 1e-4:
            print "trivial solution converged"
            break
    
    return W,y,TPD,tpd,tm            

def StabilityLiquid(R,P,T,z,Pc,Tc,omega,kij):
    # 1) Calculate fugacity coefficients for overall composition z and calculate di
    Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
    fug_z = np.exp(lf_z)*y*P
    d = np.log(z) + lf_z
    #print "d_i: ",d
    # 2) Generate initial estimates for vapor like and liquid like trial phases
    K = WilsonEqn(Pc,Tc,omega,T,P)

    # 3) Run Stability Test
    while True:
        W = z/K    
        # Normalize
        x = W/np.sum(W)
        Zc,rho,lf = SRK(R,P,T,x,Pc,Tc,omega,kij)    
        fug = np.exp(lf)*x*P
        # Update W
        Res = (fug/fug_z) * np.sum(W)
        K = K*Res
      # Reduced TPD
        tpd = np.sum(x*(np.log(x) + lf - d))
        # TPD 
        TPD = tpd * R * T
        # Modified Tangent Plane Function
        tm = 1 + np.sum(W*(np.log(W) + lf - d - 1) )
        r = np.sum(Res - 1)**2
        trivial = np.sum(K)**2
        #print "r: ",r, trivial
        if r < 1e-10:
            print "test phase converged"
            break
        if trivial < 1e-4:
            print "trivial solution converged"
            break           
    
    return W,x,TPD,tpd,tm            



# Test Parameters

compstr = ["CH4","nC4","nC10"]
omega = np.array([0.008,0.199,0.489])
Pc = np.array([454.,380.,212.]) * 1E4 * 9.8692327E-06
Tc = np.array([190.6,425.2,617.7])
Mw = np.array([16.04,58.12,142.29])
z = np.array([0.35,0.45,0.20])      # Fraction of oil input stream
kij = np.array([ [0.0000, 0.0133, 0.0422],
                     [0.0133, 0.0000, 0.0078],
                     [0.0422, 0.0078, 0.0000]  ])
# Ambient Conditions
n_mols = 1.0                        # # mols entering system
feed_rate = 100.                    # kg/hr
P         = 50                      # atm
T         = 80 + 273.15              # K
R  = 8.2057338*(10**-5)	    # m3 atm / K *mol

VF = 0.9
K = WilsonEqn(Pc,Tc,omega,T,P)
x = z / (1 + VF*(K - 1))
y = K*x
numiter = 5
numitera = 1
numiterb = 1



# Accelerated SS/DS

VF,x,y = RR_Vapor_Calc(z,VF,K,1e-10)
v = y*z

for i in np.arange(0,numiter):
    Zc_z,rho_x,lf_x = SRK(R,P,T,x,Pc,Tc,omega,kij)
    Zc_y,rho_y,lf_y = SRK(R,P,T,y,Pc,Tc,omega,kij)
    Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)           
    
                
    g = np.log(np.exp(lf_y)*y*P) - np.log(np.exp(lf_x)*x*P) 
    J = np.eye(Pc.size) * (z / x*y - 1) * (1/(VF * (1-VF)))
    dv = np.linalg.solve(J,-1.*g)
    v = v + dv
    y = v/np.sum(v)
    VF = np.sum(v)
    print v




## RR Solution
 

## Check if solution exists
#RRZero = RR_Vapor(z,0.0,K)
#RROne  = RR_Vapor(z,1.0,K)

#if np.sum(RRZero) > 0 and np.sum(RROne) < 0:
#    print("Solution Exists")
#    asym = 1/(1-K)
#    W_v,v,TPD_v,tpd_v,tm_v = StabilityVapor(R,P,T,z,Pc,Tc,omega,kij)
#    W_l,l,TPD_l,tpd_l,tm_l = StabilityLiquid(R,P,T,z,Pc,Tc,omega,kij)
#    l = x
#    y = v
#    for i in np.arange(0,numiter):
#        for j in np.arange(0,numitera):
#            f  = RR_Vapor(z,VF,K)
#            df = RRPrime_Vapor(z,VF,K)
#            VF = VF - f/df
#            print "Iteration ",i, "VF: ,",VF, "err: ",f            
#            # Stabilize Result
#            if VF > 1:
#                VF = (0.75 + 1)/2 - f/df
#            elif VF < 0:
#                VF = 0.25/2 - f/df
#            
#            if f**2 < 1e-12:
#                print "Converged"
#                break                
#        for j in np.arange(0,numiterb): 
#            # Update y,x
#            x = z / (1 + VF*(K - 1) )
#            y = K*z / (1 + VF*(K-1) )
#               
#            
#            # Solve EOS        
#            Zc_z,rho_x,lf_x = SRK(R,P,T,x,Pc,Tc,omega,kij)
#            Zc_y,rho_y,lf_y = SRK(R,P,T,y,Pc,Tc,omega,kij)
#            Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)           
#        
#            # Generate Fugacities
#            fug_x = np.exp(lf_x)*x*P
#            fug_y = np.exp(lf_y)*y*P
#            fug_z = np.exp(lf_z)*z*P

##            K = np.exp(lf_x)/np.exp(lf_y)
#            # Check Error
#            err = np.sum(fug_x/fug_y - 1)**2
#            print "Iteration ",i, "K: ,",np.sum(K), "err: ",err
#            if err < 1e-12:
#                print("Converged")
#                break
#            K = K * (fug_x/fug_y)         

#    
#elif np.sum(RRZero) < 0:
#    print("Subcooled Liquid")
#elif np.sum(RROne) > 0:
#    print("Superheated Vapor")
#    






#################################################################################
## Stability Testing
##-------------------------------------------------------------------------------

#W_v,v,TPD_v,tpd_v,tm_v = StabilityVapor(R,P,T,z,Pc,Tc,omega,kij)
#W_l,l,TPD_l,tpd_l,tm_l = StabilityLiquid(R,P,T,z,Pc,Tc,omega,kij)
##Sv = np.sum(W_v)
##Sl = np.sum(W_l)

##y = v
##x = l

#P,y = P_Bubble_Ideal(Pc,Tc,omega,T,z)
#VF = 0.0
#y = K*z / (1 + VF*(K - 1) )
#x =   z / (1 + VF*(K - 1) )

##   print P
#f = 10
#df = 10
#eptol = 1e-6
#maxiter = 100
#for i in np.arange(0,maxiter):
#    Zc_z,rho_x,lf_x = SRK(R,P,T,x,Pc,Tc,omega,kij)
#    Zc_y,rho_y,lf_y = SRK(R,P,T,y,Pc,Tc,omega,kij)
#    Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)    

#    dlf_x = dlfdp(R,P,T,x,Pc,Tc,omega,kij,eptol)
#    dlf_y = dlfdp(R,P,T,y,Pc,Tc,omega,kij,eptol)
#    dlf_z = dlfdp(R,P,T,z,Pc,Tc,omega,kij,eptol)

#    fug_x = lf_x*x*P
#    fug_y = lf_y*y*P
#    fug_z = lf_z*P        
#    err = np.sum(fug_x/fug_y -1)

#    K = np.exp(lf_x - lf_y)    
#        
#    f =  np.sum( z*(K-1) / (1 + VF*(K-1)) )
#    df = np.sum(z * K * (dlf_x - dlf_y))

#    P = P - f/df
#    y = K*z / (1 + VF*(K - 1) )
#    x =   z / (1 + VF*(K - 1) )




#    print P, f, df, err
#    if f**2 < 1e-12:
#        print("solution converged")
#        print P, f, df            
#        break
#	


