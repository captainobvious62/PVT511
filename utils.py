#!/usr/bin/python
###################
# Function Module
#
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

##############################################################################80
# Functions for Cn+ Characterization
##############################################################################80


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

    # Set Coefficients
    
    if np.size(Pc) > 1:
        AC = am * P / R**2 / T**2
        BC = bm * P / R / T
    elif np.size(Pc) == 1:
        AC = a * P / R**2 / T**2
        BC = b * P / R / T  

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
    cubic = lambda Z:   Z**3 -                                                    \
                        Z**2 * (1.0 + BC - u*BC) +                                \
                        Z    * (AC + w*BC**2 - u*BC - u*BC**2) -                  \
                                AC*BC - w*BC**2 - w*BC**3
    z = np.roots(coeff).real #op.newton(g, 1.0)  # compressibility factor

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


def Ped_Tc(coeff,rho,MW):
    Tc = coeff[0]*rho + coeff[1]*np.log(MW) + coeff[2]*MW + coeff[3]/MW
    return Tc

def Ped_Pc(coeff,rho,MW):
    Pc = coeff[4] + coeff[5]*(rho**coeff[8]) + coeff[6]/MW + coeff[7]/(MW**2)
    return np.exp(Pc)
    
def Ped_omega(coeff,rho,MW):
    m  = coeff[9] + coeff[10]*MW + coeff[11]*rho + coeff[12]*(MW**2)
    root_coeff = np.zeros([m.size,2])
    root_coeff[:,0] = 1.574 / -0.176
    root_coeff[:,1] = (0.480 / (-0.176) )
    root_coeff[:,1] = root_coeff[:,1] - m / (-0.176) 
    omega = np.zeros([m.size])
    for i in range(0,m.size):
        omega[i] = np.roots(root_coeff[i,:])
    return omega

def Ped_z(A,B,Ci):
    z = A + B*Ci
    return np.exp(z)    

def Ped_MW(Ci):
    MW = 14.*Ci - 4.
    return MW

def Ped_rho(C,D,Ci):
    rho = C + D*np.log(Ci)
    return rho

def Ped_R1(z_cnplus,A,B,SCN):
    R1 = z_cnplus - np.sum(Ped_z(A,B,SCN))
    return R1

def Ped_R2(z_cnplus,MW_cnplus,A,B,SCN):
    R2 = z_cnplus*MW_cnplus - np.sum(Ped_z(A,B,SCN)*Ped_MW(SCN))
    return R2    

def Ped_R3(rho_cnplus,A,B,C,D,SCN):
    R3 = rho_cnplus - np.sum(Ped_z(A,B,SCN)*Ped_MW(SCN)) /                     \
                    np.sum( (Ped_z(A,B,SCN)*Ped_MW(SCN))/Ped_rho(C,D,SCN))
    return R3

def Ped_R4(rho_i,C,D,Ci):
    R4 = rho_i - Ped_rho(C,D,Ci)
    return R4

def PS_Tc(z,MW,Tc):
    PS_Tc = np.sum(z*MW*Tc) / np.sum(z*MW)
    return PS_Tc

def PS_Pc(z,MW,Pc):
    PS_Pc = np.sum(z*MW*Pc) / np.sum(z*MW)
    return PS_Pc    

def PS_omega(z,MW,omega):
    PS_omega = np.sum(z*MW*omega) / np.sum(z*MW)
    return PS_omega    
    
def PS_MW(z,MW):
    PS_MW = np.sum(z*MW) / np.sum(z)
    return PS_MW    


    
def fug(coeff,y,P):
    fug = coeff*y*P
    return fug    
    
def WilsonEqn(Pc,Tc,omega,T,P):
    K = (Pc/P) * np.exp(5.373 * (1. + omega)*(1. - Tc/T))
    return K

def P_Bubble_Ideal(Pc,Tc,omega,T,z):
    A = 5.377*(1 + omega) * (1 - Tc / T)
    B = Pc * np.exp(A)
    P = np.sum(B*z)
    K = WilsonEqn(Pc,Tc,omega,T,P)    
    y = K*z    
    return P,y
    
def RR_Vapor(z,VF,K):
    f_beta = (z * (K - 1) ) / (VF*(K-1) + 1 )
    return f_beta
    
def RRPrime_Vapor(z,VF,K):
    fprime_beta = (z * (K - 1)**2 ) / (VF*(K-1) + 1 )**2
    return -1.0*fprime_beta

def RR_Vapor_Calc(z,VF0,K,eptol):
    eps = 1
    VFn = VF0
    while True: 
        VF = VFn - np.sum(RR_Vapor(z,VF0,K)) / np.sum(RRPrime_Vapor(z,VF0,K) )
        eps  = np.abs(VF - VF0)
        # Stabilize
        VFn = VF
        if np.abs(VF - VF0) < eptol:
            break            
        if VF > 1:
            VF = ((0.75+1)/2 + VF) / 2
        elif VF < 0:
            VF = 0.25/2    
        VF0 = VF            
        print eps,VFn
    x = z / (1 + VF*(K - 1))
    y = K*x        
    return VF,x,y

def RR_Liquid(z,VL,K):
    f_alfa = (z * (1 - K) ) / (K + (1 - K)*VL)
    return np.sum(f_alfa)

def RR_Liquid_Calc(z,LF0,K,eptol):
    eps = 1
    while eps > eptol:        
        LF = LF0 - RR_Liquid(z,LF0,K) /                                        \
            ( (RR_Liquid(z,LF0+eptol,K) - RR_Liquid(z,LF0-eptol,K)) / 2*eptol ) 
        eps  = np.abs(LF - LF0)
        LF0 = LF
    return LF                

def vaporStability(R,P,T,z,Pc,Tc,omega,kij):

    # Calculate mixture fugacity
    Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)   
    # Derive fugacity
    fug_z = P * z * np.exp(lf_z)
    # Vapor like phase
    K_y = WilsonEqn(Pc,Tc,omega,T,P)
    while True:    
        Y = z*K_y

        S_v = np.sum(Y)
        y_i = Y / S_v

        # Calculate mixture fugacity
        Zc_y,rho_y,lf_y = SRK(R,P,T,y_i,Pc,Tc,omega,kij)   
        # Derive fugacity
        fug_y = P * y_i * np.exp(lf_y)


        R_y = (fug_z/fug_y) * (1 / S_v)
        K_y = K_y * R_y

        # Check for trivial solution/convergence

        eps_c = np.sum(R_y - 1)**2
        eps_t = np.sum(np.log(K_y)**2)
        print eps_c, eps_t
        if eps_c < 1E-10:
            print 'Converged!'
            y_trivial = False
            break

        if eps_t < 1E-4:
            print 'Trivial Solution'
            y_trivial = True
            break        
    d = np.log(z) + lf_z
    tpd = np.sum(y_i*(np.log(y_i) + lf_y - d))
#    tm = 1 + np.sum( np.log(K_y + lf_y)
    return S_v,y_trivial,tpd,K_y

def liquidStability(R,P,T,z,Pc,Tc,omega,kij):

    # Calculate mixture fugacity
    Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)   
    # Derive fugacity
    fug_z = P * z * np.exp(lf_z)
    # Vapor like phase
    K_x = WilsonEqn(Pc,Tc,omega,T,P)
    while True:    
        X = z/K_x

        S_l = np.sum(X)
        x_i = X / S_l

        # Calculate mixture fugacity
        Zc_x,rho_x,lf_x = SRK(R,P,T,x_i,Pc,Tc,omega,kij)   
        # Derive fugacity
        fug_x = P * x_i * np.exp(lf_x)

    
        R_x = (fug_x/fug_z) * S_l
        K_x = K_x * R_x

        # Check for trivial solution/convergence

        eps_c = np.sum(R_x - 1)**2
        eps_t = np.sum(np.log(K_x)**2)
        print eps_c, eps_t
        if eps_c < 1E-10:
            print 'Converged!'
            x_trivial = False
            break

        if eps_t < 1E-4:
            print 'Trivial Solution'
            x_trivial = True
            break        

    d = np.log(z) + lf_z
    tpd = np.sum(x_i*(np.log(x_i) + lf_x - d))
    return S_l,x_trivial,tpd,K_x

def StabilityTest(R,P,T,z,Pc,Tc,omega,kij):
    S_l, x_trivial,tpd_l,K_x = liquidStability(R,P,T,z,Pc,Tc,omega,kij)
    S_v, y_trivial,tpd_v,K_y = vaporStability(R,P,T,z,Pc,Tc,omega,kij)        
    if (S_l < 1 and S_v < 1):
        singlephase = True
    elif(x_trivial == True and y_trivial == True):
        singlephase = True        
    elif(x_trivial == True and S_v <= 1.0):
        singlephase = True        
    elif(S_l <= 1.0 and y_trivial == True):
        singlephase = True                
        
    if (S_l > 1.0 or S_v > 1.0):
        singlephase = False
    return singlephase,tpd_l,tpd_v,K_x,K_y      
    


def P_Bubble_Calc(R,T,VF,z,Pc,Tc,omega,kij,eptol):
    P,y = P_Bubble_Ideal(Pc,Tc,omega,T,z)
    f = 10
    while f**2 > eptol**2:
        Zc_y,rho_y,lf_y = SRK(R,P,T,y,Pc,Tc,omega,kij)
        Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
        K = np.exp(lf_z - lf_y)
        # Bubble Point / Psat
        dlf_z = dlfdp(R,P,T,z,Pc,Tc,omega,kij,eptol)
        dlf_y = dlfdp(R,P,T,y,Pc,Tc,omega,kij,eptol)
        f = np.sum( z*(K - 1) / (1 + VF*(K - 1)) ) 
        df = np.sum(z*K * (dlf_z - dlf_y))
        P = P - f/df
        y = K*z
        print f, P
    return P,y    

def RR_F(z,K,VF):
    f  = np.sum( (z*(1-K)) / (1 - VF*(1 - K) ) )
    return f

def StabilityVapor(R,P,T,z,Pc,Tc,omega,kij):

    # Calculate mixture fugacity
    Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
    # Derive fugacity
    fug_z = P * z * np.exp(lf_z)
    # Vapor like phase
    K_y = WilsonEqn(Pc,Tc,omega,T,P)
    while True:    
        Y = z*K_y

        S_v = np.sum(Y)
        y_i = Y / S_v

        # Calculate mixture fugacity
        Zc_y,rho_y,lf_y = SRK(R,P,T,y_i,Pc,Tc,omega,kij)
        # Derive fugacity
        fug_y = P * y_i * np.exp(lf_y)


        R_y = (fug_z/fug_y) * (1 / S_v)
        K_y = K_y * R_y

        # Check for trivial solution/convergence

        eps_c = np.sum(R_y - 1)**2
        eps_t = np.sum(np.log(K_y)**2)
        print eps_c, eps_t
        if eps_c < 1E-10:
            print 'Converged!'
            y_trivial = False
            break

        if eps_t < 1E-4:
            print 'Trivial Solution'
            y_trivial = True
            break        
    d = np.log(z) + lf_z
    tpd = np.sum(y_i*(np.log(y_i) + lf_y - d))
#    tm = 1 + np.sum( np.log(K_y + lf_y)
    return S_v,y_trivial,tpd,K_y,Y

def dlfdp(R,P,T,x,Pc,Tc,omega,kij,eptol):
    Zc_x_p,rho_x_p,lf_x_p = SRK(R,P+eptol,T,x,Pc,Tc,omega,kij)
    Zc_x_m,rho_x_m,lf_x_m = SRK(R,P-eptol,T,x,Pc,Tc,omega,kij)
    dlfdp_x = (lf_x_p - lf_x_m) / 2*eptol
    return dlfdp_x
    
def PsatIdeal(z,T,Pc,Tc,omega,eptol):
    A = 5.377*(1 + omega)*(1 - Tc/T)
    B = Pc*np.exp(A)
    P =  np.sum(B)
    K = WilsonEqn(Pc,Tc,omega,T,P)    
    y = K*z
    return P,y

def Psat(R,P,T,z,y,Pc,Tc,omega,kij,eptol):
    f = 10
    while f**2 > eptol**2:
        Zc_y,rho_y,lf_y = SRK(R,P,T,y,Pc,Tc,omega,kij)
        Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
        K = np.exp(lf_z - lf_y)
        # Bubble Point / Psat
        dlf_z = dlfdp(P,T,z,Pc,Tc,omega,kij,eptol)
        dlf_y = dlfdp(P,T,y,Pc,Tc,omega,kij,eptol)
        f = np.sum(z*K - 1)
        df = np.sum(z*K * (dlf_z - dlf_y))
        P = P - f/df
        y = K*z
    return P
    
    
    
    
    












    
       
       
       
       
       
