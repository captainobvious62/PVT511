#!/usr/bin/python

##############################################################################80
# [Second Attempt at] PTE 511 Final Project
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
from utils import *
from math import *
from matplotlib.pyplot import cm
np.set_printoptions(suppress=True)
from scipy.optimize import fsolve, root, brentq, newton, bisect
# Verification Module
import thermo as th
import CoolProp.CoolProp as CP
import json
from CoolProp.Plots import PropertyPlot
##############################################################################80
# Plot Parameters
##############################################################################80

def SRK_DF(R,P,T,F8Comp):
    x = F8Comp.z.as_matrix()
    Pc = F8Comp.Pc.as_matrix()
    Tc = F8Comp.Tc.as_matrix()    
    omega = F8Comp.omega.as_matrix()    
    kij = F8Comp.kij        
    dG,Zc,rho,lnfug = SRK(R,P,T,x,Pc,Tc,omega,kij)
    return dG,Zc,rho,lnfug

def SRK(R,P,T,x,Pc,Tc,omega,kij):
    """
    Returns Gibbs Energy change, two compressibility factors (vapor and liquid),
    two density estimates, and two ln(phi) values.
    
    For the variables that return two values, the one determined by the minimum
    Gibbs energy change is placed first.    
    
    """
    # Psuedoreduced values
    Tr = T/Tc
    Pr = P/Pc

    # SRK Parameters
    u  = 1.
    w  = 0.
    oa1 = 0.42748
    oa2 = 0.08664
    fw = 0.480 + 1.574*omega - 0.176*omega**2
    alpha = (1.0 + fw * (1 - np.sqrt(Tr)))**2
    a  = (0.42748 *  (R * Tc)**2  / Pc ) * alpha
    b  = 0.08664 * R * Tc / Pc

    A  = oa1 * (Pr/Tr**2) * alpha #a * ( P / ( R * T )**2 )
    B  = oa2 * (Pr/Tr)  #b * ( P / ( R * T ) ) 

    # Incorporate Mixture 
    if np.size(Pc) > 1:
        AIJ = np.sqrt( A.reshape( [A.size,1] ) * A) * (1 - kij)
        AM  = np.sum( ( x.reshape([x.size,1]) * x ) *               \
              np.sqrt( A.reshape( [A.size,1] ) * A) *               \
              (1.0 - kij) )
        BM  = np.sum(x * B)
    elif np.size(Pc) == 1:
        AM = A
        BM = B
    coeff = np.zeros(4)
             
    coeff[0] =  1.0
    coeff[1] = -1.0         - BM      + u*BM
    coeff[2] =  AM          + w*BM**2 - u*BM - u*BM**2
    coeff[3] = -1.0*AM*BM   - w*BM**2 - w*BM**3

    z = np.real( np.roots(coeff)[np.flatnonzero(np.iscomplex(np.roots(coeff))==False)])
    z = np.array([z.min(),z.max()])
  
    if np.size(Pc) > 1:
        lf =  B.reshape([B.size,1])/BM * (z-1) -                              \
              np.log(z-BM) +                                                  \
              AM/BM*( B.reshape([B.size,1])/BM -                              \
              2./AM *                                                         \
              np.sum(x*AIJ,axis=0).reshape([x.size,1]) ) * np.log(1. + BM/z)
        
    else:
        lf = z - 1. - np.log(z - B) - A/B*np.log(1.+B/z)
    gibbs = lf*x.reshape([x.size,1])
    g_star = np.sum(lf*x.reshape([x.size,1]),axis=0)
    index = np.flatnonzero(g_star == g_star.min())
    index = index[0]
    Zc    = z[index]
    lnfug = lf[:,index]
    dG    = g_star[index]        
        
#    print('z: ',z)
#    print('gstar: ',g_star)    
#    print('====================================================')
#    print('----------------------------------------------------')
    SC_SRK = {'N2':-0.0079,'CO2':0.0833,'H2S':0.0466,'C1':0.0234,'C2':0.0605,
              'C3':0.0825,'iC4':0.0830,'nC4':0.0975,'iC5':0.1022,'nC5':0.1209,
              'nC6':0.1467,'nC7':0.1554,'nC8':0.1794,'nC9':0.1868,'nC10':0.2080}
    
    Zra = 0.2906 - 0.0878 * omega
    ci = 0.4077 * ( (R*Tc)/Pc ) * (0.2944 - Zra) 
    lnfug = lnfug - (ci*P)/(R*T)
    rho = P / (R * T * Zc)  # density
    return dG,Zc,rho,lnfug

def RR_Vapor(VF, *data):
    z, K = data
    h = ( z * (K - 1) ) / ( 1 + VF*(K - 1) )
    f_beta = np.sum(h)
    return f_beta
    
def RRPrime_Vapor(VF, *data):
    z, K = data
    h_prime = ( z*(K - 1)**2 ) / ( VF*(K - 1) + 1 )**2
    fprime_beta = np.sum(h_prime)
    return fprime_beta

def RR_Vapor_Calc(R,Pc,Tc,omega,kij,VF,z,K,eptol,numiter):
    f  = 0
    df = 0
    x = 10
    y = 10
    for i in np.arange(0,numiter):
        f = f + RR_Vapor(VF,z,K) 
        df = df - RRPrime_Vapor(VF,z,K)
        VF = VF - f/df
        # Stabilize
        if np.abs(df) < eptol:
            error3 = f
        else:
            error3 = np.abs(f/df)
        if i == 0:
            x = z / (VF * (K-1) + 1)
            y = K*x                                  
        xerror = np.abs(np.sum(x))
        yerror = np.abs(np.sum(y))
        print i,VF,f    
        print K
        print ('x: ',x)
        print ('y: ',y)  
        ################
        if VF > 1:
            VF = 1
            break
        elif VF < 0:
            VF = 0
            break
        if f**2 < eptol:
            break                        
        if error3 < eptol:
            break           
        if np.abs(np.sum(x) - 1) < eptol:
            break
        if np.abs(np.sum(y) - 1) < eptol:
            break            
        x = z / (VF * (K-1) + 1)
        y = K*x                                       
        x = x/np.sum(x)
        y = y/np.sum(y)
        
        # Generate new K Values

        dG_y,Zc_x,rho_x,lf_x = SRK(R, P, T, x, Pc, Tc, omega, kij)
        dG_y,Zc_y,rho_y,lf_y = SRK(R, P, T, y, Pc, Tc, omega, kij)
        # Calculate K from EOS
        K = np.nan_to_num(np.exp(lf_x)/np.exp(lf_y))             
        print('================================================================')
        print('================================================================')
        print(' End of iteration ',i,', Error Below' )
        print('################################################################')
        print('x error: ',xerror,'y error: ',yerror,' abs error: ',np.abs(f),'rel error: ',error3)
        print('-----------------------------------------------------------------------------------')
    if np.sum(x) != 0:
        x = x/np.sum(x)
    if np.sum(y) != 0:
        y = y/np.sum(x)
    if VF == 0:
        x = z.copy()    
    if VF == 1:
        y = z.copy()    
    print VF,f            
    return VF,x,y


def StabilityVapor2(R,P,T,z,Pc,Tc,omega,kij):
    # 1) Calculate fugacity coefficients for overall composition z and calculate di
    dG,Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
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
        dG,Zc,rho,lf = SRK(R,P,T,y,Pc,Tc,omega,kij)    
        fug = np.exp(lf)*y*P
        Res = (fug_z/fug) * (1/np.sum(W))
        K = K*Res

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
    
    return W,y,TPD,tpd,tm,Zc         

def RetardedDewPointPress(P,*data):
    R,T,z,Pc,Tc,omega,kij,maxiter,gdem_interval = data
    #######################################
    # Vapor Like Trial Phase
    #--------------------------------------
    # Initial Condition K
    K_v = WilsonEqn(Pc,Tc,omega,T,P)
    # Call EOS to get fugacity coefficients for liquid feed phase
    dG_z,Zc_z,rho_z,lf_z = SRK(R, P, T, z, Pc, Tc, omega, kij)
    fug_z = np.exp(lf_z)*z*P
    d = np.log(z) + lf_z
    print('Vapor Like Trial Phase')
    for i in np.arange(0,maxiter):
        # Assume feed z is liquid, generate vapor like trial phase w
        Y = z*K_v
        S_v = np.sum(Y)
        y = Y/S_v
        # Call EOS to get fugacity coefficients for vapor trial phase 
        dG_y,Zc_y,rho_y,lf_y = SRK(R, P, T, y, Pc, Tc, omega, kij)
        fug_y = np.exp(lf_y)*y*P
        # Calculate Fugacity Ratio Corrections for successive substitution
        # Account for GDEM Update
        if np.mod(i+1,gdem_interval) == 0:
            Res_m1 = Res_v
            b11 = np.sum(np.log(Res_m1)*np.log(Res_m1))
        
            # K Value Updates
            Res_v = (fug_z/fug_y)*(1./S_v)
            b01 = np.sum(np.log(Res_v)*np.log(Res_m1))
            lmda_v = np.abs( b11 / (b11 - b01) )
        else:
            Res_v = (fug_z/fug_y)*(1./S_v)

        # K Value Updates

        tpdv = np.sum(y*(np.log(y) + lf_y - d))
        tmv = 1 + np.sum(np.log(K_v) + lf_y)
        # Check for convergence
        trivial_v = np.sum(np.log(K_v))**2
        err_v = np.sum(Res_v - 1)**2
        if err_v < 1e-12:
            print('Vapor Like Test Phase Converged')
            break
        elif trivial_v < 1e-4:
            v_trivial = True
            print('Vapor Like Test Phase Arrived at Trivial Solution')
            break            
        else:
            if np.mod(i+1,gdem_interval) == 0:
                print('GDEM Update')
                K_v = K_v * (Res_v**lmda_v)
            else:                 
                K_v = K_v * Res_v 
        print 'err: ',err_v,'triv: ',trivial_v,'tpd: ',tpdv,'tm: ',tmv
        print S_v
        print fug_y
        print y
    print("=============================================================")        

    Sopt = S_v - 1
    return Sopt

def RetardedBubblePointPress(P, *data):
    R,T,z,Pc,Tc,omega,kij,maxiter,gdem_interval = data
    #===========================================================================
    # Liquid Like Trial Phase
    # Initial Condition K
    K_l = WilsonEqn(Pc,Tc,omega,T,P)
    print('Liquid Like Trial Phase')
    for i in np.arange(0,maxiter):
        print('Iteration ',i)
        # Assume feed z is liquid, generate liquid like trial phase w
        X = z/K_l
        S_l = np.sum(X)
        x = X/S_l
        # Call EOS to get fugacity coefficients for liquid trial phase 
        dG_x,Zc_x_stb,rho_x,lf = SRK(R, P, T, x, Pc, Tc, omega, kij)
        fug_x = np.exp(lf_x)*x*P
        # Calculate Fugacity Ratio Corrections for successive substitution
        
        # Account for GDEM Update
        if np.mod(i+1,gdem_interval) == 0:
            Res_m1 = Res_l
            b11 = np.sum(np.log(Res_m1)*np.log(Res_m1))
        
            # K Value Updates
            Res_l = (fug_x/fug_z)*(S_l)
            b01 = np.sum(np.log(Res_l)*np.log(Res_m1))
            lmda_l = np.abs( b11 / (b11 - b01) )
        else:
            Res_l = (fug_x/fug_z)*(S_l)            
        tpdl = np.sum(x*(np.log(x) + lf_x - d))
        tml = 1 + np.sum(np.log(K_l) + lf_x)
        # Check for convergence
        trivial_l = np.sum(np.log(K_l))**2
        err_l = np.sum(Res_l - 1)**2
        if err_l < 1e-10:
            print('Liquid Like Test Phase Converged')
            break
        elif trivial_l < 1e-4:
            l_trivial = True
            print('Liquid Like Test Phase Arrived at Trivial Solution')
            break            
        else:
            if np.mod(i+1,gdem_interval) == 0:
                print('GDEM Update')
                K_l = K_l * Res_l**lmda_l
            else:                 
                K_l = K_l * Res_l
        print 'err: ',err_l,'triv: ',trivial_l,'tpd: ',tpdl,'tm: ',tml
        print fug_l
        Sopt = S_l - 1
    return Sopt

def StabilityLiquid2(R,P,T,z,Pc,Tc,omega,kij):
    # 1) Calculate fugacity coefficients for overall composition z and calculate di
    dG,Zc_z,rho_z,lf_z = SRK(R,P,T,z,Pc,Tc,omega,kij)
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
        dG,Zc,rho,lf = SRK(R,P,T,x,Pc,Tc,omega,kij)    
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
        S = np.sum(W)
        #print "r: ",r, trivial
        if r < 1e-10:
            print "test phase converged"
            break
        if trivial < 1e-4:
            print "trivial solution converged"
            break           
    
    return W,x,TPD,tpd,tm,S,Zc


def WilsonEqn(Pc,Tc,omega,T,P):
    Tr = T/Tc
    Pr = P/Pc
    K = np.exp(5.373 * (1.0 + omega) * (1.0 - Tr**(-1.0) ) ) / Pr
    return K

def StabTest(R,P,T,z,Pc,Tc,omega,kij,maxiter,gdem_interval):
    # Perform Stability Analysis
    v_trivial = False
    l_trivial = False
    W_vapor = np.zeros([maxiter+1,z.size])
    #######################################
    # Vapor Like Trial Phase
    #--------------------------------------
    # Initial Condition K
    K_v = WilsonEqn(Pc,Tc,omega,T,P)
    # Call EOS to get fugacity coefficients for liquid feed phase
    dG_z,Zc_z,rho_z,lf_z = SRK(R, P, T, z, Pc, Tc, omega, kij)
    fug_z = np.exp(lf_z)*z*P
    d = np.log(z) + lf_z
    print('Vapor Like Trial Phase')
    for i in np.arange(0,maxiter):
        print('================================================================')
        print('Iteration ',i)
        print('================================================================')
        # Assume feed z is liquid, generate vapor like trial phase w
        Y = z*K_v
        print Y
        S_v = np.sum(Y)
        y = Y/S_v
        if np.isnan(np.min(y)) == True:
            print ('Convergence Failure')
            break
        # Call EOS to get fugacity coefficients for vapor trial phase 
        dG_y,Zc_y,rho_y,lf_y = SRK(R, P, T, y, Pc, Tc, omega, kij)
        fug_y = np.exp(lf_y)*y*P
        # Calculate Fugacity Ratio Corrections for successive substitution
        # Account for GDEM Update
        if np.mod(i+1,gdem_interval) == 0:
            Res_m1 = Res_v
            b11 = np.sum(np.log(Res_m1)*np.log(Res_m1))
        
            # K Value Updates
            Res_v = (fug_z/fug_y)*(1./S_v)
            b01 = np.sum(np.log(Res_v)*np.log(Res_m1))
            lmda_v = np.abs( b11 / (b11 - b01) )
        else:
            Res_v = (fug_z/fug_y)*(1./S_v)

        # K Value Updates

        tpdv = np.sum(y*(np.log(y) + lf_y - d))
        tmv = 1 + np.sum(np.log(K_v) + lf_y)
        # Check for convergence
        trivial_v = np.sum(np.log(K_v))**2
        err_v = np.sum(Res_v - 1)**2
        if err_v < 1e-12:
            print('Vapor Like Test Phase Converged')
            break
        elif trivial_v < 1e-4:
            v_trivial = True
            print('Vapor Like Test Phase Arrived at Trivial Solution')
            break            
        else:
            if np.mod(i+1,gdem_interval) == 0:
                print('GDEM Update')
                K_v = K_v * (Res_v**lmda_v)
            else:                 
                K_v = K_v * Res_v 
        print 'err: ',err_v,'triv: ',trivial_v,'tpd: ',tpdv,'tm: ',tmv
    print("=============================================================")        
    #===========================================================================
    # Liquid Like Trial Phase
    # Initial Condition K
    K_l = WilsonEqn(Pc,Tc,omega,T,P)
    print('Liquid Like Trial Phase')
    for i in np.arange(0,maxiter):
        # Assume feed z is liquid, generate liquid like trial phase w
        X = z/K_l
        S_l = np.sum(X)
        x = X/S_l
        if np.isnan(np.min(x)) == True:
            print ('Convergence Failure')
            break        
        # Call EOS to get fugacity coefficients for liquid trial phase 
        dG_x,Zc_x,rho_x,lf_x = SRK(R, P, T, x, Pc, Tc, omega, kij)
        fug_x = np.exp(lf_x)*x*P
        # Calculate Fugacity Ratio Corrections for successive substitution
        
        # Account for GDEM Update
        if np.mod(i+1,gdem_interval) == 0:
            Res_m1 = Res_l
            b11 = np.sum(np.log(Res_m1)*np.log(Res_m1))
        
            # K Value Updates
            Res_l = (fug_x/fug_z)*(S_l)
            b01 = np.sum(np.log(Res_l)*np.log(Res_m1))
            lmda_l = np.abs( b11 / (b11 - b01) )
        else:
            Res_l = (fug_x/fug_z)*(S_l)            
        tpdl = np.sum(x*(np.log(x) + lf_x - d))
        tml = 1 + np.sum(np.log(K_l) + lf_x)
        # Check for convergence
        trivial_l = np.sum(np.log(K_l))**2
        err_l = np.sum(Res_l - 1)**2
        if err_l < 1e-10:
            print('Liquid Like Test Phase Converged')
            break
        elif trivial_l < 1e-4:
            l_trivial = True
            print('Liquid Like Test Phase Arrived at Trivial Solution')
            break            
        else:
            if np.mod(i+1,gdem_interval) == 0:
                print('GDEM Update')
                K_l = K_l * Res_l**lmda_l
            else:                 
                K_l = K_l * Res_l
        print 'err: ',err_l,'triv: ',trivial_l,'tpd: ',tpdl,'tm: ',tml        
    print('========================================')
    print('Trivial Solutions: ')
    print('========================================')
    print('Vapor: ',v_trivial,' Liquid: ',l_trivial)   
    print('Sv: ',S_v,' Sl: ',S_l)
    return Zc_x,Zc_y,K_l,K_v,S_l,S_v,v_trivial,l_trivial,x,y
        
################################################################################
def FlashCalc(R,T,P,x,y,z,Tc,Pc,omega,kij,maxiter):
    for j in np.arange(0,mainiter):
        print('===================================================================')
        print('Iteration #: ',j)
        print('-------------------------------------------------------------------')
        # Initial Rachford Rice Stuff
    ################################################################################
    # STEP 2) CALCULATE KMIN AND KMAX    
    ################################################################################
        VFmin  = 1.0 / (1.0 - K.max())
        VFmax  = 1.0 / (1.0 - K.min())
        RRZero = RR_Vapor(z,0.0,K)
        RROne  = RR_Vapor(z,1.0,K)
        RRHalf = RR_Vapor(z,0.5,K)
        print('===========================================')
        print('RR(0): ',RRZero, ' RR(1): ',RROne)
        print('===========================================')    
    ################################################################################
    # STEP 3) SOLVE RACHFORD RICE SPLIT EQUATION
    ################################################################################
        if VFmin <  0.0 and VFmax > 0.0:
            print('Rachford Rice Solution Exists')
            data = z, K
            VF = op.fsolve(RR_Vapor,VFmin+VFmax/2,args=data)
    ################################################################################
    # STEP 4) CALCULATE PHASE COMPOSITIONS X AND Y
    ################################################################################        
            x = z / (VF * (K - 1) + 1)
            y = (z * K) / (VF * (K - 1) + 1)  
            print('VF: ',VF)
            print('x: ',x)
            print('y: ',y)   

        else:
            print('Solving for VF')
            VF,x,y = RR_Vapor_Calc(R,Pc,Tc,omega,kij,VF,z,K,eptol,maxiter)
            print('VF: ',VF)
            print('x: ',x)
            print('y: ',y)   
          
        # Calculate Fugacities
        print('=================================================')
        print(' EOS Calc ')
    ################################################################################
    # STEP 4) CALCULATE PHASE Z FACTORS ZL AND ZV FROM THE EOS
    ################################################################################    
        print('x / liquid')
        dG_x,Zc_x,rho_x,lf_x = SRK(R, P, T, x, Pc, Tc, omega, kij)
        print('y / vapor')
        dG_y,Zc_y,rho_y,lf_y = SRK(R, P, T, y, Pc, Tc, omega, kij)
        print('z / combined')
        dG_z,Zc_z,rho_z,lf_z = SRK(R, P, T, z, Pc, Tc, omega, kij)
        fug_l = np.exp(lf_x) * x * P
        fug_v = np.exp(lf_y) * y * P
        fug_z = np.exp(lf_z) * z * P
        #fug_l = lf_x * x * P
        #fug_v = lf_y * y * P
        #fug_z = lf_z * z * P    
        d = np.log(z) + lf_z
        DU[np.mod(j-1,KGDEM_iter),:] = np.log(fug_l/fug_v)

        # Calculate Tangent Plane
        tpdy = np.sum(y*(np.log(y) + lf_y - d))
        tpdx = np.sum(x*(np.log(x) + lf_x - d))

        # Calculate Change in Gibbs Energy
        dG = R*T*(VF*tpdy + (1-VF)*tpdx)

        # Calculate Normalized Gibbs Energy Functions
        g_star_l = np.sum(x*lf_x)
        g_star_v = np.sum(y*lf_y)
        g_star_mix = VF*g_star_v + (1 - VF)*g_star_l

        # Calculate Equal Fugacity Constraint
        efug_err = np.sum(fug_l/fug_v - 1)**2

        if efug_err < 1e-15:
            print('****************************************************************')
            print('Process Converged')
            print('****************************************************************')
            print('====================================================================')
            print('Final Output')
            print('K: ',K)        
            print('VF: ',VF)
            print('x: ',x)
            print('y: ',y)
            print('dG: ',dG)        
            print('tpdx: ',tpdx)
            print('tpdy: ',tpdy)
            print('g*mix: ',g_star_mix)       
            print('--------------------------------------------------------------------')
            print('====================================================================')
            print('P,atm: ',P,'T, K: ',T,' err: ',efug_err,' Trivial err: ',Ktriv_err)   
            print('====================================================================')
            break

        # Check for trivial solution
        Ktriv_err = np.sum(np.log(K))**2
        # Update K Values

        # GDEM Method
        if np.mod(j+1,KGDEM_iter) == 0:
            b01 = np.sum(DU[0-(0+1),:] * DU[0-(1+1)])
            b02 = np.sum(DU[0-(0+1),:] * DU[0-(2+1)])
            b11 = np.sum(DU[0-(1+1),:] * DU[0-(1+1)])    
            b12 = np.sum(DU[0-(1+1),:] * DU[0-(2+1)])
            b22 = np.sum(DU[0-(2+1),:] * DU[0-(2+1)])

            mu1 = (b02*b12 - b01*b22) / (b11*b22 - b12*b12)
            mu2 = (b01*b12 - b02*b11) / (b11*b22 - b12*b12)
            
            KGDEM = (DU[0-(0+1),:] - mu2 * DU[0-(1+1),:]) / (1 + mu1 + mu2)
            K = np.exp(K + KGDEM)
        else:
            K = K * (fug_l / fug_v)
        print('====================================================================')
        print('Iteration Output')
        print('K: ',K)        
        print('VF: ',VF)
        print('x: ',x)
        print('y: ',y)
        print('dG: ',dG)        
        print('tpdx: ',tpdx)
        print('tpdy: ',tpdy)
        print('g*mix: ',g_star_mix)       
        print('--------------------------------------------------------------------')
        print('====================================================================')
        print('P,atm: ',P,'T, K: ',T,' err: ',efug_err,' Trivial err: ',Ktriv_err)   
        print('====================================================================')
        
def VLE(R,T,P,Pc,Tc,z,omega,kij,eptol,maxiter):
    VF = 0.5
    K = WilsonEqn(Pc,Tc,omega,T,P)
    VF,x,y = RR_Vapor_Calc(R,Pc,Tc,omega,kij,VF,z,K,eptol,maxiter)
    return VF,x,y

def RR_NR(z,VF,maxiter,eptol):
    f = 0
    dfdv = 0
    for i in np.arange(0,maxiter):
        f = f + RR_Vapor(VF,z,K) 
        df = df - RRPrime_Vapor(VF,z,K)
        VF = VF - f/df
        if np.abs(f/df) < eptol:
            break
    return VF            

def checkphasesplit(v,z,P,T,Pc,Tc,omega,kij,tol,maxiter):
    eps = 1
    dG_z,Zc_z,rho_z,lf_z =SRK(R,P,T,z,Pc,Tc,omega,kij)
    for i in np.arange(0,maxiter):
        Sv = np.sum(v)
        x = v/Sv
        if np.max(np.abs(x - z)) < tol:
            phasesplit = False
            return phasesplit
        if eps < tol:
            break
        dG_x,Zc_x,rho_x,lf_x =SRK(R,P,T,x,Pc,Tc,omega,kij)
        eps = np.sum(np.log(x/z*np.exp(lf_x)/np.exp(lf_z)))**2
        x = z*np.exp(lf_z)/np.exp(lf_x)
    if i >= maxiter:
        print('Phase Split Check Did Not Converge.')
    if Sv > 1:
        phasesplit = True
    else:
        phasesplit = False
    return phasesplit     

def tpdss(z,P,T,Pc,Tc,omega,kij,tol,maxiter):
    K = WilsonEqn(Pc,Tc,omega,T,P)
    vaporphase = checkphasesplit(K*z,z,P,T,Pc,Tc,omega,kij,tol,maxiter)
    liquidphase = checkphasesplit(K/z,z,P,T,Pc,Tc,omega,kij,tol,maxiter)
    if vaporphase == True and liquidphase == True:
        twophase = True
    else:
        twophase = False
    return twophase                


def PbEstimate(z,T,Pc,Tc,omega):
    Pb = np.sum(z*Pc*np.exp(5.373*(1 + omega)*(1 - Tc/T)))
    return Pb

def PdewEstimate(v,T,Pc,Tc,omega):
    Pdew = np.sum(v / np.exp(5.373*(1 + omega * (1 - Tc/T))))
    return Pdew            
################################################################################
# BLATANTLY STOLEN FUNCTIONS
################################################################################     
# Mathematical special functions
def root_poly3(a1, a2, a3):
    """Roots for a cubic polinomy, x^3 + a1*x^2 + a2*x + a3"""
    Q = (3*a2-a1**2)/9.
    L = (9*a1*a2-27*a3-2*a1**3)/54.
    D = Q**3+L**2
    if D < 0:
        tita = acos(L/(-Q**3)**0.5)
        z1 = 2*(-Q)**0.5*cos(tita/3.+2.*pi/3)-a1/3.
        z2 = 2*(-Q)**0.5*cos(tita/3.+4.*pi/3)-a1/3.
        z3 = 2*(-Q)**0.5*cos(tita/3.)-a1/3.
        z = [z1, z2, z3]
    else:
        S1 = cbrt((L+D**0.5))
        S2 = cbrt((L-D**0.5))
        if D > 0:
            z = [S1+S2-a1/3.]
        else:
            z1 = cbrt(L)+cbrt(L)-a1/3.
            z2 = -cbrt(L)-a1/3.
            z = [z1, z2]
#    print z[0]+z[1]+z[2], -a1
#    print z[0]*z[1]+z[1]*z[2]+z[2]*z[0], a2
#    print z[0]*z[1]*z[2], -a3
    z.sort()
    return z

def ComputeAverageMW(zc,MW):
    AMW = np.sum(z*MW)
    return AMW


def cbrt(x):
    ''' 
    This method calculates the cubic root of a negative value.
    
    '''
    if x >= 0.0:                
        return x ** (1.0/3.0)
    else:
        return -(abs(x) ** (1.0/3.0))
    
def find_correct_root_of_cubic_eos(p0, p1, p2, p3, fluid_type):
    coef_a = (3.0 * p2 - (p1 ** 2)) / 3.0        
    coef_b = (2.0 * (p1 ** 3) - 9.0 * p1 * p2 + 27.0 * p3) / 27.0        
    delta = 0.25 * (coef_b ** 2) + (coef_a ** 3) / 27.0     

    if delta > 0.0:
        # 1 real root, 2 imaginary                 
        const_A =  cbrt(-0.5 * coef_b + sqrt(delta)) 
        const_B =  cbrt(-0.5 * coef_b - sqrt(delta))

        correct_root = const_A + const_B - p1 / 3.0 
    else:
        # 3 real roots
        phi = acos(-0.5 * coef_b / sqrt(-(coef_a ** 3) / 27.0))
        root_1 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0) - p1 / 3.0
        root_2 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0 + 2.0 * np.pi / 3.0) - p1 / 3.0
        root_3 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0 + 4.0 * np.pi / 3.0) - p1 / 3.0

        smallest_root = min(min(root_1,root_2), root_3)
        largest_root = max(max(root_1,root_2), root_3)

        if fluid_type is 'liquid':        
            correct_root = smallest_root
        else:
            assert fluid_type is 'vapor', 'Wrong fluid type! ' + fluid_type
            correct_root = largest_root
    
    assert correct_root > 0.0, fluid_type + ' Z-factor < 0.0! Delta is %f, %f' % (delta, correct_root)
    
    return correct_root

def find_all_roots_of_cubic_eos(p0, p1, p2, p3, fluid_type):
    coef_a = (3.0 * p2 - (p1 ** 2)) / 3.0        
    coef_b = (2.0 * (p1 ** 3) - 9.0 * p1 * p2 + 27.0 * p3) / 27.0        
    delta = 0.25 * (coef_b ** 2) + (coef_a ** 3) / 27.0     

    if delta > 0.0:
        # 1 real root, 2 imaginary                 
        const_A =  cbrt(-0.5 * coef_b + sqrt(delta)) 
        const_B =  cbrt(-0.5 * coef_b - sqrt(delta))

        correct_root = const_A + const_B - p1 / 3.0 
    else:
        # 3 real roots
        phi = acos(-0.5 * coef_b / sqrt(-(coef_a ** 3) / 27.0))
        root_1 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0) - p1 / 3.0
        root_2 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0 + 2.0 * np.pi / 3.0) - p1 / 3.0
        root_3 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0 + 4.0 * np.pi / 3.0) - p1 / 3.0

        smallest_root = min(min(root_1,root_2), root_3)
        largest_root = max(max(root_1,root_2), root_3)

        if fluid_type is 'liquid':        
            correct_root = smallest_root
        else:
            assert fluid_type is 'vapor', 'Wrong fluid type! ' + fluid_type
            correct_root = largest_root
    
    assert correct_root > 0.0, fluid_type + ' Z-factor < 0.0! Delta is %f, %f' % (delta, correct_root)
    
    return smallest_root,largest_root


def calculate_roots_of_cubic_eos(p0, p1, p2, p3):
    coef_a = (3.0 * p2 - (p1 ** 2)) / 3.0        
    coef_b = (2.0 * (p1 ** 3) - 9.0 * p1 * p2 + 27.0 * p3) / 27.0        
    delta = 0.25 * (coef_b ** 2) + (coef_a ** 3) / 27.0     

    roots = []
    if delta > 0.0:
        # 1 real root, 2 imaginary                 
        const_A =  cbrt(-0.5 * coef_b + sqrt(delta)) 
        const_B =  cbrt(-0.5 * coef_b - sqrt(delta))
        
        single_root = const_A + const_B - p1 / 3.0
        
        assert single_root > 0.0, 'Z-factor < 0.0! Delta is %f, %f' % (delta, single_root)
        
        roots.append(single_root) 
    else:
        # 3 real roots
        phi = acos(-0.5 * coef_b / sqrt(-(coef_a ** 3) / 27.0))
        root_1 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0) - p1 / 3.0
        root_2 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0 + 2.0 * np.pi / 3.0) - p1 / 3.0
        root_3 = 2.0 * sqrt(-coef_a / 3.0) * cos(phi / 3.0 + 4.0 * np.pi / 3.0) - p1 / 3.0

        smallest_root = min(min(root_1,root_2), root_3)
        assert smallest_root > 0.0, 'Z-factor < 0.0! Delta is %f, %f' % (delta, smallest_root)
        
        largest_root = max(max(root_1,root_2), root_3)
        assert largest_root > 0.0, 'Z-factor < 0.0! Delta is %f, %f' % (delta, largest_root)
        
        roots.append(smallest_root)
        roots.append(largest_root)
    
    return roots

def flash_residual_function(x, temperature, pressure, eos, global_molar_fractions):
    T = temperature
    P = pressure
    z = global_molar_fractions
    
    size = x.shape[0]
    
    K = x[0:size-1] # K-values
    K_minus_one = (K - 1.0)
    
    F_V = x[size-1]

    if F_V < 0.0:
        F_V = 0.0
    if F_V > 1.0:
        F_V = 1.0
        
    x_L = z / (F_V * K_minus_one + 1.0)
    x_V = K * x_L
    
    # Vapor
    eos.set_molar_fractions(x_V)
    eos.update_eos_coefficients(P, T)
    z_factor = eos.calculate_eos_roots('vapor')    
    f_V = eos.calculate_fugacities(P, T, z_factor)
    
    # Liquid
    eos.set_molar_fractions(x_L)
    eos.update_eos_coefficients(P, T)
    z_factor = eos.calculate_eos_roots('liquid')        
    f_L = eos.calculate_fugacities(P, T, z_factor)
    
    residual_fugacity = f_L - f_V 
    residual_mass = np.sum(z * K_minus_one / (1.0 + F_V * K_minus_one))
    residual = np.append(residual_fugacity, residual_mass)

    return residual        
    
def func_rachford_rice(x, global_molar_fractions, K_values):
    z = global_molar_fractions
    c = 1.0 / (K_values - 1.0)
    return np.sum(z / (c + x))
    
def deriv_rachford_rice(x, global_molar_fractions, K_values):
    z = global_molar_fractions
    c = 1.0 / (K_values - 1.0)
    return - np.sum(z / ((c + x) ** 2))
    
def calculate_rachford_rice(global_molar_fractions, K_values):
    min_K = np.min(K_values)
    max_K = np.max(K_values)

    min_val = 0.999 / (1.0 - max_K)
    max_val = 0.999 / (1.0 - min_K) 
  
    F_V = brentq(func_rachford_rice, min_val, max_val, args=(global_molar_fractions, K_values))
    #F_V = newton(func=func_rachford_rice, x0=0.5, fprime=deriv_rachford_rice, args=(global_molar_fractions, K_values))
    #F_V = bisect(func_rachford_rice, min_val, max_val, args=(global_molar_fractions, K_values))

    return F_V

def stability_test_residual_function(x, temperature, pressure, eos, global_molar_fractions, test_type):
    T = temperature
    P = pressure
    z = global_molar_fractions    
    u = x # Getting unknowns

    eos.set_molar_fractions(z)
    eos.update_eos_coefficients(P, T)
    z_factor = eos.calculate_eos_roots(test_type)    
    f_ref = eos.calculate_fugacities(P, T, z_factor)
    
    if test_type is 'vapor':
        other_type = 'liquid'
        K_values = u / z
        x_u = z * K_values
        
    else:
        assert test_type is 'liquid', 'Non existing test_type! ' + test_type
        other_type = 'vapor'
        K_values = z / u
        x_u = z / K_values
    
    x_u_normalized = x_u / np.sum(x_u)
    
    eos.set_molar_fractions(x_u_normalized)
    eos.update_eos_coefficients(P, T)
    z_factor = eos.calculate_eos_roots(other_type)        
    f_u = eos.calculate_fugacities(P, T, z_factor)
    
    residual = f_ref - f_u * np.sum(x_u)

    return residual

def calculate_K_values_wilson(
    pressure,
    temperature,
    critical_pressure,
    critical_temperature,
    acentric_factor
):
    P = pressure
    T = temperature
    P_c = critical_pressure
    T_c = critical_temperature
    omega = acentric_factor
    
    return (P_c / P) * np.exp(5.37 * (1.0 + omega) * (1.0 - (T_c / T)))
    
class EquationOfState():
    def __init__(
        self,
        critical_pressure, 
        critical_temperature,
        acentric_factor,
        omega_a,
        omega_b,
        binary_interaction
    ):
        self.critical_pressure = critical_pressure
        self.critical_temperature = critical_temperature
        self.acentric_factor = acentric_factor
        self.omega_a = omega_a
        self.omega_b = omega_b
        self.binary_interaction = binary_interaction
    
    def set_molar_fractions(self, molar_fractions):
        self.molar_fractions = molar_fractions
    
    def alpha_function(self, temperature, critical_temperature, acentric_factor):
        return 1.0
    
    def calculate_eos_coefficients(self, pressure, temperature):
        alpha = self.alpha_function(temperature, self.critical_temperature, self.acentric_factor)
        
        a = (self.omega_a * alpha * (R * self.critical_temperature) ** 2) / self.critical_pressure
        b = (self.omega_b * R * self.critical_temperature) / self.critical_pressure
        a *= pressure / (R * temperature) ** 2
        b *= pressure / (R * temperature)

        return a, b

    def update_eos_coefficients(self, pressure, temperature):
        x = self.molar_fractions
        
        BI = self.binary_interaction

        self.a, self.b = self.calculate_eos_coefficients(pressure, temperature)
        
        #AIJ = (1.0 - BI) * np.sqrt( self.a.reshape( [self.a.size,1] ) * self.a )        
        AIJ = (1.0 - BI) * np.sqrt(np.einsum('i,j', self.a[:], self.a[:]))

        # This variables will be used in the fugacity expression
        self.numerator_coef = np.einsum('ij,j', AIJ, x)
        
        self.mixture_a = np.dot(np.dot(x, AIJ), x)
        self.mixture_b = np.sum(x * self.b)
        
    def calculate_normalized_gibbs_energy(fugacities):        
        normalized_gibbs_energy = 0.0
        
        for x, fugacity in izip(self.molar_fractions, fugacities):
            normalized_gibbs_energy += x[i]*log(fugacity)
            
        return normalized_gibbs_energy
    
#    def calculate_fugacities_with_minimum_gibbs_energy(self, pressure, temperature):
#        # TODO: Work in progress, calculate fugacities by 
#        # calculating all roots and if it has two possible roots
#        # calculate both minimim gibbs energy and choose the 
#        # group of fugacities with minimum gibbs energy
#               
#        self.update_eos_coefficients(pressure, temperature)
#        #for Z in z_factors:
        
class SoaveRedlichKwongEos(EquationOfState):      
    def alpha_function(self, temperature, critical_temperature, acentric_factor):
        m = 0.480 + 1.574 * acentric_factor - 0.176 * (acentric_factor ** 2)
        return (1.0 + m * (1.0 - np.sqrt(temperature / critical_temperature))) ** 2
        
    def calculate_eos_roots(self, fluid_type):
        A_mix = self.mixture_a
        B_mix = self.mixture_b        
        
        p0 = 1.0
        p1 = -1.0
        p2 = A_mix - B_mix - (B_mix ** 2) 
        p3 = -(A_mix * B_mix)     

        return find_correct_root_of_cubic_eos(p0, p1, p2, p3, fluid_type)       

    def calculate_all_eos_roots(self, fluid_type):
        A_mix = self.mixture_a
        B_mix = self.mixture_b        
        
        p0 = 1.0
        p1 = -1.0
        p2 = A_mix - B_mix - (B_mix ** 2) 
        p3 = -(A_mix * B_mix)     

        return find_all_roots_of_cubic_eos(p0, p1, p2, p3, fluid_type)       


        
    def calculate_fugacities(self,pressure,temperature,z_factor):
        P = pressure
        T = temperature
        Z = z_factor
        x = self.molar_fractions        
        a = self.a
        a_mix = self.mixture_a
        b = self.b
        b_mix = self.mixture_b        
        sum_x_j_A_ij = self.numerator_coef
        SQRT_2 = sqrt(2.0)
        
        ln_f = (b / b_mix)*(Z - 1.0) - np.log( Z - b_mix ) \
             + (a_mix / b_mix) \
             * ( (b / b_mix) - 2.0 * sum_x_j_A_ij / a_mix ) \
             * np.log(1.0 + (b_mix / Z))
    
        return (x * P) * np.exp(ln_f) # [Pa]            

    def calculate_fugacities_with_minimum_gibbs_energy(self, pressure, temperature,fluid_type):
        # TODO: Work in progress, calculate fugacities by 
        # calculating all roots and if it has two possible roots
        # calculate both minimim gibbs energy and choose the 
        # group of fugacities with minimum gibbs energy
               
        self.update_eos_coefficients(pressure, temperature)
        z_factors = self.calculate_all_eos_roots(self, fluid_type)
        fugacities = np.zeros(np.size(z_factors))
        for Z in z_factors:
            calculate_fugacities(self,pressure,temperature,z_factor[Z])
        normalized_gibbs_energy = self.calculate_normalized_gibbs_energy(fugacities)                    
        index = np.argmin(normalized_gibbs_energy)
        return fugacities[index]

    
def ss_flash(
    eos,
    pressure, 
    temperature, 
    global_molar_fractions, 
    initial_K_values, 
    max_iter = 250,
    tolerance = 1.0e-3,
    print_statistics=False
):
    P = pressure
    T = temperature
    z = global_molar_fractions
    K = np.copy(initial_K_values)

    # Initialize error with some value
    error = 100.0

    counter = 0
    while error > tolerance and counter < max_iter:
        K_minus_one = K - 1.0

        F_V = calculate_rachford_rice(global_molar_fractions, K)
        if F_V > 1:
            F_V = 1
        if F_V < 0:
            F_V = 0            
        x_L = z / (F_V * K_minus_one + 1.0)
        x_V = K * x_L

        # Vapor
        eos.set_molar_fractions(x_V)
        eos.update_eos_coefficients(P, T)
        z_factor = eos.calculate_eos_roots('vapor')    
        f_V = eos.calculate_fugacities(P, T, z_factor)

        # Liquid
        eos.set_molar_fractions(x_L)
        eos.update_eos_coefficients(P, T)
        z_factor = eos.calculate_eos_roots('liquid')
        f_L = eos.calculate_fugacities(P, T, z_factor)

        f_ratio = f_L / f_V

        K *= f_ratio.ravel()

        error = np.linalg.norm(f_ratio - 1)
        counter += 1
    
    if print_statistics:
        print ('SS Flash: %d iterations, error is %g.' %(counter, error))
        
    return K, F_V, f_L
    
def ss_stability_test(eos,pressure,temperature,global_molar_fractions,test_type,
    initial_K_values, max_iter = 100,tolerance = 1.0e-5):
    P = pressure
    T = temperature
    z = global_molar_fractions
    K = np.copy(initial_K_values)

    error = 100.0

    eos.set_molar_fractions(z)
    eos.update_eos_coefficients(P, T)
    z_factor = eos.calculate_eos_roots(test_type)    
    f_ref = eos.calculate_fugacities(P, T, z_factor)
        
    counter = 0
    while error > tolerance and counter < max_iter:
        if test_type is 'vapor':
            other_type = 'liquid'
            x_u = z * K
        else:
            assert test_type is 'liquid', 'Non existing test_type! ' + test_type
            other_type = 'vapor'
            x_u = z / K        
        
        sum_x_u = np.sum(x_u)
        x_u_normalized = x_u / sum_x_u
        
        eos.set_molar_fractions(x_u_normalized)
        eos.update_eos_coefficients(P, T)
        z_factor = eos.calculate_eos_roots(other_type)        
        f_u = eos.calculate_fugacities(P, T, z_factor)        
        
        if test_type is 'vapor':            
            correction = f_ref / (f_u * sum_x_u)
        else:
            assert test_type is 'liquid', 'Non existing test_type! ' + test_type
            correction = (f_u * sum_x_u) / f_ref      
        K *= correction.ravel()    
        error = np.linalg.norm(correction - 1.0)
        counter += 1

    return sum_x_u, K
################################################################################
def calculate_stability_test(
    eos,
    pressure, 
    temperature, 
    global_molar_fractions,
    initial_K_values
):    
    sum_vapor, K_values_vapor = ss_stability_test(
        eos,
        pressure, 
        temperature, 
        global_molar_fractions,
        'vapor',
        initial_K_values
    )

    sum_liquid, K_values_liquid = ss_stability_test(
        eos,
        pressure, 
        temperature, 
        global_molar_fractions,
        'liquid',
        initial_K_values
    )
    
    sum_ln_K_vapor = np.linalg.norm(np.log(K_values_vapor)) ** 2
    sum_ln_K_liquid = np.linalg.norm(np.log(K_values_liquid)) ** 2
    sum_tol = 1.0e-8

    if sum_ln_K_vapor < 1.0e-4 and sum_ln_K_liquid < 1.0e-4:
        is_stable = True
    elif (sum_vapor-1.0) <= sum_tol and sum_ln_K_liquid < 1.0e-4:
        is_stable = True
    elif (sum_liquid-1.0) <= sum_tol and sum_ln_K_vapor < 1.0e-4:
        is_stable = True
    elif (sum_vapor-1.0) <= sum_tol and (sum_liquid-1.0) <= sum_tol:
        is_stable = True
    elif (sum_vapor-1.0) > sum_tol and sum_ln_K_liquid < 1.0e-4:
        is_stable = False
    elif (sum_liquid-1.0) > sum_tol and sum_ln_K_vapor < 1.0e-4:
        is_stable = False
    elif (sum_vapor-1.0) > sum_tol and (sum_liquid-1.0) > sum_tol:
        is_stable = False
    elif (sum_vapor-1.0) > sum_tol and (sum_liquid-1.0) <= sum_tol:
        is_stable = False
    elif (sum_vapor-1.0) <= sum_tol and (sum_liquid-1.0) > sum_tol:
        is_stable = False
    else:
        assert False, 'ERROR: No stability condition found...'

    if not is_stable:
        K_values_estimates = K_values_vapor * K_values_liquid
    else:
        K_values_estimates = np.copy(initial_K_values)

    return is_stable, K_values_estimates        
    
################################################################################
# Output Calc
################################################################################


def calculate_vapor_liquid_equilibrium(
    eos, pressure, temperature, 
    global_molar_fractions, K_values_estimates, print_statistics=False):
    size = global_molar_fractions.shape[0]

    is_stable, K_values_est = calculate_stability_test(
        eos,
        pressure, 
        temperature, 
        global_molar_fractions,
        K_values_estimates
    )    
    if not is_stable:
        K_values_from_ss_flash, F_V, f_L = ss_flash(eos, pressure, temperature, 
                                                    global_molar_fractions, K_values_est, 
                                                    tolerance = 1.0e-5, 
                                                    print_statistics=print_statistics)
        
        x0 = np.append(K_values_from_ss_flash, F_V)
        result, infodict, ier, mesg = fsolve(
            func=flash_residual_function,
            x0=x0,
            args=(temperature, pressure, eos, global_molar_fractions),
            full_output=True,
        )        
        if print_statistics:
            print ('Newton flash converged? %d, %s' %(ier, mesg))            
        K_values = result[0:size]
        F_V = result[size]
    else:
        if pressure < 50.0e5:
            F_V = 1.0
        else:
            F_V = 0.0            
        K_values = np.ones(size)        
    return F_V, K_values_estimates
    

def calculate_molar_fraction_curve(
    eos, pressure, temperature, global_molar_fractions, 
    print_statistics=False
):
    size = global_molar_fractions.shape[0]
    T = temperature
    res = []
    K_values = np.ones(size)
    
    for iteration, P in enumerate(pressure): 
        if np.linalg.norm(K_values - 1.0) < 1.0e-3:
            # Estimate initial K-values
            K_values_estimates = calculate_K_values_wilson(
                P,
                T,
                critical_pressure,
                critical_temperature,
                acentric_factor
            )
        else:
            K_values_estimates = np.copy(K_values)

        if print_statistics:
            print ('Pressure: %g bar' %(P/1.0e5))
        F_V, K_values = calculate_vapor_liquid_equilibrium(
            eos, P, T, 
            global_molar_fractions, 
            K_values_estimates, print_statistics)
        #print P/1.0e5, F_V

        res.append(F_V)

    return np.array(res)        
        

    
################################################################################
# Test Data
################################################################################


def input_properties_case_whitson_problem_18_SRK():
    '''
    TEST PROBLEM PHASE BEHAVIOUR WHITSON PROBLEM 18 APPENDIX
    
    Methane, Butane and Decane (C1, C4 and C10).
    
    Properties for the Soave-Redlich-Kwong Equation of State.
    
    '''
    temperature = (280.0 + 459.67) * 5.0 / 9.0
    pressure = 500.0 * 6894.75729

    critical_pressure = 6894.75729 * np.array([667.8, 550.7, 304.0]) # [atm]
    critical_temperature = (5.0 / 9.0) * np.array([343.0, 765.3, 1111.8]) # [K]
    acentric_factor = np.array([0.011500, 0.192800, 0.490200]) # [-]
    molar_mass = 0.001 * np.array([16.04, 58.12, 142.29]) # [g/mol]
    omega_a = 0.42748 * np.array([1.0, 1.0, 1.0]) # [-]
    omega_b = 0.08664 * np.array([1.0, 1.0, 1.0]) # [-]

    binary_interaction = np.array(
    [[0.000000,  0.000000, 0.000000],
     [0.000000,  0.000000, 0.000000],
     [0.000000,  0.000000, 0.000000]]
    )

    global_molar_fractions = np.array([0.5, 0.42, 0.08])

    return (pressure, temperature, global_molar_fractions, 
        critical_pressure, critical_temperature, acentric_factor,
        molar_mass, omega_a, omega_b, binary_interaction)    
################################################################################

################################################################################
# [kinda] Stolen Classes
################################################################################
class SRKMIX(object):
    '''
    '''    
    def __init__(self, Tcs, Pcs, omegas, zs, kijs=None, T=None, P=None, V=None,R=82.057338	):
        # Read in component parameters
        self.N = len(Tcs)
        self.cmps = range(self.N)
        self.Tc = Tcs
        self.Pc = Pcs
        self.z = zs
        self.omega = omegas
        self.kij = kijs
        # Ambient Parameters
        self.T = T
        self.P = P
        self.R = R
        self.V = V
        # SRK Parameters
        self.Tc = Tcs
        self.Pc = Pcs
        self.omega = omegas
        self.zs = zs
        if kijs is None:
            kijs = [[0]*self.N for i in range(self.N)]
        self.kijs = kijs
        
    def ComputeMixParameters(self,zc,P,T):
        
        self.b      = 0.08664 * self.R * self.Tc / self.Pc
        self.bm     = np.sum(zc*self.b)
        self.a      = 0.42747 * self.R * self.R * self.Tc * self.Tc / self.Pc
        self.alpha  = ( (1.0 + 0.480 + 1.574*self.omega - 0.176*self.omega*self.omega) * (1.0 - np.sqrt(T/self.Tc)) )**2
        self.A  = self.alpha*self.a
        self.B  = self.b
        self.am     = np.sum( zc.reshape([np.size(zc),1])*zc *                 \
                      np.sqrt(self.A.reshape([np.size(self.A),1]) *    \
                      self.A) * (1.0 - self.kij) )

        self.AIJ = np.sqrt( A.reshape( [A.size,1] ) * A) * (1 - kij)
       
       # Reduced Values
        self.AM     = self.am * P / (self.R * self.R * T * T)
        self.BM     = self.bm * P / (self.R * T)
        
        # Cubic EOS coefficients
        self.c0     =  1.0
        self.c1     = -1.0
        self.c2     = self.AM - self.BM - self.BM*self.BM
        self.c3     = -1.0*self.AM*self.BM
    
        
        
    def SolveCubicEOS(self):
        
        coeffs = np.zeros(4)
        coeffs[0] = self.c0
        coeffs[1] = self.c1
        coeffs[2] = self.c2
        coeffs[3] = self.c3
        z = np.real( np.roots(coeff)[np.flatnonzero(np.iscomplex(np.roots(coeff))==False)])
        ZPair = np.array([z.min(),z.max()])
        return ZPair        
    
    
    def ComputeDensity(self,zc,P,T):
        
        # Update mix parameters
        self.ComputeMixParameters(zc,P,T)
        
        # Compute Average AMU
        self.AMU = np.sum(zc*self.MW)
        
        
        
    def ComputelnPhi(self,ZPair,x,P,T):       
        
        A = self.A
        B = self.B
        if self.N > 1:
            AM  = self.AM
            BM  = self.BM
            AIJ = self.AIJ
            lf  =  B.reshape([B.size,1])/BM * (ZPair-1) -                     \
              np.log(ZPair-BM) +                                              \
              AM/BM*( B.reshape([B.size,1])/BM -                              \
              2./AM *                                                         \
              np.sum(x*AIJ,axis=0).reshape([x.size,1]) ) * np.log(1. + BM/ZPair)
        
        else:
            lf = ZPair - 1. - np.log(ZPair - B) - A/B*np.log(1.+B/ZPair)            
        return lf   
    
    def ComputeFugacity(self,lf,x,P):
        fug = np.exp(lf)*x*P
        return fug
        
    def SRKRun(self,x,P,T):
        dG,Zc,rho,lnfug = SRK(self.R,P,T,x,self.Pc,self.Tc,self.omega,self.kij)            
        return dG,Zc,rho,lnfug 
##############################################################################80
# "Read In" Parameters
##############################################################################80
c_min = 11      # Only considered for C7 and above
c_max = 200
eptol = 1.0E-8

T_res = 394.25  # K
P_sat = 145.8*0.986923  # atm
Vol_Feed    = 1.0   # cm**3

# Make inputted data useful
TBP_ind   =   11             # Index for first TBD value on input
exp_ind   =   c_min - 6 + 11 # Index for Cn+ value on input table

R  = 82.057338	    # cm3 atm / K *mol
        
##############################################################################80
# "Read In" Data
##############################################################################80
# From Compositional Analysis
F8Comp = pd.read_csv('F8.dat',sep='\s+',header=0)
# Convert from g/cc to kg/m**3
F8Comp['rho'] = F8Comp['rho']
F8Comp['MW'] = F8Comp['MW']
F8Comp.Tres = 394.25 # K
F8Comp.Psat = 145.8*0.98692327 # bar ==> atm
##############################################################################80
# SRK Coefficients
##############################################################################80
coeff       =   np.zeros(13)
coeff[0]    =   304.143
coeff[1]    =   48.4052
coeff[2]    =   0.710774
coeff[3]    =   3800.73
coeff[4]    =   3.05081
coeff[5]    =   -0.903352
coeff[6]    =   233.768
coeff[7]    =   -12715.4
coeff[8]    =   0.25
coeff[9]    =   0.496902
coeff[10]   =   0.00558442
coeff[11]   =   0.010564
coeff[12]   =   -5.243E-6
##############################################################################80
# SRK Shift Parameters
##############################################################################80
shift = np.zeros(11)
shift[0]    =   -0.0079	# Nitrogen
shift[1]    =    0.0833	# CO2
shift[2]    =    0.0466	# H2S
shift[3]    =    0.0234	# Methane
shift[4]    =    0.0605	# Ethane
shift[5]    =    0.0825	# Propane
shift[6]    =    0.083	# i-Butane
shift[7]    =    0.0975	# n-Butane
shift[8]    =    0.1022	# i-Pentane
shift[9]    =    0.1209	# n-Pentane
shift[10]   =    0.1467	# n-Hexane
################################################################################
# Binary Interaction Coefficients (to nC5)
################################################################################
kij = np.zeros([10,10])

kij_part = np.array ([  [0,  	0,  	0.13,	0.025,	0.01,	0.09,	0.095,	0.095,	0.1,	0.11],
                        [0,  	0,	    0.135,	0.105,	0.13,	0.125,	0.12,	0.115,	0.115,	0.115],	
                        [0.13,	0.135,	0,	    0.07,	0.085,	0.08,	0.075,	0.075,	0.07,	0.07]])
kij[:3,:] = kij_part
kij[:,:3] = np.transpose(kij_part)
kij_base = kij.copy()


##############################################################################80
# Try to solve for A,B,C,D, *AGAIN*
##############################################################################80
# Initial Characterization Constant Estimates
A = -3.0
B = -0.1
C =  0.5
D =  0.1



# misc
c6offset = c_min - 7
numcomps = c_max - c_min + 1


# Data Structures
SCN     = np.arange(c_min,c_max+1)
Vec_R   = np.zeros(2)
Vec_R0  = np.zeros(2)
Vec_V   = np.zeros(2)
Mat_J   = np.zeros([2,2])



# Plus Frac z
z_cnplus = np.sum(F8Comp['z'][10+c6offset:].as_matrix())

# Plus Frac MW
MW_cnplus = np.sum( F8Comp['z'][10+c6offset:].as_matrix() *                    \
                    F8Comp['MW'][10+c6offset:].as_matrix() ) /                 \
            np.sum( F8Comp['z'][10+c6offset:].as_matrix() )

# Plus Frac rho
rho_cnplus = np.sum( F8Comp['z'][10+c6offset:].as_matrix() *                   \
                    F8Comp['MW'][10+c6offset:].as_matrix() ) /                 \
                    np.sum( (F8Comp['z'][10+c6offset:].as_matrix() *           \
                    F8Comp['MW'][10+c6offset:].as_matrix() ) / F8Comp['rho'][10+c6offset:].as_matrix() )

# Last TBP Frac rho
rho_cnminus = F8Comp['rho'][10+c6offset-1]



# Initialize solution vector
Vec_V[0]    = A
Vec_V[1]    = B


# Solve for A and B
for i in range(0,1000):
# Initial residual function
    Vec_R0[:]   = Vec_R[:] 
    Vec_R[0]    = Ped_R1(z_cnplus,A,B,SCN)
    Vec_R[1]    = Ped_R2(z_cnplus,MW_cnplus,A,B,SCN)

    Mat_J[0,0] =  (Ped_R1(z_cnplus,A+eptol,B,SCN) - Ped_R1(z_cnplus,A-eptol,B,SCN)) / (2*eptol)
    Mat_J[0,1] =  (Ped_R1(z_cnplus,A,B+eptol,SCN) - Ped_R1(z_cnplus,A,B-eptol,SCN)) / (2*eptol)
    Mat_J[1,0] =  (Ped_R2(z_cnplus,MW_cnplus,A+eptol,B,SCN) - Ped_R2(z_cnplus,MW_cnplus,A-eptol,B,SCN)) / (2*eptol)
    Mat_J[1,1] =  (Ped_R2(z_cnplus,MW_cnplus,A,B+eptol,SCN) - Ped_R2(z_cnplus,MW_cnplus,A,B-eptol,SCN)) / (2*eptol)            
        
    # Get Delta V result
    Vec_dV = np.linalg.solve(Mat_J,-1.*Vec_R)
    # Update V
    Vec_V = Vec_V + Vec_dV
    A = Vec_V[0]
    B = Vec_V[1]
    if np.all(np.abs(Vec_R0 - Vec_R) < eptol*2):
        break
    
    #print(Vec_V, Vec_R)
    
# Construct MF array for expansion
z_exp   = Ped_z(A,B,SCN)
MW_exp  = Ped_MW(SCN)


# Solve for C and D

Vec_V[0]    = C
Vec_V[1]    = D
Vec_R0[:]   = 0.0
Vec_R[:]    = 0.0

for i in range(0,1000):
    Vec_R0[:]   = Vec_R[:] 
    Vec_R[0]    = Ped_R3(rho_cnplus,A,B,C,D,SCN)
    Vec_R[1]    = Ped_R4(rho_cnminus,C,D,c_min-1)

    Mat_J[0,0] =  ( Ped_R3(rho_cnplus,A,B,C+eptol,D,SCN)  - Ped_R3(rho_cnplus,A,B,C-eptol,D,SCN)  ) / (2*eptol)
    Mat_J[0,1] =  ( Ped_R3(rho_cnplus,A,B,C,D+eptol,SCN)  - Ped_R3(rho_cnplus,A,B,C,D-eptol,SCN)  ) / (2*eptol)
    Mat_J[1,0] =  ( Ped_R4(rho_cnminus,C+eptol,D,c_min-1) - Ped_R4(rho_cnminus,C-eptol,D,c_min-1) ) / (2*eptol)
    Mat_J[1,1] =  ( Ped_R4(rho_cnminus,C,D+eptol,c_min-1) - Ped_R4(rho_cnminus,C,D-eptol,c_min-1) ) / (2*eptol)
        
    # Get Delta V result
    Vec_dV = np.linalg.solve(Mat_J,-1.*Vec_R)
    # Update V
    Vec_V = Vec_V + Vec_dV
    C = Vec_V[0]
    D = Vec_V[1]
    if np.all(np.abs(Vec_R0 - Vec_R) < eptol*2):
        break
    
    #print(Vec_V, Vec_R)


##############################################################################80
# Cn+ Property Estimates (c_min to c_max)
##############################################################################80
z_cnexp         = Ped_z(A,B,SCN)
rho_cnexp       = Ped_rho(C,D,SCN)
MW_cnexp        = Ped_MW(SCN)
Tc_cnexp        = Ped_Tc(coeff,rho_cnexp,MW_cnexp)
Pc_cnexp        = Ped_Pc(coeff,rho_cnexp,MW_cnexp)
omega_cnexp     = Ped_omega(coeff,rho_cnexp,MW_cnexp)

##############################################################################80
# Plotting / Overall Data
##############################################################################80
z           =   np.zeros(c_max - 7 + 1)
Pc          =   np.zeros(c_max - 7 + 1)
Tc          =   np.zeros(c_max - 7 + 1)
MW          =   np.zeros(c_max - 7 + 1)
omega       =   np.zeros(c_max - 7 + 1)
rho         =   np.zeros(c_max - 7 + 1)
C_num       =   np.arange(7,c_max + 1)

z_exp           =   np.zeros(c_max - 7 + 12)
y_exp           =   np.zeros(c_max - 7 + 12)
x_exp           =   np.zeros(c_max - 7 + 12)
Pc_exp          =   np.zeros(c_max - 7 + 12)
Tc_exp          =   np.zeros(c_max - 7 + 12)
MW_exp          =   np.zeros(c_max - 7 + 12)
omega_exp       =   np.zeros(c_max - 7 + 12)
rho_exp         =   np.zeros(c_max - 7 + 12)
comps_exp       =   np.ndarray(c_max - 7 + 12)
comps_exp       =   comps_exp.astype('object')
om_a_exp = 0.42748 * np.ones(c_max - 7 + 12) # [-]
om_b_exp = 0.08664 * np.ones(c_max - 7 + 12) # [-]


z_in        =   F8Comp['z'].as_matrix()
MW_in       =   F8Comp['MW'].as_matrix()
rho_in      =   F8Comp['rho'].as_matrix()
y_in        =   F8Comp['y'].as_matrix()
Pc_in       =   Ped_Pc(coeff,rho_in,MW_in)
Tc_in       =   Ped_Tc(coeff,rho_in,MW_in)
omega_in    =   Ped_omega(coeff,rho_in,MW_in)
om_a_in = 0.42748 * np.ones(z_in.size) # [-]
om_b_in = 0.08664 * np.ones(z_in.size) # [-]

z[-1*z_cnexp.size:] = z_cnexp[:]
Pc[-1*z_cnexp.size:] = Pc_cnexp[:]
Tc[-1*z_cnexp.size:] = Tc_cnexp[:]
MW[-1*z_cnexp.size:] = MW_cnexp[:]
omega[-1*z_cnexp.size:] = omega_cnexp[:]
rho[-1*z_cnexp.size:] = rho_cnexp[:]

# MF Gasses:
z_c6minus = np.sum(z_in[0:11])


if c_min > 7:
    i           = c_min - 7
    z[:i]       = z_in[11:i     + 11]
    MW[:i]      = MW_in[11:i    + 11]
    rho[:i]     = rho_in[11:i   + 11]    
    Pc[:i]      = Pc_in[11:i    + 11]    
    Tc[:i]      = Tc_in[11:i    + 11]        
    omega[:i]   = omega_in[11:i + 11]        

# Derived Values

MassFrac  = (z*MW) / np.sum(z*MW) * 100
cummass   = np.cumsum(z) / np.sum(z)
#cummass[c_min:]   = z[c_min:]
#cummass[:c_min]   =   np.cumsum(cummass[:c_min])
#cummass[c_min:]   =   np.cumsum(cummass[c_min:])  + (1 - z_cnplus)
    
    
##############################################################################80
# Generate Plots of Extrapolated Values
##############################################################################80
# MF, rho, cmin to cmax
fig, ax_MF  = plt.subplots(figsize=(10, 6))
ax_rho       = ax_MF.twinx()
ax_MF.plot(C_num,z,'r--')   
ax_rho.plot(C_num,rho,'b')
plt.tight_layout()

ax_MF.set_xlim(7,SCN[-1])
ax_MF.set_ylabel('Mole Fraction')
ax_MF.set_xlabel('Carbon Number')
ax_rho.set_ylabel('Density, g/cc')
plt.savefig('1_rho_MF.png',frameon=True,bbox_inches='tight',padinches=0.1)
plt.clf()

# Pc, Tc, cmin to cmax
fig, ax_Pc  = plt.subplots(figsize=(10, 5))
ax_Tc       = ax_Pc.twinx()
ax_Tc.plot(C_num,Tc,'r--')   
ax_Pc.plot(C_num,Pc,'b')
ax_Pc.set_xlim(7,SCN[-1])
plt.tight_layout()

ax_Tc.set_ylabel('Critical Temperature, Deg. Kelvin')
ax_Pc.set_xlabel('Carbon Number')
ax_Pc.set_ylabel('Critical Pressure, atm')
plt.savefig('1_Pc_Tc.png',frameon=True,bbox_inches='tight',padinches=0.1)
plt.clf()

# Mass Frac, % of total cmin to cmax
fig, ax_MW  = plt.subplots(figsize=(10, 6))
ax_sum      = ax_MW.twinx()
ax_MW.plot(C_num,MassFrac,'r--')   
ax_sum.plot(C_num,cummass,'b')
plt.tight_layout()

ax_sum.set_ylabel('Cumulative Mass, %')
ax_sum.set_ylim(0,1)
ax_MW.set_xlim(7,SCN[-1])
ax_MW.set_xlabel('Carbon Number')
ax_MW.set_ylabel('Mass Fraction of Plus')
plt.savefig('1_MW_fracsum.png',frameon=True,bbox_inches='tight',padinches=0.1)
plt.clf()

##############################################################################80
# Lump into groups
##############################################################################80
# Initial Guess
C_G1 = np.arange(7,c_min) 
C_G2 = np.arange(c_min,20)
C_G3 = np.arange(20,30)
C_G4 = np.arange(30,200+1)

# Adjust for index commencing with 7
C_G1 = C_G1 - 7
C_G2 = C_G2 - 7
C_G3 = C_G3 - 7
C_G4 = C_G4 - 7

# 14 base components defined + 4 PCs
ncomps_PS   =   11+4
z_PS        =   np.zeros(ncomps_PS)
Tc_PS       =   np.zeros(ncomps_PS)
Pc_PS       =   np.zeros(ncomps_PS)
MW_PS       =   np.zeros(ncomps_PS)
omega_PS    =   np.zeros(ncomps_PS)
x_PS        =   np.zeros(ncomps_PS)
y_PS        =   np.zeros(ncomps_PS)
om_a_PS = 0.42748 * np.ones(ncomps_PS) # [-]
om_b_PS = 0.08664 * np.ones(ncomps_PS) # [-]

# PS BIPs

kij_PS      = np.zeros([ncomps_PS,ncomps_PS])
kij_PS[:10,:10] = kij
kij_PS[10:,:3]  = kij[-1,:3]
kij_PS[:3,10:]  = kij[:3,-1:]
#  N2
MW_PS[0]      = 28.016
Pc_PS[0]      = 33.6
Tc_PS[0]      = 126.2
omega_PS[0]   = 0.040

# CO2
MW_PS[1]    = 44.01
Pc_PS[1]    = 72.9
Tc_PS[1]    = 304.2
omega_PS[1] = 0.228 

# H2S
MW_PS[2]    = 34.08
Pc_PS[2]    = 88.9
Tc_PS[2]    = 373.2
omega_PS[2] = 0.081

# C1
MW_PS[3]  = 16.043
Pc_PS[3]  = 45.4
Tc_PS[3]  = 190.4
omega_PS[3]  = 0.008

# C2
MW_PS[4]  = 30.069
Pc_PS[4]  = 48.2
Tc_PS[4]  = 305.4
omega_PS[4]  = 0.098

# C3
MW_PS[5]  = 44.096
Pc_PS[5]  = 41.9
Tc_PS[5]  = 369.8
omega_PS[5]  = 0.152

# iC4
MW_PS[6]  = 58.123
Pc_PS[6]  = 36.0
Tc_PS[6]  = 408.2
omega_PS[6]  = 0.176

# nC4
MW_PS[7]  = 58.123
Pc_PS[7]  = 37.5
Tc_PS[7]  = 425.2
omega_PS[7]  = 0.193

# iC5
MW_PS[8]  = 72.15
Pc_PS[8]  = 33.4
Tc_PS[8]  = 460.4
omega_PS[8]  = 0.227 

# nC5
MW_PS[9]  = 72.15
Pc_PS[9]  = 33.3
Tc_PS[9]  = 469.7
omega_PS[9]  = 0.251

# C6
MW_PS[10] = 86.177
Pc_PS[10] = 29.3
Tc_PS[10] = 507.5
omega_PS[10] = 0.296

# Convert to kg/m**3
#MW_PS[0:10] = MW_PS[0:10]/1000

# Manually input gas cut
y_PS[0] = 0
y_PS[1] = 0
y_PS[2] = 0
y_PS[3] = 0.88
y_PS[4] = 0.07
y_PS[5] = 0.03
y_PS[6] = 0.005
y_PS[7] = 0.005
y_PS[8] = 0.005
y_PS[9] = 0.005
y_PS[10] = 0
#x_PS = x_PS - y_PS


MW_in[0:11] = MW_PS[0:11]
Pc_in[0:11] = Pc_PS[0:11]
Tc_in[0:11] = Tc_PS[0:11]
omega_in[0:11] = omega_PS[0:11]
y_in[0:11] = y_PS[0:11]
#x_in = 1.0 - y_in


# MF for above 11

z_PS[:11] = z_in[:11]
#===============================================================================
# Update Pandas dataframe for input case for use as a baseline
#F8Comp['x'] = x_in
F8Comp['Pc'] = Pc_in
F8Comp['Tc'] = Tc_in
F8Comp['omega'] = omega_in
kij = np.zeros([Pc_in.size,Pc_in.size])
kij[:10,:10] = kij_base
kij[10:,0] = kij_base[9,0]
kij[10:,1] = kij_base[9,1]
kij[0,9:] = kij_base[0,9]
kij[1,9:] = kij_base[1,9]
F8Comp.kij = kij.copy()
F8Comp['om_a'] = om_a_in
F8Comp['om_b'] = om_b_in
# Calculate Liquid Components at given conditions

K_in = calculate_K_values_wilson(F8Comp.Psat,F8Comp.Tres,Pc_in,Tc_in,omega_in)
eos = SoaveRedlichKwongEos(Pc_in,Tc_in,omega_in, om_a_in, om_b_in, kij)

#===============================================================================
# Calculate parameters for rest to C_min    

#Pc_PS[TBP_ind:c_min-1] = Pc[TBP_ind:c_min]
#Tc_PS[TBP_ind:c_min-1] = Tc[TBP_ind:c_min-1]
#MW_PS[TBP_ind:c_min-1] = MW[TBP_ind:c_min-1]

# C_G1
Tc_PS[-4]    = PS_Tc(z[C_G1],MW[C_G1],Tc[C_G1])
Pc_PS[-4]    = PS_Pc(z[C_G1],MW[C_G1],Pc[C_G1])
omega_PS[-4] = PS_omega(z[C_G1],MW[C_G1],omega[C_G1])
z_PS[-4]     = np.sum(z[C_G1])
MW_PS[-4]     = PS_MW(z[C_G1],MW[C_G1])

# C_G2
Tc_PS[-3]    = PS_Tc(z[C_G2],MW[C_G2],Tc[C_G2])
Pc_PS[-3]    = PS_Pc(z[C_G2],MW[C_G2],Pc[C_G2])
omega_PS[-3] = PS_omega(z[C_G2],MW[C_G2],omega[C_G2])
z_PS[-3]     = np.sum(z[C_G2])
MW_PS[-3]     = PS_MW(z[C_G2],MW[C_G2])

# C_G3
Tc_PS[-2]    = PS_Tc(z[C_G3],MW[C_G3],Tc[C_G3])
Pc_PS[-2]    = PS_Pc(z[C_G3],MW[C_G3],Pc[C_G3])
omega_PS[-2] = PS_omega(z[C_G3],MW[C_G3],omega[C_G3])
z_PS[-2]     = np.sum(z[C_G3])
MW_PS[-2]     = PS_MW(z[C_G3],MW[C_G3])

# C_G4
Tc_PS[-1]    = PS_Tc(z[C_G4],MW[C_G4],Tc[C_G4])
Pc_PS[-1]    = PS_Pc(z[C_G4],MW[C_G4],Pc[C_G4])
omega_PS[-1] = PS_omega(z[C_G4],MW[C_G4],omega[C_G4])
z_PS[-1]     = np.sum(z[C_G4])
MW_PS[-1]     = PS_MW(z[C_G4],MW[C_G4])

# Fix excess mass discrepancy
z_PS[np.argmax(z_PS[-4:])] = z_PS[np.argmax(z_PS[-4:])] + (1.000000000001 - np.sum(z_PS))
#===============================================================================
# Create dataframe for expanded set
MW_exp[0:ncomps_PS] = MW_PS
Pc_exp[0:ncomps_PS] = Pc_PS
Tc_exp[0:ncomps_PS] = Tc_PS
omega_exp[0:ncomps_PS] = omega_PS
y_exp[0:ncomps_PS] = y_PS
z_exp[0:ncomps_PS] = z_PS

z_exp[-1*z_cnexp.size:] = z_cnexp
MW_exp[-1*MW_cnexp.size:] = MW_cnexp
rho_exp[-1*rho_cnexp.size:] = rho_cnexp
omega_exp[-1*omega_cnexp.size:] = omega_cnexp
Pc_exp[-1*Pc_cnexp.size:] = Pc_cnexp
Tc_exp[-1*Tc_cnexp.size:] = Tc_cnexp

comps_exp[0:ncomps_PS] = F8Comp['comp'].as_matrix()[0:ncomps_PS]
for i in np.arange(ncomps_PS,comps_exp.size - 1):
    name = 'C' + np.str(i-ncomps_PS + c_min)
    comps_exp[i] = name
name = 'C' + np.str(i-ncomps_PS + 1 + c_min) + '+'
comps_exp[-1] = name

#x_exp = 1.0 - y_exp
#F8Exp = pd.DataFrame({'comp':comps_exp,'z':z_exp,'MW':MW_exp,'rho':rho_exp,'y':y_exp,'x':x_exp,
#                      'Pc':Pc_exp,'Tc':Tc_exp,'omega':omega_exp})
F8Exp = pd.DataFrame({'comp':comps_exp,'z':z_exp,'MW':MW_exp,'rho':rho_exp,'y':y_exp,
                      'Pc':Pc_exp,'Tc':Tc_exp,'omega':omega_exp,'om_a':om_a_exp,'om_b':om_b_exp})
kij = np.zeros([Pc_exp.size,Pc_exp.size])
kij[:10,:10] = kij_base
kij[10:,0] = kij_base[9,0]
kij[10:,1] = kij_base[9,1]
kij[0,9:] = kij_base[0,9]
kij[1,9:] = kij_base[1,9]
F8Exp.kij = kij.copy()
F8Exp.Tres = F8Comp.Tres
F8Exp.Pres = F8Comp.Psat
#===============================================================================
# DataFrame for reduced set

comps_PS = comps_exp[:ncomps_PS]
comps_PS[-4] = 'CL1'
comps_PS[-3] = 'CL2'
comps_PS[-2] = 'CL3'
comps_PS[-1] = 'CL4'
#F8Reduced = pd.DataFrame({'comp':comps_PS,'z':z_PS,'MW':MW_PS,'y':y_PS,'x':x_PS,
#                      'Pc':Pc_PS,'Tc':Tc_PS,'omega':omega_PS})
F8Reduced = pd.DataFrame({'comp':comps_PS,'z':z_PS,'MW':MW_PS,'y':y_PS,
                      'Pc':Pc_PS,'Tc':Tc_PS,'omega':omega_PS,'om_a':om_a_PS,'om_b':om_b_PS})
kij = np.zeros([Pc_PS.size,Pc_PS.size])
kij[:10,:10] = kij_base
kij[10:,0] = kij_base[9,0]
kij[10:,1] = kij_base[9,1]
kij[0,9:] = kij_base[0,9]
kij[1,9:] = kij_base[1,9]
F8Reduced.kij = kij.copy()
F8Reduced.Tres = F8Comp.Tres
F8Reduced.Pres = F8Comp.Psat                     
################################################################################
# Simple Verification Set
#-------------------------------------------------------------------------------
compstrs = ["Methane","n-Butane","n-Decane"]
omegas = np.array([0.008,0.199,0.489])
Pcs = np.array([44.806316,37.503084,20.922773])
Tcs = np.array([190.6,425.2,617.7])
Mws = np.array([16.04,58.12,142.29])
zs = np.array([0.35,0.45,0.20])      # Fraction of oil input stream
kijs = np.array([ [0.0000, 0.0133, 0.0422],
                     [0.0133, 0.0000, 0.0078],
                     [0.0422, 0.0078, 0.0000]  ])
xs = np.array([3.18e-1,4.72e-1,2.10e-1])
ys = np.array([9.89e-1,1.13e-2,3.e-6])                    
ls = (ys - zs)/(ys-xs)
SimpleSet = pd.DataFrame({'comp':compstrs,'z':zs,'MW':Mws,'Pc':Pcs,'Tc':Tcs,
                          'omega':omegas,'x':xs,'y':ys,'K':ys/xs,'l':ls,'v':1-ls})
SimpleSet.kij = kijs.copy()
SimpleSet.VF = np.sum(ls*zs)
SimpleSet.LF = 1-np.sum(ls*zs)
SimpleSet.Tres = 75 + 273.15 # K
SimpleSet.Pres = 2.9607698 # atm
################################################################################
# Basic Verification Set
#-------------------------------------------------------------------------------
compstrs = ["CO2","C3"]
omegas = np.array([0.225,0.152])
Pcs = np.array([72.83,41.94]) # atm
Tcs = np.array([304.2,369.8])
Mws = np.array([44.01,44.096])
zs = np.array([0.5,0.5])      # Fraction of oil input stream
kijs = np.array([0.125])
BasicSet = pd.DataFrame({'comp':compstrs,'z':zs,'MW':Mws,'Pc':Pcs,'Tc':Tcs,
                          'omega':omegas})
BasicSet.kij = kijs.copy()

################################################################################
# Single Component Verification Set
#-------------------------------------------------------------------------------
compstrs = ["CO2"]
omegas = np.array([0.225])
Pcs = np.array([72.83]) # atm
Tcs = np.array([304.2]) # K
Mws = np.array([44.01]) # g/mol
zs = np.array([1.0])      # Fraction of oil input stream
#kijs = np.array([0.125])
SingleCompSet = pd.DataFrame({'comp':compstrs,'z':zs,'MW':Mws,'Pc':Pcs,'Tc':Tcs,
                          'omega':omegas})
SingleCompSet.kij = np.array([0])
################################################################################
# Experimental Verification Data
###############################################################################

################################################################################
# Bubble Point Estimates
################################################################################
# CME at Reservoir Temperature
##############################
# P, bar
P           = np.array([500.5,450.5,400,351,300.5,250.5,200,160,145.8,144,138,132,123,112,94,78,62,50])
# relative volume
vrel        = np.array([0.946,0.951,0.957,0.964,0.972,0.98,0.989,0.997,1,1.003,1.016,1.032,1.059,1.103,1.211,1.371,1.652,2.023])
# cell density
celldens    = np.array([762.78,758.73,753.58,748.5,742.39,735.84,729.39,723.59,721.5,719.42,709.72,699.3,681.2,654.02,595.95,526.32,436.87,356.63])
T = 394.2   # K
CME_T_res = pd.DataFrame({'P':P*0.98692327,'vrel':vrel,'rho':celldens})
CME_T_res.Tres = T
# DDE at Reservoir Temperature
##############################
# P, bar==>atm, vrel, liquid density, kg/m**3
DDE = np.array([    [351,	0.9644,	748.3],
                    [300.5,	0.9716,	742.4],
                    [250.5,	0.98,	736.2],
                    [200,	0.989,	729.3],
                    [160,	0.9974,	723.6],
                    [145.8,	1,	    721.5],
                    [114.8,	0.9638,	732],
                    [86.5,	0.9282,	743.8],
                    [60,	0.8927,	757.4],
                    [34.8,	0.8565,	769.8],
                    [16,	0.8222,	780.4] ] )
P = DDE[:,0]                    
vrel = DDE[:,1]
rho = DDE[:,2]
DDE_T_res = pd.DataFrame({'P':P*0.98692327,'vrel':vrel,'rho':rho})
T = 394.2   # K
DDE_T_res.Tres = T

# Swelling Test at Reservoir Temperature
#---------------------------------------
# x gas, psat, bar==>atm, liquid density, kg/m**3

SWELL = np.array([  [0,	        145.8,	721.6],
                    [0.2308,	225.6,	672.6],
                    [0.3333,	276.5,	649.8],
                    [0.3939,	316.4,	636.6],
                    [0.4286,	335.1,	627.3],
                    [0.4595,	355,	620.1]])
x_gas = SWELL[:,0]
psat = SWELL[:,1]
rho = SWELL[:,2]
T = 394.2   # K
SWELL_T_RES = pd.DataFrame({'x_gas':x_gas,'Psat':psat*0.98692327,'rho':rho})
SWELL_T_RES.Tres = T


# CME during swelling test (0)
#--------------------------
x_gas = 0.2308
# P,bar==>atm,  vrel,   cell density, kg/m**3
CME0    = np.array([    [501,	0.9434,	712.76],
                        [451,	0.9508,	707.21],
                        [401,	0.9595,	700.77],
                        [351,	0.9695,	693.48],
                        [301,	0.9807,	685.87],
                        [276,	0.9868,	681.66],
                        [251,	0.9932,	677.05],
                        [225.6,	1,     	672.5],
                        [219.5,	1.0087,	666.67],
                        [203,	1.0363,	648.93],
                        [181,	1.0846,	619.96],
                        [157,	1.1593,	580.05],
                        [126,	1.3132,	512.03],
                        [115,	1.3933,	482.63],
                        [91,	1.6505,	407.5],
                        [73.5,	1.9639,	342.47] ] )
P=CME0[:,0]
vrel=CME0[:,1]
rho=CME0[:,2]
CME0_ST = pd.DataFrame({'P':P*0.98692327,'vrel':vrel,'rho':rho})
CME0_ST.x_gas = x_gas




# CME during swelling test (1)
#--------------------------
x_gas = 0.3333
# P,bar==>atm,  vrel,   cell density, kg/m**3
CME1    =   np.array([  [501,	0.9476,	685.87],
                        [451,	0.9565,	679.35],
                        [401,	0.9669,	672.04],
                        [351,	0.979,	663.57],
                        [326,	0.9856,	659.2],
                        [301,	0.9927,	654.45],
                        [276.5,	1,  	649.77],
                        [261,	1.0174,	638.57],
                        [246,	1.0374,	626.17],
                        [221,	1.0799,	601.68],
                        [202.5,	1.1212,	579.37],
                        [188.5,	1.1599,	560.22],
                        [161,	1.2627,	514.67],
                        [138,	1.3901,	467.51],
                        [101,	1.75,	371.33],
                        [84,	2.0427,	318.07]])

P=CME1[:,0]
vrel=CME1[:,1]
rho=CME1[:,2]
CME1_ST = pd.DataFrame({'P':P*0.98692327,'vrel':vrel,'rho':rho})
CME1_ST.x_gas = x_gas

# CME during swelling test (2)
#--------------------------
x_gas = 0.4286
# P,bar==>atm,  vrel,   cell density, kg/m**3
CME2    =   np.array([  [501,	0.9566,	655.74],
                        [451,	0.9673,	648.51],
                        [426,	0.9734,	644.33],
                        [401,	0.98,	640.21],
                        [376,	0.9872,	635.32],
                        [351,	0.9948,	630.52],
                        [335.1,	1,  	627.35],
                        [326.5,	1.007,	623.05],
                        [312,	1.0202,	615.01],
                        [295,	1.0383,	604.23],
                        [276,	1.0624,	590.32],
                        [249.5,	1.105,	567.54],
                        [228,	1.1499,	545.55],
                        [202.5,	1.2201,	514.14],
                        [170.5,	1.3479,	465.33],
                        [143.5,	1.5125,	414.77],
                        [105.5,	1.9201,	326.69]     ])
 
P=CME2[:,0]
vrel=CME2[:,1]
rho=CME2[:,2]
CME2_ST = pd.DataFrame({'P':P*0.98692327,'vrel':vrel,'rho':rho})
CME2_ST.x_gas = x_gas


################################################################################
# Flash Calculation, given fluid - Saturation Pressure Calculation
# Or, Finding the Bubble Point
#-------------------------------------------------------------------------------
# Loop Flags
mainiter = 15
maxiter = 50
gdem_interval = 10
stability_interval = 5
KGDEM_iter = 10
errtol = 1e-10
Tref = 100 + 273.15 # K
Pref = 1 # atm
eptol = 1e-6

# Initial Values
VF = 0.5


##########################
# Apply reference conditions
#  Temperature
SimpleSet.Tref = Tref
# Pressure
SimpleSet.Pref = Pref

############################################
# Read values from dataframe
############################################
#Pc = SimpleSet.Pc.as_matrix()
#Tc = SimpleSet.Tc.as_matrix()
#z = SimpleSet.z.as_matrix()
#omega = SimpleSet.omega.as_matrix()
#kij = SimpleSet.kij
Pc = F8Comp.Pc.as_matrix()
Tc = F8Comp.Tc.as_matrix()
z = F8Comp.z.as_matrix()
omega = F8Comp.omega.as_matrix()
kij = F8Comp.kij
################################################################################
# Verification using Thermo
# Following CME test from bar to pascal
P_range = np.linspace(500,1,201)*100000. # Pa
T = 394.2        # K

# Baseline Volume
#eos = SRKMIX(T=T,P=P_range[0],Tcs=Tc,Pcs=Pc*101325.,omegas=omega,zs=z,kijs=kij)
#Vs_fast = eos.volume_solutions(T, P_range[0], eos.b, eos.delta, eos.epsilon, eos.a_alpha) 
#Vs = np.real(Vs_fast[np.flatnonzero(np.iscomplex(Vs_fast) == False)])


T_range = np.zeros(P_range.size)
#for i in np.arange(0,P_range.size):
#    eos = th.eos_mix.SRKMIX(T=T,P=P_range[i],Tcs=Tc,Pcs=Pc*101325.,omegas=omega,zs=z,kijs=kij)
#    print eos.phase
    #Vs_fast = eos.volume_solutions(T, P_range[i], eos.b, eos.delta, eos.epsilon, eos.a_alpha) 
    #Vs = np.real(Vs_fast[np.flatnonzero(np.iscomplex(Vs_fast) == False)])    
    #T_range[i] = eos.solve_T(P_range[i],Vs)
#############################################
# VLE Flash for surface conditions
#############################################
P = 106  # SimpleSet.Pref
T = 394.25 # SimpleSet.Tref

######################################
# Generate Bubble Point
#######################################

def gate(m, *data):
    R,x,T,Tc,Pc,omega,kij = data
    K = m[:-1]
    P = m[-1]
    f = objfun(K, P, R,x,T,Tc,Pc,omega,kij) 
    return f
    
def objfun(K, P, R,x,T,Tc,Pc,omega,kij):
    """
    Objective function for Bubble Point Pressure Newton Solver.
    """
    
    y = K*x
    dG_l,Zc_l,rho_l,lf_l =SRK(R,P,T,x,Pc,Tc,omega,kij)
    dG_v,Zc_v,rho_v,lf_v =SRK(R,P,T,y,Pc,Tc,omega,kij)

    f = np.zeros(np.size(K) + 1)
    f[:-1] = K - np.exp(lf_l)/np.exp(lf_v)
    f[-1] = np.sum(y)
    return f

def jacobfun(f0, x, fun, *data):
    """
    Jacobian matrix creator for newton solver. x should be an array
    with dimensions n by 1.
    """
    J = np.zeros([x.size,x.size])
    eps = 1e-6
    for i in np.arange(0,x.size):
        dx = np.zeros(x.size) 
        dx[i] = eps
#        print('====================================================')
        #f0 = fun(x-dx, *data)
        f1 = fun(x+dx, *data)
#        print x
#        print dx
#        print x+dx
#        print f0
#        print f1
#        print('====================================================')
        J[:,i] = (f1 - f0) / (eps)
    return J                 

def updatekss(K,x,P,T,Pc,Tc,omega,kij):
    y = K*x
    dG_l,Zc_l,rho_l,lf_l =SRK(R,P,T,x,Pc,Tc,omega,kij)
    dG_v,Zc_v,rho_v,lf_v =SRK(R,P,T,y,Pc,Tc,omega,kij)    
    K = np.exp(lf_l)/np.exp(lf_v)
    return K


def RRNo14(K,z):
    VFmin  = 1.0 / (1.0 - K.max())
    VFmax  = 1.0 / (1.0 - K.min())
    RRZero = RR_Vapor(z,0.0,K)
    RROne  = RR_Vapor(z,1.0,K)
    RRHalf = RR_Vapor(z,0.5,K)
    print('===========================================')
    print('RR(0): ',RRZero, ' RR(1): ',RROne)
    print('===========================================')    
################################################################################
# STEP 3) SOLVE RACHFORD RICE SPLIT EQUATION
################################################################################
    if VFmin <  0.0 and VFmax > 0.0:
        print('Rachford Rice Solution Exists')
        data = z, K
        VF = op.fsolve(RR_Vapor,VFmin+VFmax/2,args=data)
    else:
        print('No solution recovered,returning a default value')
    return VF
        
def Pb_multi_Newton(K,x,P,T,Pc,Tc,omega,kij,eptol,maxiter):
    K = updatekss(K,x,P,T,Pc,Tc,omega,kij)

    # Normalize Result

    # Define initial conditions
    x_opt = np.zeros(K.size + 1)
    x_opt[:-1] = K
    x_opt[-1] = P
    for i in np.arange(0,maxiter):
        data = (R,x,T,Tc,Pc,omega,kij)
        # Update using Newton Raphson
        f = gate(x_opt,*data)
        J = jacobfun(f,x_opt,gate,*data)
        dx = np.linalg.solve(-1.0*J,f)
        # Update x
        x_opt = x_opt + dx
        eps = np.abs(np.max(f))
        if eps < eptol:
            print('Converged after ',i, 'iterations.')
            break        
        print('Newton Solver Iteration ',i,' Close')
        print('============================================================')
        print('dx: ',dx)
        print('x: ',x_opt)
        print('Residual: ',f)            
    # Update K and Bubble Point Pressure
    K = x_opt[:-1]
    P = x_opt[-1]        
    y = K*x   
    if i+1 >= maxiter:
        print('Convergence Failure')
        failflag = True
    else:
        failflag = False        
    return P,y,failflag

def PbMultiSS(K,x,P,T,Pc,Tc,omega,kij,eptol,maxiter):
    for i in np.arange(0,maxiter):
        f,K,P = PbSSObj(K,x,P,T,Pc,Tc,omega,kij)
        if np.any(np.abs(f) < eptol):
            break
    if i + 1 >= maxiter:
        print('Convergence Failure')
        failflag = True
    else:
        failflag = False
    y = K*x
    return P,y,failflag
    
def PbSSObj(K,x,P,T,Pc,Tc,omega,kij):
    y = K*x
    dG_l,Zc_l,rho_l,lf_l =SRK(R,P,T,x,Pc,Tc,omega,kij)
    dG_v,Zc_v,rho_v,lf_v =SRK(R,P,T,y,Pc,Tc,omega,kij)
    fug_v = y*np.exp(lf_v)*P    
    fug_l = x*np.exp(lf_l)*P    
    K = fug_l / fug_v
    P = P + fug_l / fug_v
    f = x*K - 1
    return f,K,P
#############################################
# Estimate initial K values
#############################################

# Read in table properties
def EOSBuild(F8Comp):
    ncomps,nprops = F8Comp.shape
    dump = ncomps*[None]
    for i in np.arange(0,ncomps):
        dump[i] = [
                    {
                        "CAS" : str("000-00-0"+str(i)),
                        "Tc"  : F8Comp.Tc[i],
                        "Tc_units": "K",
                        "acentric": F8Comp.omega[i],
                        "aliases": [
                        ],
                        "molemass": F8Comp.MW[i] / 1000,
                        "molemass_units": "kg/mol",
                        "name": F8Comp.comp.as_matrix()[i],
                        "pc": F8Comp.Pc[i]*101325,
                        "pc_units": "Pa"
                    }
                ]            
        CP.add_fluids_as_JSON("SRK",json.dumps(dump[i]))             
    # Define mixture
    mix_string = str(F8Comp.comp.as_matrix()[0])
    for i in np.arange(1,ncomps):
        mix_string = mix_string + '&' + str(F8Comp.comp.as_matrix()[i])
 
    # Add Mixture
    SRK = CP.AbstractState("SRK",mix_string)
    SRK.set_mole_fractions(F8Comp.z.as_matrix())

    # Add BIPs
    for i in np.arange(0,ncomps):
        for j in np.arange(0,ncomps):
           SRK.set_binary_interaction_double(i,j,"kij",F8Comp.kij[i,j])

    return SRK

#SRK.build_phase_envelope("") # SLOOOOWWWW

                
def PropPlot(F8Comp):
    ncomps,nprops = F8Comp.shape
    z = F8Comp.z.as_matrix()
    dump = ncomps*[None]
    for i in np.arange(0,ncomps):
        dump[i] = [
                    {
                        "CAS" : str("000-00-0"+str(i)),
                        "Tc"  : F8Comp.Tc[i],
                        "Tc_units": "K",
                        "acentric": F8Comp.omega[i],
                        "aliases": [
                        ],
                        "molemass": F8Comp.MW[i] / 1000,
                        "molemass_units": "kg/mol",
                        "name": F8Comp.comp.as_matrix()[i],
                        "pc": F8Comp.Pc[i]*101325,
                        "pc_units": "Pa"
                    }
                ]            
        CP.add_fluids_as_JSON("SRK",json.dumps(dump[i]))             
    # Define mixture
    mix_string = str(F8Comp.comp.as_matrix()[0])
    plot_string = 'HEOS::' + str(F8Comp.comp.as_matrix()[0]) + '[' + str(z[0]) + ']'
    for i in np.arange(1,ncomps):
        mix_string = mix_string + '&' + str(F8Comp.comp.as_matrix()[i])
        plot_string = plot_string + '&' + str(F8Comp.comp.as_matrix()[i]) + '[' + str(z[i]) + ']'

    # Add Mixture
    SRK = CP.AbstractState("SRK",mix_string)
    SRK.set_mole_fractions(F8Comp.z.as_matrix())
    
    # Generate Plot Object
    plot = PropertyPlot(plot_string,'PT',unit_system='EUR')
    return SRK,plot


def VLEProcess(P,T,F8Comp):
    """
    INPUT PRESSURE IS IN ATM!!!
    """
#    F8 = EOSBuild(F8Comp) 
    P       = P*101325
    ncomps  = F8Comp.z.size
    z       = F8Comp.z.as_matrix()
    Pc      = F8Comp.Pc.as_matrix()*101325
    Tc      = F8Comp.Tc.as_matrix()
    omega   = F8Comp.omega.as_matrix()
    MW      = F8Comp.MW.as_matrix() / 1000.
    KIJ     = F8Comp.kij

    # SRK Coefficients
    om_a = 0.42748 * np.ones(ncomps) # [-]
    om_b = 0.08664 * np.ones(ncomps) # [-]

    #(P_sat*101325, T_res,z,Pc,Tc,omega,MW,om_a,om_b,KIJ) = props

    # Estimate initial K-values
    initial_K_values = calculate_K_values_wilson(
        P,T,Pc,Tc,omega )
        
    # Create EoS object and set properties

    eos = SoaveRedlichKwongEos(Pc,Tc,omega, om_a, om_b, KIJ)

    is_stable, K_values_est = calculate_stability_test(eos,P,T,z,initial_K_values)

#    print ('System is stable?', is_stable)
#    print ('K_values estimates:', K_values_est)


    K_values_from_ss_flash, F_V, f_L = ss_flash(eos,P,T,z,K_values_est, tolerance = 1.0e-3)
#    print ('K_values Successive Subst:', K_values_from_ss_flash)
#    print ('Vapor molar fraction:', F_V)
#    print ('\n-----\nFugacities obtained:', f_L)
    return F_V    

def FullVLEProcess(P,T,F8Comp):
    """
    INPUT PRESSURE IS IN ATM!!
    """
#    F8 = EOSBuild(F8Comp) 
    
    ncomps  = F8Comp.z.size
    z       = F8Comp.z.as_matrix()
    Pc      = F8Comp.Pc.as_matrix()*101325
    Tc      = F8Comp.Tc.as_matrix()
    omega   = F8Comp.omega.as_matrix()
    MW      = F8Comp.MW.as_matrix() / 1000.
    KIJ     = F8Comp.kij
    P       = P*101325
    # SRK Coefficients
    om_a = 0.42748 * np.ones(ncomps) # [-]
    om_b = 0.08664 * np.ones(ncomps) # [-]

    #(P_sat*101325, T_res,z,Pc,Tc,omega,MW,om_a,om_b,KIJ) = props

    # Estimate initial K-values
    initial_K_values = calculate_K_values_wilson(
        P,T,Pc,Tc,omega )
        
    # Create EoS object and set properties

    eos = SoaveRedlichKwongEos(Pc,Tc,omega, om_a, om_b, KIJ)

    is_stable, K_values_est = calculate_stability_test(eos,P,T,z,initial_K_values)

#    print ('System is stable?', is_stable)
#    print ('K_values estimates:', K_values_est)


    K_values_from_ss_flash, F_V, f_L = ss_flash(eos,P,T,z,K_values_est, tolerance = 1.0e-3)
#    print ('K_values Successive Subst:', K_values_from_ss_flash)
#    print ('Vapor molar fraction:', F_V)
#    print ('\n-----\nFugacities obtained:', f_L)
    return F_V,K_values_from_ss_flash,f_L


def FullVLETuning(P,T,F8Comp):
    """
    INPUT PRESSURE IS IN ATM!!
    """
#    F8 = EOSBuild(F8Comp) 
    
    ncomps  = F8Comp.z.size
    z       = F8Comp.z.as_matrix()
    Pc      = F8Comp.Pc.as_matrix()*101325
    Tc      = F8Comp.Tc.as_matrix()
    omega   = F8Comp.omega.as_matrix()
    MW      = F8Comp.MW.as_matrix() / 1000.
    KIJ     = F8Comp.kij
    P       = P*101325
    # SRK Coefficients
    om_a    = F8Comp.om_a.as_matrix()
    om_b    = F8Comp.om_b.as_matrix()
#    om_a = 0.42748 * np.ones(ncomps) # [-]
#   om_b = 0.08664 * np.ones(ncomps) # [-]

    #(P_sat*101325, T_res,z,Pc,Tc,omega,MW,om_a,om_b,KIJ) = props

    # Estimate initial K-values
    initial_K_values = calculate_K_values_wilson(
        P,T,Pc,Tc,omega )
        
    # Create EoS object and set properties

    eos = SoaveRedlichKwongEos(Pc,Tc,omega, om_a, om_b, KIJ)

    is_stable, K_values_est = calculate_stability_test(eos,P,T,z,initial_K_values)

#    print ('System is stable?', is_stable)
#    print ('K_values estimates:', K_values_est)


    K_values_from_ss_flash, F_V, f_L = ss_flash(eos,P,T,z,K_values_est, tolerance = 1.0e-3)
#    print ('K_values Successive Subst:', K_values_from_ss_flash)
#    print ('Vapor molar fraction:', F_V)
#    print ('\n-----\nFugacities obtained:', f_L)
    return F_V,K_values_from_ss_flash,f_L
                

                
def EnvGrid(F8Comp,Pmin,Pmax,Tmin,Tmax,nsamps):
    Pgrid = np.linspace(Pmin,Pmax,nsamps)
    Tgrid = np.linspace(Tmin,Tmax,nsamps)

    VFGrid = np.zeros([nsamps,nsamps])

    for i in np.arange(0,nsamps):
        for j in np.arange(0,nsamps):
            try:
                VFGrid[i,j] = VLEProcess(Pgrid[i],Tgrid[j],F8Comp)
            except:
                VFGrid[i,j] = -999
    plotGrid = np.clip(VFGrid,0.0,1.0)            
    fig = plt.figure()
    im = plt.imshow(plotGrid,interpolation='bilinear',origin='lower')       
    
    plt.tight_layout()
    fig.savefig('VF_Envelope.png')   
    return VFGrid            
            
def PhaseLine(P,T,VF,Comp):
    out = VLEProcess(P,T,Comp)
    return out - VF                    

def PhaseOpt(f, *params):
    (VF,Comp) = params
    try:
        out = VLEProcess(f[0],f[1],Comp)
    except:
        out = 999
    return np.clip((out-VF),-0.1,1.1) 

def PhaseOptFixedP(T, *params):
    (VF,P,Comp) = params
    try:
        out = VLEProcess(P,T,Comp)
    except:
        out = 999
    return np.clip((out-VF),-0.1,1.1) 

def PhaseOptFixedT(P, *params):
    (VF,T,Comp) = params
    try:
        out = VLEProcess(P,T,Comp)
    except:
        out = 999
    return np.clip((out-VF),-0.1,1.1) 

def FindTforVF(VF,P,Comp):
    data = (VF,P,Comp)
    T = op.fsolve(PhaseOptFixedP,500,args=data)
    return T

def FindPforVF(VF,T,Comp):
    data = (VF,T,Comp)
    P = op.fsolve(PhaseOptFixedT,100,args=data)
    return P   

################################################################################
# Baseline (Given) Fluid Model
#-------------------------------------------------------------------------------
# EOS Object


# Reservoir Temperature CME
CME_F8_rho = np.zeros(CME_T_res.rho.size)

# Conditions for initial estimate
T = 350.        # K
P = 50.*101325  # Pa





nsamps = 101

Pmax = 200
Pmin = 1

Tmax = 1000
Tmin = 1


#grid = EnvGrid(F8Reduced,Pmin,Pmax,Tmin,Tmax,nsamps)
def BruteForce(F8Reduced,Pmin,Pmax,Pstep,Tmin,Tmax,Tstep,VF):
    """
    Function to brute force collect fluid PT / Vapor Fraction data.
    """
    prange = slice(Pmin,Pmax,Pstep)
    trange = slice(Tmin,Tmax,Tstep)
    rrange = (prange,trange)
    rb_x0,rb_fval,rb_grid,rb_Jout = op.brute(PhaseOpt,rrange,args=(VF,F8Reduced),full_output=True,finish=op.fmin)
    resbrute = (rb_x0,rb_fval,rb_grid,rb_Jout)
    return resbrute
    
# Actual Phase Envelope Points
#PE_upper = np.flatnonzero( (-0.02 < rb_Jout) & (rb_Jout < 0.02) == True)
#PE_lower = np.flatnonzero( (1-0.02 < rb_Jout) & (rb_Jout < 1+0.02) == True)
#rescon_prange = np.linspace(501,11,4901)*0.986923


R = 0.082057338 # m**3 * atm * K**-1 * mol**-1
T = 394.2 # K
n = 1.0 # # mols in

ncomps = F8Reduced.Tc.size
om_a = 0.42748 * np.ones(ncomps) # [-]
om_b = 0.08664 * np.ones(ncomps) # [-]
sat_index = 8
run_num = 0

def rescon(Tc,Pc,z,omega,ncomps,MW,kij,om_a,om_b,CME_T_res,DDE_T_res,tol,R,T,n,run_num):
    K = calculate_K_values_wilson(P,T,Pc,Tc,omega)
    eos = SoaveRedlichKwongEos(Pc,Tc,omega, om_a, om_b, kij)
    rescon_prange = np.unique(np.concatenate( (CME_T_res.P.as_matrix(),DDE_T_res.P.as_matrix()) ))
    rescon_prange[::-1].sort()
    rescon_steps = rescon_prange.size



    ################################################################################
    # CME/DDE at Reservoir Temperature
    ################################################################################

    # CME at 394.2K
    CME_res_dG_z  = np.zeros(rescon_steps)
    CME_res_dG_y  = np.zeros(rescon_steps)
    CME_res_dG_x  = np.zeros(rescon_steps)
    CME_res_Zc_z  = np.zeros(rescon_steps)
    CME_res_Zc_y  = np.zeros(rescon_steps)
    CME_res_Zc_x  = np.zeros(rescon_steps)
    CME_res_rho_z  = np.zeros(rescon_steps)
    CME_res_rho_y  = np.zeros(rescon_steps)
    CME_res_rho_x  = np.zeros(rescon_steps)
    CME_res_rho_mix  = np.zeros(rescon_steps)
    CME_res_lf_z  = np.zeros([rescon_steps,z.size])
    CME_res_lf_y  = np.zeros([rescon_steps,z.size])
    CME_res_lf_x  = np.zeros([rescon_steps,z.size])

    CME_res_VF  = np.zeros(rescon_steps)
    CME_res_lf  = np.zeros([rescon_steps,z.size])
    CME_res_K  = np.zeros([rescon_steps,z.size])
    CME_res_x  = np.zeros([rescon_steps,z.size])
    CME_res_y  = np.zeros([rescon_steps,z.size])
    CME_res_fl  = np.zeros([rescon_steps,z.size])

    CME_res_Vl  = np.zeros(rescon_steps)
    CME_res_Vg  = np.zeros(rescon_steps)
    CME_res_V  = np.zeros(rescon_steps)
    CME_res_Vsum  = np.zeros(rescon_steps)

    CME_res_AMWl  = np.zeros(rescon_steps)
    CME_res_AMWg  = np.zeros(rescon_steps)

    CME_res_rho_g_eos  = np.zeros(rescon_steps)
    CME_res_rho_l_eos  = np.zeros(rescon_steps)
    CME_res_rho_eos  = np.zeros(rescon_steps)

    AMW = np.sum(MW * z)

    for i in np.arange(0,rescon_steps):
        # This function input requires pressure in atm
        CME_res_dG_z[i],CME_res_Zc_z[i],CME_res_rho_z[i],CME_res_lf_z[i,:] = SRK(R,rescon_prange[i],T,z,Pc,Tc,omega,kij)
        CME_res_VF[i],CME_res_K[i,:],CME_res_fl[i,:] = FullVLEProcess(rescon_prange[i],T,F8Reduced)
        
        
        
        CME_res_VF = np.clip(CME_res_VF,-10,10) 
        CME_res_x[i,:] = z / (1.0 + CME_res_VF[i]*(CME_res_K[i,:] - 1.0) )
        CME_res_y[i,:] = ( CME_res_K[i,:] * z ) / (1.0 + CME_res_VF[i]*(CME_res_K[i,:] - 1.0) )
        # Normalize AGAIN
        CME_res_y[i,:] = CME_res_y[i,:] / np.sum( CME_res_y[i,:] )
        CME_res_x[i,:] = CME_res_x[i,:] / np.sum( CME_res_x[i,:] )
        CME_res_dG_y[i],CME_res_Zc_y[i],CME_res_rho_y[i],CME_res_lf_y[i,:] = SRK(R,rescon_prange[i],T,CME_res_y[i,:],Pc,Tc,omega,kij)    
        CME_res_dG_x[i],CME_res_Zc_x[i],CME_res_rho_x[i],CME_res_lf_x[i,:] = SRK(R,rescon_prange[i],T,CME_res_x[i,:],Pc,Tc,omega,kij)
      
        # Solve for liquid and vapor volumes
        CME_res_V[i]  = (CME_res_Zc_z[i] * n * R * T) / rescon_prange[i]
        CME_res_Vl[i] = (CME_res_Zc_x[i] * n * (1-CME_res_VF[i]) * R * T) / rescon_prange[i]
        CME_res_Vg[i] = (CME_res_Zc_y[i] * n * CME_res_VF[i] * R * T) / rescon_prange[i]
        CME_res_Vsum[i] = CME_res_Vl[i] + CME_res_Vg[i]
        
        CME_res_AMWg[i] = np.sum(MW * CME_res_y[i,:])# * CME_res_VF[i])
        CME_res_AMWl[i] = np.sum(MW * CME_res_x[i,:])# * ( 1 - CME_res_VF[i]) )
        
        # Density from my EOS solver
        CME_res_rho_z[i] = CME_res_rho_z[i] * (CME_res_AMWg[i] * CME_res_VF[i] + CME_res_AMWl[i]*(1-CME_res_VF[i]))
        CME_res_rho_y[i] = CME_res_rho_y[i]*CME_res_AMWg[i]    
        CME_res_rho_x[i] = CME_res_rho_x[i]*CME_res_AMWl[i]
        CME_res_rho_mix[i]  = CME_res_rho_y[i]*CME_res_VF[i] + CME_res_rho_x[i] * (1-CME_res_VF[i])


        
        CME_res_rho_g_eos[i]  = np.nan_to_num(CME_res_AMWg[i] / CME_res_Vg[i])
        CME_res_rho_l_eos[i]  = np.nan_to_num(CME_res_AMWl[i] / CME_res_Vl[i])
        CME_res_rho_eos[i]    = CME_res_rho_g_eos[i]*CME_res_VF[i] + CME_res_rho_l_eos[i] * (1-CME_res_VF[i])
        #(CME_res_AMWg[i]*CME_res_VF[i] + CME_res_AMWl[i]*(1-CME_res_VF[i])) / CME_res_Vsum[i]

    # Averaged / Sanitized Values
        CME_V = (np.clip(CME_res_Vg,0,1e6) + CME_res_V)
        

    ################################################################################
    # Differential Depletion Test
    ################################################################################

    #for i in np.arange(0,DDE_T_res.P.size):

    ################################################################################
    # Plotting
    ################################################################################

    # CME @ Reservoir Conditions

    CMEindex = np.flatnonzero(np.in1d(rescon_prange,CME_T_res.P.as_matrix()))
    DDEindex = np.flatnonzero(np.in1d(rescon_prange,DDE_T_res.P.as_matrix()))

    VrelCME = CME_res_Vsum / CME_res_Vsum[8]
    VrelCME = VrelCME[CMEindex]
    #===============================================================================
    # CME Plot
    fig, ax_CME = plt.subplots(figsize=(10, 6))
    plt.scatter(CME_T_res.P,CME_T_res.vrel)
    #plt.plot(DDE_T_res.P,DDE_T_res.vrel,'y')
    #plt.scatter(CME_T_res.P,CME_res_V / CME_res_V[sat_index])
    #plt.scatter(CME_T_res.P,(np.clip(CME_res_Vg,0,1e6) + CME_res_V)/ CME_res_V[sat_index],c='r' )
    plt.plot(rescon_prange[CMEindex],VrelCME,'r')
    #plt.scatter(rescon_prange,CME_res_Vl / CME_res_Vsum[8],marker='*' )
    #plt.show()
    plt.tight_layout()

    ax_CME.set_ylabel('Relative Volume, m**3/m**3')
    ax_CME.set_xlabel('Pressure, atm')
    #ax_Pc.set_xlabel('Carbon Number')
    #ax_Pc.set_ylabel('Critical Pressure, atm')
    #plt.savefig('1_Pc_Tc.png',frameon=True,bbox_inches='tight',padinches=0.1)
    cmetitle = 'CME at Reservoir Temperature, Iteration # ' + str(run_num)
    ax_CME.set_title(cmetitle)
    cmefilename = 'CME_at_Tres_'+ str(run_num) + '.png'
    fig.savefig(cmefilename,frameon=True,bbox_inches='tight',padinches=1.0)
    plt.clf()

    #===============================================================================
    # CME Cell Density Plot
    fig, ax_CME_rho = plt.subplots(figsize=(10, 6))
    plt.scatter(CME_T_res.P,CME_T_res.rho)
    #plt.plot(DDE_T_res.P,DDE_T_res.vrel,'y')
    #plt.scatter(CME_T_res.P,CME_res_V / CME_res_V[sat_index])
    #plt.scatter(CME_T_res.P,(np.clip(CME_res_Vg,0,1e6) + CME_res_V)/ CME_res_V[sat_index],c='r' )
    plt.plot(rescon_prange[CMEindex],CME_res_rho_mix[CMEindex],'r')
    #plt.scatter(rescon_prange,CME_res_Vl / CME_res_Vsum[8],marker='*' )
    #plt.show()
    plt.tight_layout()

    ax_CME_rho.set_ylabel('Cell Density, kg/m**3')
    ax_CME_rho.set_xlabel('Pressure, atm')
    #ax_Pc.set_xlabel('Carbon Number')
    #ax_Pc.set_ylabel('Critical Pressure, atm')
    #plt.savefig('1_Pc_Tc.png',frameon=True,bbox_inches='tight',padinches=0.1)
    cmerhotitle = 'CME Cell Density at Reservoir Temperature, Iteration # ' + str(run_num)
    ax_CME_rho.set_title(cmerhotitle)
    cmerhofilename = 'CME_rho_at_Tres_'+ str(run_num) + '.png'
    fig.savefig(cmerhofilename,frameon=True,bbox_inches='tight',padinches=1.0)
    plt.clf()

    #===============================================================================
    # DDE Plot
    fig, ax_DDE = plt.subplots(figsize=(10, 6))
    #plt.plot(CME_T_res.P,CME_T_res.vrel,'r')
    plt.scatter(DDE_T_res.P,DDE_T_res.vrel)
    plt.plot(rescon_prange[DDEindex],CME_res_Vl[DDEindex] / CME_res_Vsum[sat_index],'r')
    plt.tight_layout()
    ax_DDE.set_ylabel('Relative Volume, m**3/m**3')
    ax_DDE.set_xlabel('Pressure, atm')
    ddetitle = 'DDE at Reservoir Temperature, Iteration # ' + str(run_num)
    ax_DDE.set_title(ddetitle)
    #plt.show()
    ddefilename = 'DDE_at_Tres_'+ str(run_num) + '.png'
    fig.savefig(ddefilename,frameon=True,bbox_inches='tight',padinches=1.0)
    plt.clf()

    #===============================================================================
    # DDE Liquid Density Plot
    fig, ax_DDE_rho = plt.subplots(figsize=(10, 6))
    plt.scatter(DDE_T_res.P,DDE_T_res.rho)
    #plt.plot(DDE_T_res.P,DDE_T_res.vrel,'y')
    #plt.scatter(CME_T_res.P,CME_res_V / CME_res_V[sat_index])
    #plt.scatter(CME_T_res.P,(np.clip(CME_res_Vg,0,1e6) + CME_res_V)/ CME_res_V[sat_index],c='r' )
    plt.plot(rescon_prange[DDEindex],CME_res_rho_x[DDEindex],'r')
    #plt.scatter(rescon_prange,CME_res_Vl / CME_res_Vsum[8],marker='*' )
    #plt.show()
    plt.tight_layout()

    ax_DDE_rho.set_ylabel('Liquid Density, kg/m**3')
    ax_DDE_rho.set_xlabel('Pressure, atm')
    #ax_Pc.set_xlabel('Carbon Number')
    #ax_Pc.set_ylabel('Critical Pressure, atm')
    #plt.savefig('1_Pc_Tc.png',frameon=True,bbox_inches='tight',padinches=0.1)
    dderhotitle = 'DDE Liquid Density at Reservoir Temperature, Iteration # ' + str(run_num)
    ax_DDE_rho.set_title(cmerhotitle)
    dderhofilename = 'DDE_rho_at_Tres_'+ str(run_num) + '.png'
    fig.savefig(dderhofilename,frameon=True,bbox_inches='tight',padinches=1.0)
    plt.clf()
    #-------------------------------------------------------------------------------
    CME_Vrel = CME_res_Vl / CME_res_Vsum[sat_index]
    return CMEindex, DDEindex, rescon_prange, CME_Vrel, CME_res_rho_x, CME_res_rho_mix


################################################################################
# Attempt at optimizing this mess
# The following variables are game:
#
# Pc, Tc, omega, om_a, om_b
# DO NOT TOUCH kij!!
tol = 1e-6



Tc = F8Reduced.Tc.as_matrix()
Pc = F8Reduced.Pc.as_matrix()
z = F8Reduced.z.as_matrix()
omega = F8Reduced.omega.as_matrix()
ncomps = F8Reduced.Tc.size
MW = F8Reduced.MW.as_matrix()
kij = F8Reduced.kij

numPC = 4
input_matrix = np.zeros([z.size,5])
input_matrix[:,0] = Pc
input_matrix[:,1] = Tc
input_matrix[:,2] = omega
input_matrix[:,3] = om_a
input_matrix[:,4] = om_b

input_matrix = input_matrix.T.flatten()


def rescon_opt(input_matrix, *data):

    (CME_T_res,DDE_T_res,z,ncomps,MW,kij,tol,R,T,n,run_num) = data

    Pc = input_matrix[z.size*0:z.size*1]
    Tc = input_matrix[z.size*1:z.size*2]
    omega = input_matrix[z.size*2:z.size*3]
    om_a = input_matrix[z.size*3:z.size*4]
    om_b = input_matrix[z.size*4:z.size*5]

    CME_check_matrix = np.zeros([CME_T_res.P.size,2])
    CME_check_matrix[:,0] = CME_T_res.rho.as_matrix()
    CME_check_matrix[:,1] = CME_T_res.vrel.as_matrix()

    DDE_check_matrix = np.zeros([DDE_T_res.P.size,2])
    DDE_check_matrix[:,0] = DDE_T_res.rho.as_matrix()
    DDE_check_matrix[:,1] = DDE_T_res.vrel.as_matrix()

    output_CME_matrix =  np.zeros([CME_T_res.P.size,2])
    output_DDE_matrix = np.zeros([DDE_T_res.P.size,2])

    CMEindex, DDEindex, rescon_prange, CME_Vrel, CME_res_rho_x, CME_res_rho_mix = rescon(Tc,Pc,z,omega,ncomps,MW,kij,om_a,om_b,CME_T_res,DDE_T_res,tol,R,T,n,run_num)

    output_CME_matrix[:,0] = CME_res_rho_mix[CMEindex]
    output_CME_matrix[:,1] = CME_Vrel[CMEindex]

    output_DDE_matrix[:,0] = CME_res_rho_x[DDEindex]
    output_DDE_matrix[:,1] = CME_Vrel[DDEindex]

    error_list = np.zeros(4)
    error_list[0] = np.abs(np.sum(CME_check_matrix[:,0] - output_CME_matrix[:,0]))
    error_list[1] = np.abs(np.sum(CME_check_matrix[:,1] - output_CME_matrix[:,1]))
    error_list[2] = np.abs(np.sum(DDE_check_matrix[:,0] - output_DDE_matrix[:,0]))
    error_list[3] = np.abs(np.sum(DDE_check_matrix[:,1] - output_DDE_matrix[:,1]))



    return error_list #np.abs(np.sum(error_list)), CME_check_matrix,DDE_check_matrix, output_CME_matrix,output_DDE_matrix

run_num = 1
data = (CME_T_res,DDE_T_res,z,ncomps,MW,kij,tol,R,T,n,run_num)

# Close remaining plots left open
plt.close('all')
#def minattempt(input_matrix,data,numiter,run_num):
    
