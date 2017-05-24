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
from matplotlib.pyplot import cm
np.set_printoptions(suppress=True)

# Verification Module
import thermo as th
import CoolProp.CoolProp as CP
import json
from CoolProp.Plots import PropertyPlot
##############################################################################80
# Plot Parameters
##############################################################################80

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
F8Comp.Psat = 145.8*98692327 # atm
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

kij_part = np.array ([[0,  	0,  	0.13,	0.025,	0.01,	0.09,	0.095,	0.095,	0.1,	0.11],
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


z_in        =   F8Comp['z'].as_matrix()
MW_in       =   F8Comp['MW'].as_matrix()
rho_in      =   F8Comp['rho'].as_matrix()
y_in        =   F8Comp['y'].as_matrix()
Pc_in       =   Ped_Pc(coeff,rho_in,MW_in)
Tc_in       =   Ped_Tc(coeff,rho_in,MW_in)
omega_in    =   Ped_omega(coeff,rho_in,MW_in)

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
#fig, ax_MF  = plt.subplots(figsize=(10, 6))
#ax_rho       = ax_MF.twinx()
#ax_MF.plot(C_num,z,'r--')   
#ax_rho.plot(C_num,rho,'b')
#plt.tight_layout()

#ax_MF.set_xlim(7,SCN[-1])
#ax_MF.set_ylabel('Mole Fraction')
#ax_MF.set_xlabel('Carbon Number')
#ax_rho.set_ylabel('Density, g/cc')
#plt.savefig('1_rho_MF.png',frameon=True,bbox_inches='tight',padinches=0.1)
#plt.clf()

## Pc, Tc, cmin to cmax
#fig, ax_Pc  = plt.subplots(figsize=(10, 5))
#ax_Tc       = ax_Pc.twinx()
#ax_Tc.plot(C_num,Tc,'r--')   
#ax_Pc.plot(C_num,Pc,'b')
#ax_Pc.set_xlim(7,SCN[-1])
#plt.tight_layout()

#ax_Tc.set_ylabel('Critical Temperature, Deg. Kelvin')
#ax_Pc.set_xlabel('Carbon Number')
#ax_Pc.set_ylabel('Critical Pressure, atm')
#plt.savefig('1_Pc_Tc.png',frameon=True,bbox_inches='tight',padinches=0.1)
#plt.clf()

## Mass Frac, % of total cmin to cmax
#fig, ax_MW  = plt.subplots(figsize=(10, 6))
#ax_sum      = ax_MW.twinx()
#ax_MW.plot(C_num,MassFrac,'r--')   
#ax_sum.plot(C_num,cummass,'b')
#plt.tight_layout()

#ax_sum.set_ylabel('Cumulative Mass, %')
#ax_sum.set_ylim(0,1)
#ax_MW.set_xlim(7,SCN[-1])
#ax_MW.set_xlabel('Carbon Number')
#ax_MW.set_ylabel('Mass Fraction of Plus')
#plt.savefig('1_MW_fracsum.png',frameon=True,bbox_inches='tight',padinches=0.1)
#plt.clf()

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
x_PS        =   np.ones(ncomps_PS)
y_PS        =   np.zeros(ncomps_PS)
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
x_PS = x_PS - y_PS


MW_in[0:11] = MW_PS[0:11]
Pc_in[0:11] = Pc_PS[0:11]
Tc_in[0:11] = Tc_PS[0:11]
omega_in[0:11] = omega_PS[0:11]
y_in[0:11] = y_PS[0:11]
x_in = 1.0 - y_in


# MF for above 11

z_PS[:11] = z_in[:11]
#===============================================================================
# Update Pandas dataframe for input case for use as a baseline
F8Comp['x'] = x_in
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

x_exp = 1.0 - y_exp
F8Exp = pd.DataFrame({'comp':comps_exp,'z':z_exp,'MW':MW_exp,'rho':rho_exp,'y':y_exp,'x':x_exp,
                      'Pc':Pc_exp,'Tc':Tc_exp,'omega':omega_exp})
#===============================================================================
# DataFrame for reduced set

comps_PS = comps_exp[:ncomps_PS]
comps_PS[-4] = 'CL1'
comps_PS[-3] = 'CL2'
comps_PS[-2] = 'CL3'
comps_PS[-1] = 'CL4'
F8Reduced = pd.DataFrame({'comp':comps_PS,'z':z_PS,'MW':MW_PS,'y':y_PS,'x':x_PS,
                      'Pc':Pc_PS,'Tc':Tc_PS,'omega':omega_PS})
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
eos = SRKMIX(T=T,P=P_range[0],Tcs=Tc,Pcs=Pc*101325.,omegas=omega,zs=z,kijs=kij)
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
    plot_string = 'SRK::' + str(F8Comp.comp.as_matrix()[0]) + '[' + str(z[0]) + ']'
    for i in np.arange(1,ncomps):
        mix_string = mix_string + '&' + str(F8Comp.comp.as_matrix()[i])
        plot_string = plot_string + '&' + str(F8Comp.comp.as_matrix()[i]) + '[' + str(z[i]) + ']'

    # Add Mixture
    SRK = CP.AbstractState("SRK",mix_string)
    SRK.set_mole_fractions(F8Comp.z.as_matrix())
    
    # Generate Plot Object
    #plot = PropertyPlot(plot_string,'PT',unit_system='EUR')
    return SRK,plot_string
                
                
                
################################################################################
# Baseline (Given) Fluid Model
#-------------------------------------------------------------------------------
# EOS Object

F8 = EOSBuild(F8Comp) 
# Reservoir Temperature CME
CME_F8_rho = np.zeros(CME_T_res.rho.size)

## Build Phase Envelope/Find Critical Points - ONLY RUN ONCE - TAKES FOREVER
#F8.build_phase_envelope("")
#PE = F8.get_phase_envelope_data()
#pts = F8.all_critical_points()


for i in np.arange(0,CME_F8_rho.size):
    F8.update(CP.PT_INPUTS,101325.*CME_T_res.P[i], F8Comp.Tres)
    CME_F8_rho[i] = F8.rhomass()

        






