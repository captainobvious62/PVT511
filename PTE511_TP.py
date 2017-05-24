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
np.set_printoptions(suppress=True)

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
        AIJ = (1 - kij) * (A * A.T)**0.5
        AM  = np.sum( ( x.reshape([x.size,1]) * x ) *               \
              np.sqrt( A.reshape( [A.size,1] ) * A) *                   \
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
        
    print('z: ',z)
    print('gstar: ',g_star)    
    print('====================================================')
    print('----------------------------------------------------')
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

def RR_Vapor_Calc(VF,z,K,eptol,numiter):
    eps = 1
    for i in np.arange(0,numiter):
        f = RR_Vapor(VF,z,K) 
        df = RRPrime_Vapor(VF,z,K)
        VF = VF - f/df
        # Stabilize
        print i,VF,f      
        if f**2 < eptol:
            break            
        if VF > 1:
            VF = 1
            break
        elif VF < 0:
            VF = 0
            break

    x = z / (VF * (K-1) + 1)
    y = K*x        
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
            VF,x,y = RR_Vapor_Calc(VF,z,K,1e-6,maxiter)
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
compstrs = ["CH4","nC4","nC10"]
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
# Bubble Point Estimates
################################################################################

################################################################################
# Initial estimates from ideal case

def EnvelopeInitialT(z,T,P,Pc,Tc,omega,eptol):
    fT = 10
    while fT**2 > eptol**2:
        K   = WilsonEqn(Pc,Tc,omega,T,P)
        fT  = np.sum(K*z) - 1
        dfT =( np.sum(WilsonEqn(Pc,Tc,omega,T+eptol,P)*z - 1)  -    \
               np.sum(WilsonEqn(Pc,Tc,omega,T-eptol,P)*z - 1) ) / (2*eptol)
        T = T - fT/dfT
    return T


#T = EnvelopeInitialT(z,T,P,Pc_PS,Tc_PS,omega_PS,eptol)

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
        print P, f
    return P

def dlfdp(R,P,T,x,Pc,Tc,omega,kij,eptol):
    Zc_x_p,rho_x_p,lf_x_p = SRK(R,P+eptol,T,x,Pc,Tc,omega,kij)
    Zc_x_m,rho_x_m,lf_x_m = SRK(R,P-eptol,T,x,Pc,Tc,omega,kij)
    dlfdp_x = (lf_x_p - lf_x_m) / 2*eptol
    return dlfdp_x    
   


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

T = 175 + 273.15 # K
P = 400 # atm

#Tc = SingleCompSet['Tc'].as_matrix()
#Pc = SingleCompSet['Pc'].as_matrix()
#z  = SingleCompSet['z'].as_matrix()
#omega = SingleCompSet['omega'].as_matrix()
#kij = SingleCompSet.kij
#Tc = SimpleSet['Tc'].as_matrix()
#Pc = SimpleSet['Pc'].as_matrix()
#z  = SimpleSet['z'].as_matrix()
#omega = SimpleSet['omega'].as_matrix()
#kij = SimpleSet.kij

Tc = F8Comp['Tc'].as_matrix()
Pc = F8Comp['Pc'].as_matrix()
z  = F8Comp['z'].as_matrix()
omega = F8Comp['omega'].as_matrix()
kij = F8Comp.kij


#GDEM Initalization

DU = np.zeros([KGDEM_iter,z.size])


# Initial Conditions
################################################################################
# STEP 1) ESTIMATE K VALUES
################################################################################
K = WilsonEqn(Pc,Tc,omega,T,P)

# Flash Calculation
j = 0
VF = SimpleSet.VF

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
        VF,x,y = RR_Vapor_Calc(VF,z,K,1e-6,maxiter)
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

    if efug_err < errtol:
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

################################################################################
# Lazy attempt at bubble point pressure through stability testing
################################################################################
T = 300 # K
data = (R,T,z,Pc,Tc,omega,kij)
Pb = op.bisect

    
         
#   if np.mod(j+1,stability_interval) == 0:
#    if Ktriv_err < 1e-4:


#VF,x,y = RR_Vapor_Calc(z,0.0,K,1e-6,100)

## Fugacities
#dG_x,Zc_x,rho_x,lf_l = SRK(R, P, T, x, Pc, Tc, omega, kij)
#dG_y,Zc_y,rho_y,lf_v = SRK(R, P, T, y, Pc, Tc, omega, kij)
#dG_z,Zc_z,rho_z,lf_z = SRK(R, P, T, z, Pc, Tc, omega, kij)
#fug_l = np.exp(lf_l) * x * P
#fug_v = np.exp(lf_v) * y * P
#fug_z = np.exp(lf_z) * z * P

## Decrease in reduced Gibbs Energy
#g_r = VF*np.sum(y*(np.log(y)+lf_v)) + (1-VF)*np.sum(x*(np.log(x)+lf_l)) - np.sum(z*(np.log(z) + lf_z))

## Tangent Plane Distance
#tpd_x = np.sum(x*(np.log(x)+lf_l - np.log(z) - lf_z))
#tpd_y = np.sum(y*(np.log(y)+lf_v - np.log(z) - lf_z))
## If g_r is negative, two phases are at equilibrium 
## enforcing a reduction of g_r in all subsequent iteration steps
## eliminates the possibility of arriving at a trivial solution
#g_r_tpd = VF*tpd_x + (1-VF)*tpd_y
#if g_r_tpd < 0:
#    K = (lf_l*x*P)*(lf_v*y*P)*K
#    VF,x,y = RR_Vapor_Calc(z,VF,K,1e-6,100)
#    tpd_x = np.sum(x*(np.log(x)+lf_l - np.log(z) - lf_z))
#    tpd_y = np.sum(y*(np.log(y)+lf_v - np.log(z) - lf_z))
#    g_r_tpd = VF*tpd_x + (1-VF)*tpd_y
#    dG_x,Zc_x,rho_x,lf_l = SRK(R, P, T, x, Pc, Tc, omega, kij)
#    dG_y,Zc_y,rho_y,lf_v = SRK(R, P, T, y, Pc, Tc, omega, kij)
#    dG_z,Zc_z,rho_z,lf_z = SRK(R, P, T, z, Pc, Tc, omega, kij)
#    fug_l = np.exp(lf_l) * x * P
#    fug_v = np.exp(lf_v) * y * P
#    fug_z = np.exp(lf_z) * z * P
#    
#    
#    # ITERATION 2 TEST
#    K = (lf_l*x*P)*(lf_v*y*P)*K
#    VF,x,y = RR_Vapor_Calc(z,VF,K,1e-6,100)
#    tpd_x = np.sum(x*(np.log(x)+lf_l - np.log(z) - lf_z))
#    tpd_y = np.sum(y*(np.log(y)+lf_v - np.log(z) - lf_z))
#    g_r_tpd = VF*tpd_x + (1-VF)*tpd_y
#    dG_x,Zc_x,rho_x,lf_l = SRK(R, P, T, x, Pc, Tc, omega, kij)
#    dG_y,Zc_y,rho_y,lf_v = SRK(R, P, T, y, Pc, Tc, omega, kij)
#    dG_z,Zc_z,rho_z,lf_z = SRK(R, P, T, z, Pc, Tc, omega, kij)
#    fug_l = np.exp(lf_l) * x * P
#    fug_v = np.exp(lf_v) * y * P
#    fug_z = np.exp(lf_z) * z * P

################################################################################

# Stability Test
#Wl,x,TPD_l,tpd_l,tm_l = StabilityLiquid2(R,P,T,z,Pc,Tc,omega,kij)
#Wv,y,TPD_v,tpd_v,tm_v = StabilityVapor2(R,P,T,z,Pc,Tc,omega,kij)

## Newton Solver
#V = np.ones(z.size) * VF*z
#g_zero = np.log(fug_v) - np.log(fug_l)
## Timestep one

#for i in np.arange(0,maxiter):
#    print("Iteration ",i)
#    print("==================================")
#    
#    g = np.log(fug_v) - np.log(fug_l)
#    J = np.eye(g.size) * (z / (x*y) - 1) * (1/ (VF * (1-VF) ) )
#    dV = np.linalg.solve(J,-1.*g)
#    V = V + dV
#    L = z - V
#    VF = np.sum(V)
#    y = V / VF
#    x = L / (1-VF)
#    dG_x,Zc_x,rho_x,lf_l = SRK(R, P, T, x, Pc, Tc, omega, kij)
#    dG_y,Zc_y,rho_y,lf_v = SRK(R, P, T, y, Pc, Tc, omega, kij)    
#    fug_l = np.exp(lf_l) * x * P
#    fug_v = np.exp(lf_v) * y * P

#    print np.sum(g - g_zero), np.sum(dV), VF
#    if np.sum(g - g_zero)**2 < 1e-12:
#        print "converged"
#        break

    
#    f0 = RR_Vapor(z,0.0,K)
#    f1 = RR_Vapor(z,1.0,K)
#    if f0 < 0 or f1 > 0:
#        Wl,x,TPD_l,tpd_l,tm_l = StabilityLiquid2(R,P,T,z,Pc,Tc,omega,kij)
#        Wv,y,TPD_v,tpd_v,tm_v = StabilityVapor2(R,P,T,z,Pc,Tc,omega,kij)
#    else:        
#        VF,x,y = RR_Vapor_Calc(z,VF,K,1e-6,200)
#        print VF
#        x = z * (1 + VF*(K - 1))
#        y = K * x
#    dG_x,Zc_x,rho_x,lf_l = SRK(R, P, T, x, Pc, Tc, omega, kij)
#    dG_y,Zc_y,rho_y,lf_v = SRK(R, P, T, y, Pc, Tc, omega, kij)
#    fug_l = np.exp(lf_l) * x * P
#    fug_v = np.exp(lf_v) * y * P
#    err = np.sum( fug_l/fug_v) - 1
#    print err
#    K = K * (fug_l/fug_v)




## Residual Vector
#def f_envelope(z,T,P,P0,K,VF,lf_y,lf_x):
#    f = np.zeros([len(K)+2])
#    f[:len(K)] = np.log(K) + lf_y - lf_x   
#    f[-2] = np.sum(z*(K-1) / (1 + VF*(K-1)))
#    f[-1] = np.log(P) - np.log(P0)
#    return f
#f_nc  = lambda K,lf_y,lf_x: np.log(K) - (lf_x - lf_y)
#f_nc1 = lambda z,K,VF: np.sum( (z*(K-1)) / (1 + VF*(K-1)) )
#f_nc2 = lambda P,P0: np.log(P) - np.log(P0)





## Solve ########################################################################s
#f = -1.*f_envelope(z,T,P,P0,K,VF,lf_y,lf_x)


#J = np.zeros([len(V),len(V)])

## Lazy Jacobian Method
#for i in np.arange(0,len(K)):
#    J[i,i] = 1.0 / K[i]

#J[:len(K),i] = np.sum( (z/(VF*(K-1) + 1)) - ( (z*VF*(K-1)) / (VF*(K-1)+1)**2))
#J[-1,-1] = 1 / P 
#J_lazy = J.copy()




## Computationally Expensive Method
##def df_envelope(z,T,P,P0,K,VF,lf_y,lf_x):
#df = np.zeros([len(K)+2,len(K)+2])

#dKplus  = K + eptol
#dKminus = K - eptol

#for i in np.arange(0,len(K)):
#    K_FDplus     = K.copy()    
#    K_FDminus    = K.copy()
#    K_FDplus[i]  = dKplus[i]
#    K_FDminus[i] = dKminus[i]
#    for j in np.arange(0,len(K)):
#        xplus  = y/K_FDplus
#        xminus = y/K_FDminus

#        yplus  = x*K_FDplus
#        yminus = x*K_FDminus
#        
##        xplus.clip(min=0)
##        xminus.clip(min=0)
##        yplus.clip(min=0)
##        yminus.clip(min=0)
##        K_FDplus.clip(min=0)
##        K_FDminus.clip(min=0)
#        
#        Zc,lf_xplus  = EOS(P,T,xplus,Pc_PS,Tc_PS,omega_PS,kij_PS)
#        Zc,lf_xminus = EOS(P,T,xminus,Pc_PS,Tc_PS,omega_PS,kij_PS)        

#        Zc,lf_yplus  = EOS(P,T,yplus,Pc_PS,Tc_PS,omega_PS,kij_PS)
#        Zc,lf_yminus = EOS(P,T,yminus,Pc_PS,Tc_PS,omega_PS,kij_PS)    
#        
#        fncplus      = f_nc(K_FDplus[i],lf_yplus[i],lf_xplus[i])
#        fncminus     = f_nc(K_FDminus[i],lf_yminus[i],lf_xminus[i])    

#        J[j,i] = (fncplus - fncminus) / (2*eptol)

#f = np.zeros([len(K+2)])
#f[:len(K)] = np.log(K) + lf_y - lf_x   
#f[-2] = np.sum(z*(K-1) / (1 + VF*(K-1)))
#f[-1] = np.log(P) - np.log(P0)




#singlephase,tpd_l,tpd_v,K_l,K_v=StabilityTest(P,T,z,Pc_PS,Tc_PS,omega_PS,kij_PS)
#K = K_l
#while True:
#    
#    # Initial RR Run
#    VF,x,y = RR_Vapor_Calc(z,VF,K,1e-8)
## Run EOS to update K, etc.
#    Zc_v,lf_v = EOS(P,T,y,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    Zc_l,lf_l = EOS(P,T,x,Pc_PS,Tc_PS,omega_PS,kij_PS)    
#    
## Generate K factors
#    fu_l = lf_l*x*P
#    fu_v = lf_v*y*P    
#    if np.sum(fu_l) / np.sum(fu_v)  < eptol:
#        break
#    else:        
#        K = K*(fu_l / fu_v)
#        print('VF, eps')
#        print VF, eps
#        
#        


#f = RR_F(z_PS,K,VF)
#df = (RR_F(z_PS,K,VF+eptol) - RR_F(z_PS,K,VF-eptol))/(2*eptol)
#VF = VF + f/df

#y = K*z_PS / (1 - VF + K*VF)
#x = z_PS / (1 - VF + K*VF)




#while True:
#    # Run EOS to update K, etc.
#    Zc_z_PS,lnfugcoeff_z_PS = EOS(Pb,T_res,z_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    Zc_v_PS,lnfugcoeff_v_PS = EOS(Pb,T_res,y_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    Zc_l_PS,lnfugcoeff_l_PS = EOS(Pb,T_res,x_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)

#    # Update K from fugacity
#    K_PS = np.exp(lnfugcoeff_z_PS - lnfugcoeff_v_PS)

#    # BP Equilibrium
#    f_bp = np.sum(z_PS*K_PS) - 1.

#    # Derivative

#    Zc_z_plus, lf_z_plus = EOS(Pb+eptol,T_res,z_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    Zc_z_minus,lf_z_minus = EOS(Pb-eptol,T_res,z_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)

#    Zc_v_plus, lf_v_plus = EOS(Pb+eptol,T_res,y_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    Zc_v_minus,lf_v_minus = EOS(Pb-eptol,T_res,y_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)

#    Zc_l_plus, lf_l_plus = EOS(Pb+eptol,T_res,x_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    Zc_l_minus,lf_l_minus = EOS(Pb-eptol,T_res,x_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)

#    dlf_dp_z = (lf_z_plus - lf_z_minus) / (2.*eptol)
#    dlf_dp_v = (lf_v_plus - lf_v_minus) / (2.*eptol)
#    dlf_dp_l = (lf_l_plus - lf_l_minus) / (2.*eptol)


#    df_dp_bp = np.sum(z_PS*K_PS * (dlf_dp_z-dlf_dp_v))
#    Pb = Pb - f_bp / df_dp_bp
#    print(Pb,f_bp)
#    y_PS = K_PS * z_PS
#    if f_bp**2 < 1E-12:
#        break

################################################################################
# Direct Solution
################################################################################
# Define # mols in
#mols_PS = 1.0                   # mol/mol
# Define total mass in
#massin_PS = np.sum(z_PS*MW_PS)  #  g
#P_test  =   Pb*1.5 #atm

# Vapor Fraction initial guess
#VF_PS =  np.sum(y_PS*z_PS)
#LF_PS = 1 - VF_PS


# Check if single phase
#singlephase_PS,tpd_l,tpd_v,K_x_stb,K_y_stb = StabilityTest(P_test,T_res,z_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)


## Begin with Wilson
##K_PS = WilsonEqn(Pc_PS,Tc_PS,omega_PS,T_res,P_test)
#K_PS = K_x_stb
#eps = 1
## Begin Iteration
#while True:

#    VF_PS = RR_Vapor_Calc(z_PS,VF_PS,K_PS,eptol)
#    LF = 1 - VF_PS
#    # Update MF Phase fractions
#    y_PS = (z_PS * K_PS) / (1 + VF_PS*(K_PS - 1))
#    x_PS = z_PS / (1 + VF_PS*(K_PS - 1))
#    # Total Mass / Phase
#    mass_v_PS  =   np.sum(y_PS*MW_PS)
#    mass_l_PS  =   np.sum(x_PS*MW_PS)
#    # EOS Compressibility Vapor 
#    Zc_v_PS,lnfugcoeff_v_PS = EOS(P_test,T_res,y_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    Zc_l_PS,lnfugcoeff_l_PS = EOS(P_test,T_res,x_PS,Pc_PS,Tc_PS,omega_PS,kij_PS)
#    # Total Vol / Phase
#    vol_v_PS = (Zc_v_PS * mols_PS*VF_PS * T_res) / P_test
#    vol_l_PS = (Zc_l_PS * mols_PS*(1-VF_PS) * T_res) / P_test 
#    # Fugacities
#    fug_v_PS = P_test * y_PS * np.exp(lnfugcoeff_v_PS)
#    fug_l_PS = P_test * x_PS * np.exp(lnfugcoeff_l_PS)

#    # Update K Vals
#    eps = np.abs(np.sum( (fug_l_PS/fug_v_PS )**2 -1))
#    print(eps)
#    if eps < eptol:
#        break
#    K_PS = (y_PS/x_PS) * (fug_l_PS/fug_v_PS)

#    
        
                
        
        
                
################################################################################
# Output Data Tables
################################################################################
# Given Data
#F8Comp.to_latex(buf='F8_props.tex',float_format=np.str)

