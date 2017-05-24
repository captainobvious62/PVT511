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



def SRK_Z(R,P,T,x,Pc,Tc,omega,kij):
    
    # Psuedoreduced values
    Tr = T/Tc
    Pr = P/Pc

    # SRK Parameters
    u  = 1.
    w  = 0.
    b  = 0.08664 * R * Tc / Pc
    fw = 0.48 + 1.574*omega - 0.176*omega**2
    a  = 0.42748 * R**2 * T**2  / Pc  
    alpha = (1.0 + fw * (1 - np.sqrt(Tr)))**2

    # Incorporate Mixture 
    if np.size(Pc) > 1:
        a = np.sum(x * x.reshape([np.size(x),1]) *                                 \
            ( np.sqrt(a * a.reshape([np.size(a),1])) * (1 - kij) )  )
    
        b = np.sum(x * b)

    # Set Coefficients
    
    AC = a * alpha * P / R**2 / T**2
    BC = b * P / R / T

    coeff = np.zeros(4)
        # build cubic polynomial
    def g(z):
        """
        Cubic polynomial in z from EOS. This should be zero.
        :param z: float compressibility factor
        """
        return z**3 - (1 - BC) * z**2 + (AC - 2*BC - 3*BC**2) * z - (
                AC * BC - BC**2 - BC**3)
                
                
    coeff[0] = 1.0
    coeff[1] = -1.0*BC - 1.0 + u*BC
    coeff[2] = AC + w*BC**2 - u*BC - u*BC**2
    coeff[3] = -1.0*AC*BC - w*BC**2 - w*BC**3
    cubic = lambda Z:   Z**3 -                                                    \
                        Z**2 * (1.0 + BC - u*BC) +                                \
                        Z    * (AC + w*BC**2 - u*BC - u*BC**2) -                  \
                                AC*BC - w*BC**2 - w*BC**3
    Z = np.roots(coeff).real
    z = op.newton(g, 1.0)  # compressibility factor
    rho = P / (R * T * z)  # density
    fugacity_coeff = np.exp(z - 1 - np.log(z - BC) - AC / np.sqrt(8) / BC * np.log(
                (z + (1 + np.sqrt(2)) * BC) / (z + (1 - np.sqrt(2)) * BC)))
    #cubic = lambda Z:   Z**3 -                                                    \
    #                    Z**2 * (1.0 + BC - u*BC) +                                \
    #                    Z    * (AC + w*BC**2 - u*BC - u*BC**2) -                  \
    #                            AC*BC - w*BC**2 - w*BC**3

    # V = maximum root,
    # l = minimum root
    # EOS is run per phase, thus not specific here
    lnfug_v = np.nan_to_num(lnfugcoef(np.max(z),x,Tc,Pc,u,w,a,AC,BC,kij))
    lnfug_l = np.nan_to_num(lnfugcoef(np.min(z),x,Tc,Pc,u,w,a,AC,BC,kij))
    print z
    if np.size(Pc) == 1:
        z.clip(min=0)
        Zn = np.zeros(2)
        lnfug = np.zeros(2)
        # Liquid
        Zn[0] = z.min()
        lnfug[0] = lnfug_l
        # Vapor
        Zn[1] = z.max()
        lnfug[1] = lnfug_v
    elif np.size(Pc) > 1:                    
        dG = np.sum(x*lnfug_l) - np.sum(x*lnfug_v)
        if dG > 0:
            Zn = np.max(z)
            lnfug = lnfug_v
        elif dG < 0:
            Zn = np.min(z)
            lnfug = lnfug_v
        elif dG == 0:
            Zn = z
            lnfug = 0.0

#    B_ratio = b/np.sum(x*b)
#    A_i = (am*P) / (R*T)**2.
#    B_i = (bm*P) / (R*T)
#    C_i = np.sqrt(U**2. - 4*W)
#   
#    gamma_i = (2.*np.sqrt(a))/np.sum(x*a) * ( np.dot((x * np.sqrt(a)),(1.-kij)))    

#    Zn = Z.min()
#    lin_fug_l = np.nan_to_num(B_ratio * (Zn - 1.) - np.log(Zn - B_i) + (A_i/(B_i*C_i))*(B_ratio - gamma_i))

#    Zn = Z.max()
#    lin_fug_v = np.nan_to_num(B_ratio * (Zn - 1.) - np.log(Zn - B_i) + (A_i/(B_i*C_i))*(B_ratio - gamma_i))
    return z,rho,fugacity_coeff

def lnfugcoef(Z,x,Tc,Pc,U,W,a,AM,BM,kij):
    # Fugacity Coefficient
    
    bi_o_b = (Tc/Pc) / (np.sum( (x*Tc)/Pc))
    gamma_i = (2.*np.sqrt(a))/np.sum(x*a) * ( np.dot((x * np.sqrt(a)),(1.-kij)))    
    alfa = np.sqrt(U**2 - 4*W)
    
    india  = bi_o_b*(Z-1) - np.log(Z-BM)
    juliet = (AM/(BM*alfa))*(bi_o_b-gamma_i)
    kilo   = np.log( (2.*Z + BM*(U+alfa)) / (2.*Z + BM*(U-alfa)))
    
    lnfug = india + juliet + kilo
    return lnfug        
           
################################################################################
R  = 8.2057338*(10**-5)	    # m3 atm / K *mol

# EOS Verification Test - CO2, 1.0 mol
###########################################
MW      = 44.01   # kg/mol
Tc      = 304.2     # K
Pc      = 72.9      # atm
omega   = 0.228     
kij     = np.array([0.0])
x       = 1.
mol     = .5 
M       = mol*MW

# Test Values

P = np.linspace(0,200,101) #np.array([1,5,10.,15])  # atm
T = 150 # K

Zc = np.zeros([np.size(P),3])

for i in np.arange(0,np.size(P)):
    Zc[i,:] =  SRK_Z(R,P[i],T,x,Pc,Tc,omega,kij)


#EOS_Z = np.vectorize(SRK_Z, excluded=['R','x','Pc','omega','kij'],
#        otypes=[np.ndarray] )
#EOS_Fugacity = np.vectorize(SRK_Fugacity, excluded=['R','x','Pc','omega','kij'],
#        otypes=[np.ndarray] )        
#bZc,lf = SRKEOS(R,P[:],T,x,Pc,Tc,omega,kij)

V = P / (Zc[:,0]*R*T)
rho = V/MW
#rho = P / (R * T * z)  # density
T = 250
for i in np.arange(0,np.size(P)):
    Zc[i,:] =  SRK_Z(R,P[i],T,x,Pc,Tc,omega,kij)
V2 = P / (Zc[:,0]*R*T)
T = 350
for i in np.arange(0,np.size(P)):
    Zc[i,:] =  SRK_Z(R,P[i],T,x,Pc,Tc,omega,kij)
V3 = P / (Zc[:,0]*R*T)



# CO2 Test Plot
fig, ax_V  = plt.subplots(figsize=(10, 6))
ax_rho       = ax_V.twinx()
ax_V.plot(V,P,'r--')   
ax_V.plot(V2,P,'b--')
ax_V.plot(V3,P,'g--')
#ax_rho.plot(P,rho,'b')
plt.tight_layout()

#ax_V.set_xlim(7,SCN[-1])
ax_V.set_xlabel('rho, mol/m**3')
ax_V.set_ylabel('Pressure, atm')
#ax_rho.set_ylabel('Density, g/cc')
plt.savefig('CO2_SRK_test.png',frameon=True,bbox_inches='tight',padinches=0.1)
plt.clf()

###############################################################################
# Test Parameters

compstr = ["CH4","nC4","nC10"]
omega = np.array([0.008,0.199,0.489])
Pc = np.array([454.,380.,212.]) * 1E4
Tc = np.array([190.6,425.2,617.7])
Mw = np.array([16.04,58.12,142.29])
z = np.array([0.35,0.45,0.20])      # Fraction of oil input stream
kij = np.array([ [0.0000, 0.0133, 0.0422],
                     [0.0133, 0.0000, 0.0078],
                     [0.0422, 0.0078, 0.0000]  ])
# Ambient Conditions
n_mols = 1.0                        # # mols entering system
feed_rate = 100.                    # kg/hr
P         = 6.89*1e06               # atm ==> Pa
T         = 137.9+273.15            # deg C ==> K


    
    
## Psuedoreduced values
#Tr = T/Tc
#Pr = P/Pc

## SRK Parameters
#u  = 1.
#w  = 0.
#b  = (0.08664*R*Tc)/Pc
#fw = 0.48 + 1.574*omega - 0.176*omega**2
#a  = ( (0.42748 * (R**2) * (Tc**2)) / Pc ) * (1.0 + fw*(1 - Tr**0.5))**2

## Incorporate Mixture 
#if np.size(Pc) > 1:
#    a = np.sum(x * x.reshape([np.size(x),1]) *                                 \
#        ( np.sqrt(a * a.reshape([np.size(a),1])) * (1 - kij) )  )

#    b = np.sum(x * b)

## Set Coefficients

#AC = a*P / (R*T)**2
#BC = b*P / R*T

#coeff = np.zeros(4)

#coeff[0] = 1.0
#coeff[1] = -1.0*(1.0 + BC - u*BC)
#coeff[2] = AC + w*BC**2 - u*BC - u*BC**2
#coeff[3] = -1.0*AC*BC - w*BC**2 - w*BC**3
#cubic = lambda Z:   Z**3 -                                                    \
#                    Z**2 * (1.0 + BC - u*BC) +                                \
#                    Z    * (AC + w*BC**2 - u*BC - u*BC**2) -                  \
#                            AC*BC - w*BC**2 - w*BC**3
#Z = np.roots(coeff).real

## V = maximum root,
## l = minimum root
## EOS is run per phase, thus not specific here
#lnfug_v = np.nan_to_num(lnfugcoef(np.max(Z),x,Tc,Pc,u,w,a,AC,BC,kij))
#lnfug_l = np.nan_to_num(lnfugcoef(np.min(Z),x,Tc,Pc,u,w,a,AC,BC,kij))

#dG = np.sum(x*lnfug_l) - np.sum(x*lnfug_v)
#if dG > 0:
#    Zn = np.max(Z)
#    lnfug = lnfug_v
#elif dG < 0:
#    Zn = np.min(Z)
#    lnfug = lnfug_v
#elif dG == 0:
#    Zn = Z
#    lnfug = 0.0
#    
