# Psuedoreduced values

import numpy as np

T = 348.15
P = 2.9607698
Tc = np.array([ 190.6,  425.2,  617.7])
Pc = np.array([ 44.806316,  37.503084,  20.922773])
omega = np.array([ 0.008,  0.199,  0.489])
z = np.array([ 0.35,  0.45,  0.2 ])
x = np.array([0.318,0.472,0.210])
R = 82.057338
kij = np.array([[ 0.    ,  0.0133,  0.0422],
       [ 0.0133,  0.    ,  0.0078],
       [ 0.0422,  0.0078,  0.    ]])


Tr = T/Tc
Pr = P/Pc

# SRK Parameters
u  = 1.
w  = 0.

fw = 0.480 + 1.574*omega - 0.176*omega**2
alpha = (1.0 + fw * (1 - Tr**0.5))**2
a  = (0.42748 *  ( (R**2 * Tc**2)  / Pc ) ) * alpha
b  = 0.08664 * ( ( R * Tc ) / Pc )

A  = a * ( P / ( R * T )**2 )
B  = b * ( P / ( R * T ) ) 

# Reshape for linear algebra

x = x.reshape([x.size,1])

A = A.reshape([A.size,1])
B = A.reshape([B.size,1])
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
    AIJ = (1 - kij) * (A * A.T)**0.5
    AM  = np.sum(x*x.T*AIJ)
    BM  = np.sum(x*B)
elif np.size(Pc) == 1:
    AM = A
    BM = B
coeff = np.zeros(4)
    # build cubic polynomial
def g(z):
    """
    Cubic polynomial in z from EOS. This should be zero.
    :param z: float compressibility factor
    """
    return z**3 - (1 + BM - u*BM)  * z**2 +           \
                  (AM + w*BM**2 - u*BM - u*BM**2) * z    -           \
                  AM*BM   - w*BM**2 - w*BM**3
            
            
coeff[0] =  1.0
coeff[1] = -1.0         - BM      + u*BM
coeff[2] =  AM          + w*BM**2 - u*BM - u*BM**2
coeff[3] = -1.0*AM*BM   - w*BM**2 - w*BM**3
z = np.real( np.roots(coeff)[np.flatnonzero(np.iscomplex(np.roots(coeff))==False)])


z = np.array([z.min(),z.max()])

if np.size(Pc) > 1:
    lf =  B/BM*(z-1) - np.log(z-BM) + AM/BM*( B/BM - 2./AM *                   \
          np.sum(x*AIJ,axis=1).reshape([x.size,1]) ) * np.log(1. + BM/z)

else:
    lf = z - 1. - np.log(z - B) - A/B*np.log(1.+B/z)

gibbs = lf*x
g_star = np.sum(lf*x,axis=0)
index = np.flatnonzero(g_star == g_star.min())
index = index[0]
Zc    = z[index]
lnfug = lf[:,index]
dG    = g_star[index] 
    
    
















# Pick root that results in lowest Gibbs Energy
# V = maximum root,
# l = minimum root
# EOS is run per phase, thus not specific here
# Vapor Like

