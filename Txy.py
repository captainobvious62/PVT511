# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 18:54:41 2016

@author: ingoldo
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 13:13:10 2016

@author: ingoldo

Program for the calulation of vapor liquid equilibria using SRK
 equation of state.
"""
import numpy as np
import math
import matplotlib.pyplot as plt
plt.close("allef ReturnPY(T, params):
 
 
 (Pgoal, x1) = params
 
 k12b = [0.0133, 0, 0] # [-]


 Tc1 = 304.2 # critical temperature CO2 [K]
 Tc2 = 469.66 # critical temperature pentane [K]

 w1 = 0.231 # asymmetric factor CO2 []
 w2 = 0.251 # asymmetric factor pentane []

 pc1 = 7.38152 * 10**6 # critical pressure CO2 [Pa]
 pc2 = 3.36906 * 10**6 # critical pressure pentane [Pa]

 # T = 188 # Temperature [K]

 # Prow = zeros(numel(x1row),1); # Initializing pressure-vector 
 # y1row = zeros(numel(x1row),1); # Initializing vapor mole fraction vector 

 R = 8.3145 # Gas constant [J mol^-1 K^-1]

 a11c = 0.45724 * R**2 * Tc1**2 / pc1 # a at critical point for CO2 [J^2 mol^-1 / Pa]
 a22c = 0.45724 * R**2 * Tc2**2 / pc2 # a at critical point for pentane [J^2 mol^-1 / Pa] 
 
 
 check = 0.2 # Initial value for criterion
 x2 = 1 - x1 # specifying x2

 Pneu = 1 * 10**5 # Initial guess for pressure [Pa] 
 y1_new = x1
 y2_new = 1 - x1


 while abs(check) > 0.005: # checking if convergence criterium fulfilled

     y1 = y1_new
     y2 = y2_new
     P = Pneu
 
 
     k12 = k12b[0] + k12b[1] * T + k12b[2] / T # Calculating mixing parameter
     k21 = k12

     kappa1 = 0.37464 + 1.54226 * w1 - 0.26992 * w1**2 # Parameter for critical point correction CO2
     kappa2 = 0.37464 + 1.54226 * w2 - 0.26992 * w2**2 # Parameter for critical point correction pentane

     alpha1 = (1 + kappa1 * (1 - math.sqrt(T / Tc1)))**2 # Parameter for critical point correction CO2
     alpha2 = (1 + kappa2 * (1 - math.sqrt(T / Tc2)))**2 # Parameter for critical point correction pentane

     a11 = a11c * alpha1 # Critical point corrected a CO2 [J^2 mol^-1 / Pa]
     a22 = a22c * alpha2 # Critical point corrected a pentane [J^2 mol^-1 / Pa]

     a12 = (1 - k12) * math.sqrt(a11 * a22) # Mixing rule [J^2 mol^-1 / Pa] 
     a21 = (1 - k21) * math.sqrt(a11 * a22) # Mixing rule [J^2 mol^-1 / Pa]

     b1 = 0.0778 * R * Tc1 / pc1 # b for CO2 [J/(mol Pa)]
     b2 = 0.0778 * R * Tc2 / pc2 # b for pentane [J/(mol Pa)]
 
     b = x1 * b1 + x2 * b2 # b for the mixture in the liquid phase [J/(mol Pa)]
     a = x1**2 * a11 + 2 * a12 * x1 * x2 + x2**2 * a22 # a for the mixture in the liquid phase

 # Coefficients of the EoS model equation
     A = a * P / (R**2 * T**2)
     B = b * P / (R * T)

 # Cubic equation for the gas phase
     c = [1, -(1-B), A - 2 * B - 3 * B**2, -(A * B - B**2 - B**3)]


 # Roots finds the three solutions for the cubic equation
     r = np.roots(c)

     Zliq = min(r) # liquid compressibility factor is the smallest solution

 # a and b for the gas phase
     bg = y1 * b1 + y2 * b2
     ag = y1**2 * a11 + 2 * a12 * y1 * y2 + y2**2 * a22

     A = ag * P / (R**2 * T**2)
     B = bg * P / (R * T)
 # Cubic equation
     c[0] = 1
     c[1] = -(1 - B)
     c[2] = A - 2 * B - 3 * B**2
     c[3] = -(A * B - B**2 - B**3)

 # Roots finds the three solutions for the cubic equation
     r = np.roots(c)

     Zgas = max(r) # gas compressibility factor is the biggest solution

 # Coefficients of the EoS model equation for the liquid phase
     A = a * P / (R**2 * T**2)
     B = b * P / (R * T)

 # Fugacity coefficients. Formula from http://kshmakov.org/fluid/note/3/
     liquid_1 = -math.log(abs(Zliq) - B) + (abs(Zliq) - 1) * b1 / b - A / (math.sqrt(8) * B) * (1 / a * (2 * x1 * a11 + 2 * x2 * a12) - b1/b) * math.log((abs(Zliq) + B * (math.sqrt(2) + 1)) / (abs(Zliq) - B * (math.sqrt(2) - 1)))
     liquid_2 = -math.log(abs(Zliq) - B) + (abs(Zliq) - 1) * b2 / b - A / (math.sqrt(8) * B) * (1 / a * (2 * x1 * a12 + 2 * x2 * a22) - b2/b) * math.log((abs(Zliq) + B * (math.sqrt(2) + 1)) / (abs(Zliq) - B * (math.sqrt(2) - 1)))
 # Values of A and B for the gas phase
     A = ag * P / (R**2 * T**2)
     B = bg * P / (R * T)
 # Fugacity coefficients gas phase
     vapor_1 = -math.log(abs(Zgas) - B) + (abs(Zgas)-1) * b1/bg - A / (math.sqrt(8) * B) * (1 / ag * (2 * y1 * a11 + 2 * y2 * a12) - b1 / bg) * math.log((abs(Zgas) + B * (math.sqrt(2) + 1)) / (abs(Zgas) - B * (math.sqrt(2) - 1)))
     vapor_2 = -math.log(abs(Zgas) - B) + (abs(Zgas)-1) * b2/bg - A / (math.sqrt(8) * B) * (1 / ag * (2 * y1 * a21 + 2 * y2 * a22) - b2 / bg) * math.log((abs(Zgas) + B * (math.sqrt(2) + 1)) / (abs(Zgas) - B * (math.sqrt(2) - 1)))

 # sum of y1 and y2
     S = x1 * math.exp(liquid_1) / math.exp(vapor_1) + (1 - x1) * math.exp(liquid_2) / math.exp(vapor_2)
     check = (1 - S)

 # New guesses for y1 and y2
     y1_new = (x1 * math.exp(liquid_1) / math.exp(vapor_1)) / (S)
     y2_new = ((1 - x1) * math.exp(liquid_2) / math.exp(vapor_2)) / (S)
     Pneu = P * S
     return Pneu, y1_new, y2_new



def TempProblem(T, params):
 (Pgoal, x1) = params
 P, y1, y2 = ReturnPY(T, params)
 return P - Pgoal





# Matlab code
y1row = list()
Trow = list()
# binary interaction parameters

x1row = np.linspace(0.0, 1.0, 100) # Molar fractions liquid phase [-]
T = 200
Pgoal = 100000

for h in range(0, len(x1row)):

     x1 = x1row[h] # Specifiying a value for x1
     params = (Pgoal, x1)
     Perr = TempProblem(T, params)
     while abs(Perr) > 100:

         if Perr > 0:
             T = T-0.1*Perr/Pgoal
         elif Perr < 0: 
             T = T-0.1*Perr/Pgoal
         Perr = TempProblem(T, params)
 
     P, y1, y2 = ReturnPY(T, params)

     y1row.append(y1) # Saving the value for y1 into the vector

     Trow.append(T)

# Plotting the data
plt.plot(x1row, Trow)
plt.plot(y1row, Trow)
ChemCADx = [0,0.04,0.09,0.14,0.19,0.24,0.29,0.34,0.39,0.44,0.49,0.54,0.59,0.64,
0.69,0.74,0.79,0.84,0.89,0.94,1]
ChemCADT = [308.983,235.627,209.686,198.868,192.719,188.78,186.121,184.299,
183.073,182.295,181.872,181.737,181.842,182.148,182.618,183.216,183.898,184.602,
185.234,185.631,185.363]
ChemCADy = [0,0.968996,0.995635,0.998421,0.999178,0.999479,0.999626,0.999707,
0.999753,0.99978,0.999795,0.999799,0.999795,0.999782,0.999758,0.99972,0.999658,
0.999561,0.999409,0.999215,1]

plt.plot(ChemCADx,ChemCADT,marker='oinestyle='')
plt.plot(ChemCADy,ChemCADT,marker='oinestyle='')
plt.plot([0.607161, 0.85661], [284.096, 262.213], marker='o', linestyle='')
plt.xlabel('Molar fraction of $CO_2$ (x/y)')
plt.ylabel('Temperature [K]')
axes = plt.gca()


plt.figure()
plt.plot(x1row, y1row)
plt.plot([0, 1], [0, 1])
ChemCADx = [0,0.04,0.09,0.14,0.19,0.24,0.29,0.34,0.39,0.44,0.49,0.54,0.59,0.64,
0.69,0.74,0.79,0.84,0.89,0.94,1]
ChemCADy = [0,0.968996,0.995635,0.998421,0.999178,0.999479,0.999626,0.999707,
0.999753,0.99978,0.999795,0.999799,0.999795,0.999782,0.999758,0.99972,0.999658,
0.999561,0.999409,0.999215,1]
plt.plot(ChemCADx,ChemCADy,marker='oinestyle='')
plt.xlabel('Molar fraction of $CO_2$ liquid phase (x)')
plt.ylabel('Molar fraction of $CO_2$ vapor phase (y)')





