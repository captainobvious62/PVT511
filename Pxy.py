# 
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
plt.close("all Matlab code")
y1row = list()
Prow = list()
# binary interaction parameters
k12b = [0.1278, 0, 0] # [-]


Tc1 = 304.2 # critical temperature CO2[K]
Tc2 = 469.66 # critical temperature pentane[K]

w1 = 0.231 # asymmetric factor CO2 []
w2 = 0.251 # asymmetric factor pentane []

pc1 = 7.38152 * 10**6 # critical pressure CO2 [Pa]
pc2 = 3.36906 * 10**6 # critical pressure pentane [Pa]

T = 304 # Temperature [K]

x1row = np.linspace(0.0, 1.0, 150) # Molar fractions liquid phase [-]
# Prow = zeros(numel(x1row),1); # Initializing pressure-vector 
# y1row = zeros(numel(x1row),1); # Initializing vapor mole fraction vector 

R = 8.3145 # Gas constant [J mol^-1 K^-1]

a11c = 0.45724 * R**2 * Tc1**2 / pc1 # a at critical point for CO2 [J^2 mol^-1 / Pa]
a22c = 0.45724 * R**2 * Tc2**2 / pc2 # a at critical point for pentane [J^2 mol^-1 / Pa]

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

for h in range(0, len(x1row)):

    x1 = x1row[h] # Specifiying a value for x1
    y1_new = 0.9 # Initial guess for y1
    y2_new = 1 - y1_new # Initial guess for y2
    check = 0.2 # Initial value for criterion
    x2 = 1 - x1 # specifying x2

    Pneu = 0.5 * 10**5 # Initial guess for pressure [Pa] 

    b = x1 * b1 + x2 * b2 # b for the mixture in the liquid phase [J/(mol Pa)]
    a = x1**2 * a11 + 2 * a12 * x1 * x2 + x2**2 * a22 # a for the mixture in the liquid phase

    while abs(check) > 0.000005: # checking if convergence criterium fulfilled

        y1 = y1_new
        y2 = y2_new
        P = Pneu

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
        liquid_2 = -math.log(abs(Zliq) - B) + (abs(Zliq) - 1) * b2 / b - A / (math.sqrt(8) * B) * (1 / a * (2 * x1 * a21 + 2 * x2 * a22) - b2/b) * math.log((abs(Zliq) + B * (math.sqrt(2) + 1)) / (abs(Zliq) - B * (math.sqrt(2) - 1)))
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

     # y1row[h] = y1 # Saving the value for y1 into the vector 
     y1row.append(y1) # Saving the value for y1 into the vector

     # Prow[h] = P # Saving the found value for the pressure
     Prow.append(P/10**5)

# Plotting the data
plt.plot(x1row, Prow)
plt.plot(y1row, Prow)
ChemCADx = [0,0.04,0.09,0.14,0.19,0.24,0.29,0.34,0.39,0.44,0.49,0.54,0.59,0.64,
0.69,0.74,0.79,0.84,0.89,0.94,1]
ChemCADP = [83924,400179,795795,1.19E+06,1.58E+06,1.97E+06,2.36E+06,2.73E+06,
3.10E+06,3.45E+06,3.79E+06,4.12E+06,4.43E+06,4.71E+06,4.98E+06,5.24E+06,5.50E+06,
5.77E+06,6.09E+06,6.51E+06,8.04E+06]
ChemCADP = [x / 10**5 for x in ChemCADP]
ChemCADx2 = [0,0.780251,0.884557,0.919441,0.936741,0.946959,0.953616,0.958226,
0.961551,0.964019,0.965891,0.967344,0.968502,0.969469,0.970337,0.97127,0.972208,
0.973682,0.97586,0.980514,1]
plt.plot(ChemCADx,ChemCADP,marker='oinestyle='')
plt.plot(ChemCADx2,ChemCADP,marker='oinestyle='')
plt.plot(0.476294, 162934/(10**5),marker='o', linestyle='')
plt.xlabel('Molar fraction of $CO_2$ (x/y)')
plt.ylabel('Pressure / bar')

axes = plt.gca()
# axes.set_xlim([0.0,1.0])




plt.figure()
plt.plot(x1row, y1row)
plt.plot([0, 1], [0, 1])
ChemCADx = [0,0.04,0.09,0.14,0.19,0.24,0.29,0.34,0.39,0.44,0.49,0.54,0.59,0.64,0.69,0.74,
0.79,0.84,0.89,0.94,1]
ChemCADy = [0,0.780251,0.884557,0.919441,0.936741,0.946959,0.953616,0.958226,
0.961551,0.964019,0.965891,0.967344,0.968502,0.969469,0.970337,0.97127,0.972208,
0.973682,0.97586,0.980514,1]
plt.plot(ChemCADx,ChemCADy,marker='oinestyle='')
plt.xlabel('Molar fraction of $CO_2$ liquid phase (x)')
plt.ylabel('Molar fraction of $CO_2$ vapor phase (y)')
