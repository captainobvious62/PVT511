function f=Pressure(v,T)

Tc = 273.15 + 32.28; 
Pc = 48.84; 
Omega= 0.09860; 
R = 83.14;
b = 0.08664*R*Tc/Pc;
m = 0.480 + 1.574*Omega - 0.176*Omega^2;
Tre = T/Tc;
a = 0.42748*(R*Tc)^2/Pc*(1 + m*(1 - sqrt(Tre)))^2;
f=R*T./(v - b) - a./(v.*(v + b));

