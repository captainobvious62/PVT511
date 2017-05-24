global Pv

clc
Tc = 273.15 + 32.28; 
Pc = 48.84; 
Omega= 0.09860; 
R = 83.14;
b = 0.08664*R*Tc/Pc;
m = 0.480 + 1.574*Omega - 0.176*Omega^2;
j=1;
for i=1:5:30
    v=0.001:1:2500;
    T(i)=273.15+i;
    Tre = T(i)/Tc;
    a = 0.42748*(R*Tc)^2/Pc*(1 + m*(1 - sqrt(Tre)))^2;
    P=R*T(i)./(v - b) - a./(v.*(v + b));
    figure(1)
    h=plot(v,P);
    set(h,'color',rand(1,3));
    hold on
    axis([0 2500 0 100])
    xlabel('Volume in cm3/mol')
    ylabel('pressure in bar')
    title('Isotherms for ethane')
    Pv=[Pv P'];
    figure(2)
    h=plot(v,P);
    set(h,'color',rand(1,3));
    hold on
    axis([0 700 0 50])
    xlabel('Volume in cm3/mol')
    ylabel('pressure in bar')
    title('Isotherms for ethane')
    clear v
    X=fsolve(@(x) equations3(x,T(i)),[80 120 600])
    plot(max(X),feval(@(v) Pressure(v,T(i)),max(X)),'bo')
    plot(min(X),feval(@(v) Pressure(v,T(i)),max(X)),'bo')
    psat(j)=feval(@(v) Pressure(v,T(i)),max(X));
    Tsat(j)=T(i);
    j=j+1;
end
    figure(3)
    plot(Tsat,psat,'bo')
    hold on
    T=273.15:1:305.15;
    Psattheo=exp(44.0103-2.56882e3./T-4.97635*log(T)+1.46447e-5*T.^2)/100;
    plot(T,Psattheo,'c')
    xlabel('Temperature in K')
    ylabel('pressure in bar')
    title('Vapour pressure versus temperature for ethane')
    
    i=1;
    v=0.001:1:2500;
    j=1;
    T(i)=273.15+i;
    Tre = T(i)/Tc;
    a = 0.42748*(R*Tc)^2/Pc*(1 + m*(1 - sqrt(Tre)))^2;
    P=R*T(i)./(v - b) - a./(v.*(v + b));
    Pv=P';
    figure(4)
    plot(v,P,'r')
    hold on
    axis([0 700 -20 50])
    clear v
    X=fsolve(@(x) equations3(x,T(i)),[80 120 600])
    plot(max(X),feval(@(v) Pressure(v,T(i)),max(X)),'bo')
    plot(min(X),feval(@(v) Pressure(v,T(i)),max(X)),'bo')
    xlabel('Volume in cm3/mol')
    ylabel('pressure in bar')
    title('Isotherm at 1°C for ethane')

    T=274.15
    z=Pv(:,1)-feval(@(v) Pressure(v,T),max(X))*ones(2500,1);
    j=1;
    for p=1:2500
        if (abs(z(p,1))<1) 
            indice(j)=p;
            j=j+1;
        end
    end
    min(indice)
    max(indice)
    
    y=Pv(min(indice)+2:max(indice)-40,1)';
    x=feval(@(v) Pressure(v,T),max(X))*ones(max(indice)-40-min(indice)-2+1,1);
    x=x';
    [ph,msg]=jbfill(min(indice)+2:max(indice)-40,y,x,rand(1,3),rand(1,3),0)

