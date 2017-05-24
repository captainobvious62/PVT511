function f=equations3(x,T)

f(1)=-quad(@(v) Pressure(v,T),x(1),x(2))+...
    feval(@(v) Pressure(v,T),x(1))*(x(2)-x(1))...
    +quad(@(v) Pressure(v,T),x(3),x(2))...
    -feval(@(v) Pressure(v,T),x(2))*(x(2)-x(3));
f(2)=feval(@(v) Pressure(v,T),x(1))-feval(@(v) Pressure(v,T),x(3));
f(3)=feval(@(v) Pressure(v,T),x(2))-feval(@(v) Pressure(v,T),x(3));

end
