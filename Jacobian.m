 % Determinant of Jacobian matrix |J|
function [J]=Jacobian(s,t,Xc,Yc)
T=[0 1-t t-s s-1; t-1 0 s+1 -s-t; s-t -s-1 0 t+1; 1-s s+t -t-1 0];
J=(1/8)*Xc*T*transpose(Yc);
end