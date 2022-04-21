% Bmatrix
function [B]=Bi(s,t,Xc,Yc)
x1=Xc(1); x2=Xc(2); x3=Xc(3); x4=Xc(4);
y1=Yc(1); y2=Yc(2); y3=Yc(3); y4=Yc(4);
a=0.25*(y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s));
b=0.25*(y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t));
c=0.25*(x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t));
d=0.25*(x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s));

N1s=0.25*(t-1); N2s=0.25*(1-t); N3s=0.25*(1+t); N4s=0.25*(-1-t);
N1t=0.25*(s-1); N2t=0.25*(-1-s); N3t=0.25*(1+s); N4t=0.25*(1-s);

B1=[a*N1s-b*N1t 0; 0 c*N1t-d*N1s; c*N1t-d*N1s a*N1s-b*N1t];
B2=[a*N2s-b*N2t 0; 0 c*N2t-d*N2s; c*N2t-d*N2s a*N2s-b*N2t];
B3=[a*N3s-b*N3t 0; 0 c*N3t-d*N3s; c*N3t-d*N3s a*N3s-b*N3t];
B4=[a*N4s-b*N4t 0; 0 c*N4t-d*N4s; c*N4t-d*N4s a*N4s-b*N4t];
J=Jacobian(s,t,Xc,Yc);
B=(1/J)*[B1 B2 B3 B4];

end
