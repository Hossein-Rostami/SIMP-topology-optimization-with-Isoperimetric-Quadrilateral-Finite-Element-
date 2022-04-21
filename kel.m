% quadrilateral elements 
function [ke]=kel(Xc,Yc)
E=200e3;
v=0.3; 
D=(E/(1-(v*v)))*[1 v 0; v 1 0; 0 0 0.5*(1-v)];
h=1; % assumed constant thickness
% Four- point Gaussian quadrature
s1=-0.5773; t1=-0.5773;
s2=-0.5773; t2=0.5773;
s3=0.5773; t3=-0.5773;
s4=0.5773; t4=0.5773;

B11=Bi(s1,t1,Xc,Yc); J11=Jacobian(s1,t1,Xc,Yc);
B22=Bi(s2,t2,Xc,Yc); J22=Jacobian(s1,t1,Xc,Yc);
B33=Bi(s3,t3,Xc,Yc); J33=Jacobian(s1,t1,Xc,Yc);
B44=Bi(s3,t3,Xc,Yc); J44=Jacobian(s1,t1,Xc,Yc);

% Weights 
W1=1; W2=1; W3=1; W4=1;% assume W1=W2=W3=W4

ke=transpose(B11)*D*B11*J11*h*W1*W1+transpose(B22)*D*B22*J22*h*W2*W2+...
    transpose(B33)*D*B33*J33*h*W3*W3+transpose(B44)*D*B44*J44*h*W4*W4;
end





