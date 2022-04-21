% Topology optimization code with SIMP
% Hossein Rostami, Avdanced Topology Optimization course ENGR 5016G
% Advisor: Professor Ahmad Barari, University of Ontario Institude of Technology
% 2D quadrilateral elements mesh, imported from ANSYS

clear
clc
% Read the mesh from text file(generated in ANSYS)
nodes=readmatrix('nd4.txt'); 
elements=readmatrix('el4.txt');
elements(:,2)=[]; % remove element type names, not necessary
elements(:,1)=[];

% Initializing parameters
nod=size(nodes,1);%number of nodes
DOF=2*nod; % degrees of freedom for 2D nodes
d=NaN(DOF,1);
f=NaN(DOF,1); % The force matrix to apply known forces 

%B.Cs.
%case 2, clamped both sides
for i=1:nod
    if nodes(i,2)==0 || nodes(i,2)==500; % clamped left edge
        d(2*i-1)=0;
        d(2*i)=0;
    else
        f(2*i-1)=0;
        f(2*i)=0;
    end
end

% force on all elements on the top
for i=1:nod
    if nodes(i,3)==100
        f(2*i)=-1e1;
    end
end

% calculation of K matrix in global coordinate 
Kg=zeros(DOF);
ele=size(elements);
el=ele(:,1); % number of elements
Loc=zeros(el,8,8); % all of the local stifness matrices 
glob=zeros(el,8);
for e=1:el
    for i=1:4; %x and y for each node
        xc(i)=nodes(elements(e,i),2);
        yc(i)=nodes(elements(e,i),3);
    end
%     patch(Xc,Yc, 'g');
%     hold on
    Xc(e,:)=xc;
    Yc(e,:)=yc;
    ke=kel(Xc(e,:),Yc(e,:));
    Loc(e,:,:)=ke;
    %nodes
    n1=elements(e,1); n2=elements(e,2); n3=elements(e,3); n4=elements(e,4);
    glob(e,:)=[2*n1-1 2*n1 2*n2-1 2*n2 2*n3-1 2*n3 2*n4-1 2*n4];
    Ke=zeros(DOF); Ke(glob(e,:),glob(e,:))=ke; % stiffness matrix of element in global coordinate
    Kg=Kg+Ke;

% Kg(glob,glob)=Kg(glob,glob)+ke;
end

free=[]; % free degree of freedoms
for i=1:DOF
    if d(i) ~= 0
        free=[free,i];
    end
end

d(free)=(Kg(free,free))\(f(free));
% 
% % calculate stress
E=200e3;
v=0.3;  
D=(E/(1-(v*v)))*[1 v 0; v 1 0; 0 0 0.5*(1-v)];
for e=1:size(elements)

    n1=elements(e,1); n2=elements(e,2); n3=elements(e,3); n4=elements(e,4);
    [B]=Bi(0,0,Xc(e,:),Yc(e,:));
    dis=[d(2*n1-1); d(2*n1); d(2*n2-1); d(2*n2); d(2*n3-1); d(2*n3); d(2*n4-1); d(2*n4)];
    sigma=D*B*dis;
    Stress(e,:)=sigma;
    disx(e,:)=[dis(1) dis(3) dis(5) dis(7)];
    disy(e,:)=[dis(2) dis(4) dis(6) dis(8)];
    Xcd(e,:)=Xc(e,:)+disx(e,:);
    Ycd(e,:)=Yc(e,:)+disy(e,:);
    Mises(e)=sqrt(sigma(1)*sigma(1)-sigma(1)*sigma(2)+sigma(2)*sigma(2)+3*sigma(3)*sigma(3));
end

%Topology Optimization
p=3; %penalization factor
elements(:,5)=.6; %initial solution
S1=Xc(:,1).*(Yc(:,2)-Yc(:,3))+Xc(:,2).*(Yc(:,3)-Yc(:,1))+Xc(:,3).*(Yc(:,1)-Yc(:,2));
S2=Xc(:,1).*(Yc(:,3)-Yc(:,4))+Xc(:,3).*(Yc(:,4)-Yc(:,1))+Xc(:,4).*(Yc(:,1)-Yc(:,3));
elements(:,6)=.5*(S1+S2); % surface area of each element to calculate volume fraction by density
elements(:,6)=elements(:,6)/sum(elements(:,6));
elem=[1 2 3 4 5 6 7 8];
vc=.9; % desired volume fraction

% Definition of passive elements(solid or void in some elements)
xh1=40; xh2=460; % center location of hole
yh1=50; yh2=50;
t=12; %thickness of passive region
rh=25+t;
for i=1:el
    xp(i)=mean(Xc(i,:)); % approximate center of element
    yp(i)=mean(Yc(i,:));
    Dist1(i)=sqrt((xp(i)-xh1)^2+(yp(i)-yh1)^2);
    Dist2(i)=sqrt((xp(i)-xh2)^2+(yp(i)-yh2)^2);
    if Dist1(i)<rh || Dist2(i)<=rh
        passive(i)=1;
        elements(i,5)=1;
    else
        passive(i)=0;
    end
end

loop=0;
change=1;
while change >.01
    loop=loop+1;
    xold=elements(:,5);
    Kg=zeros(DOF);
    for e=1:el
    kke(e,:,:)=(elements(e,5)^p)*Loc(e,:,:);
    Ke=zeros(DOF); Ke(glob(e,:),glob(e,:))=kke(e,:,:); % stiffness matrix of element in global coordinate and penalized
    Kg=Kg+Ke;
    end
    d(free)=(Kg(free,free))\(f(free));
    c=0;
    for i=1:el
        n1=elements(i,1); n2=elements(i,2); n3=elements(i,3); n4=elements(i,4); % connected nodes
        Ue=[d(2*n1-1); d(2*n1); d(2*n2-1); d(2*n2); d(2*n3-1); d(2*n3); d(2*n4-1); d(2*n4)]; %displacements of nodes
        ki=zeros(8); ki(elem,elem)=Loc(i,:,:);
        c=c+(elements(i,5)^p)*(Ue'*ki*Ue); % objective of compliance
        dc(i)=-p*(elements(i,5)^(p-1))*(Ue'*ki*Ue); % sensitivity
    end
    [dcn]=check(elements,el,dc,xp,yp);
    dc=dcn';
    [xnew]=OC(elements,nodes,vc,dc,passive);
      elements(:,5)=xnew;
      comp(loop)=c;
    change=max(abs(elements(:,5)-xold));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',(sum(elements(:,5).*elements(:,6)))) ...
        ' ch.: ' sprintf('%6.3f',change )])
end


for i=1:el
%     if elements(i,5)==1
    patch(Xc(i,:), Yc(i,:),-elements(i,5));
    hold on
%     end
end
colormap(gray)
axis equal
% for i=1:el
% patch(Xcd(i,:), Ycd(i,:), Ycd(i,:)-Yc(i,:)+Xcd(i,:)-Xc(i,:));
% colorbar
% end