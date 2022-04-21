function [dcn]=check3(elements,el,dc,xp,yp)
dcn=zeros(el,1);
r=2; % radius of sensitivity filter
distance=zeros(el,el);
sum=zeros(el,1);
for i = 1:el
    for j=1:el
        distance(i,j)=(xp(i)-xp(j))^2+(yp(i)-yp(j))^2;
        if distance(i,j)<r*r
            H=r-sqrt(distance(i,j));
            sum(i)=sum(i)+H;
            dcn(i)=dcn(i)+H*elements(j,5)*dc(j);
        end
    end
    dcn(i)=dcn(i)/(elements(i,5)*sum(i));
end

end