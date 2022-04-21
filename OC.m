%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(elements,nodes,vc,dc,passive);
xold=elements(:,5);
l1 = 0; l2 = 100000; move = .2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(xold-move,min(1.,min(xold+move,xold.*sqrt(-dc'./lmid)))));
  xnew(find(passive)) = 1;
  if sum(xnew.*elements(:,6)) - sum(vc*elements(:,6))> 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end