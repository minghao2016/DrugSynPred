function [minX,res,ii] = conjGrad(initialX,b,A,maxIt)
% function minX = conjGrad(initialX,b,A)
% conjugate gradient method to find the minimum of vector X in the system 
% b = A*x where A is the nxn Data matrix (kernel matrix here) and b is
% class label vector of nx1
% maxIt is the maximum iteration number


% required parameters:
r = b - A*initialX;
d = r;
x = initialX;

for ii = 1:maxIt
   
    step = (r'*r)/((d'*A)*d);
    x = x + step*d;
    rn = r - step*A*d;
    
    Beta = (rn'*rn)/(r'*r);
    dn = rn+Beta*d;

    d = dn;
    r = rn;
    if(sum(r.^2) < 1e-16)
        break
    end

end

res = r;
minX = x;

