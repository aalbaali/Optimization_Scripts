% CG method
function [argmin iterations minval] = cgFR(func,x0, d0, t0, sigma, rho,epsilon, gamma)

if nargin == 2
    t0 = 1;
    
    sigma = 1e-4;
    rho =0.9;
    epsilon = 1e-6;
    gamma = 2;
    d0 = NaN;
end
    

xOld = x0;
[fval fgradOld] = func(xOld);
if isnan(d0)
    d0 = -fgradOld;
end

d = d0;

k = 0;
while(norm(fgradOld)>epsilon)
%     if (rem(k,100)==0)
%         k;
%     end
%     t = Wolfe_Powll_rule(func, xOld, d, t0, gamma, rho, sigma); 
    
    l =0;
    while(func(xOld+(0.5^l)*d)>func(xOld)+(0.5^l)*sigma*fgradOld'*d)
        l = l+1;
     end
    t = 0.5^l;


    xNew = xOld+t*d;
    [~, fgradNew] = func(xNew);
%     beta =  (fgradNew'*fgradNew)/(fgradOld'*fgradOld); % FR
%     beta =  (fgradNew'*(fgradNew-fgradOld))/(fgradOld'*fgradOld); % PR
    beta =  (fgradNew'*(fgradNew-fgradOld))/((fgradNew-fgradOld)'*d); % HS
    
    d = -fgradNew + beta*d;
    
    fgradOld = fgradNew;
    xOld = xNew;
    
    k = k+1;
    
end

iterations = k;
argmin = xOld;
minval = fval;

end