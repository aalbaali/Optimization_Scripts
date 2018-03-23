%Globalized BFGS
function [argmin iterations minval] = BFGSglobal(func, x0, H0, t0, sigma, rho,epsilon, gamma)

if nargin == 2
    H0 = eye(size(x0,1));
    t0 = 1;
    sigma = 1e-4;
    rho =0.9;
    epsilon = 1e-6;
    gamma = 2;
end
    

xOld = x0;
[fval fgradOld] = func(xOld);
H = H0;
k = 0;
while(norm(fgradOld)>epsilon)
    if (rem(k,100)==0)
        k;
    end
    
    d = -H\fgradOld; %solving the linear system
    
    t = Wolfe_Powll_rule(func, xOld, d, t0, gamma, rho, sigma);
%     t = mb_nocLineSearch(func, xOld, d, t0, gamma, rho, sigma); %wolfe from online source
    xNew = xOld + t*d;
    s = xNew - xOld;    
   [fval fgradNew] = func(xNew);
    y = fgradNew - fgradOld;
    H = H + (y*y')/(y'*s)-(H*s*s'*H)/(s'*H*s);
    xOld = xNew;
    fgradOld = fgradNew;
    k = k+1;
    
end

iterations = k;
argmin = xOld;
minval = fval;

end