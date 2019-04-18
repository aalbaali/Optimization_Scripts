%Globalized BFGS
function [argmin iterations minval fgrad xList] = BFGSglobal(func, x0, H0, t0, sigma, rho,epsilon, gamma)
fgrad = [];
xList = {};
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
    fgrad = [fgrad;norm(fgradOld)];
    xList{end+1} = xOld;
    if (rem(k,100)==0)
        k;
    end
    
    d = -H\fgradOld; %solving the linear system
    
    t = Wolfe_Powll_rule(func, xOld, d, t0, gamma, rho, sigma);
%     t = mb_nocLineSearch(func, xOld, d, t0, gamma, rho, sigma); %wolfe from online source
%     t = 1e-3;
    xNew = xOld + t*d;    
   [fval fgradNew] = func(xNew);
   
%    % Recently added for the brachistochrone problem
%    while norm(imag(fgradNew)) > epsilon
%        t = 0.5*t;
%        xNew = xOld + t*d;    
%        [fval fgradNew] = func(xNew);
%    end
    s = xNew - xOld;    
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