function [argmin, iterations, minVal] = globalInexactNewton(func, x0, epsilon, rho, p, beta, sigma, c1, c2, maxIterations)
if nargin == 2
    epsilon = 1e-6;
    rho = 1e-8;
    p = 2.1;
    beta = 0.5;
    sigma = 1e-4;
    c1 = 1e-2;
    c2 = 1;
    maxIterations = 1000;
end

if(epsilon<0 || rho <= 0 || p <=2 || beta <= 0 || beta >= 1 || sigma<=0 || sigma >=0.5)
    disp("Error with paramaters");
    return;
end
argmin = NaN;
iterations = NaN;
minVal = NaN;

x = x0;
k = 0;

if (nargin < 10)
    maxIterations = 10000;
end


[fVal fGrad fHess] = func(x);




while (norm(fGrad) > epsilon)
    %use conjugate gradient method to solve for d
    %if A is not symmetric positive definite then d = NaN;    
    eta = min(c1/(k+1), c2*norm(fGrad));
    d = CGforLinEqns(fHess, -fGrad, -fGrad, eta*norm(fGrad));
    if(size(d,1) ~= size(x0,1))
%         fprintf("CG not compatible for step k = %i. => Gradient method used\n",k);
%         return;
    else
%         fprintf("CG Compatible at k = %i\n",k);
    end
        
    if ( any(isnan(d)) || fGrad'*d > -rho*norm(d).^p)
        d = -fGrad;
    end
       
    %determining tk
    l = 0; %lower case L    
    while(func(x+(beta^l)*d) > func(x) + sigma*beta^l*fGrad'*d)
        l = l+1;
    end
    t = beta^l;    
    
    x = x + t*d;
    
    k = k + 1;
    
    if k == maxIterations
        disp("max iterations reached");
        break;
    end      
    
    [fVal fGrad fHess] = func(x);
end

argmin = x;
iterations = k;
minVal = func(x);