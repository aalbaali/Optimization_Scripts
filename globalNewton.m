% Implementing the globalized Newton method

function [argmin iterations minVal]= globalNewton (func, x0,rho,p,beta,sigma,eps, maxIterations)

if nargin == 2
    rho = 1e-8;
    p = 2.1;
    beta = 0.5;
    sigma = 1e-4;
    eps = 1e-6;
    maxIterations = 200;
end

if (rho <=0 || p<=2 || beta>=1 || beta<=0 || sigma >= 0.5 || sigma <=0 ||eps<0)
    disp('One of the paramaters is not in the proper range');
    return;    
end

x = x0;
k = 0;

[fVal fGrad fHess] = func(x);


while norm(fGrad) > eps
    [fVal fGrad fHess] = func(x);
    
    %if eqn solvable, use Newton, otherwise use gradient
    if (1/cond(fHess) > 1e-12)
        d = -fHess\fGrad;
        if  (fGrad'*d>-rho*norm(d)^p)
            d = -fGrad;
        end
    else
        d = -fGrad;
    end
    
    l = 0; %lower case L

    t = 0;
    while(func(x+beta.^l*d)>func(x)+beta^l*sigma*fGrad'*d)
        l = l+1;
    end

    t = beta^l;
    k = k+1;
    x = x+t*d;
    
    if k == maxIterations
        break;
    end      

end

argmin = x;
iterations = k;
minVal = func(x);
