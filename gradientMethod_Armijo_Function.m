%Amro Al Baali
%Feb 9, 2018
%Consule Algorithm 3.2.1 from ECSE 507 Notes: Gradient method with Armijo
%rule
%The function must return the function value and a gradient.

function [argmin, iterations] = gradientMethod_Armijo_Function(func,x0,sigma,beta, eps)
if (nargin == 2)
    sigma = 1e-5;
    beta = 0.5;
    eps = 1e-5;
elseif (beta>=1 || beta<=0 || eps <=0)
    disp('Wrong paramaters used beta should be Beta and sigma in (0,1), eps >0');
end

x = x0;
k = 0;
maxIterations = 1e8;

[fval fgrad] = func(x);

while norm(fgrad) > eps
    [fval fgrad] = func(x);
    
    d = -fgrad;
    l = 0; %lower case L

    t = 0;
    while(func(x+(beta^l)*d)>func(x)+(beta^l)*sigma*fgrad'*d)
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
