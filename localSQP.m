% func, x0, H0, t0, sigma, rho,epsilon,
% subject to equality constraints only h = h(x)
function [argmin, iterations, minval, xArr, fArr, mu] = localSQP(func,h,x0,mu0,epsilon, maxIterations)
if nargin <= 5
    maxIterations = 1e5;
    if nargin <= 4
        epsilon = 1e-5;
    end
end

x = x0;
mu = mu0;
[fVal fGrad fHess] = func(x);
[hVal hGrad hHess] = h(x);
lagGradx = fGrad+mu*hGrad;  % gradient of Lagrangian w.r.t. x
lagHessxx = fHess+mu*hHess;

xArr = [x'];
fArr = [fVal'];
    
k = 0;
while ~(norm(lagGradx)< epsilon && norm(hVal) < epsilon) && k < maxIterations
    [qpSol, d, mu] = quadraticProgram(lagHessxx,fGrad, hGrad', -hVal);
    
    x = x + d;    
    [fVal fGrad fHess] = func(x);
    [hVal hGrad hHess] = h(x);
    lagGradx = fGrad+mu*hGrad;  % gradient of Lagrangian w.r.t. x
    lagHessxx = fHess+mu*hHess;
    
    xArr = [xArr; x'];
    fArr = [fArr; fVal'];

    k = k + 1;
end
mu = -mu;
argmin = x;
iterations = k;
minval = fVal;
end
