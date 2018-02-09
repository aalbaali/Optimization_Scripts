function [argmin, iterations] = acceleratedGradientMethod(func,x0,eps, L)

y0 = x0;
y = y0;
alphaOld = 1;
xOld = x0;
[fval fgrad] = func(y);
k = 0;
maxIterations= 2000;
while (norm(fgrad) > eps )    
   [fval, fgradY] = objFun(y);
    xNew = y -1/L*fgradY;
    alphaNew = 0.5*(1+sqrt(1+4*alphaOld^2));
    y = xNew+(alphaOld-1)/(alphaNew)*(xNew-xOld);
    [fval fgradX] = objFun(xNew);
    
    xOld = xNew;
    alphaOld = alphaNew;
    k = k+1;
    
end
argmin = y;
iterations = k;


