%Amro Al Baali
%Feb 9, 2018

close all;
clear all;
home;

iterationsTot = zeros(2,6);
value1 = zeros(4,6);
value2 = zeros(4,6);
%% The function to be minimized
for i=1:1:6;
    delta = 10^(-i);
    Q = eye(4);
    Q(4,4) = delta;

    c = ones(4,1);
    gamma = 0;


    objFun = @(x) quadraticFunc(x,Q,c,gamma);
    %% optimization criteria
    eps = 1e-5; %termination criteria

    x0 = 10*rand(4,1);

    [argmin iterations]= gradientMethod_Armijo_Function(objFun,x0);
    % [argmin, iterations]= acceleratedGradientMethod(objFun,x0, delta, 1);
    % disp(iterations);
    % x0 = 10*rand(4,1);
    
    
    iterationsTot(1,i) = iterations; 
    value1(:,i) = argmin;
    
    [argmin iterations]= gradientMethod_Armijo_Function(objFun,x0,1e-5,0.5,1e-1);
    
   
    x0 = argmin;

    y0 = x0;
    y = y0;
    alphaOld = 1;
    k = 0;
    L = norm(Q);
    xOld = x0;
    maxIterations = 10000;

    [fval, fgradX] = objFun(xOld);
    disp('gradX');
    disp(norm(fgradX));
    while (norm(fgradX)>1e-5)
        [fval, fgradY] = objFun(y);
        xNew = y -1/L*fgradY;
        alphaNew = 0.5*(1+sqrt(1+4*alphaOld^2));
        y = xNew+(alphaOld-1)/(alphaNew)*(xNew-xOld);
        [fval fgradX] = objFun(xNew);

        xOld = xNew;
        alphaOld = alphaNew;
        k = k+1;
        if k == maxIterations
            disp('max iterations reached');
        end
    end

    disp(k);
    disp(xNew);
    

   
   
   iterationsTot(2,i) = k;
   
   
   value2(:,i) = xNew;
   
end

%Amro AL Baali
%Feb 9, 2018
% f(x) = 1/2*x'Qx+c'x+gamma
%fval is the function value while 
%fgrad is the gradient value
function [fval, fgrad] = quadraticFunc (x,Q,c,gamma)
fval = 1/2*x'*Q*x+c'*x+gamma;

fgrad = 1/2*Q*x+Q'*x+c;




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
maxIteratins = 1e4;

[fval fgrad] = func(x);

while norm(fgrad) > eps
    [fval fgrad] = func(x);
    
    d = -fgrad;
    l = 0; %lower case L

    t = 0;
    while(func(x+beta.^l*d)>func(x)+beta^l*sigma*fgrad'*d)
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

    
    




