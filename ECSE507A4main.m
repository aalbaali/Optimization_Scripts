%% The main script 
%Amro Al Baali
%Feb 9, 2018

close all;
clear all;
home;

iterationsTot = zeros(2,6);
value1 = zeros(4,6);
value2 = zeros(4,6);
% The function to be minimized
for i=1:1:6;
    delta = 10^(-i);
    Q = eye(4);
    Q(4,4) = delta;

    c = ones(4,1);
    gamma = 0;


    objFun = @(x) quadraticFunc(x,Q,c,gamma);
    % optimization criteria
    eps = 1e-5; %termination criteria

    x0 = 10*rand(4,1);

    [argmin iterations]= gradientMethod_Armijo_Function(objFun,x0);
    % [argmin, iterations]= acceleratedGradientMethod(objFun,x0, delta, 1);
    % disp(iterations);
    % x0 = 10*rand(4,1);
    
    
    iterationsTot(1,i) = iterations; 
    value1(:,i) = argmin;
    
%     [argmin iterations]= gradientMethod_Armijo_Function(objFun,x0,1e-5,0.5,1e-1);
    
   
%     x0 = argmin;

    y0 = x0;
    y = y0;
    alphaOld = 1;
    k = 0;
    L = norm(Q);
    xOld = x0;
    maxIterations = 10000;

    [fval, fgradX] = objFun(xOld);
%     disp('gradX');
%     disp(norm(fgradX));
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

%     disp(k);
%     disp(xNew);
%     

   
   
   iterationsTot(2,i) = k;
   
   
   value2(:,i) = xNew;
   
end


for i=1:6
    delta=10^(-i);
%     Q(4,4)=delta;
%       [xopt_grad,it_grad]=gradq(Q,c,x0);
%       [xopt_accgrad,it_accgrad]=accgradq(Q,c,x0);
      fprintf('%4.6f || %9.0f | %1.0f  %1.0f  %1.0f  %14.6f|| %7.0f | %1.0f  %1.0f  %1.0f   %15.6f||\n', delta, iterationsTot(1,i), value1(1,i), value1(2,i), value1(3,i), value1(4,i), iterationsTot(2,i), value2(1,i), value2(2,i), value2(3,i), value2(4,i));

end

