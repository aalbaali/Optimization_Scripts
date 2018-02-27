% Amro Al Baali
%Feb 26, 2018
%ECSE507 A6

clear all;
close all;
home;

epsilon = 1e-6;
sigma = 1e-1;
rho = 0.9;

gamma = 2;
t0 = 1;


%Rosenbrock
H0 = eye(2);
x0 = [-1.2 1]';


% c = [1 1]';
% Q=eye(2); Q(2,2)=0;
func = @Rosenbrock;
% func = @(x) quadraticFunc(x,Q,c,0);


% [argmin, iterations, minval] = BFGSglobal(func, x0, H0, t0, sigma, rho,epsilon, gamma);


t = [15:10:85, 105 185 245 305 365];
y = [0.038, 0.085, 0.1, 0.103, 0.093, 0.095, 0.088, 0.08, 0.073, 0.05, 0.038, 0.028, 0.02];
m = size(y,2);


% 
% yFunc = @(t,x) Bateman(t,x);
% objFunc = @(x) f(x,yFunc, y, t);
x0 = [0.05, 0.1, 0.4]';
H0 = eye(3);
% 
% interval = 0:0.1:365;
% [y2] = yFunc(interval,x0);
% plot(interval,y2);
% hold on;
% plot(t,y,'x');


funcModel = @(t,x) Bateman(t,x);
objFunc = @(x) leastSquares(x, funcModel, y, t);


[argmin, iterations, minval] = BFGSglobal(objFunc, x0, H0, t0, sigma, rho,epsilon, gamma);


% 
% 
% 
% 
% function [fVal fgrad] = Bateman(t,x)    
%     x1 = x(1);
%     x2 = x(2);
%     x3 = x(3);
%     
%     fVal = x3*(exp(-x1*t)-exp(-x2*t));
%     
%     fgrad = [-x3*exp(-x1*t);
%             x3*exp(-x2*t);
%             exp(-x1*t)-exp(-x2*t)];
% end
% 
% 
% function [fVal fgrad] = f(x, yFunc, y, t)
% m = size(y,2);
% sum1 = 0;
% sum2 = 0; %for the derivative
%     for i=1:m
%         [yi ygrad] = yFunc(t(i),x);
%         sum1 = sum1+(yi-y(i))^2;
%         sum2 = sum2+(yi-y(i))'*ygrad;
%     end
% fVal = 0.5*sum1;
% fgrad = sum2;
% 
% end
