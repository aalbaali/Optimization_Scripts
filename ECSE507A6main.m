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


% [argmin, iterations, minval] = BFGSglobal(func, x0, H0, t0, sigma, rho,epsilon, gamma)


t = [15:10:85, 105 185 245 305 365];
y = [0.038, 0.085, 0.1, 0.103, 0.093, 0.095, 0.088, 0.08, 0.073, 0.05, 0.038, 0.028, 0.02];
% t = [245 305 365];
% y = [0.038, 0.028, 0.02];
m = length(y);


% 
% yFunc = @(t,x) Bateman(t,x);
% objFunc = @(x) f(x,yFunc, y, t);
x0 = [0.05, 0.1, 0.4]';
H0 = eye(3);
% 




objFunc = @(x) leastSquares(x, @Bateman, y, t);


[argmin, iterations, minval] = BFGSglobal(objFunc, x0, H0, t0, sigma, rho, epsilon, gamma);

% %fmincon
% [x1,fval1,exitflag1,output1,lambda1,grad1,hessian1] = fmincon(objFunc,x0);
% options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(objFunc, x0, [],[],[],[],[],[],[],options);

interval = 0:0.1:365;
for i=1:length(interval)
    y2(i) = Bateman(interval(i),x);
    y3(i) = Bateman(interval(i),argmin);
end


% plot(interval,y2);
hold on;
plot(t,y,'x');
plot(interval,y3);


