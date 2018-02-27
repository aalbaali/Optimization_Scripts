%Amro Al Baali
%Feb 15, 2018
%ECSE 507 Assignment 5

%Finding the minimum of the Rosenbrock function and the Brown function
%using the globalized Newton method
clear all;
close all;
clc;

rho = 1e-8;
p=2.1;
beta = 0.5;
sigma = 1e-4;
eps = 1e-6;
kMax = 200;


rosenbrock = @(x) Rosenbrock(x);

x0 = [1.2; 1];
[argmin_r, k_r] = globalNewton(rosenbrock,x0,rho,p,beta,sigma,eps,kMax);


brown = @(x) Brown(x);

x0 = [1; 1];
[argmin_b k_b] = globalNewton(brown,x0,rho,p,beta,sigma,eps,kMax);

disp('Rosenbrock minimizer: ');
%%
% 

% 
disp(argmin_r);
disp('Rosenbrock iterations: ');
disp(k_r);;
disp('Brown minimizer: ');
disp(argmin_b);
disp('Brown iterations: ');
disp(k_b);


%$$k_b$