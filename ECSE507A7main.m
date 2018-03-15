% Amro Al Baali
%Mar 15, 2018
%ECSE507 A6

clear all;
close all;
home;
format long;

beta = 0.5;
epsilon = 1e-6;
sigma = 1e-4;
rho = 1e-8;
p = 2.1;

c1 = 1e-2;
c2 = 1;


%% Rosenbrock
x0 = [-1.2 1]';

func = @Rosenbrock;


[R_argmin, R_iterations, R_minval] = globalInexactNewton(func, x0, epsilon, rho, p, beta, sigma, c1, c2);



%% Function 2
n = 10;

