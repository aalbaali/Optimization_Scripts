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
for i=1:n
    x0(i) = 1-i/n;
end

func = @funcA7;

%f1 for function 1 (n=10)
[f1_argmin, f1_iterations, f1_minval] = globalInexactNewton(func, x0, epsilon, rho, p, beta, sigma, c1, c2);

n = 100;
for i=1:n
    x0(i) = 1-i/n;
end
[f2_argmin, f2_iterations, f2_minval] = globalInexactNewton(func, x0, epsilon, rho, p, beta, sigma, c1, c2);


%% Results
fprintf('Rosenbrock argmin:\n     %0.13f\n     %0.13f\n',R_argmin(1),R_argmin(2));
fprintf('Rosenbrock iterations:\n     %i\n',R_iterations);

fprintf('Example 2 n = 10:\n');
for i=1:10
    fprintf("%0.13f\n",f1_argmin(i));
end
fprintf('Example 2 n = 10: iterations: %i\n',f1_iterations);

fprintf('Example 2 n = 100:\n');

for i=1:9
    r = rem(i,10);
    fprintf("%0.13f     %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f\n",...
        f2_argmin(r), f2_argmin(r+10), f2_argmin(r+20), f2_argmin(r+20), f2_argmin(r+30), f2_argmin(r+40), f2_argmin(r+50), f2_argmin(r+60),...
            f2_argmin(r+70), f2_argmin(r+80), f2_argmin(r+80), f2_argmin(r+90));
end
r = 10;
fprintf("%0.13f     %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f      %0.13f\n",...
        f2_argmin(r), f2_argmin(r+10), f2_argmin(r+20), f2_argmin(r+20), f2_argmin(r+30), f2_argmin(r+40), f2_argmin(r+50), f2_argmin(r+60),...
            f2_argmin(r+70), f2_argmin(r+80), f2_argmin(r+80), f2_argmin(r+90));
fprintf('Example 2 n = 100: iterations: %i\n',f2_iterations);