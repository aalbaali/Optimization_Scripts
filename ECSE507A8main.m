% Amro Al Baali
% March 22, 2018
% For A8 of ECSE 507
clear all;
close all;
home;


% in = input('1: Gradient method, 2: Globalized Newton, 3: Globalized BFGS, 4: Globalized inexact Newton method, 5: All methods\n');
in = 5;


eps = 1e-6;
x0 = [-1.2;1];
if in == 5
    first = 1;
    last = 4;
else
    first = in;
    last = in;
end

for in = first:last
    
    fprintf("Method: ");
    if in == 1
        fprintf("Gradient Method\n");

        [R_argmin, R_iteratinos] = gradientMethod_Armijo_Function(@Rosenbrock,x0);    
        [H_argmin, H_iteratinos] = gradientMethod_Armijo_Function(@Himmelblau, x0);

    elseif in == 2
        fprintf("Newton Method\n");
        [R_argmin, R_iteratinos] = globalNewton(@Rosenbrock,x0);    
        [H_argmin, H_iteratinos] = globalNewton(@Himmelblau, x0);    
    elseif in == 3
        fprintf("Globalized BFGS\n");
        [R_argmin, R_iteratinos] = BFGSglobal(@Rosenbrock,x0);    
        [H_argmin, H_iteratinos] = BFGSglobal(@Himmelblau, x0);    
    elseif in == 4
        fprintf("Globalized inexact Newton\n");
        [R_argmin, R_iteratinos] = globalInexactNewton(@Rosenbrock,x0);    
        [H_argmin, H_iteratinos] = globalInexactNewton(@Himmelblau, x0);    
    end 


    fprintf("Function               argmin                  Iterations\n");
    fprintf("Rosenbrock             %.13f                    %i\n",R_argmin(1),R_iteratinos);
    fprintf("                       %.13f                     \n\n",R_argmin(2));

    
    fprintf("Himmelblau             %.13f                   %i\n",H_argmin(1),H_iteratinos);
    fprintf("                        %.13f                     \n",H_argmin(2));
end 


