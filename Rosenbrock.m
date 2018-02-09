function [val grad] = Rosenbrock(x)
x1=x(1);
x2=x(2);

val = 100*(x2-x1^2)^2+(1-x1)^2;

grad  =[-400*x1*(x2-x1^2)-2*(1-x1);
    200*(x2-x1^2)];