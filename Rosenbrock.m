function [val grad hess] = Rosenbrock(x)
x1=x(1);
x2=x(2);

val = 100*(x2-x1^2)^2+(1-x1)^2;

grad  =[-400*x1*(x2-x1^2)-2*(1-x1);
    200*(x2-x1^2)];

hess = [-400*(x2-x1^2)+800*x1^2+2, -400*x1;
        -400*x1, 200];
    

end