function [val grad] = Bateman(t,x)
x1 = x(1);
x2 = x(2);
x3 = x(3);

val = x3*(exp(-x1*t)-exp(-x2*t));

grad = [-x3*exp(-x1*t);
        x3*exp(-x2*t);
        exp(-x1*t)-exp(x2*t)];
end
