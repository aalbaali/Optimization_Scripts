function [val grad hess] = Himmelblau(x)
x1 = x(1);
x2 = x(2);
val = (x1^2+x2-11)^2+(x1+x2^2-7)^2;

grad = [4*x1*(x1^2+x2-11)+2*(x1+x2^2-7);
    2*(x1^2+x2-11)+2*(x1+x2^2-7)*2*x2];

hess = [4*(x1^2+x2-11)+8*x1^2+2,     4*x1+4*x2;
        4*x1+4*x2,      2+4*(x1+x2^2-7)+8*x2];
    
end
