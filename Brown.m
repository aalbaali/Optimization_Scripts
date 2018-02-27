%% This is the Brow function

function [val grad hess] = Brown(x)
x1 = x(1);
x2 = x(2);


val = (x1-1e6)^2+(x2-2e-6)^2+(x1*x2-2)^2;

grad = [2*(x1-1e6)+2*x2*(x1*x2-2);
        2*(x2-2e-6)+2*x1*(x1*x2-2)];
    
hess = [2+2*x2^2, 2*(x1*x2-2)+2*x1*x2;
        2*(x1*x2-2)+2*x1*x2, 2*x2+2*x1^2];

    
    
    