%% Least Squares
function [val grad] = leastSquares(x, func, yList, tList)

t = tList;
y = yList;

m = length(yList);

x1 = x(1);
x2 = x(2);
x3 = x(3);

fVal = 0;
fGrad1 = 0;
fGrad2 = 0;
fGrad3 = 0;
fGrad = 0;

for i=1:m
    [val grad] = func(t(i),x);

    fVal = fVal + (val-y(i))^2;
    
    fGrad1 = fGrad1 + (val-y(i))*grad(1);
    fGrad2 = fGrad2 + (val-y(i))*grad(2);
    fGrad3 = fGrad3 + (val-y(i))*grad(3);
end

val = 0.5*fVal;
grad = [fGrad1; fGrad2; fGrad3];

end
