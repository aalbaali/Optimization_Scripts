t = [245 305 365];
y = [0.038, 0.028, 0.02];
m = length(y);
fVal = 0;
x1 =  0.005293431637108;
x2 = 0.099998667023739;
x3 = 0.139358214006687;
x =[x1;x2;x3];

[fvalFunc fgradFunc] = leastSquares(x, @Bateman, y, t);

fVal = 0;
fGrad1 = 0;
fGrad2 = 0;
fGrad3 = 0;
for i=1:m
    fVal = fVal +(x3*(exp(-x1*t(i))-exp(-x2*t(i)))-y(i))^2;
    fGrad1 = fGrad1 + (x3*(exp(-x1*t(i))-exp(-x2*t(i)))-y(i))*(-t(i)*x3*exp(-x1*t(i)));
    fGrad2 = fGrad2 + (x3*(exp(-x1*t(i))-exp(-x2*t(i)))-y(i))*(t(i)*x3*exp(-x2*t(i)));
    fGrad3 = fGrad3 + (x3*(exp(-x1*t(i))-exp(-x2*t(i)))-y(i))*(exp(-x1*t(i))-exp(-x2*t(i)));
end

fVal = fVal*0.5;
fGrad = [fGrad1; fGrad2; fGrad3];

normValDiff = norm(fvalFunc - fVal)
normGradDiff = norm(fgradFunc-fGrad)