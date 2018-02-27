function [val grad] = leastSquares(x, func, yList, tList)
m = size(yList,2);
sum1 = 0;
sum2 = 0;

for i=1:13
    [fval fgrad] = func(tList(i),x);
    sum1 = sum1 + (fval-yList(i))^2;
    sum2 = sum2 + (fval-yList(i))'*fgrad;
end

val = 0.5*sum1;
grad = sum2;


end
