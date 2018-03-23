%% Example 2
function [val, grad, hess] = funcA7(x)

val = NaN;
grad = NaN;
hess = NaN;

n = length(x);
m = n+2;

Fn1 = 0;
Fn2 = 0;

for i = 1:n        
    Fi = x(i)-1; %Fi
    F(i,1) = Fi;    
%     fSum = fSum + Fi^2;    
    Fn1 =Fn1+i*(Fi);  %F = x(i)-1 (where j = i in this case);   
end

Fn2 = Fn1^2;

F(n+1,1) = Fn1;
F(n+2,1) = Fn2;

%f(x)
val = 0.5*F'*F;

Fn1Grad = [1:n]';

DF = [eye(n), Fn1Grad, 2*Fn1*Fn1Grad]';

grad = DF'*F;

Fn2Hess = 2*(Fn1Grad)*(Fn1Grad)';

hess = DF'*DF + Fn2*Fn2Hess;

end        