%% CG for Linear Eqns
function [x, iterations] = CGforLinEqns(A, b, x0, epsilon)
% This only solves equations with A being symmetric and positive definite

[m,n] = size(A);

x = NaN;    
iterations = NaN;

if(m~=n || norm(A-A')> 1e-4 || any(eig(A)<= ones(size(A,1),1)*1e-4))
%     disp("A is not symmetric positive definite");   
    return;
end

x = x0;
k = 0;

gOld = A*x-b;
d = -gOld;

while (norm(gOld) > epsilon)
    t = (norm(gOld)^2)/(d'*A*d);
    
    x = x + t*d;
    gNew = gOld + t*A*d;
    
    beta = (norm(gNew)^2)/(norm(gOld)^2);
    
    d = -gNew+beta*d;
    
    gOld = gNew;
    
    
    k = k+1;
    if (k == n+1)
        disp("something's wrong: k > n");        
        return;
    end
end
iterations = k;   
