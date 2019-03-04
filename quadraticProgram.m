%%%
% Solving 1/2*x'Qx + c'x s.t. Ax=b.
%     Where Q=Q'.
%     Uses linear solvers of choice.
% 
%     Basically solves [Q -A'; A  0] = [-c; b];
function x = quadraticProgram(Q,c,A,b)
if nargin == 2
    A = [];
    b = [];
end

xHat = [Q -A'; A zeros(size(A,1),size(A,1))]\[-c;b];
x = xHat(1:size(Q,1));
end

    