%%%
% Solving 1/2*x'Qx + c'x s.t. Ax=b.
%     Where Q=Q'.
%     Uses linear solvers of choice.
% 
%     Basically solves [Q -A'; A  0] = [-c; b];
function [qpSol, x, mu] = quadraticProgram(Q,c,A,b)
if nargin == 2
    A = [];
    b = [];
end

qpSol = [Q A'; A zeros(size(A,1),size(A,1))]\[-c;b];
x = qpSol(1:size(Q,1));
mu = qpSol(size(Q,1)+1,end);
end

    