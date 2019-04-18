%%%
% Solving 1/2*x'Qx + c'x s.t. Ax=b.
%     Where Q=Q'.
%     Uses linear solvers of choice.
% 
%     Basically solves [Q -A'; A  0] = [-c; b];
function [qpSol, x, mu] = quadraticProgram(Q,c,A,b)
% function [qpSol, x, mu] = quadraticProgram(Q,c,A,b,muOld)
% function [fVal, fGrad, fHess] = quadraticProgram(x, Q,c,A,b)
if nargin == 2
    A = [];
    b = [];
end
p = length(b);

qpSol = [Q A'; A zeros(size(A,1))]\[-c;b];
% qpSol = [Q A'; A zeros(size(A,1))]\[-c-A'*muOld;b];

x = qpSol(1:size(Q,1));
if (p > 0)
    mu = qpSol(size(Q,1)+1:end);
else
    mu = NaN;
end

    