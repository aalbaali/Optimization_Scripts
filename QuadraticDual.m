% This script is to test my quadratic dual funciton

clear all;
close all;

x0p = [1 1]';

Q = [3 1;1 5];
c = [1 2]';
gamma = 0;

A = [1 2];
b = 2;


quadFunc = @(x,Q,c,gamma) 1/2*x'*Q*x+c'*x+gamma;

f1 = @(x) quadFunc(x,Q,c,gamma);

argminAnalyt = -inv(Q)*c;

[argminPrimal, minPrimal] = fmincon(f1,x0p,[],[],A,b,zeros(size(Q,1),1),inf(size(Q,1),1));

Qinv = inv(Q);
Lagrange = @(x,lam,mu) 1/2*x'*Q*x+c'*x-lam'*x+mu'*(b-A*x);
gPos = @(lam,mu) g(lam,mu,Lagrange,true,Qinv,c,A,b);
gNeg = @(lam,mu) g(lam,mu,Lagrange,false,Qinv,c,A,b);
 


n = size(Q,1);
m = size(A,1);

Anew = [eye(n), zeros(n,m)];
Bnew = [eye(n), A'];
bNew = [zeros(n,1);b];


Qinv = inv(Q);

Qnew = -Bnew'*Qinv*Bnew;
cNew = (c'*Qinv*Bnew+bNew')';
gamma = -1/2*c'*Qinv*c;

x0d = ones(n+m,1);

fDual = @(y) -(1/2*y'*Qnew*y+cNew'*y+gamma);

[argminDual, minDual] = fmincon(fDual, x0d, -Anew,zeros(n,1));

% [argminDual, minDual] = fmincon(@(y) gNeg(y(1:n),y(n+1:n+m)), x0d,-Anew,zeros(n,1));

minDual = -minDual;


% 
% p = size(A,1); % mu
% m = size(Q,1); % lambda
% x0d = ones(m+p,1);
% 
% [argminDual, valDual] = fmincon(@(y) gNeg(y(1:m),y(m+1:m+p)), x0d,[],[],[],[],[zeros(m,1);-inf(p,1)],inf(m+p,1));
% 
% valDual = -valDual;
% 
% 
% % dual = fmincon(g,x0d);
% 
% Anew = [eye(m), zeros(m,p)];
% 
% Bnew = [-eye(m), -A'];
% 
% Qinv = inv(Q);
% 
% % bNew = [zeros(m,1);b];
% bNew = b;
% 
% Qnew = -Bnew'*Qinv*Bnew;
% cNew = (-c'*Qinv*Bnew+bNew')';
% gammaNew = -1/2*c'*Qinv*c;
% 
% f2 = @(x) -quadFunc(x,Qnew,cNew,gammaNew);
% 
% 
% 
% [argminDual, minDual] = fmincon(f2,x0d,Anew,zeros(m,1));
% minDual = -minDual;

