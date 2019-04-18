% func, x0, H0, t0, sigma, rho,epsilon,
% subject to equality constraints only h = h(x)
% UPDATE_2019-04-15:
%   accepts multiple equality constraints. Insert them as a "p x 1" CELL array.
%       e.g., H = {h1; h2; h3}, 
%            where h1 = h(x) and returns [value, gradient, hessian]
%   mu as an array of size "p x 1".
%
% Working example:
%   localSQP(@(x) quadraticFunc(x,eye(2),[1;1],0),{@(x) quadraticFunc(x,zeros(2),[1;0],-1); @(x) quadraticFunc(x,zeros(2),[0;1],-1)},[0.5;0.5],[0.1;0.2])
function [argmin, iterations, minval, xArr, fArr, mu, fGradArr] = localSQP(func,H,x0,mu0,epsilon, maxIterations)
if nargin <= 5
    maxIterations = 1e4;
    if nargin <= 4
        epsilon = 1e-5;
    end
end

% if it's a single constraint, convert AND it's not a cell constraint,
% convert it to cell
if (length(H) == 1 && ~iscell(H))
    H = {H};
elseif (length(H) == 0)
    ME = MException('MyComponent:noSuchVariable', ...
        "No constraints are entered! You DO NOT need SQP if it's an unconstrained problem!!");
    throw(ME);    
elseif (length(H) ~= length(mu0))
    ME = MException('MyComponent:noSuchVariable', ...
        "mu0 should have the same length as H.");
    throw(ME);
end

n = length(x0);
x = x0;
mu = mu0;   % lagrange multipliers. Array of size p by 1.
p = length(H); % number of constraints
H_cell = cell(p,1);

[fVal, fGrad, fHess] = func(x);
lagGradx = fGrad;
lagHessxx = fHess;
Hval = NaN(1,p);    % matrix of constraints values
Hgrad = NaN(n,p);   % matrix of constraints' gradients
for i = 1:p
    h = H{i};
    [hVal, hGrad, hHess] = h(x);
    H_cell{i}.hVal = hVal;
    H_cell{i}.hGrad = hGrad;
    H_cell{i}.hHess = hHess;
    
    Hval(i) = hVal;
    Hgrad(:,i) = hGrad;
    
    lagGradx = lagGradx + hGrad*mu(i);  % gradient of Lagrangian w.r.t. x
    lagHessxx = lagHessxx + hHess*mu(i);    
end



xArr = [x'];
fArr = [fVal'];
fGradArr = [norm(fGrad)];
k = 0;
alphaConst = 0.5;
% alphaConst = 1e-6;
beta = 0.9;
muOld = nan;
lag_func = @(x) lagFunc(x(1:n),x(n+1:end),func,H);
sigma = 1e-4;
quadPen = @(x) quadraticPenaltyFun(x,1e-1,func,H);
dx = Inf;
dMu = Inf;
[qpSol, dx, muOld] = quadraticProgram(lagHessxx,fGrad, Hgrad', -Hval');
while (norm(lagGradx)> epsilon || norm(Hval) > epsilon) && k < maxIterations
% while (norm(lagGradx)> epsilon || norm(Hval) > epsilon || norm(fGrad(1:2))>1e-1) && k < maxIterations
% while ((norm(fGrad(1:2))>epsilon || norm(Hval)>epsilon) && k < maxIterations)
    [qpSol, dx, muNew] = quadraticProgram(lagHessxx,fGrad, Hgrad', -Hval');
%     [qpSol, dx, muNew] = quadraticProgram(lagHessxx,fGrad, Hgrad', -Hval', muOld);
%     dx = dx/norm(dx);
%     qpSol = qpSol/norm(qpSol);
%     dx = qpSol(1:n);
%     muNew = qpSol(n+1:end);
    
%     if norm(dx) < epsilon
%         x = x+dx;
%         mu = muNew;
%         break;
%     end
    dMu = muNew-muOld;
%     dMu = muNew;
    
    l = 0;
    t = 0;
%     [~,fgrad] = lag_func([x;muOld]);
    [~,fgrad] = quadPen(x);
    
%     while(lag_func([x;muOld]+(beta^l)*[dx;dMu])>lag_func([x;muOld])+(beta^l)*sigma*fgrad'*[dx;dMu])
    
    while(quadPen(x+alphaConst*(beta^l)*dx)>quadPen(x)+alphaConst*(beta^l)*sigma*fgrad'*dx)
        l = l+1;
        if l > 50
            break
        end
    end
%     t = alphaConst*beta^l;
    t = beta^l;
    
%     alpha = t;
%     alpha = Wolfe_Powll_rule(quadPen,x,dx);
    alpha = 1;
    
    
%     alpha = 
    mu = muOld+alpha*dMu;
    % muOld = mu;
    x = x + alpha*dx;    
    % x = x + dx;    
    [fVal, fGrad, fHess] = func(x);
    lagGradx = fGrad;
    lagHessxx = fHess;
    Hval = NaN(1,p);    % matrix of constraints values
    Hgrad = NaN(n,p);   % matrix of constraints' gradients
    for i = 1:p
        h = H{i};
        [hVal, hGrad, hHess] = h(x);
        H_cell{i}.hVal = hVal;
        H_cell{i}.hGrad = hGrad;
        H_cell{i}.hHess = hHess;

        Hval(i) = hVal;
        Hgrad(:,i) = hGrad;

        lagGradx = lagGradx + hGrad*mu(i);  % gradient of Lagrangian w.r.t. x
        lagHessxx = lagHessxx + hHess*mu(i);    
    end
    
    xArr = [xArr; x'];
    fArr = [fArr; fVal'];
    fGradArr = [fGradArr;norm(fGrad)];

    k = k + 1;
    if k == maxIterations-1
        disp('Max iterations reached');
    end
end
% [a,b] = lagFunc(argmin,mu,func,H);
mu = -mu;
argmin = x;
iterations = k;
minval = fVal;

    
end

function [lagVal, lagGrad] = lagFunc(x, mu, func, H)        
        [fVal,fGrad] = func(x);
        hValSum = 0;
        hGradSum = 0;
        hValArr = NaN(length(H),1);
        for i=1:length(H)
            h_i = H{i};
            [hVal,hGrad] = h_i(x);
            hValSum = hValSum + hVal*mu(i);
            hGradSum= hGradSum + hGrad*mu(i);
            hValArr(i) = hVal;
        end
        
        lagVal = fVal + hValSum;
        lagGrad = [fGrad + hGradSum; hValArr];
    end

