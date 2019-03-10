% eqConstCellArr is a cell array of equality constraints.
function [val, grad, hess] = quadraticPenaltyFun(x,alpha,objFunc, eqConstCellArr)

%% evaluate function
m = length(eqConstCellArr);
% constCellMatrix = cell{3,m};    % first row is val, second is gradient, third is hessian
constValSquareSum = 0;
constValGradSum = 0;
constGradOuterProductSum = 0;
constValHessSum = 0;
for i=1:m
    [constVal, constGrad, constHess] = eqConstCellArr{i}(x);
%     constCellMatrix{1,i} = constVal;
%     constCellMatrix{2,i} = constGrad;
%     constCellMatrix{3,i} = constHess;     
    
    constValSquareSum = constValSquareSum + constVal^2;
    constValGradSum = constValGradSum + constVal*constGrad;
    constGradOuterProductSum = constGradOuterProductSum + constGrad*constGrad';
    constValHessSum = constValHessSum + constVal*constHess;
end

[objFuncVal, objFuncGrad, objFuncHess] = objFunc(x);

val = objFuncVal + (alpha/2)*constValSquareSum;
grad = objFuncGrad + alpha*constValGradSum;
hess = objFuncHess + alpha*(constGradOuterProductSum + constValHessSum);
end

