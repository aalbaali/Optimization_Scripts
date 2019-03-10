% alphaLst: list of alphas to iterate over. Ex: [1,0.1,0.001,etc];
function [argmin, iterations, minVal, xCell] = quadraticPenalty(objFunc,eqConstCellArr,x0, alphaLst,epsilon, maxIterations)

xCell = {};
iterations = 0;


penaltyFunc = @(x,alpha) quadraticPenaltyFun(x,alpha,objFunc,eqConstCellArr);
    for alpha = alphaLst
        [x0, penIterations, minVal, xPenCell] = globalNewton (@(x) penaltyFunc(x,alpha), x0);%,rho,p,beta,sigma,eps, maxIterations)
        iterations = iterations + penIterations;
        xCell = [xCell,xPenCell];        
    end
    
    argmin = xCell{end};
end
