% There's no MATLAB programming question in A9. 
% Hovever, this script is to test the values (minimum) of Q2

function sol = ECSE507A9main(gamma, x0)
  if nargin == 0  
    gamma = -sqrt(2);
    x0 = [0.5;0.5];
  elseif nargin == 1
      x0 = [0.5;0.5];
  end

f = @(x) -(x(1)+1).^2-(x(2)+1).^2;
% g1 = @(x) [[];[x(1)^2+x(2)^2-2]];
    function [nonleq, nonlineq] = g1(x)
        nonleq = [];
        nonlineq = x(1)^2+x(2)^2-2;
    end
g2 = @(x) x(1)-gamma;

% g2 can be re-written as an affine constraint
A = [1 0];
b = gamma;



sol = fmincon(f,x0,A,b,[],[],-Inf(2,1),[gamma;inf],@g1);


end