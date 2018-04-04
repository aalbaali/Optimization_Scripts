function val = g(lam,mu,Lagrange,positive,Qinv,c,A,b)
    
x0p = [1 1]';

if positive    
    [~, val] = fmincon(@(x) Lagrange(x,lam,mu),x0p);
    val = -1/2*c'*Qinv*c-1/2*(-lam+A'*mu)'*Qinv*(-lam+A'*mu)-c'*Qinv*(-lam+A'*mu)+mu'*b;
else
    
%     [~, val] = fmincon(@(x) Lagrange(x,lam,mu),x0p);
%     val = -val;
%     Qinv = inv(Q);
    
%     val = -1/2*(c-lam-A'*mu)'*Qinv*(c-lam-A'*mu)+mu'*b;
    val = -1/2*c'*Qinv*c-1/2*(-lam-A'*mu)'*Qinv*(-lam-A'*mu)-c'*Qinv*(-lam-A'*mu)+mu'*b;
    val = -val;
    
end