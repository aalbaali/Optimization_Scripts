%Implementation of Wolfe-Powell rule from Alg1 in A6.
function [ti] = Wolfe_Powll_rule(func, x, d, t0, gamma, rho, sigma)

% phi = @(t) func(x+t*d);
% 
% psi = @(t) phi(t)-phi(0)- sigma*t*phi_dot(0);
% 
% function phi_dot_t = phi_dot(t)    
%     [fval fgrad] = func(x+t*d);
%     phi_dot_t = fgrad'*d;
% end
% 
% ti = t0;
% t = phaseA(t0);
% return;
% 
%     function t = phaseA(t)
%         if (psi(ti) >=0)
%             t = phaseB(0,ti, ti);
%             return;
%         else
%             if (phi_dot(ti) >= rho*phi_dot(0))
%                 t = ti;
%                 return;
%             else
%                 ti = gamma*ti;
%                 t = phaseA(ti);
%                 return;
%             end
%         end
%     end
% 
% 
%     function t = phaseB(a,b,t)
%     aj = a;
%     bj = b;
%     tj = aj+(bj-aj)/2;
% 
%         if (psi(tj) < 0)
%             if (phi_dot(tj) >= rho*phi_dot(0))
%                 t = tj;%STOP 2
%                 return;
%             else
%                 aj = tj;
%                 t = phaseB(aj,bj,tj);
%                 return;
%             end
%         else
%             aj = a;
%             bj = tj;
%             t = phaseB(aj,bj,tj);
%             return;
%         end
% 
%     end
% 
% end
if (t0<=0 || gamma <=1)
    disp('error with the parameters');
    return;
end

[fval1 fgrad1]=func(x);
if(fgrad1'*d>=0)
    disp('error!!!');
    return;
end

f = func;
phi = @(t) f(x+t*d);

psi = @(t) phi(t)-phi(0)- sigma*t*phi_dot(0);



ti = t0;
% psi_i = psi(ti); %psi at t_i; this is to reduce the number of computations (it's computed 3 times in the loop)
% phi_dot_i = phi_dot(ti); %phi dot at t_i
phi_dot0 = phi_dot(0); %phi dot at t = 0;
i = 0;
while(1)
    
     if(rem(i,100)==0)
            i;
     end
        
     
    psi(ti);    
    if ( psi(ti) >=0)
        a = 0;
        b = ti;        
        break;
    else
        if (phi_dot(ti)>=rho*phi_dot0)
            t = ti;
            return;     %STOP 1
        else
            ti = gamma*ti;
        end
    end
    
%     psi_i = psi(ti); %psi at t_i; this is to reduce the number of computations (it's computed 3 times in the loop)    
%     phi_dot_i = phi_dot(ti); 
i = i+1;

end
aj = a;
bj = b;
% tj = aj+(bj-aj)/2;


% disp(tj);

% psi_j = psi(tj); %psi at t_i; this is to reduce the number of computations (it's computed 3 times in the loop)
% phi_dot_j = phi_dot(tj); %phi dot at t_i
j = 0;
if (psi(ti) >=0)
    while(j<1000)
        if(rem(j,100)==0)
            j;
        end
        tj = aj+(bj-aj)/2;   %B1
        
        
        if (psi(tj) >= 0)             
            bj = tj;
%             tj = aj+(bj-aj)/2;
%             [fval fgrad] = f(x+tj*d);
%             phi_dot_t = fgrad'*d;
%             disp(bj)s;             
            continue; %go to B1
        else  %psi(tj)<0
            if(phi_dot(tj) >= rho*phi_dot0)
                t = tj;
                return;     %STOP 2
            else %psi(tj)<0 and phi_dot(tj)< 0
                aj = tj;
%                 tj = aj+(bj-aj)/2;
                continue; %go to B1
            end
        end    
        j = j+1;
    end
else
    disp('something wrong');
end


%phi doet used in the psi function
function phi_dot_t = phi_dot(t)    
    [fval fgrad] = f(x+t*d);
    phi_dot_t = fgrad'*d;    
end

j = j+1;


end



