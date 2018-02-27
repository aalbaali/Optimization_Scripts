% function alpha = mb_nocLineSearch(f,gradF,x,dir,slope0,of)
%This function originally by Mark Bangert from https://www.mathworks.com/matlabcentral/fileexchange/34835-optimization-tutorial?focused=5215429&tab=function
%Modified by Amro Al Baali

function alpha = mb_nocLineSearch(f, x, dir, t0, gamma, rho, sigma)

% implementation of strong wolfe line search
% 
% original c++ code by m.siggel@dkfz.de
% import to matlab by m.bangert@dkfz.de
% reference: nocedal: numerical optimization 3.4 line search methods
%
% f       - function handle of objective function
% gradF   - function handle of gradient
% x       - current iterate
% dir     - search direction
% slope0  - gradient of function at x in direction dir
% alphaLo - alpha at lower objective function value for zoom
% alphaHi - alpha at higher objectvie function value for zoom
% of_0
% oflo

% parameter for sufficient decrease condition
c1 = sigma;
% parameter for curvature condition
c2 = rho;

if c1 > c2
    error('c1 > c2\n');
end

alphaMax = 100;
alpha    = t0;

alpha_0  = 0;       %a
alpha_1  = alpha;   %b


phi = @(t) f(x+t*dir);

psi = @(t) phi(t)-phi(0)- sigma*t*phi_dot(0);

phi_dot = @(t) phi_dotXD(t, f, x, dir);
slope0 = phi_dot(0); 
of_x = phi(0);
of_0 = phi(0);
iter = 0;

while 1

    xc     = x+alpha_1*dir;
    [of gradF] = f(xc);
    slopec = gradF'*dir;
    
    % check if current iterate violates sufficient decrease
    if  (of > of_0 + slope0*c1*alpha_1) || ((of >= of_x ) && (iter > 0))
        % there has to be an acceptable point between alpha_0 and alpha_1
        % (because c1 > c2)
        alpha = nocZoom(f, x, dir, slope0, alpha_0, alpha_1, of_0, of_x , c1, c2);
        break;
    end
    % current iterate has sufficient decrease, but are we too close?
    if(abs(slopec) <= -c2*slope0)
        % strong wolfe fullfilled, quit
        alpha = alpha_1;
        break;
    end
    % are we behind the minimum?
    if (slopec >= 0)
        % there has to be an acceptable point between alpha_0 and alpha_1
        alpha = nocZoom(f, x, dir, slope0, alpha_1, alpha_0, of_0, of, c1, c2);
        break;
    end

    alpha_0 = alpha_1;
%     alpha_1 = min(alphaMax, alpha_1*3);
    alpha_1 = min(alphaMax, alpha_1*gamma);
    of_x = of;

    iter = iter + 1;
end

end

function alpha = nocZoom(f,x,dir,slope0,alphaLo,alphaHi,of_0,ofLo,c1,c2)
% this function is only called by mb_nocLineSearch - everything else does
% not make sense!

while 1
    
    alpha = (alphaLo+alphaHi)/2;
    xc    = x + alpha*dir;
    [of gradF] = f(xc);
    
    if of > of_0 + c1*alpha*slope0 || of >= ofLo
        % if we do not observe sufficient decrease in point alpha, we set
        % the maximum of the feasible interval to alpha
        alphaHi = alpha;       
    else
        slopec = gradF'*dir;
        % strong wolfe fullfilled?
        if abs(slopec) <= -c2*slope0
            return;
        end
        if slopec*(alphaHi-alphaLo) >= 0 % if slope positive and alphaHi > alphaLo  
            alphaHi = alphaLo;
            alphaLo = alpha;
            ofLo    = of;
        end
    end
    
end

end

function phi_dot_t = phi_dotXD(t, f, x, d)    
    [fval fgrad] = f(x+t*d);
    phi_dot_t = fgrad'*d;    
end
