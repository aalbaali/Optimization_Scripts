%Amro Al Baali
%ECSE 507, Winter 2018, Assignment 3 Question 3.6
%January 31, 2018
%Implelemnting the gradient method on the Rosenbrock function

clear all;

home;

[X Y] = meshgrid(-1.5:0.01:1.5,0.6:0.01:1.3);

Z = 100*(Y-X.^2).^2+(1-X.^2).^2;
N = 30; %number of level sets

figure
contour(X,Y,Z,N);
xlabel('$x_1$','Interpreter','latex','fontsize',14);
ylabel('$x_2$','Interpreter','latex','fontsize',14);
% title('Contour plots for the Rosenbrock function','interpreter','latex');
grid on;
hold on;

x0 = [-1.2;1];
Beta = 0.5;
sigma = 1e-4;
eps = 1e-4;

[val, grad]=Rosenbrock(x0);
x=x0;

k=0;
scatter(x(1),x(2),'r');
grid on;
while norm(grad)>eps
    [val, grad]=Rosenbrock(x);
    d = -grad;
    l=0; %L
    while Rosenbrock(x+Beta^l*d)>Rosenbrock(x)+Beta^l*sigma*grad'*d
        l=l+1;
    end
    t = Beta^l;
    
    x = x+t*d;
    k=k+1;
    
    line(x(1),x(2));
    scatter(x(1),x(2),'b','x');    
end
disp('x: ');    disp(x);
disp('Number of iterations: ');     disp(k);