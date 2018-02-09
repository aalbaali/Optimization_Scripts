% close all;
% clear all; 
home;

[X Y] = meshgrid(-1.5:0.01:1.1,-0.1:0.01:1.1);

Z = 100*(Y-X.^2).^2+(1-X.^2).^2;
N = 30; %number of level sets

figure
contour(X,Y,Z,N);
xlabel('$x_1$','Interpreter','latex','fontsize',14);
ylabel('$x_2$','Interpreter','latex','fontsize',14);
% title('Contour plots for the Rosenbrock function','interpreter','latex');
grid on;