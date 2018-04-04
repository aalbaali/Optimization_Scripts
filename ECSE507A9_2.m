close all;
clear all; 
home;

gamma = 1.5;



xmin = -1;
xmax = 2;
dx = 0.1;


ymin = -1;
ymax = 2;
dy = 0.1;
func = @(x) -(x(1)+1).^2-(x(2)+1).^2;

[x1 x2] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);


obj = -(x1+1).^2-(x2+1).^2;
con1 = x1.^2+x2.^2;
con2 = x1;
tan = x1+x2;

N = 30; %number of level sets

fig  = figure;
fig.Position = [640 518 580 460];
contour(x1,x2,obj,N); hold on;

[C h] = contour(x1,x2,con1,[2 2],'b');
h.LevelList=[0:0.1:2];
h.LevelStep = 0.1;

% 
% [C h] = contour(x1,x2,con2,[gamma gamma],'r');
% h.LevelList=[xmin:0.05:gamma];


[C h] = contour(x1,x2,tan,[(2) (2)],'r');
% h.LevelList=[xmin:0.05:gamma];


[C]=plot(1,1,'xk');
C.LineWidth = 6;

% contour(x1,x2,con2,[-sqrt(2)]); hold on;

xlabel('$x_1$','Interpreter','latex','fontsize',14);
ylabel('$x_2$','Interpreter','latex','fontsize',14);
% title('Contour plots for the Rosenbrock function','interpreter','latex');
grid on;
%  leg = legend('Obj func','$x_1^2+x_2^2-2\leq0$','$x_1\leq\gamma$');
 leg = legend('Obj func','$x_1^2+x_2^2-2\leq0$','$T_X(\bar{x})$');
 leg.Interpreter = 'latex';
 leg.FontSize = 12;
 leg.Position = [0.6376 0.7901 0.2557 0.1252];
 
 tit = title('Contour plot with $\gamma = 1.5$');
 tit.Interpreter = 'latex';


