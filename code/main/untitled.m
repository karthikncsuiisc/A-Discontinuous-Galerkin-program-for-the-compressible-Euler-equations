clear all
clc


format long e
fileID = fopen('results_1d_burger.dat','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

x=A(:,1);
velp0=A(:,2)

sizex=size(x);
Npoint=sizex(1);


format long e

fig=figure(1);
title('Comparision of solution at t=2 with DDG and N=64');
xlabel('x')
ylabel('u')
hold on;
plot(x,velp0,'color','k','Linewidth',2)
hold on

%axis([min(x) max(x) -1.1 1.1])
grid on

%saveas(fig,['solution_ddg_t_2_diff.bmp'])
