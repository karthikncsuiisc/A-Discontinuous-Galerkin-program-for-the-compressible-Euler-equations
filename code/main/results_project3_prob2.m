close all
clear all
clc

format long e
fileID = fopen('results_1d_burger_initial_condition.dat','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

xini=A(:,1);
velpini=A(:,2)

fileID = fopen('results_1d_burger.dat','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

x=A(:,1);
velp0=A(:,2)

sizex=size(x);
Npoint=sizex(1);
N=size(x);

flag=1;
for i=1:N(1)
if(x(i)>=29.0&&x(i)<=31.0)
xreq(flag)=x(i);
velreq(flag)=(velp0(i));
flag=flag+1;
end
end



fig=figure(1);
title('Solution of visocous burger equaiton at t=30 sec using DDG scheme');
xlabel('x')
ylabel('u')
hold on

plot(xini,velpini,'Linewidth',2)
hold on;
plot(x,velp0,'Linewidth',2)
axis([min(x) max(x) 1.1*min(velp0) 1.1*max(velp0)])
grid on

% hleg1 = legend('Initial','current');
% set(hleg1,'Location','NorthWest')
saveas(fig,['solution_ddg_t_30_prob2_br2.png'])


fig=figure(2);
title('Solution of visocous burger equaiton at t=30 sec using DDG scheme');
xlabel('x')
ylabel('u')
hold on
plot(xreq,velreq,'color','k','Linewidth',2)
grid on
saveas(fig,['solution_ddg_t_30_prob2_2_br2.png'])

