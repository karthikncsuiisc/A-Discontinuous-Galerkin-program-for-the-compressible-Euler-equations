clear all
close all
clc
fileID = fopen('L2norm_function_1d_burger.dat','r');

format long e

formatSpec = '%lf %lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

%N=A(:,1);
h=A(1:4,1)
%h=[1 2 4 8]';
Er=A(1:4,2)

f=fit(log10(h),log10(Er),'poly1')
p1=f.p1;
p2=f.p2;
fig=figure(1)
title(['L2norm vs grid size']);
xlabel('log10(havg)')
ylabel('log10(L2norm)')
hold on
plot(log10(h),log10(Er),'*','Markersize',10)
grid on
%grid minor

h=A(5:8,1);
Er=A(5:8,2);

f=fit(log10(h),log10(Er),'poly1')
p11=f.p1;
p22=f.p2;

hold on
plot(log10(h),log10(Er),'o','Markersize',10)
hold on
plot(log10(h),p11*log10(h)+p22,'linewidth',2)
hold on

plot(log10(h),p1*log10(h)+p2,'linewidth',2)
hold on

grid on

hleg1 = legend('DDG','BR2');
%hleg1 = legend('Initial',['DGp0 L2norm=',num2str(p0norm)],['DGP1 L2norm=',num2str(p1norm)],['Exact at t=',num2str(t)]);
set(hleg1,'Location','NorthWest')
saveas(fig,['Errorfunction_diff_1d.png'])