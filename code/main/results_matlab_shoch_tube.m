clear all
%close all
clc
gama=1.4;

fileID = fopen('pres.txt','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

xpres=A(:,1);
presexact=A(:,2);

fileID = fopen('rho.txt','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

xrho=A(:,1);
rhoexact=A(:,2);

fileID = fopen('vel.txt','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

xrho=A(:,1);
velexact=A(:,2)./sqrt(gama*presexact./rhoexact);


fileID = fopen('results.dat','r');
formatSpec = '%lf %lf %lf %lf';
sizeA = [4 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

x=A(:,1);
rho=A(:,2);
pres=A(:,4);
vel=A(:,3)./sqrt(gama*pres./rho);

figure(1)
%subplot(2,2,1)
xlabel('x')
ylabel('rho')
hold on
plot(x,rho,'-o','Linewidth',2);
hold on
plot(xrho,rhoexact,'Linewidth',2);
grid on
hold on
hleg1 = legend('DGP1','Exact');
set(hleg1,'Location','SouthWest')





%subplot(2,2,2)
figure(2)
xlabel('x')
ylabel('Mach Number')
hold on
plot(x,vel,'-o','Linewidth',2);
hold on
plot(xrho,velexact,'Linewidth',2);
grid on
hold on
% hleg1 = legend('DGP1','Exact');
% set(hleg1,'Location','NorthWest')


% subplot(2,2,3)
figure(3)
xlabel('x')
ylabel('pressure')
hold on
plot(x,pres,'-o','Linewidth',2);
hold on
plot(xrho,presexact,'Linewidth',2);
grid on
hold on
% hleg1 = legend('DGP1','Exact');
% set(hleg1,'Location','NorthWest')
