close all
clear all
clc


format long e
fileID = fopen('results_1d_burger_p0_sin_test_just_burger.dat','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

x=A(:,1);
velp0=A(:,2);

fileID = fopen('results_1d_burger_p1_sin_test_just_burger.dat','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

x=A(:,1);
velp1=A(:,2);

fileID = fopen('results_1d_burger_p2_sin_test_just_burger.dat','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
A=A';

x2=A(:,1);
velp2=A(:,2);

velini=sin(x);

format long e

%--------------------------------------------------------------
% exact solution using newton raphson method

sizex=size(x);
Npoint=sizex(1)
t=2;
Niter=100000;
tol=1e-8;
fx=0.1;
fdashx=0.2;
sizeexact=1;
flag=0;

p0norm=0;
p1norm=0;
p2norm=0;


for ipoint=1:Npoint
   
    zeta=asin(velp0(ipoint));
    
    for i=1:Niter
    fx=(sin(zeta)*t+zeta-x(ipoint));
    fdashx=cos(zeta)*t+1;
    
    zeta=zeta-fx/fdashx;
    
format long e
    if(abs(fx/fdashx)<tol&&fx<tol)
        fprintf('Ipoint %d iter %d xval %f residue %f \n',ipoint,i,x(ipoint),fx/fdashx);
        break;
    end
    
    if(i==Niter)
       fprintf('solution did not converge for x= %f \n',x(ipoint)); 
       flag=1;
    break
    end
    
    end
    if(flag==0)
        if(x(ipoint)>pi&&sin(zeta)>0)
            continue;
        end
        xexact(sizeexact)=x(ipoint);
        velexact(sizeexact)=sin(zeta);
        p0norm=p0norm+(velp0(ipoint)-velexact(sizeexact))^2;
        p1norm=p1norm+(velp1(ipoint)-velexact(sizeexact))^2;
        p2norm=p2norm+(velp2(ipoint)-velexact(sizeexact))^2;
        sizeexact=sizeexact+1;
    end
    flag=0;
end

p0norm=sqrt(p0norm/(sizeexact-1))
p1norm=sqrt(p1norm/(sizeexact-1))
p2norm=sqrt(p2norm/(sizeexact-1))


fig=figure('Name',['Solution of burger equation at t=',num2str(t),' with N=32']);
title(['Solution of burger equation at t=',num2str(t),' with N=32']);
xlabel('x')
ylabel('u')
hold on;
%-----------------------
plot(x,velini,'color','k')
hold on
plot(x,velp0,'-gs','color','k','Linewidth',1,'MarkerSize',10)
hold on
plot(x,velp1,'-o','color','k','Linewidth',1,'MarkerSize',10)
hold on
plot(x2,velp2,'-x','color','k','Linewidth',1,'MarkerSize',10)
hold on
plot(xexact,velexact,'*','color','k','Linewidth',1,'MarkerSize',10)
hold on
%---------------------
% plot(x,velini,'Linewidth',2)
% hold on
% plot(x,velp0,'Linewidth',2)
% hold on
% plot(x,velp1,'Linewidth',2)
% hold on
% plot(x2,velp2,'Linewidth',2)
% hold on
% plot(xexact,velexact,'*','Linewidth',1)
% hold on
%-------------------------------

axis([min(x) max(x) -1.2 1.5])
grid on

hleg1 = legend('Initial','DGp0','DGP1','DGP2',['Exact at t=',num2str(t)]);
%hleg1 = legend('Initial',['DGp0 L2norm=',num2str(p0norm)],['DGP1 L2norm=',num2str(p1norm)],['Exact at t=',num2str(t)]);
set(hleg1,'Location','NorthEast')

saveas(fig,'burger_sin_t_2.bmp')

