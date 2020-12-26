//Program to solve godunov_flux
//--------------------------------------------------------------------------------------------------
//	  Ghost								domain cells									  Ghost cells	
//	  0	  1	   2   3     4    5     6    ..................          nelem  .   nelem+1		  nelem+2         
//	|---|---|____|____|____|_____|_____|_____________________|______|_____|_____|_____|-----|-----|	
//	0   1   2	 3    4    5     6     7 .................        nelem      nelem+1     nelem+4	
//	Points location
//--------------------------------------------------------------------------------------------------

#include "functions.h"
#include "oned_header.h"

void godunov_method_1d(double tt)
{
double pressol,velsol,rhosol;
double gama=1.4;
int max_iterations=0;
	
for(int node=2;node<=nelem+2;node++)
{
	double u1,u2,u3;
	double rhol,vell,presl,rhor,velr,presr;
	
	
	//Fi+ calcualtion
    if(disc_order==0)
    {
	u1=unkel[node-1][0][0];
	u2=unkel[node-1][1][0];
	u3=unkel[node-1][2][0];
    }
    else if(disc_order==1)
    {
	u1=unkel[node-1][0][0]+unkel[node-1][0][1];
	u2=unkel[node-1][1][0]+unkel[node-1][1][1];
	u3=unkel[node-1][2][0]+unkel[node-1][2][1];
    }
    else if(disc_order==2)
    {
	u1=unkel[node-1][0][0]+unkel[node-1][0][1]+unkel[node-1][0][2]/3.0;
	u2=unkel[node-1][1][0]+unkel[node-1][1][1]+unkel[node-1][1][2]/3.0;
	u3=unkel[node-1][2][0]+unkel[node-1][2][1]+unkel[node-1][2][2]/3.0;
    }

	rhol=u1;
	vell=u2/u1;
//	presl=max((gama-1)*(u3-0.5*pow(u2,2)/u1),1e-8);
	presl=(gama-1)*(u3-0.5*pow(u2,2)/u1);
	
	if(disc_order==0)
	{
	u1=unkel[node][0][0];
	u2=unkel[node][1][0];
	u3=unkel[node][2][0];
    }
    else if(disc_order==1)
	{
	u1=unkel[node][0][0]-unkel[node][0][1];
	u2=unkel[node][1][0]-unkel[node][1][1];
	u3=unkel[node][2][0]-unkel[node][2][1];
	
    }
    else if(disc_order==2)
	{
	u1=unkel[node][0][0]-unkel[node][0][1]+unkel[node][0][2]/3.0;
	u2=unkel[node][1][0]-unkel[node][1][1]+unkel[node][1][2]/3.0;
	u3=unkel[node][2][0]-unkel[node][2][1]+unkel[node][2][2]/3.0;
    }

	rhor=u1;
	velr=u2/u1;
//	presr=max((gama-1)*(u3-0.5*pow(u2,2)/u1),1e-8);
	presr=(gama-1)*(u3-0.5*pow(u2,2)/u1);
	
//-----------------------------------------------------------------------------
//	Riemann solver
//-----------------------------------------------------------------------------

double cl,cr;
double Ksol,Kold;
double a1,a2,a3,a4,a5;
double tol;
double fdashx,diffint,dummy,fx;
	
//cl=sqrt(gama*presl/rhol);
//cr=sqrt(gama*presr/rhor);

cl=sqrt(gama*max(presl,hrd_lmtr_cin)/rhol);
cr=sqrt(gama*max(presr,hrd_lmtr_cin)/rhor);

Ksol=(presl+presr)/(2*presr); 					//P2/Pr
Kold=Ksol;
long long int Niter=20000000;
int iter;
tol=1e-6;
a2=(gama-1)/(2*cl);
a1=presl/presr;
a3=cr/gama;
a4=(gama+1)/(2*gama);
a5=(-2*gama)/(gama-1);

for(iter=0;iter<Niter;iter++)
{
fdashx=0;
diffint=0;

dummy=1+a2*(vell-velr-a3*(Ksol-1)/sqrt((a4*(Ksol-1)+1)));
fx=Ksol*pow(dummy,a5)-a1;
diffint=1/sqrt(a4*(Ksol-1)+1);
diffint=diffint-0.5*(Ksol-1)*a4/(pow(((a4*(Ksol-1)+1)),1.5));
fdashx=Ksol*a5*pow(dummy,(a5-1));
fdashx=fdashx*(-a3)*(a2)*diffint;
fdashx=fdashx+pow(dummy,a5);
Ksol=Kold-fx/fdashx;
if(abs(Ksol-Kold)<tol) break;
Kold=Ksol;

if(iter==Niter-1)
{
	cout<<"Maximum Number of iterations reached for reimann solver"<<endl;
	cout<<"Number of iterations: "<<iter<<endl;
	cout<<"Node number and coordinate"<<node<<" "<<Ksol-Kold<<endl;
	cout<<"Density left :"<<rhol<<endl;
	cout<<"Density right :"<<rhor<<endl;
	cout<<"pressure left :"<<presl<<endl;
	cout<<"pressure right :"<<presr<<endl;
	cout<<"velocity left :"<<vell<<endl;
	cout<<"velocity right :"<<velr<<endl;
	
onedresults();	
exit(1);
}

}

	max_iterations=max(max_iterations,iter);

double p2,rho2,v2,p3,v3,rho3,vs,c3;
double xcenter=0,alpha;
double xh,xt,xc,xs;

alpha=(gama+1)/(gama-1);
p2=Ksol*presr;
v2=vell+2*cl*(1-pow(p2/presl,(gama-1)/(2*gama)))/(gama-1);
v3=v2;
p3=p2;
vs=velr+(Ksol-1)*cr*cr/(gama*v2-gama*velr+1e-10);
rho3=rhol*pow((p3/presl),(1/gama));
c3=sqrt(gama*p3/rho3);

xh=(vell-cl)*tt;
xt=(v3-c3)*tt;
xc=v2*tt;
xs=vs*tt;

if(xcenter<xh)
{
	velsol=vell;
	pressol=presl;
	rhosol=rhol;
}
else if(xcenter<=xt&&xcenter>=xh)
{
	velsol=(2.0)/(gama-1)*(cl+(gama-1)/2.0*vell);
	pressol=presl*pow(velsol/cl,2*gama/(gama-1));
	rhosol=rhol*pow(velsol/cl,2.0/(gama-1));
}
else if(xcenter<=xc && xcenter>=xt)
{
	velsol=v2;
	pressol=p2;
	rhosol=rhol*pow(p3/presl,1.0/gama);
} 
else if(xcenter>=xc&&xcenter<=xs)
{
	velsol=v2;
	pressol=p2;
	rhosol=rhor*(1+alpha*Ksol)/(alpha+Ksol);
}
else if(xcenter>xs)
{
	velsol=velr;
	pressol=presr;
	rhosol=rhor;
}


//--------------------------------------------------------------------------
// End of Riemann solver
//---------------------------------------------------------------------------

	double rhoE3;
	rhoE3=pressol/(gama-1)+0.5*rhosol*velsol*velsol;
	
	flux[node][0]=rhosol*velsol;		
	flux[node][1]=rhosol*velsol*velsol+pressol;
	flux[node][2]=velsol*(rhoE3+pressol);


}

//cout<<"Maximum Iterations:"<<max_iterations<<endl;

	
	
	
}
