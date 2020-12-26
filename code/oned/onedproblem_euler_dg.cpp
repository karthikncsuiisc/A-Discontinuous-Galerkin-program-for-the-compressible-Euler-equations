//Function to Solve 1d Euler equation
// The cells are of equal length. The function gives the results as an output file
//
//-------------------------------------------------------------------------------------------------
// Cells discrption
//--------------------------------------------------------------------------------------------------
//	  Ghost								domain cells									  Ghost cells	
//	  0	  1	   2   3     4    5     6    ..................          nelem  .   nelem+1		  nelem+2         
//	|---|---|____|____|____|_____|_____|_____________________|______|_____|_____|_____|-----|-----|	
//	0   1   2	 3    4    5     6     7 .................        nelem      nelem+1     nelem+4	
//	Points location
//--------------------------------------------------------------------------------------------------
#include "functions.h"
#include "oned_header.h"

void onedproblem_euler_dg()
{
double deltax=1.0/(double)nelem;

cout<<"solving euler equation"<<endl;
//Adding ghost cells after calcualting lenght of element. 
int nelemg=nelem+4;
npoint_1d=nelem+5;
	
coord_1d=new double[npoint_1d];

// Defining the mesh for 1d problem
for(int i=0; i<npoint_1d;i++)
{
coord_1d[i]=leftcoord+(i-2)*deltax;
//cout<<i<<" "<<coord_1d[i]<<endl;
}

double deltaT;

double phytime=t0;

flux=new double* [npoint_1d];
for(int i=0;i<npoint_1d;i++) flux[i]=new double[3];

rhspo= new double** [nelemg];
for(int i=0;i<nelemg; i++) rhspo[i]=new double* [3];
for(int i=0;i<nelemg; i++)
for(int j=0;j<3;j++)
rhspo[i][j]=new double[ndegr];

unkel= new double** [nelemg];
for(int i=0;i<nelemg; i++) unkel[i]=new double* [3];
for(int i=0;i<nelemg; i++)
for(int j=0;j<3;j++)
unkel[i][j]=new double[ndegr];

double ***unkeln;
unkeln= new double** [nelemg];
for(int i=0;i<nelemg; i++) unkeln[i]=new double* [3];
for(int i=0;i<nelemg; i++)
for(int j=0;j<3;j++)
unkeln[i][j]=new double[ndegr];

if(recons==1)
{
	switch(disc_order)
	{
		case 0: 
				switch(flux_type)
				{
					case 0: cout<<"reconstruction is not implemented for solving euler equations, P0 method using stegar warming flux"<<endl; exit(1);
					case 1: if(limiter==0){ cout<<"It is required to use limiter to solve reconstructed euler equation"<<endl; exit(2);}
							cout<<"Solving Euler equations using P0 reconstruction to p1 and with Vanleer flux"<<endl; break; 
					case 2: cout<<"reconstruction is not implemented for solving euler equations, P0 method using Gudnov flux"<<endl; exit(3);
				}
		break;
		case 1: cout<<"reconstruction is not implemented for solving euler equations using P1 method"<<endl; exit(4);
		case 2: cout<<"reconstruction is not implemented for solving euler equations using P2 method"<<endl; exit(5);
		default: cout<<"wrong discretization order"<<endl; exit(1);
	}
}

double M[3];

M[0]=deltax;
M[1]=deltax/3.0;
M[2]=deltax/45.0;

initial_condition_1d();

while(phytime<tf)
{

double maxval=-1e10;
double maxdt=1e10;
double rho,vel,pres,c0;
double u1,u2,u3;
for(int i=2;i<=nelem+1;i++) 
{
	maxval=max(maxval,abs(unkel[i][0][0]));
	
	u1=unkel[i][0][0];
	u2=unkel[i][1][0];
	u3=unkel[i][2][0];

	rho=u1;
	vel=u2/u1;
	pres=(gama-1)*(u3-0.5*pow(u2,2)/u1);
	c0=sqrt(gama*max(pres,hrd_lmtr_cin)/rho);
	maxdt=min(maxdt,deltax/(abs(vel)+c0));
}
if(CFL_use==0) deltaT=dt;
else deltaT=CFL*maxdt;

if(phytime+deltaT>tf) deltaT=tf-phytime;

phytime=phytime+deltaT;

cout<<"-------------------------------------------"<<endl;
cout<<"Time: "<<phytime<<" Time step: "<<deltaT<<endl;
cout<<"Max value density: "<<maxval<<endl;


if(phytype==0) periodicboundary_1d();
else if(phytype==1) noflux_1d();
else {cout<<"wrong physics choice"<<endl; exit(1);}

for(int i=2;i<=nelem+1;i++)
for(int ivar=0;ivar<3;ivar++)
for(int degree=0;degree<ndegr; degree++)
	unkeln[i][ivar][degree]=unkel[i][ivar][degree];

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage1
//-----------------------------------------------------------------------
if(phytype==0) periodicboundary_1d();
else if(phytype==1) noflux_1d();

if(flux_type==0) stegar_warming_flux_1d();
else if(flux_type==1) van_leer_FVS_1d();
else if(flux_type==2) godunov_method_1d(phytime);
else {cout<<"wrong chhoice of flux method"<<endl; exit(1);}

if(ndegr==1) dgp0_fun_1d();
else if(ndegr==2) dgp1_fun_1d();
else if(ndegr==3) dgp2_fun_1d();

for(int i=2;i<=nelem+1;i++)
for(int nvar=0;nvar<3;nvar++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][nvar][degree]=unkeln[i][nvar][degree]+deltaT/M[degree]*rhspo[i][nvar][degree];

if(timedep_method==2)
{
	cout<<"stage 2"<<endl;
//-----------------------------------------------------------------------
//							Two Stage RK Method: Stage2
//-----------------------------------------------------------------------
if(phytype==0) periodicboundary_1d();
else if(phytype==1) noflux_1d();

if(flux_type==0) stegar_warming_flux_1d();
else if(flux_type==1) van_leer_FVS_1d();
else if(flux_type==2) godunov_method_1d(phytime);

if(ndegr==1) dgp0_fun_1d();
else if(ndegr==2) dgp1_fun_1d();
else if(ndegr==3) dgp2_fun_1d();

for(int i=2;i<=nelem+1;i++)
for(int nvar=0;nvar<3;nvar++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][nvar][degree]=0.5*unkeln[i][nvar][degree]+0.5*(unkel[i][nvar][degree]+deltaT/M[degree]*rhspo[i][nvar][degree]);
}

else if(timedep_method==3)
{	
//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage2
//-----------------------------------------------------------------------
if(phytype==0) periodicboundary_1d();
else if(phytype==1) noflux_1d();

if(flux_type==0) stegar_warming_flux_1d();
else if(flux_type==1) van_leer_FVS_1d();
else if(flux_type==2) godunov_method_1d(phytime);

if(ndegr==1) dgp0_fun_1d();
else if(ndegr==2) dgp1_fun_1d();
else if(ndegr==3) dgp2_fun_1d();

for(int i=2;i<=nelem+1;i++)
for(int nvar=0;nvar<3;nvar++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][nvar][degree]=0.75*unkeln[i][nvar][degree]+0.25*(unkel[i][nvar][degree]+deltaT/M[degree]*rhspo[i][nvar][degree]);
	
//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage3
//-----------------------------------------------------------------------
if(phytype==0) periodicboundary_1d();
else if(phytype==1) noflux_1d();

if(flux_type==0) stegar_warming_flux_1d();
else if(flux_type==1) van_leer_FVS_1d();
else if(flux_type==2) godunov_method_1d(phytime);

if(ndegr==1) dgp0_fun_1d();
else if(ndegr==2) dgp1_fun_1d();
else if(ndegr==3) dgp2_fun_1d();

for(int i=2;i<=nelem+1;i++)
for(int nvar=0;nvar<3;nvar++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][nvar][degree]=(1.0/3.0)*unkeln[i][nvar][degree]+(2/3.0)*(unkel[i][nvar][degree]+deltaT/M[degree]*rhspo[i][nvar][degree]);

//-----------------------------------------------------------------------
//							END of Three Stage RK Method
//-----------------------------------------------------------------------
}

}

onedresults();

}
