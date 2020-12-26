//Function to Solve 1d advection equation
// The cells are of equal length. The function gives the results as an output file: results_1d_wave.dat
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

void onedproblem_wave_dg()
{
double deltax=1.0/(double)nelem;

//----------------------------------------------------
//Adding ghost cells after calcualting lenght of element. 
//---------------------------------------------------
int nelemg=nelem+4;
npoint_1d=nelem+5;
	
coord_1d=new double[npoint_1d];

//------------------------------------------------------
// Defining the mesh for 1d problem
//------------------------------------------------------
for(int i=0; i<npoint_1d;i++)
{
coord_1d[i]=leftcoord+(i-2)*deltax;
cout<<i<<" "<<coord_1d[i]<<endl;
}

//------------------------------------------------------
// Defninig the time step
//------------------------------------------------------
double deltaT;

if(CFL_use==0) deltaT=dt;
else deltaT=CFL*deltax/abs(a0x);

double phytime=t0;

//-----------------------------------------------------
// Initiating the variables
//-----------------------------------------------------
flux=new double* [npoint_1d];
for(int i=0;i<npoint_1d;i++) flux[i]=new double[1];

rhspo= new double** [nelemg];
for(int i=0;i<nelemg; i++) rhspo[i]=new double* [1];
for(int i=0;i<nelemg; i++) rhspo[i][0]=new double[ndegr];

unkel= new double** [nelemg];
for(int i=0;i<nelemg; i++) unkel[i]=new double* [1];
for(int i=0;i<nelemg; i++) unkel[i][0]=new double[ndegr];

double ***unkeln;
unkeln= new double** [nelemg];
for(int i=0;i<nelemg; i++) unkeln[i]=new double* [1];
for(int i=0;i<nelemg; i++) unkeln[i][0]=new double[ndegr];

double M[3];

M[0]=deltax;
M[1]=deltax/3.0;
M[2]=deltax/45.0;

initial_condition_1d_wave();

while(phytime<=tf)
{

phytime=phytime+deltaT;

double maxval=-1e10;
for(int i=2;i<=nelem+1;i++) maxval=max(maxval,abs(unkel[i][0][0]));
cout<<"-------------------------------------------"<<endl;
cout<<"Time: "<<phytime<<endl;
cout<<"Max value density: "<<maxval<<endl;

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

for(int i=0;i<=nelem+3;i++)
for(int degree=0;degree<ndegr; degree++)
unkeln[i][0][degree]=unkel[i][0][degree];

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage1
//-----------------------------------------------------------------------
cout<<"stage1"<<endl;
if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

rieman_solver_oned_wave();

	for(int i=2;i<=nelem+1;i++)
	{
	if(ndegr>=1) rhspo[i][0][0]=-(flux[i+1][0]-flux[i][0]);	
    if(ndegr>=2) rhspo[i][0][1]=2.0*a0x*unkel[i][0][0]-(flux[i+1][0]+flux[i][0]);	
    if(ndegr==3) rhspo[i][0][2]=2.0*a0x*unkel[i][0][1]/3.0-(flux[i+1][0]-flux[i][0])/3.0;	
	}

for(int i=2;i<=nelem+1;i++)
for(int degree=0;degree<3; degree++)
	unkel[i][0][degree]=unkeln[i][0][degree]+deltaT/M[degree]*rhspo[i][0][degree];

//-----------------------------------------------------------------------
//							Two Stage RK Method: Stage2
//-----------------------------------------------------------------------
if(timedep_method==2)
{
cout<<"stage2"<<endl;
if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

rieman_solver_oned_wave();

	for(int i=2;i<=nelem+1;i++)
	{
	if(ndegr>=1) rhspo[i][0][0]=-(flux[i+1][0]-flux[i][0]);	
    if(ndegr>=2) rhspo[i][0][1]=2.0*a0x*unkel[i][0][0]-(flux[i+1][0]+flux[i][0]);	
    if(ndegr==3) rhspo[i][0][2]=2.0*a0x*unkel[i][0][1]/3.0-(flux[i+1][0]-flux[i][0])/3.0;	
	}

for(int i=2;i<=nelem+1;i++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][0][degree]=0.5*unkeln[i][0][degree]+0.5*(unkel[i][0][degree]+deltaT/M[degree]*rhspo[i][0][degree]);
}

else if(timedep_method==3)
{
//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage2
//-----------------------------------------------------------------------
cout<<"stage2"<<endl;

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

rieman_solver_oned_wave();

	for(int i=2;i<=nelem+1;i++)
	{
	if(ndegr>=1) rhspo[i][0][0]=-(flux[i+1][0]-flux[i][0]);	
    if(ndegr>=2) rhspo[i][0][1]=2.0*a0x*unkel[i][0][0]-(flux[i+1][0]+flux[i][0]);	
    if(ndegr==3) rhspo[i][0][2]=2.0*a0x*unkel[i][0][1]/3.0-(flux[i+1][0]-flux[i][0])/3.0;	
	}

for(int i=2;i<=nelem+1;i++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][0][degree]=0.75*unkeln[i][0][degree]+0.25*(unkel[i][0][degree]+deltaT/M[degree]*rhspo[i][0][degree]);

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage3
//-----------------------------------------------------------------------
cout<<"stage3"<<endl;

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

rieman_solver_oned_wave();

	for(int i=2;i<=nelem+1;i++)
	{
	if(ndegr>=1) rhspo[i][0][0]=-(flux[i+1][0]-flux[i][0]);	
    if(ndegr>=2) rhspo[i][0][1]=2.0*a0x*unkel[i][0][0]-(flux[i+1][0]+flux[i][0]);	
    if(ndegr==3) rhspo[i][0][2]=2.0*a0x*unkel[i][0][1]/3.0-(flux[i+1][0]-flux[i][0])/3.0;	
	}

for(int i=2;i<=nelem+1;i++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][0][degree]=(1.0/3.0)*unkeln[i][0][degree]+(2/3.0)*(unkel[i][0][degree]+deltaT/M[degree]*rhspo[i][0][degree]);

//-----------------------------------------------------------------------
//							END of Three Stage RK Method
//-----------------------------------------------------------------------
}

}


ofstream fout;
fout.open("results_1d_wave.dat", ios::out | ios::trunc);

for(int i=2;i<=nelem+1;i++)
{
double xval=0.5*(coord_1d[i+1]+coord_1d[i]);

fout<<xval<<" "<<unkel[i][0][0]<<endl;
}

}
