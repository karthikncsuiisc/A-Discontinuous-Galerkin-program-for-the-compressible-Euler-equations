//Function to Solve the equaiton of type ut+(a*u+b*u^2/2)x=c*uxx
// The cells are of equal length. The function gives the results as an output file: results_1d_wave.dat
//
//-------------------------------------------------------------------------------------------------
// Cells discrption
//--------------------------------------------------------------------------------------------------
//	  Ghost								domain cells									  Ghost cells	
//	  0	  1	   2   3     4    5     6    ..................          nelem  .   nelem+1		  nelem+2         
//	|---|---|____|____|____|_____|_____|_____________________|______|_____|_____|_____|-----|-----|	
//	0   1   2	 3    4    5     6     7 .................        nelem      nelem+1     nelem+3	
//	Points location
//--------------------------------------------------------------------------------------------------
#include "functions.h"
#include "oned_header.h"

void onedproblem_burger_dg()
{
double deltax=(rightcoord-leftcoord)/(double)nelem;

	
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

if(diff_type==1)
{
	liftopr=new double** [nelemg];
	for(int i=0;i<nelemg; i++) liftopr[i]=new double* [2];
	for(int i=0;i<nelemg; i++) {liftopr[i][0]=new double[ndegr]; liftopr[i][1]=new double[ndegr];}
}	

double M[3];

M[0]=deltax;
M[1]=deltax/3.0;
M[2]=deltax/45.0;

initial_condition_1d_burger();

phytime=0;

//if(disc_order!=1&&coefc!=0) 
if(disc_order==0&&coefc!=0) 
{ 
	cout<<"diffusion equation is coded only for P1"<<endl; 
	cout<<" the value of c is "<<coefc<<endl;
	exit(1);
}

while(phytime<tf)
{

double maxval=0;

deltaT=dt;

if(CFL_use==1)
{
deltaT=1e10;

for(int i=2;i<=nelem+1;i++)
{
deltaT=min(deltaT,deltax/(abs(coefa+coefb*unkel[i][0][0])+1e-16));
deltaT=min(deltaT,deltax*deltax/(abs(coefc)+1e-16));
}
deltaT=CFL*deltaT;

}

for(int i=2;i<=nelem+1;i++)	maxval=maxval+unkel[i][0][0]*deltax;

if(phytime+deltaT>tf) deltaT=tf-phytime;

phytime=phytime+deltaT;

cout<<"-------------------------------------------"<<endl;
cout<<"Time: "<<phytime<<" Time step "<<deltaT<<endl;
cout<<"Total : "<<maxval<<endl;

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

for(int i=0;i<=nelem+3;i++)
for(int degree=0;degree<ndegr; degree++)
unkeln[i][0][degree]=unkel[i][0][degree];

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage1
//-----------------------------------------------------------------------

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

for(int i=0;i<=nelem+3;i++) for(int degree=0;degree<ndegr; degree++) rhspo[i][0][degree]=0.0;
for(int i=0;i<=nelem+4;i++) for(int degree=0;degree<ndegr; degree++) flux[i][degree]=0.0;


if(coefa!=0||coefb!=0) rhs_flux_contrb_burger();

if(diff_type==0&&coefc!=0) rhs_flux_contrb_diff_ddg();

		switch(disc_order)
		{
			case 0: 
					for(int i=2;i<=nelem+1;i++)
					rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
			break;
			
			case 1: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p1();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p1();
					if(diff_type==1&&coefc!=0) rhs_contrb_diff_br2_p1();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);					
					}
			break;

			case 2: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p2();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p2();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);
					 rhspo[i][0][2]=rhspo[i][0][2]-(flux[i+1][0]-flux[i][0])/3.0;				
					}
			break;
			
			default: cout<<"wrong choice of discretization order"<<endl; exit(1);
		}	
	
for(int i=2;i<=nelem+1;i++)
for(int degree=0;degree<ndegr; degree++)
{
	unkel[i][0][degree]=unkeln[i][0][degree]+deltaT/M[degree]*rhspo[i][0][degree];
}
//-----------------------------------------------------------------------
//							Two Stage RK Method: Stage2
//-----------------------------------------------------------------------
if(timedep_method==2)
{

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

for(int i=0;i<=nelem+3;i++) for(int degree=0;degree<ndegr; degree++) rhspo[i][0][degree]=0.0;
for(int i=0;i<=nelem+4;i++) for(int degree=0;degree<ndegr; degree++) flux[i][degree]=0.0;

if(coefa!=0||coefb!=0) rhs_flux_contrb_burger();
if(diff_type==0&&coefc!=0) rhs_flux_contrb_diff_ddg();

		switch(disc_order)
		{
			case 0: 
					for(int i=2;i<=nelem+1;i++)
					rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
			break;
			
			case 1: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p1();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p1();
					if(diff_type==1&&coefc!=0) rhs_contrb_diff_br2_p1();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);					
					}
			break;

			case 2: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p2();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p2();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);
					 rhspo[i][0][2]=rhspo[i][0][2]-(flux[i+1][0]-flux[i][0])/3.0;				
					}
			break;
			
			default: cout<<"wrong choice of discretization order"<<endl; exit(1);
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

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

for(int i=0;i<=nelem+3;i++) for(int degree=0;degree<ndegr; degree++) rhspo[i][0][degree]=0.0;
for(int i=0;i<=nelem+4;i++) for(int degree=0;degree<ndegr; degree++) flux[i][degree]=0.0;


if(coefa!=0||coefb!=0) rhs_flux_contrb_burger();
if(diff_type==0&&coefc!=0) rhs_flux_contrb_diff_ddg();

		switch(disc_order)
		{
			case 0: 
					for(int i=2;i<=nelem+1;i++)
					rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
			break;
			
			case 1: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p1();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p1();
					if(diff_type==1&&coefc!=0) rhs_contrb_diff_br2_p1();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);					
					}
			break;

			case 2: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p2();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p2();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);
					 rhspo[i][0][2]=rhspo[i][0][2]-(flux[i+1][0]-flux[i][0])/3.0;				
					}
			break;
			
			default: cout<<"wrong choice of discretization order"<<endl; exit(1);
		}	


for(int i=2;i<=nelem+1;i++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][0][degree]=0.75*unkeln[i][0][degree]+0.25*(unkel[i][0][degree]+deltaT/M[degree]*rhspo[i][0][degree]);

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage3
//-----------------------------------------------------------------------

if(bcl==2) periodicboundary_1d_wave();
else if(bcl==1) noflux_1d_wave();

for(int i=0;i<=nelem+3;i++) for(int degree=0;degree<ndegr; degree++) rhspo[i][0][degree]=0.0;
for(int i=0;i<=nelem+4;i++) for(int degree=0;degree<ndegr; degree++) flux[i][degree]=0.0;


if(coefa!=0||coefb!=0) rhs_flux_contrb_burger();
if(diff_type==0&&coefc!=0) rhs_flux_contrb_diff_ddg();

		switch(disc_order)
		{
			case 0: 
					for(int i=2;i<=nelem+1;i++)
					rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
			break;
			
			case 1: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p1();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p1();
					if(diff_type==1&&coefc!=0) rhs_contrb_diff_br2_p1();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);					
					}
			break;

			case 2: 
					if(coefa!=0||coefb!=0) rhs_contrb_burger_p2();
					if(diff_type==0&&coefc!=0) rhs_contrb_ddg_p2();
					for(int i=2;i<=nelem+1;i++)
					{	
					 rhspo[i][0][0]=rhspo[i][0][0]-(flux[i+1][0]-flux[i][0]);
					 rhspo[i][0][1]=rhspo[i][0][1]-(flux[i+1][0]+flux[i][0]);
					 rhspo[i][0][2]=rhspo[i][0][2]-(flux[i+1][0]-flux[i][0])/3.0;				
					}
			break;
			
			default: cout<<"wrong choice of discretization order"<<endl; exit(1);
		}	

for(int i=2;i<=nelem+1;i++)
for(int degree=0;degree<ndegr; degree++)
	unkel[i][0][degree]=(1.0/3.0)*unkeln[i][0][degree]+(2/3.0)*(unkel[i][0][degree]+deltaT/M[degree]*rhspo[i][0][degree]);

//-----------------------------------------------------------------------
//							END of Three Stage RK Method
//-----------------------------------------------------------------------
}

}


fstream fout;

/*
string phifile;
ostringstream convert1;
convert1 <<nelem;
phifile="results_1d_burger_"+convert1.str()+".dat";
fout.open(phifile.c_str(), ios::out | ios::trunc);
*/

fout.open("results_1d_burger.dat", ios::out | ios::trunc);

//double uleft,uright;
double unode,uexact;
double L2norm=0,utot=0;

for(int inode=2;inode<=nelem+2;inode++)
{
double xval=0.5*(coord_1d[inode]+coord_1d[inode+1]);
//uleft=unkel[inode-1][0][0]+unkel[inode-1][0][1];
//uright=unkel[inode][0][0]-unkel[inode][0][1];

unode=unkel[inode][0][0];

/*//sin wave pure diffusion problem output: comment if you are not solving this problem
uexact=exp(-coefc*phytime)*sin(xval);
L2norm=L2norm+pow(unode-uexact,2);
utot=utot+uexact*uexact;
fout<<xval<<" "<<unode<<" "<<uexact<<endl;
//sin wave diffusion problem output*/
fout<<xval<<" "<<unode<<endl;
}
fout.close();

/*//sin wave pure diffusion problem output: comment if you are not solving this problem
L2norm=sqrt(L2norm/(nelem+1));
utot=sqrt(utot/(nelem+1));
cout<<"L2norm :"<<L2norm<<" Relative norm :"<<L2norm/utot<<endl;
fout.open("L2norm_function_1d_burger", ios::app);
fout<<abs(leftcoord-rightcoord)/nelem<<" "<<L2norm<<endl;
fout.close();
//sin wave pure diffusion problem output*/


}
