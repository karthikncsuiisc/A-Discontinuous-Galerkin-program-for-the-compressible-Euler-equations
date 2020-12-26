//FUnction solves advection method sing P0 method

#include "functions.h"
#include "twod_header.h"

void tvd_rk3_advection()
{
double phytime=t0;

cout<<"----------------------------------------------------------"<<endl;
cout<<"                  Solving using TVD RK3                   "<<endl;

//-----------------------------------------------------
// Initiating the variables
//-----------------------------------------------------

	
	rhspo= new double** [nelem+1];
	for(int ielem=1;ielem<=nelem; ielem++) rhspo[ielem]=new double* [nvar];
	for(int ielem=1;ielem<=nelem; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) rhspo[ielem][ivar]=new double[ndegr];
	
	unkel= new double** [nelem+1+nface];
	for(int ielem=1;ielem<=nelem+nface; ielem++) unkel[ielem]=new double* [nvar];
	for(int ielem=1;ielem<=nelem+nface; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) unkel[ielem][ivar]=new double[ndegr];

	double ***unkeln;
	unkeln= new double** [nelem+1];
	for(int ielem=1;ielem<=nelem; ielem++) unkeln[ielem]=new double* [nvar];
	for(int ielem=1;ielem<=nelem; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) unkeln[ielem][ivar]=new double[ndegr];
    
	double ***M;
	M=new double** [nelem];
	for(int ielem=1;ielem<=nelem;ielem++)
	M[ielem]=new double* [1];
	for(int ielem=1;ielem<=nelem;ielem++)
	for(int j=0;j<1;j++)
	M[ielem][j]=new double[1];

	double x1,x2,x3,y1,y2,y3,xcmax,ycmax;
	int ip1,ip2,ip3;
	for(int ielem=1; ielem<=nelem;ielem++)
	{
    ip1=inpoel[ielem][1];
    ip2=inpoel[ielem][2];
    ip3=inpoel[ielem][3];
    
    x1=coord[ip1][1];
    y1=coord[ip1][2];
    x2=coord[ip2][1];
    y2=coord[ip2][2];
    x3=coord[ip3][1];
    y3=coord[ip3][2];
    M[ielem][0][0]=triangle_area(x1,y1,x2,y2,x3,y3);
	}	 

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
unkel[ielem][ivar][degree]=0.0;

initial_conditions_2d();

//cout<<a0x<<" "<<a0y<<endl;
double deltaT=dt;
double utotal;

//-----------------------------------------------------
// Starting the time stepping
//-----------------------------------------------------


while(phytime<=tf)
{

utotal=0;

cout<<"Time: "<<phytime<<endl;

double maxval=-1e10;

for(int ielem=1;ielem<=nelem;ielem++) 
{
	if(maxval<unkel[ielem][0][0])
	{
	maxval=unkel[ielem][0][0];

	ip1=inpoel[ielem][1];
	ip2=inpoel[ielem][2];
	ip3=inpoel[ielem][3];

	x1=coord[ip1][1];
    y1=coord[ip1][2];
    x2=coord[ip2][1];
    y2=coord[ip2][2];
    x3=coord[ip3][1];
    y3=coord[ip3][2];
	xcmax=(x1+x2+x3)/3.0;
	ycmax=(y1+y2+y3)/3.0;
    }
	utotal=utotal+unkel[ielem][0][0];
}

cout<<"Max value of u1: "<<maxval<<endl;
cout<<"Total value of u1: "<<utotal<<endl;
cout<<"x(umax) :"<<xcmax<<" y(umax) :"<<ycmax<<endl;

phytime=phytime+deltaT;

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
unkeln[ielem][ivar][degree]=unkel[ielem][ivar][degree];

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage1
//-----------------------------------------------------------------------

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
rhspo[ielem][ivar][degree]=0.0;

boundary_condition_adv_2d();
flux2D_adv_calc();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	for(int degree=0;degree<ndegr; degree++)
	unkel[ielem][ivar][degree]=unkeln[ielem][ivar][degree]
							  +deltaT/M[ielem][0][0]*rhspo[ielem][ivar][degree];
/*maxval=-1e10;
for(int i=1;i<=nelem;i++) maxval=max(maxval,rhspo[i][0][0]);
cout<<"Max value of rhspo: "<<maxval<<endl;						  */
							  
	
//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage2
//-----------------------------------------------------------------------

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
rhspo[ielem][ivar][degree]=0.0;

boundary_condition_adv_2d();
flux2D_adv_calc();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	for(int degree=0;degree<ndegr; degree++)
	unkel[ielem][ivar][degree]=0.75*unkeln[ielem][ivar][degree]
							  +0.25*(unkel[ielem][ivar][degree]
							  +deltaT/M[ielem][0][0]*rhspo[ielem][ivar][degree]);

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage3
//-----------------------------------------------------------------------

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
rhspo[ielem][ivar][degree]=0.0;

boundary_condition_adv_2d();
flux2D_adv_calc();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	for(int degree=0;degree<ndegr; degree++)
	unkel[ielem][ivar][degree]=(1.0/3.0)*unkeln[ielem][ivar][degree]
							  +(2/3.0)*(unkel[ielem][ivar][degree]
							  +deltaT/M[ielem][0][0]*rhspo[ielem][ivar][degree]);

//-----------------------------------------------------------------------
//							END of Three Stage RK Method
//-----------------------------------------------------------------------

cout<<"-------------------------------------------"<<endl;

}


}
