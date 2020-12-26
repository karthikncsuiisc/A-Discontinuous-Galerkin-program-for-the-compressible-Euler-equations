//Function which solves the advection equaiton using P1 mehod for 2D problems


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
	M[ielem]=new double* [3];
	for(int ielem=1;ielem<=nelem;ielem++)
	for(int j=0;j<3;j++)
	M[ielem][j]=new double[3];

	if(ndegr>=2) gauss_point_values(); //Calling function to calculate values of shape function at gauss points

for(int ielem=1; ielem<=nelem;ielem++)
{
    M[ielem][0][0]=area[ielem];

//-------------------Calculating M11,M12(M12=M21),M22 for P1 method----------------------------
	if(ndegr>=2)
	{
		M[ielem][1][1]=0.0;
		M[ielem][1][2]=0.0;	
		M[ielem][2][2]=0.0;	
		for(int igp=1;igp<=3;igp++)
		{
			M[ielem][1][1]=M[ielem][1][1]+Bxquad[ielem][igp]*Bxquad[ielem][igp]*Wequad[igp]*area[ielem];
			M[ielem][1][2]=M[ielem][1][2]+Bxquad[ielem][igp]*Byquad[ielem][igp]*Wequad[igp]*area[ielem];
			M[ielem][2][2]=M[ielem][2][2]+Byquad[ielem][igp]*Byquad[ielem][igp]*Wequad[igp]*area[ielem];
		}
	
		M[ielem][2][1]=M[ielem][1][2];
//		cout<<"ielem :"<<ielem<<" "<<M[ielem][0][0]<<" "<<M[ielem][1][1]<<" "<<M[ielem][1][2]<<" "<<M[ielem][2][1]<<" "<<M[ielem][2][2]<<endl;
	}
//---------------End of calculating M11,M12(M12=M21),M22----------------------------
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
	utotal=utotal+unkel[ielem][0][0]*area[ielem];
}
cout<<"Total value of u1: "<<utotal<<endl;

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
