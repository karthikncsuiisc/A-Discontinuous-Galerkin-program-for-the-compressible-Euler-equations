//FUnction solve 2D euler equaitons using P0 methods and reconstructed P0P1 methods

#include "functions.h"
#include "twod_header.h"

void tvd_rk3_method()
{

ofstream fout;

string filename;
filename=grid_name+"_residul.dat";
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"Time	mass relative_mass x-momentum relative-x y-momentum	relative-y Energy relative-energy"<<endl;
	
double phytime;
double presmin;

cout<<"----------------------------------------------------------"<<endl;
cout<<"    Solving using TVD RK3   the Euler Equaitons               "<<endl;

//-----------------------------------------------------
// Initiating the variables
//-----------------------------------------------------

	
	rhspo= new double** [nelem+1];
	for(int ielem=0;ielem<=nelem; ielem++) rhspo[ielem]=new double* [nvar];
	for(int ielem=0;ielem<=nelem; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) rhspo[ielem][ivar]=new double[ndegr];
	
	cout<<"check 1"<<endl;
	
	unkel= new double** [nelem+1+nface];
	for(int ielem=0;ielem<=nelem+nface; ielem++) unkel[ielem]=new double* [nvar];
	for(int ielem=0;ielem<=nelem+nface; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) unkel[ielem][ivar]=new double[ndegr];
	cout<<"check 2"<<endl;
	
	for(int ielem=0;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)
	for(int degree=0;degree<ndegr; degree++)
	unkel[ielem][ivar][degree]=0.0;
	cout<<"check 3"<<endl;
	
	double ***unkeln;
	unkeln= new double** [nelem+1];
	for(int ielem=0;ielem<=nelem; ielem++) unkeln[ielem]=new double* [nvar];
	for(int ielem=0;ielem<=nelem; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) unkeln[ielem][ivar]=new double[ndegr];
    cout<<"check 4"<<endl;
	
	double ***M;
	M=new double** [nelem+1];
	for(int ielem=0;ielem<=nelem;ielem++)
	M[ielem]=new double* [1];
	for(int ielem=0;ielem<=nelem;ielem++)
	for(int j=0;j<1;j++)
	M[ielem][j]=new double[1];
	
	cout<<"check 5"<<endl;
	
//-----------------------------------------------------------------------------
//					Calculating Mass matrix(using three gauss points)
//-----------------------------------------------------------------------------
	
	for(int ielem=1; ielem<=nelem;ielem++)   M[ielem][0][0]=area[ielem];

//-----------------------------------------------------------------------------
//					End of Calculating Mass matrix
//-----------------------------------------------------------------------------

initial_conditions_2d();

int ileft,iright;
int freqoutput=0,freqmax=100,flag_steady=0;
double deltaT;
double *deltaTloc;
double *veltot;
double velxi,velyi,ci;
double velxj,velyj,cj;
double u1,u2x,u2y,u3,pres0,velx0,vely0,rho0,Mach;
double utotal1rel,utotal2rel,utotal3rel,utotal4rel,utotalrel;
double utotal1ref,utotal2ref,utotal3ref,utotal4ref;
double utotal1,utotal2,utotal3,utotal4;

deltaTloc=new double[nelem+1];
veltot=new double[nelemtot+1];

if(studytype==0) 
{
	tf=(double)Niter;
	t0=0;
}

phytime=t0;

cout<<ndegr<<endl;

if(recons==1)
{
	recons_coef_calculation();

	gradu= new double** [nelem+1+nbface];
	for(int ielem=0;ielem<=nelem+nbface; ielem++) gradu[ielem]=new double* [nvar];
	
	for(int ielem=0;ielem<=nelem+nbface; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) gradu[ielem][ivar]=new double[2];
}

//-----------------------------------------------------------------------------
// 					Starting the time stepping
//-----------------------------------------------------------------------------

while(phytime<tf)
{
	
presmin=1e10;

boundary_condition_euler_2d();

if(recons==0) flux2D_eul_calc();
else if(recons==1) flux2D_eul_calc_recons();

					//-------------------------------------------------------
					//	Time step calculation for euler equation
					//-------------------------------------------------------

if(CFL_use==1)
{

deltaT=1e10;

for(int ielem=1; ielem<=nelemtot; ielem++) veltot[ielem]=0;
for(int iface=1;iface<=naface;iface++) 
{
	
	ileft=intface[iface][1];
	iright=intface[iface][2];
	
	u1=unkel[ileft][0][0];
	u2x=unkel[ileft][1][0];
	u2y=unkel[ileft][2][0];
	u3=unkel[ileft][3][0];
	
	rho0=u1;
	velxi=u2x/u1;
	velyi=u2y/u1;
	pres0=(gama-1)*(u3-0.5*(velxi*velxi+velyi*velyi)*rho0);
	ci=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);

	u1=unkel[iright][0][0];
	u2x=unkel[iright][1][0];
	u2y=unkel[iright][2][0];
	u3=unkel[iright][3][0];
	
	rho0=u1;
	velxj=u2x/u1;
	velyj=u2y/u1;
	pres0=(gama-1)*(u3-0.5*(velxj*velxj+velyj*velyj)*rho0);
	cj=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	
	presmin=min(pres0,presmin);
	
	veltot[ileft]=veltot[ileft]+0.5*(ci+cj+abs((velxi+velxj)*nx[iface]+(velyi+velyj)*ny[iface]))*length[iface];
	veltot[iright]=veltot[iright]+0.5*(ci+cj+abs((velxi+velxj)*nx[iface]+(velyi+velyj)*ny[iface]))*length[iface];
	
}

if(presmin<=0)
{
	 cout<<"Minimum Pressure is negative :"<<presmin<<endl;
	 cout<<"Use a small positive value of hard limiter(hrdlmtr)"<<endl;
	 cout<<"in the input file to avoid NAN values of sound speed(C)"<<endl;
	 cout<<"----------------------"<<endl; 
if(hrd_lmtr_cin<0) exit(1);
}


for(int ielem=1; ielem<=nelem; ielem++)
{
	deltaTloc[ielem]=CFL*area[ielem]/veltot[ielem];
	deltaT=min(deltaTloc[ielem],deltaT);
}
}
else deltaT=dt;

					//-------------------------------------------------------
					//	End of Time step calculation for euler equation
					//-------------------------------------------------------
if(studytype==1&&phytime+deltaT>=tf) deltaT=tf-phytime;

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

boundary_condition_euler_2d();
if(recons==0) flux2D_eul_calc();
else if(recons==1) flux2D_eul_calc_recons();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	for(int degree=0;degree<ndegr; degree++)
	{
		
	if(studytype==0&&timstp_option==0)
		unkel[ielem][ivar][degree]=unkeln[ielem][ivar][degree]+deltaTloc[ielem]/M[ielem][0][0]*rhspo[ielem][ivar][degree];
	else
		unkel[ielem][ivar][degree]=unkeln[ielem][ivar][degree]+deltaT/M[ielem][0][0]*rhspo[ielem][ivar][degree];
		
	}

if(timedep_method==3)
{
	
//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage2
//-----------------------------------------------------------------------

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
rhspo[ielem][ivar][degree]=0.0;

boundary_condition_euler_2d();
if(recons==0) flux2D_eul_calc();
else if(recons==1) flux2D_eul_calc_recons();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	for(int degree=0;degree<ndegr; degree++)
	{
	if(studytype==0&&timstp_option==0)
		unkel[ielem][ivar][degree]=0.75*unkeln[ielem][ivar][degree]
							  +0.25*(unkel[ielem][ivar][degree]+deltaTloc[ielem]/M[ielem][0][0]*rhspo[ielem][ivar][degree]);
	else
		unkel[ielem][ivar][degree]=0.75*unkeln[ielem][ivar][degree]
							  +0.25*(unkel[ielem][ivar][degree]+deltaT/M[ielem][0][0]*rhspo[ielem][ivar][degree]);
	}

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage3
//-----------------------------------------------------------------------

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
rhspo[ielem][ivar][degree]=0.0;

boundary_condition_euler_2d();

if(recons==0) flux2D_eul_calc();
else if(recons==1) flux2D_eul_calc_recons();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	for(int degree=0;degree<ndegr; degree++)
	{
	if(studytype==0&&timstp_option==0)
		unkel[ielem][ivar][degree]=(1.0/3.0)*unkeln[ielem][ivar][degree]
							  +(2/3.0)*(unkel[ielem][ivar][degree]+deltaTloc[ielem]/M[ielem][0][0]*rhspo[ielem][ivar][degree]);
	else
		unkel[ielem][ivar][degree]=(1.0/3.0)*unkeln[ielem][ivar][degree]
							  +(2/3.0)*(unkel[ielem][ivar][degree]+deltaT/M[ielem][0][0]*rhspo[ielem][ivar][degree]);
	}

//-----------------------------------------------------------------------
//							END of Three Stage RK Method
//-----------------------------------------------------------------------
}


	utotal1=0; utotal2=0.0; utotal3=0; utotal4=0;
					//--------------------------------------------------
					//			Printing output for euler equations
					//--------------------------------------------------
		for(int ielem=1;ielem<=nelem;ielem++)
		{
		utotal1=utotal1+pow(unkel[ielem][0][0]-unkeln[ielem][0][0],2)*area[ielem];
		utotal2=utotal2+pow(unkel[ielem][1][0]-unkeln[ielem][1][0],2)*area[ielem];
		utotal3=utotal3+pow(unkel[ielem][2][0]-unkeln[ielem][2][0],2)*area[ielem];
		utotal4=utotal4+pow(unkel[ielem][3][0]-unkeln[ielem][3][0],2)*area[ielem];
		}

		utotal1=sqrt(utotal1);
		utotal2=sqrt(utotal2);
		utotal3=sqrt(utotal3);
		utotal4=sqrt(utotal4);

		if(phytime==t0)
		{
		utotal1ref=utotal1;
		utotal2ref=utotal2;
		utotal3ref=utotal3;
		utotal4ref=utotal4;
		}

		utotal1rel=utotal1/(utotal1ref+1e-16);
		utotal2rel=utotal2/(utotal2ref+1e-16);
		utotal3rel=utotal3/(utotal3ref+1e-16);
		utotal4rel=utotal4/(utotal4ref+1e-16);
		
		freqoutput++;
		if(freqoutput==freqmax||phytime==t0||phytime+deltaT==tf)
		{
		cout<<"------------------------------------------------------------------------------"<<endl;
		
		cout<<"Grid Name :"<<grid_name<<endl;
		if(studytype==0)
		cout<<"Iteration Number: "<<phytime<<" Minimum psuedo time step: "<<deltaT<<endl;
		else if(studytype==1)
		cout<<"Time: "<<phytime<<" Time step: "<<deltaT<<endl;

		cout<<"Mass Residue"<<endl;
		cout<<"Absolute : "<<utotal1<<endl;
		cout<<"Relative : "<<utotal1rel<<endl;
		cout<<"x-Momentum Residue"<<endl;
		cout<<"Aboslute :"<<utotal2<<endl;
		cout<<"Relative :"<<utotal2rel<<endl;
		cout<<"y-Momentum Residue "<<endl; ;
		cout<<"Absolute :"<<utotal3<<endl;
		cout<<"Relative :"<<utotal3rel<<endl;
		cout<<"Energy Residue "<<endl;
		cout<<"Absolute :"<<utotal4<<endl;
		cout<<"Relative :"<<utotal4rel<<endl;
		
		fout<<phytime;
		fout<<" "<<utotal1;
		fout<<" "<<utotal1rel;
		fout<<" "<<utotal2;
		fout<<" "<<utotal2rel;
		fout<<" "<<utotal3;
		fout<<" "<<utotal3rel;
		fout<<" "<<utotal4;
		fout<<" "<<utotal4rel<<endl;
		freqoutput=0;
		}
		
					//--------------------------------------------------
					//			Checking the termination condition
					//--------------------------------------------------
	utotalrel=max(utotal1rel,utotal2rel);
	utotalrel=max(utotalrel,utotal3rel);
	utotalrel=max(utotalrel,utotal4rel);

	if(studytype==0&&utotalrel<tol)
	{
		cout<<"------------------------------------------------------------------------------"<<endl;
		cout<<"Iteration Number: "<<phytime<<" Minimum psuedo time step: "<<deltaT<<endl;
		cout<<"Mass Residue"<<endl;
		cout<<"Absolute : "<<utotal1<<endl;
		cout<<"Relative : "<<utotal1rel<<endl;
		cout<<"x-Momentum Residue"<<endl;
		cout<<"Aboslute :"<<utotal2<<endl;
		cout<<"Relative :"<<utotal2rel<<endl;
		cout<<"y-Momentum Residue "<<endl; ;
		cout<<"Absolute :"<<utotal3<<endl;
		cout<<"Relative :"<<utotal3rel<<endl;
		cout<<"Energy Residue "<<endl;
		cout<<"Absolute :"<<utotal4<<endl;
		cout<<"Relative :"<<utotal4rel<<endl;
		
		fout<<phytime;
		fout<<" "<<utotal1;
		fout<<" "<<utotal1rel;
		fout<<" "<<utotal2;
		fout<<" "<<utotal2rel;
		fout<<" "<<utotal3;
		fout<<" "<<utotal3rel;
		fout<<" "<<utotal4;
		fout<<" "<<utotal4rel<<endl;

		flag_steady=1;
		break;
	}
					//--------------------------------------------------
					//			End of Checking the termination condition
					//--------------------------------------------------

	
    if(studytype==0) phytime=phytime+1;
	else phytime=phytime+deltaT;

}
	
	fout.close();
	

	
	cout<<"-----------------------------------------------------------------------------"<<endl;

	if(studytype==0&&flag_steady==0)
	{
		cout<<"Steady state solution did not converged for grid :"<<grid_name<<endl;
		cout<<"Maximum Number of iterations reached :"<<phytime<<endl;
		results2D();
		exit(1);

	}
	else if(studytype==0&&flag_steady==1)
	{
		cout<<"Solution converged"<<endl;
		cout<<"Number of iterations :"<<phytime<<endl;
	}


	for(int ielem=0;ielem<=nelem; ielem++) for(int ivar=0;ivar<nvar;ivar++) delete [] unkeln[ielem][ivar];
	for(int ielem=0;ielem<=nelem; ielem++) delete [] unkeln[ielem];
	delete [] unkeln;
	
	for(int ielem=0;ielem<=nelem; ielem++)
	for(int ivar=0;ivar<nvar;ivar++) delete [] rhspo[ielem][ivar];
	for(int ielem=0;ielem<=nelem; ielem++) delete [] rhspo[ielem];
	delete [] rhspo;

	for(int ielem=0;ielem<=nelem;ielem++)	for(int j=0;j<1;j++) delete []	M[ielem][j];
	for(int ielem=0;ielem<=nelem;ielem++)	delete [] M[ielem];
	delete [] M;
	
	if(phytype==1&&recons==1)
	{
	for(int ielem=0;ielem<=nelem+nbface;ielem++) delete [] rhsel[ielem];
	delete [] rhsel;
	
	for(int ielem=0;ielem<=nelem+nbface; ielem++)	for(int ivar=0;ivar<nvar;ivar++) delete [] gradu[ielem][ivar];
	for(int ielem=0;ielem<=nelem+nbface; ielem++) delete [] gradu[ielem];
	delete [] gradu;
	}

}
