//Function to Solve Euler equaiton using P1 method


#include "functions.h"
#include "twod_header.h"

void tvd_RK_euler_p1()
{

ofstream fout;

string filename;
filename=grid_name+"_residul.dat";
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"Time	mass relative_mass x-momentum relative-x y-momentum	relative-y Energy relative-energy"<<endl;
	
double phytime;
double presmin;

int ileft,iright;
int freqoutput=0,freqmax=200,flag_steady=0;
double deltaT;
double *deltaTloc;
double *veltot;
double velxi,velyi,ci;
double velxj,velyj,cj;
double u1,u2x,u2y,u3,velx0,vely0,pres0,rho0;
double utotal1rel,utotal2rel,utotal3rel,utotal4rel,utotalrel;
double utotal1ref,utotal2ref,utotal3ref,utotal4ref;
double utotal1,utotal2,utotal3,utotal4;
double det,R1,R2;
double Machn,rhog,veltg,cg,velng,r1,r2,r3,r4;
double vdotn,velxref,velyref,velnref,veltref,velnb,velxin,velyin,veltb,cref;

cout<<"----------------------------------------------------------"<<endl;
cout<<"    Solving using TVD RK3   the Euler Equaitons               "<<endl;

//-----------------------------------------------------
// Initiating the variables
//-----------------------------------------------------
	ndegr=3;
	
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
	
	double *M00,*M11,*M22,*M12;
	M00=new double[nelem+1];
	M11=new double[nelem+1];
	M22=new double[nelem+1];
	M12=new double[nelem+1];

	cout<<"check 5"<<endl;

	gauss_point_values();

	cout<<"check 6"<<endl;
	
//-----------------------------------------------------------------------------
//					Calculating Mass matrix(using three gauss points)
//-----------------------------------------------------------------------------
	for(int ielem=1; ielem<=nelem;ielem++)
	{
		M00[ielem]=area[ielem];

		M11[ielem]=0.0;
		M12[ielem]=0.0;	
		M22[ielem]=0.0;	
		for(int igp=1;igp<=3;igp++)
		{
			M11[ielem]=M11[ielem]+Bxquad[ielem][igp]*Bxquad[ielem][igp]*Wequad[igp]*area[ielem];
			M12[ielem]=M12[ielem]+Bxquad[ielem][igp]*Byquad[ielem][igp]*Wequad[igp]*area[ielem];
			M22[ielem]=M22[ielem]+Byquad[ielem][igp]*Byquad[ielem][igp]*Wequad[igp]*area[ielem];
		}
		det=M11[ielem]*M22[ielem]-M12[ielem]*M12[ielem];
		if(det==0){ cout<<ielem<<" det :"<<det<<endl; exit(1);}

//		cout<<"ielem :"<<ielem<<" "<<M00[ielem]<<" "<<M11[ielem]<<" "<<M12[ielem]<<" "<<M22[ielem]<<endl;
//		cout<<"ielem :"<<ielem<<" "<<deltax[ielem]<<" "<<deltay[ielem]<<endl;
		
	}
	
/*	cout<<"shape function"<<endl;

for(int iface=1;iface<=naface;iface++)
{
	cout<<"----------------iface------------"<<endl;
	cout<<Bxline_quad[iface][1][1]<<" "<<Bxline_quad[iface][1][2]<<endl;
	cout<<Bxline_quad[iface][2][1]<<" "<<Bxline_quad[iface][2][2]<<endl;
	cout<<Byline_quad[iface][1][1]<<" "<<Byline_quad[iface][1][2]<<endl;
	cout<<Byline_quad[iface][2][1]<<" "<<Byline_quad[iface][2][2]<<endl;
}
exit(1);*/
		cout<<"check 7"<<endl;
//-----------------------------------------------------------------------------
//					End of Calculating Mass matrix
//-----------------------------------------------------------------------------

initial_conditions_2d();

deltaTloc=new double[nelem+1];
veltot=new double[nelem+1];


if(studytype==0) 
{
	tf=(double)Niter;
	t0=0;
}

phytime=t0;

cout<<"ndegr :"<<ndegr<<endl;

if(recons==1)
{
cout<<"reconstruction is not implemented for P1"<<endl;
exit(1);
}

//-----------------------------------------------------------------------------
// 					Starting the time stepping
//-----------------------------------------------------------------------------

while(phytime<tf)
{
	
presmin=1e10;
					 //-------------------------------------------------------
					//	Time step calculation for euler equation
					//-------------------------------------------------------

if(CFL_use==1)
{

deltaT=1e10;

for(int ielem=1; ielem<=nelem; ielem++) veltot[ielem]=0;
for(int iface=1;iface<=naface;iface++) 
{
	ileft=intface[iface][1];
	iright=intface[iface][2];

	
	u1=unkel[ileft][0][0]+Bxline_lin[iface][1]*unkel[ileft][0][1]+Byline_lin[iface][1]*unkel[ileft][0][2];
   u2x=unkel[ileft][1][0]+Bxline_lin[iface][1]*unkel[ileft][1][1]+Byline_lin[iface][1]*unkel[ileft][1][2];
   u2y=unkel[ileft][2][0]+Bxline_lin[iface][1]*unkel[ileft][2][1]+Byline_lin[iface][1]*unkel[ileft][2][2];
	u3=unkel[ileft][3][0]+Bxline_lin[iface][1]*unkel[ileft][3][1]+Byline_lin[iface][1]*unkel[ileft][3][2];
	rho0=u1;
	velxi=u2x/u1;
	velyi=u2y/u1;
	pres0=(gama-1)*(u3-0.5*(velxi*velxi+velyi*velyi)*rho0);
	ci=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	presmin=min(pres0,presmin);
	
	if(iright<=nelem)
	{
	u1=unkel[iright][0][0]+Bxline_lin[iface][2]*unkel[iright][0][1]+Byline_lin[iface][2]*unkel[iright][0][2];
   u2x=unkel[iright][1][0]+Bxline_lin[iface][2]*unkel[iright][1][1]+Byline_lin[iface][2]*unkel[iright][1][2];
   u2y=unkel[iright][2][0]+Bxline_lin[iface][2]*unkel[iright][2][1]+Byline_lin[iface][2]*unkel[iright][2][2];
	u3=unkel[iright][3][0]+Bxline_lin[iface][2]*unkel[iright][3][1]+Byline_lin[iface][2]*unkel[iright][3][2];
	
	rho0=u1;
	velxj=u2x/u1;
	velyj=u2y/u1;
	pres0=(gama-1)*(u3-0.5*(velxj*velxj+velyj*velyj)*rho0);
	cj=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	presmin=min(pres0,presmin);
	}
	else
	{ 
		vdotn=(velxi*nx[iface]+velyi*ny[iface]);
		if(bound_cond[iface][1]==2.0) //Wall boundary condition
		{
				velxj=velxi-2*vdotn*nx[iface];
				velyj=velyi-2*vdotn*ny[iface];
				
				pres0=(gama-1)*(u3-0.5*(velxj*velxj+velyj*velyj)*rho0);
				cj=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
		}
		else if(bound_cond[iface][1]==4.0) //Far field boundary condition
		{
				Machn=vdotn/ci;
			
			    velxref=velref*cos(alpha*3.14159265/180.0);
				velyref=velref*sin(alpha*3.14159265/180.0);
				if(Machn<=1.0)
				{
						velxj=velxref;
						velyj=velyref;
		
						rho0=rhoref;
						pres0=presref;
						cj=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
				}
				else if(Machn>=1.0)
				{
						velxj=velxi;
						velyj=velyi;
						cj=ci;
				}
				
				else if(Machn<=0&&Machn>-1.0)
				{
					velnref=velxref*(nx[iface])+velyref*(ny[iface]);
					veltref=velxref*(-ny[iface])+velyref*(nx[iface]);
			
					velnb=velxin*nx[iface]+velyin*ny[iface];
					veltb=velxin*(-ny[iface])+velyin*(nx[iface]);
					cref=sqrt(gama*presref/rhoref);
					
					r1=velnref-2.0*cref/(gama-1);
					r2=veltref;
					r3=presref/pow(rhoref,gama);
					r4=velnb+2.0*ci/(gama-1);
			
					velng=0.5*(r4+r1);
					cg=(gama-1)*0.25*(r4-r1);
					veltg=r2;
			
					velxj=velng*nx[iface]-veltg*ny[iface];
					velyj=velng*ny[iface]+veltg*nx[iface];
					
					cj=cg;
				}
				else if(Machn<1&&Machn>0.0)
				{
					
					velnref=velxref*(nx[iface])+velyref*(ny[iface]);
					veltref=velxref*(-ny[iface])+velyref*(nx[iface]);
			
					velnb=velxi*nx[iface]+velyi*ny[iface];
					veltb=velxi*(-ny[iface])+velyi*(nx[iface]);
			
					cref=sqrt(gama*presref/rhoref);
					
					r1=velnref-2.0*cref/(gama-1);
					r2=veltb;
					r3=pres0/pow(rho0,gama);
					r4=velnb+2*ci/(gama-1);
			
					velng=0.5*(r4+r1);
					cg=(gama-1)*0.25*(r4-r1);
					veltg=r2;
			
					velxj=velng*nx[iface]-veltg*ny[iface];
					velyj=velng*ny[iface]+veltg*nx[iface];
					
					cj=cg;
				 }
		}
	}
	
	veltot[ileft]=veltot[ileft]+0.5*(ci+cj+abs((velxi+velxj)*nx[iface]+(velyi+velyj)*ny[iface]))*length[iface];
	if(iright<=nelem)
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

flux2D_eul_calc_p1();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	{
		unkel[ielem][ivar][0]=unkeln[ielem][ivar][0]+deltaT/M00[ielem]*rhspo[ielem][ivar][0];

		det=M11[ielem]*M22[ielem]-M12[ielem]*M12[ielem]+1e-16;
		
		R1=M22[ielem]*rhspo[ielem][ivar][1]-M12[ielem]*rhspo[ielem][ivar][2];
		R2=M11[ielem]*rhspo[ielem][ivar][2]-M12[ielem]*rhspo[ielem][ivar][1];

		unkel[ielem][ivar][1]=unkeln[ielem][ivar][1]+R1*deltaT/det;
		unkel[ielem][ivar][2]=unkeln[ielem][ivar][2]+R2*deltaT/det;
		
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

flux2D_eul_calc_p1();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	{
		
		unkel[ielem][ivar][0]=0.75*unkeln[ielem][ivar][0]
    						 +0.25*(unkel[ielem][ivar][0]+deltaT/M00[ielem]*rhspo[ielem][ivar][0]);
		
		det=M11[ielem]*M22[ielem]-M12[ielem]*M12[ielem]+1e-16;
		
		R1=M22[ielem]*rhspo[ielem][ivar][1]-M12[ielem]*rhspo[ielem][ivar][2];
		R2=M11[ielem]*rhspo[ielem][ivar][2]-M12[ielem]*rhspo[ielem][ivar][1];

		unkel[ielem][ivar][1]=0.75*unkeln[ielem][ivar][1]
    						 +0.25*(unkel[ielem][ivar][1]+R1*deltaT/det);

		unkel[ielem][ivar][2]=0.75*unkeln[ielem][ivar][2]
    						 +0.25*(unkel[ielem][ivar][2]+R2*deltaT/det);
    						     						 
	}

//-----------------------------------------------------------------------
//							Three Stage RK Method: Stage3
//-----------------------------------------------------------------------

for(int ielem=1;ielem<=nelem;ielem++)
for(int ivar=0;ivar<nvar;ivar++)
for(int degree=0;degree<ndegr; degree++)
rhspo[ielem][ivar][degree]=0.0;

flux2D_eul_calc_p1();

	for(int ielem=1;ielem<=nelem;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)	
	{
		
		unkel[ielem][ivar][0]=(1.0/3.0)*unkeln[ielem][ivar][0]
    						 +(2.0/3.0)*(unkel[ielem][ivar][0]+deltaT/M00[ielem]*rhspo[ielem][ivar][0]);
		
		det=M11[ielem]*M22[ielem]-M12[ielem]*M12[ielem]+1e-16;
		
		R1=M22[ielem]*rhspo[ielem][ivar][1]-M12[ielem]*rhspo[ielem][ivar][2];
		R2=M11[ielem]*rhspo[ielem][ivar][2]-M12[ielem]*rhspo[ielem][ivar][1];

		unkel[ielem][ivar][1]=(1.0/3.0)*unkeln[ielem][ivar][1]
    						 +(2.0/3.0)*(unkel[ielem][ivar][1]+R1*deltaT/det);

		unkel[ielem][ivar][2]=(1.0/3.0)*unkeln[ielem][ivar][2]
    						 +(2.0/3.0)*(unkel[ielem][ivar][2]+R2*deltaT/det);
	}

//-----------------------------------------------------------------------
//							END of Three Stage RK Method
//-----------------------------------------------------------------------
}

//results2D();
//exit(1);

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

	delete [] M11;
	delete [] M12;
	delete [] M22;
	
}
