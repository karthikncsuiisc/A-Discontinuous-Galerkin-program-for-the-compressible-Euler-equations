//Subroutine for calculaitng the error function(entropy)

#include "functions.h"
#include "twod_header.h"

void error_function()
{
	ofstream fout;
	double velxj,velyj;
	double u1,u2x,u2y,u3,velx0,vely0,pres0,rho0;

					//--------------------------------------------------
					//  	Entropy Calculation for euler equations
					//--------------------------------------------------
	entropy=new double[nelem+1];
	double errorfun=0.0;
	for(int ielem=1;ielem<=nelem;ielem++)
	{
		entropy[ielem]=0;
		
		u1=unkel[ielem][0][0];
       u2x=unkel[ielem][1][0];
	   u2y=unkel[ielem][2][0];
  		u3=unkel[ielem][3][0];
  		
		rho0=u1;
		velxj=u2x/u1;
		velyj=u2y/u1;
		pres0=(gama-1)*(u3-0.5*(velxj*velxj+velyj*velyj)*rho0);
		entropy[ielem]=pres0/pow(rho0,gama);
	}
	
					//--------------------------------------------------
					//  	Error function calculaiton
					//--------------------------------------------------
	
	if(disc_order==1)
	{
		for(int ielem=1;ielem<=nelem;ielem++)
		{	
		for(int igp=1; igp<4;igp++)
		{
		
		u1=unkel[ielem][0][0]+unkel[ielem][0][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][0][2]*Byquad[ielem][igp];

       u2x=unkel[ielem][1][0]+unkel[ielem][1][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][1][2]*Byquad[ielem][igp];

	   u2y=unkel[ielem][2][0]+unkel[ielem][2][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][2][2]*Byquad[ielem][igp];
       					
  		u3=unkel[ielem][3][0]+unkel[ielem][3][1]*Bxquad[ielem][igp]
       						 +unkel[ielem][3][2]*Byquad[ielem][igp];
       						   
       	rho0=u1;
		velx0=u2x/u1;
		vely0=u2y/u1;
		pres0=(gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0);
		
		errorfun=errorfun+Wequad[igp]*pow(pres0/pow(rho0,gama)-presref/pow(rhoref,gama),2.0)*area[ielem]/pow(presref/pow(rhoref,gama),2.0);
		}
		}
	}
	else if(disc_order==0)
	{
	for(int ielem=1;ielem<=nelem;ielem++)
	errorfun=errorfun+pow(entropy[ielem]-presref/pow(rhoref,gama),2.0)*area[ielem]/pow(presref/pow(rhoref,gama),2.0);
	}

	
	fout.open("Error_function", ios::app);
	cout<<"-----------------------------------------------------------------------------"<<endl;
	cout<<"Error Function :"<<sqrt(errorfun)<<endl;
//	fout<<endl;
	fout<<grid_name<<" "<<1.0/sqrt((double)nelem)<<"  "<<sqrt(errorfun)<<" Dicreti P"<<disc_order;
	fout<<" Recons "<<recons<<" Lim "<<limiter<<" Rel_tol "<<tol<<endl;
	fout.close();


}


