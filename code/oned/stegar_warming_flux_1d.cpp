//Subroutine to calculate stegar_warming_flux for 1d euler euqations
// The cells are of equal length. The function gives the results as an output file: results_1d_wave.dat
//
//-------------------------------------------------------------------------------------------------
// Cells discrption
//--------------------------------------------------------------------------------------------------
//	  Ghost								domain cells									  Ghost cells	
//	  0	  1	   2   3     4    5     6    ..................          nelem  .   nelem+1		  nelem+3         
//	|---|---|____|____|____|_____|_____|_____________________|______|_____|_____|_____|-----|-----|	
//	0   1   2	 3    4    5     6     7 .................        nelem      nelem+1     nelem+3	
//	Points location
//--------------------------------------------------------------------------------------------------

#include "functions.h"
#include "oned_header.h"

void stegar_warming_flux_1d()
{
	
//Calculating the Steger-warning FVS flux calculation
for(int node=2;node<=nelem+2;node++)
{
	double u1,u2,u3;
	double rho_dum,vel_dum,pres_dum,c0;
	double lambda1,lambda2,lambda3;
	double fiplus0,fiplus1,fiplus2;
	
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


	rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
	c0=sqrt(gama*max(pres_dum,hrd_lmtr_cin)/rho_dum);
	lambda1=0.5*(vel_dum+abs(vel_dum));
	lambda2=0.5*(vel_dum+c0+abs(vel_dum+c0));
	lambda3=0.5*(vel_dum-c0+abs(vel_dum-c0));
	
	fiplus0=2*(gama-1)*lambda1+lambda2+lambda3;
	fiplus0=fiplus0*(0.5*rho_dum/gama);
	
	fiplus1=2*(gama-1)*lambda1*vel_dum+lambda2*(vel_dum+c0)+lambda3*(vel_dum-c0);
	fiplus1=fiplus1*(0.5*rho_dum/gama);
	
	fiplus2=(gama-1)*lambda1*pow(vel_dum,2)+0.5*lambda2*pow(vel_dum+c0,2)+0.5*lambda3*pow(vel_dum-c0,2);
	fiplus2=fiplus2+(3-gama)/(2*gama-2)*(lambda2+lambda3)*pow(c0,2);
	fiplus2=0.5*fiplus2*rho_dum/gama;

	//F(i+1)-calculation
	double fiminus0,fiminus1,fiminus2;
	
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

	rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
	c0=sqrt(gama*max(pres_dum,hrd_lmtr_cin)/rho_dum);
	lambda1=0.5*(vel_dum-abs(vel_dum));
	lambda2=0.5*(vel_dum+c0-abs(vel_dum+c0));
	lambda3=0.5*(vel_dum-c0-abs(vel_dum-c0));
	
	fiminus0=2*(gama-1)*lambda1+lambda2+lambda3;
	fiminus0=fiminus0*(0.5*rho_dum/gama);
	
	fiminus1=2*(gama-1)*lambda1*vel_dum+lambda2*(vel_dum+c0)+lambda3*(vel_dum-c0);
	fiminus1=fiminus1*(0.5*rho_dum/gama);
	
	fiminus2=(gama-1)*lambda1*pow(vel_dum,2)+0.5*lambda2*pow(vel_dum+c0,2)+0.5*lambda3*pow(vel_dum-c0,2);
	fiminus2=fiminus2+(3-gama)/(2*gama-2)*(lambda2+lambda3)*pow(c0,2);
	fiminus2=fiminus2*0.5*rho_dum/gama;

	flux[node][0]=fiplus0+fiminus0;
	flux[node][1]=fiplus1+fiminus1;
	flux[node][2]=fiplus2+fiminus2;
}
//End of Steger-warning FVS flux method

	
	
	
}
