//Program to calculate Van leer flux vector splitting
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

void van_leer_FVS_1d()
{

	double u1,u2,u3;
	double fiplus0,fiplus1,fiplus2;
	double fiminus0,fiminus1,fiminus2;
	double Mach,dx;

	double u1_0,u2_0,u3_0;
	double u1_1,u2_1,u3_1;
	double u1_2,u2_2,u3_2;
	double rho_dum,vel_dum,pres_dum,c0;
	double rho0,vel0,pres0;
	double rho1,vel1,pres1;
	double rho2,vel2,pres2;


for(int node=2;node<=nelem+2;node++)
{
	
//-----------------------------------------------------------------------------------------------------------------------
//	Fi+ calculation
//-----------------------------------------------------------------------------------------------------------------------
switch(disc_order)
{
case 0:
	switch(recons)
	{
	case 0:
		u1=unkel[node-1][0][0];
		u2=unkel[node-1][1][0];
		u3=unkel[node-1][2][0];

		rho_dum=u1;
		vel_dum=u2/u1;
		pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);

	break;
	case 1:
		switch(reconsvariable)
		{
						//---------------------------------------------------
						//  	Reconstruction on conserved variables
						//---------------------------------------------------
		case 0:

			dx=coord_1d[node]-coord_1d[node-1];
		
			u1_1=unkel[node-1][0][0];
			u2_1=unkel[node-1][1][0];
			u3_1=unkel[node-1][2][0];
		
			u1_0=unkel[node-2][0][0];
			u2_0=unkel[node-2][1][0];
			u3_0=unkel[node-2][2][0];
		
			u1_2=unkel[node][0][0];
			u2_2=unkel[node][1][0];
			u3_2=unkel[node][2][0];
		
			switch(limitertype)
			{
			case 0:
				u1=u1_1+dx/2.0*minmod((u1_1-u1_0)/dx,(u1_2-u1_1)/dx);
				u2=u2_1+dx/2.0*minmod((u2_1-u2_0)/dx,(u2_2-u2_1)/dx);
				u3=u3_1+dx/2.0*minmod((u3_1-u3_0)/dx,(u3_2-u3_1)/dx);
			break;
			case 1:
				u1=u1_1+dx/2.0*MClimiter((u1_1-u1_0)/dx,(u1_2-u1_0)/(2.0*dx),(u1_2-u1_1)/dx);
				u2=u2_1+dx/2.0*MClimiter((u2_1-u2_0)/dx,(u2_2-u2_0)/(2.0*dx),(u2_2-u2_1)/dx);
				u3=u3_1+dx/2.0*MClimiter((u3_1-u3_0)/dx,(u3_2-u3_0)/(2.0*dx),(u3_2-u3_1)/dx);
			break;
			default: cout<<"wrong choice of limiter type"<<endl; exit(1);
			}
			
			rho_dum=u1;
			vel_dum=u2/u1;
			pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
		break;
						//---------------------------------------------------
						//  	End of Reconstruction on conserved variables
						//---------------------------------------------------
		
						//---------------------------------------------------
						//  	Reconstruction on primitive variables
						//---------------------------------------------------
		case 1:
			
			dx=coord_1d[node]-coord_1d[node-1];
		
			u1=unkel[node-2][0][0];
			u2=unkel[node-2][1][0];
			u3=unkel[node-2][2][0];
		
			rho0=u1;
			vel0=u2/u1;
			pres0=(gama-1)*(u3-0.5*pow(u2,2)/u1);
			
			u1=unkel[node-1][0][0];
			u2=unkel[node-1][1][0];
			u3=unkel[node-1][2][0];
		
			rho1=u1;
			vel1=u2/u1;
			pres1=(gama-1)*(u3-0.5*pow(u2,2)/u1);
		
			u1=unkel[node][0][0];
			u2=unkel[node][1][0];
			u3=unkel[node][2][0];
		
			rho2=u1;
			vel2=u2/u1;
			pres2=(gama-1)*(u3-0.5*pow(u2,2)/u1);

			switch(limitertype)
			{
			case 0:
				rho_dum=rho1+dx/2.0*minmod((rho1-rho0)/dx,(rho2-rho1)/dx);
				vel_dum=vel1+dx/2.0*minmod((vel1-vel0)/dx,(vel2-vel1)/dx);
				pres_dum=pres1+dx/2.0*minmod((pres1-pres0),(pres2-pres1)/dx);
			break;
			case 1:
				rho_dum=rho1+dx/2.0*MClimiter((rho1-rho0)/dx,(rho2-rho0)/(2.0*dx),(rho2-rho1)/dx);
				vel_dum=vel1+dx/2.0*MClimiter((vel1-vel0)/dx,(vel2-vel0)/(2.0*dx),(vel2-vel1)/dx);
				pres_dum=pres1+dx/2.0*MClimiter((pres1-pres0)/dx,(pres2-pres0)/(2.0*dx),(pres2-pres1)/dx);
			break;
			default: cout<<"wrong choice of limiter type"<<endl; exit(1);
			}
		break;
						//---------------------------------------------------
						//  	End of Reconstruction on primitive variables
						//---------------------------------------------------
		default: cout<<"wrong choice of reconstruction variable : "<<reconsvariable<<endl; exit(1);
		}
	break;
	default: cout<<"wrong choice of reconstruction  "<<recons<<endl; exit(1);
	}
break;
case 1:
	u1=unkel[node-1][0][0]+unkel[node-1][0][1];
	u2=unkel[node-1][1][0]+unkel[node-1][1][1];
	u3=unkel[node-1][2][0]+unkel[node-1][2][1];
	
	rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
break;
case 2:
	u1=unkel[node-1][0][0]+unkel[node-1][0][1]+unkel[node-1][0][2]/3.0;
	u2=unkel[node-1][1][0]+unkel[node-1][1][1]+unkel[node-1][1][2]/3.0;
	u3=unkel[node-1][2][0]+unkel[node-1][2][1]+unkel[node-1][2][2]/3.0;
	
	rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
break;
default: cout<<"wrong choice of discretization order"<<endl; exit(1);
}

//cout<<pres_dum<<" "<<rho_dum<<" "<<vel_dum<<endl;



	c0=sqrt(gama*max(pres_dum,hrd_lmtr_cin)/rho_dum);
	Mach=vel_dum/c0;
	
	if(Mach<-1.0)
	{
	fiplus0=0;
	fiplus1=0;
	fiplus2=0;
    }
	else if(Mach>=-1.0&&Mach<=1.0)
	{
	fiplus0=rho_dum*c0*0.25*pow((Mach+1),2);
	fiplus1=fiplus0*((gama-1)*vel_dum+2.0*c0)/gama;

	double factor;
	factor=0.5*abs(vel_dum)*abs(vel_dum)-0.5*vel_dum*vel_dum;
	factor=factor+((gama-1)*vel_dum+2.0*c0)*((gama-1)*vel_dum+2.0*c0)/(2*gama*gama-2.0);
	fiplus2=fiplus0*factor;
    }
	else if(Mach>1.0)
	{
	fiplus0=rho_dum*vel_dum;
	fiplus1=rho_dum*vel_dum*vel_dum+pres_dum;
	fiplus2=vel_dum*(u3+pres_dum);
    }
    
//-----------------------------------------------------------------------------------------------------------------------
//	F(i+1)-calculation
//-----------------------------------------------------------------------------------------------------------------------

switch(disc_order)
{
case 0:
	switch(recons)
	{
	case 0:
		u1=unkel[node][0][0];
		u2=unkel[node][1][0];
		u3=unkel[node][2][0];

		rho_dum=u1;
		vel_dum=u2/u1;
		pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
		
	break;
	case 1:
		switch(reconsvariable)
		{
						//---------------------------------------------------
						//  	Reconstruction on conserved variables
						//---------------------------------------------------
		case 0:
			dx=coord_1d[node+1]-coord_1d[node];
		
			u1_1=unkel[node][0][0];
			u2_1=unkel[node][1][0];
			u3_1=unkel[node][2][0];
		
			u1_0=unkel[node-1][0][0];
			u2_0=unkel[node-1][1][0];
			u3_0=unkel[node-1][2][0];
		
			u1_2=unkel[node+1][0][0];
			u2_2=unkel[node+1][1][0];
			u3_2=unkel[node+1][2][0];
		
			switch(limitertype)
			{
			case 0:
				u1=u1_1-dx/2.0*minmod((u1_1-u1_0)/dx,(u1_2-u1_1)/dx);
				u2=u2_1-dx/2.0*minmod((u2_1-u2_0)/dx,(u2_2-u2_1)/dx);
				u3=u3_1-dx/2.0*minmod((u3_1-u3_0)/dx,(u3_2-u3_1)/dx);
			break;
			case 1:
				u1=u1_1-dx/2.0*MClimiter((u1_1-u1_0)/dx,(u1_2-u1_0)/(2.0*dx),(u1_2-u1_1)/dx);
				u2=u2_1-dx/2.0*MClimiter((u2_1-u2_0)/dx,(u2_2-u2_0)/(2.0*dx),(u2_2-u2_1)/dx);
				u3=u3_1-dx/2.0*MClimiter((u3_1-u3_0)/dx,(u3_2-u3_0)/(2.0*dx),(u3_2-u3_1)/dx);
			break;
			default: cout<<"wrong choice of limiter type"<<endl; exit(1);
			}
			
			rho_dum=u1;
			vel_dum=u2/u1;
			pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
		break;
						//---------------------------------------------------
						//  	End of Reconstruction on conserved variables
						//---------------------------------------------------
		
						//---------------------------------------------------
						//  	Reconstruction on primitive variables
						//---------------------------------------------------
		case 1:
			
			dx=coord_1d[node+1]-coord_1d[node];
		
			u1=unkel[node-1][0][0];
			u2=unkel[node-1][1][0];
			u3=unkel[node-1][2][0];
		
			rho0=u1;
			vel0=u2/u1;
			pres0=(gama-1)*(u3-0.5*pow(u2,2)/u1);
		
			u1=unkel[node][0][0];
			u2=unkel[node][1][0];
			u3=unkel[node][2][0];
		
			rho1=u1;
			vel1=u2/u1;
			pres1=(gama-1)*(u3-0.5*pow(u2,2)/u1);
		
			u1=unkel[node+1][0][0];
			u2=unkel[node+1][1][0];
			u3=unkel[node+1][2][0];
		
			rho2=u1;
			vel2=u2/u1;
			pres2=(gama-1)*(u3-0.5*pow(u2,2)/u1);

			switch(limitertype)
			{
			case 0:
				rho_dum=rho1-dx/2.0*minmod((rho1-rho0)/dx,(rho2-rho1)/dx);
				vel_dum=vel1-dx/2.0*minmod((vel1-vel0)/dx,(vel2-vel1)/dx);
				pres_dum=pres1-dx/2.0*minmod((pres1-pres0)/dx,(pres2-pres1)/dx);
			break;
			case 1:
				rho_dum=rho1-dx/2.0*MClimiter((rho1-rho0)/dx,(rho2-rho0)/(2.0*dx),(rho2-rho1)/dx);
				vel_dum=vel1-dx/2.0*MClimiter((vel1-vel0)/dx,(vel2-vel0)/(2.0*dx),(vel2-vel1)/dx);
				pres_dum=pres1-dx/2.0*MClimiter((pres1-pres0)/dx,(pres2-pres0)/(2.0*dx),(pres2-pres1)/dx);
			break;
			default: cout<<"wrong choice of limiter type"<<endl; exit(1);
			}
		break;
		default: cout<<"wrong choice of reconstruction variable : "<<reconsvariable<<endl; exit(1);
		}
	break;
	default: cout<<"wrong choice of reconstruction"<<endl; exit(1);
	}
break;

case 1:
	u1=unkel[node][0][0]-unkel[node][0][1];
	u2=unkel[node][1][0]-unkel[node][1][1];
	u3=unkel[node][2][0]-unkel[node][2][1];
	
	rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
break;	

case 2:
	u1=unkel[node][0][0]-unkel[node][0][1]+unkel[node][0][2]/3.0;
	u2=unkel[node][1][0]-unkel[node][1][1]+unkel[node][1][2]/3.0;
	u3=unkel[node][2][0]-unkel[node][2][1]+unkel[node][2][2]/3.0;
	
	rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=(gama-1)*(u3-0.5*pow(u2,2)/u1);
break;

default: cout<<"wrong choice of discretization order"<<endl; exit(1);
}

//cout<<"Fi+1"<<endl;
//cout<<pres_dum<<" "<<rho_dum<<" "<<vel_dum<<endl;

	c0=sqrt(gama*max(pres_dum,hrd_lmtr_cin)/rho_dum);
	Mach=vel_dum/c0;
	
	if(Mach>1.0)
	{
	fiminus0=0;
	fiminus1=0;
	fiminus2=0;
    }
	else if(Mach>=-1.0&&Mach<=1.0)
	{
	fiminus0=-rho_dum*c0*0.25*pow((Mach-1),2);
	fiminus1=fiminus0*((gama-1)*vel_dum-2.0*c0)/gama;

	double factor;
	factor=0.5*abs(vel_dum)*abs(vel_dum)-0.5*vel_dum*vel_dum;
	factor=factor+((gama-1)*vel_dum-2*c0)*((gama-1)*vel_dum-2*c0)/(2*gama*gama-2.0);
	fiminus2=fiminus0*factor;
    }
	else if(Mach<-1.0)
	{
	fiminus0=rho_dum*vel_dum;
	fiminus1=rho_dum*vel_dum*vel_dum+pres_dum;
	fiminus2=vel_dum*(u3+pres_dum);
    }    
    
	flux[node][0]=fiplus0+fiminus0;
	flux[node][1]=fiplus1+fiminus1;
	flux[node][2]=fiplus2+fiminus2;
}
	
	
}
