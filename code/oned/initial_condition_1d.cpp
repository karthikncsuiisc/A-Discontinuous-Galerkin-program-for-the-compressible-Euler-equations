//Function to prescribe initial conditions
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

void initial_condition_1d()
{

double initvel=1,initpres=1,initrho=1;

for(int i=0;i<=nelem+3;i++)
{
double xval=0.5*(coord_1d[i]+coord_1d[i+1]);

if(phytype==0)
{
	if(xval>=0&&xval<0.6) unkel[i][0][0]=1+exp(-200*pow(xval-0.3,2));
	else if(xval>=0.6&&xval<=0.8) unkel[i][0][0]=2;
	else unkel[i][0][0]=1;
	
	unkel[i][1][0]=unkel[i][0][0]*initvel;
	unkel[i][2][0]=initpres/(gama-1)+0.5*unkel[i][0][0]*initvel*initvel;
}	
else if(phytype==1)
{
	if(xval<0.5)
	{
	initvel=0;
	initpres=1;
	
	unkel[i][0][0]=1.0;
	unkel[i][1][0]=unkel[i][0][0]*initvel;
	unkel[i][2][0]=initpres/(gama-1)+0.5*unkel[i][0][0]*initvel*initvel;
	}

	if(xval>0.5)
	{
	initvel=0;
	initpres=0.1;
	
	unkel[i][0][0]=0.125;
	unkel[i][1][0]=unkel[i][0][0]*initvel;
	unkel[i][2][0]=initpres/(gama-1)+0.5*unkel[i][0][0]*initvel*initvel;
	}
}

if(ndegr>=2)
{
	unkel[i][0][1]=0;
	unkel[i][1][1]=0;
	unkel[i][2][1]=0;
}

if(ndegr==3)
{
	unkel[i][0][2]=0;
	unkel[i][1][2]=0;
	unkel[i][2][2]=0;
}

}



}
