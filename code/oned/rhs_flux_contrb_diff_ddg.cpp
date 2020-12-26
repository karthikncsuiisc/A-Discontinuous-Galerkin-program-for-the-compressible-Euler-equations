//Program contributes diffusive terms to rhspo using DDG method

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

void rhs_flux_contrb_diff_ddg()
{

double ul,ur,uxl,uxr;
double fluxdiff;


double beta0,beta1;

if(disc_order==0)
{
	cout<<"DDG diffusive scheme is not capplicable for P0 methods"<<endl;  exit(1);
}

for(int node=2;node<=nelem+2;node++)
{
	//Fi+ calcualtion
switch(disc_order)
{
    case 1: ul=unkel[node-1][0][0]+unkel[node-1][0][1]; break;
 	
 	case 2: ul=unkel[node-1][0][0]+unkel[node-1][0][1]+unkel[node-1][0][2]/3.0; 
			uxl=2.0*(unkel[node-1][0][1]+unkel[node-1][0][2])/(coord_1d[node+1]-coord_1d[node]); break;
			
    case 0: cout<<"DDG diffusive scheme is not implemented for P0 methods"<<endl;  exit(1);
	
    default: cout<<"wrong choice of discretization order"<<endl; exit(1);
}

	//F(i+1)-calculation
switch(disc_order)
{
    case 1: ur=unkel[node][0][0]-unkel[node][0][1];break;

    case 2: ur=unkel[node][0][0]-unkel[node][0][1]+unkel[node][0][2]/3.0;
   			uxr=2.0*(unkel[node][0][1]-unkel[node][0][2])/(coord_1d[node+1]-coord_1d[node]);
	break;

}	


		double deltaX,deltaxi,deltaxj;
		deltaxi=(coord_1d[node]-coord_1d[node-1]);
		deltaxj=(coord_1d[node+1]-coord_1d[node]);
		if(disc_order==1)
		{
		beta0=1.0;

		deltaX=0.5*(deltaxi+deltaxj);
		fluxdiff=beta0*(ur-ul)/deltaX
				+(unkel[node-1][0][1]/deltaxi+unkel[node][0][1]/deltaxj);
		}
		else if(disc_order==2)
		{
		beta0=2.0;
		beta1=1.0/12.0;
		
		deltaX=0.5*(deltaxi+deltaxj);
		fluxdiff=beta0*(ur-ul)/deltaX+(uxr+uxl)/2.0
				+beta1*(unkel[node][0][2]/deltaxj-unkel[node-1][0][2]/deltaxi)*2.0*deltaX;
		}


	flux[node][0]=flux[node][0]-coefc*fluxdiff;
	
}

	
	
}
