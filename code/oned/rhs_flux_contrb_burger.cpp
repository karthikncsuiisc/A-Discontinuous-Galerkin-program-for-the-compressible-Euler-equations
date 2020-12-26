//Program to calculte rieman solution for advectione equaiton using limiters and reconstruction
//Implemented for P0,P1,P2.. Minmod and MC limiter as been used

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

void rhs_flux_contrb_burger()
{

double ul,ur,fluxleft,fluxright;
double fluxadve,fluxburg;

for(int node=2;node<=nelem+2;node++)
{
	//Fi+ calcualtion
switch(disc_order)
{
	case 0: ul=unkel[node-1][0][0];	break;
    case 1: ul=unkel[node-1][0][0]+unkel[node-1][0][1]; break;
  	case 2: ul=unkel[node-1][0][0]+unkel[node-1][0][1]+unkel[node-1][0][2]/3.0; break;
    default: cout<<"wrong choice of discretization order"<<endl; exit(1);
}

	//F(i+1)-calculation
switch(disc_order)
{
	case 0: ur=unkel[node][0][0]; break;
    case 1: ur=unkel[node][0][0]-unkel[node][0][1];break;
    case 2: ur=unkel[node][0][0]-unkel[node][0][1]+unkel[node][0][2]/3.0; break;

}	
//	advection contribution
	if(coefa>=0) fluxadve=coefa*ul;
	else         fluxadve=coefa*ur;

//	burger term contribution
    fluxright=0.5*coefb*pow(0.5*(ur-abs(ur)),2);
    fluxleft=0.5*coefb*pow(0.5*(ul+abs(ul)),2);

    fluxburg=fluxleft+fluxright;

	flux[node][0]=flux[node][0]+fluxadve+fluxburg;
	
}

	
	
}
