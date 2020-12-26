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

void rieman_solver_oned_wave()
{

double ul,ur,K0,uxx,dx;
double a,b,c;

K0=0.0; //Fromms method

for(int node=2;node<=nelem+2;node++)
{
	//Fi+ calcualtion
switch(disc_order)
{
	case 0:
	  switch(recons)
	  {
 	       case 0: ul=unkel[node-1][0][0]; break;
	       case 1: 
					switch(limiter)
					{
					 case 0:
						ul=unkel[node-1][0][0]+0.25*(1-K0)*(unkel[node-1][0][0]-unkel[node-2][0][0]);
						ul=ul				 +0.25*(1+K0)*(unkel[node][0][0]-unkel[node-1][0][0]);
					 break;
					 case 1:
						if(limitertype==0)
						{
							dx=coord_1d[node]-coord_1d[node-1];
							a=(unkel[node-1][0][0]-unkel[node-2][0][0])/dx;
							b=(unkel[node][0][0]-unkel[node-1][0][0])/dx;
							ul=unkel[node-1][0][0]+minmod(a,b)*dx/2.0;
						}
						else if(limitertype==1)
						{
							dx=coord_1d[node]-coord_1d[node-1];
							a=(unkel[node-1][0][0]-unkel[node-2][0][0])/dx;
							b=(unkel[node][0][0]-unkel[node-2][0][0])/(2.0*dx);
							c=(unkel[node][0][0]-unkel[node-1][0][0])/dx;
							ul=unkel[node-1][0][0]+MClimiter(a,b,c)*dx/2;
						}
					break;
					}
			break;
      }
    break;
    case 1:
    	  switch(recons)
	  {
 	       case 0: ul=unkel[node-1][0][0]+unkel[node-1][0][1]; break;
	       case 1: 
				   switch(limiter)
				   {
					case 0:
						dx=coord_1d[node]-coord_1d[node-1];
						uxx=unkel[node-1][0][1]+0.25*(1-K0)*(unkel[node-1][0][1]-unkel[node-2][0][1]);
						uxx=uxx				  +0.25*(1+K0)*(unkel[node][0][1]-unkel[node-1][0][1]); 
						uxx=2*uxx/dx;
						ul=unkel[node-1][0][0]+unkel[node-1][0][1]+dx*dx/4.0*uxx; break;
					case 1:
						cout<<"Limiter is not implemented for P1"<<endl; exit(1);
				   }
	  }
    break;
    case 2:
    	  switch(recons)
	  {
 	       case 0: ul=unkel[node-1][0][0]+unkel[node-1][0][1]+unkel[node-1][0][2]/3.0; break;
	       case 1: cout<<"Reconstruction is not implemented for P2 advection equation, implemented only for p0"<<endl; exit(2);
      }
     break;
}

	//F(i+1)-calculation
switch(disc_order)
{
	case 0:
	  switch(recons)
	  {
 	       case 0: ur=unkel[node][0][0]; break;
	       case 1: 
					switch(limiter)
					{
						case 0:
							ur=unkel[node][0][0]-0.25*(1+K0)*(unkel[node][0][0]-unkel[node-1][0][0]);
							ur=ur  			    -0.25*(1-K0)*(unkel[node+1][0][0]-unkel[node][0][0]); break;
						case 1:
							if(limitertype==0)
							{
								dx=coord_1d[node]-coord_1d[node-1];
								a=(unkel[node][0][0]-unkel[node-1][0][0])/dx;
								b=(unkel[node+1][0][0]-unkel[node][0][0])/dx;
								ur=unkel[node][0][0]-minmod(a,b)*dx/2.0;
							}
							else if(limitertype==1)
							{
								dx=coord_1d[node+1]-coord_1d[node];
								a=(unkel[node][0][0]-unkel[node-1][0][0])/dx;
								b=(unkel[node+1][0][0]-unkel[node-1][0][0])/(2.0*dx);
								c=(unkel[node+1][0][0]-unkel[node][0][0])/dx;
								ur=unkel[node][0][0]-MClimiter(a,b,c)*dx/2;
							}
						break;
					}							
      }
    break;
    case 1:
    	  switch(recons)
	  {
 	       case 0: ur=unkel[node][0][0]-unkel[node][0][1];break;
	       case 1: 
				   switch(limiter)
				   {
					   case 0:
							dx=coord_1d[node]-coord_1d[node-1];
							uxx=unkel[node][0][1]-0.25*(1+K0)*(unkel[node][0][1]-unkel[node-1][0][1]);
							uxx=uxx			    -0.25*(1-K0)*(unkel[node+1][0][1]-unkel[node][0][1]); 
							uxx=2*uxx/dx;
							ur=unkel[node][0][0]-unkel[node][0][1]+dx*dx/4.0*uxx;
					   break;
					   case 1:
							cout<<"Limiter is not implemented for P1"<<endl; exit(3);
					   break;
					}
			break;
      }
    break;
    case 2:
    	  switch(recons)
	  {
 	       case 0: ur=unkel[node][0][0]-unkel[node][0][1]+unkel[node][0][2]/3.0; break;
	       case 1: cout<<"Reconstruction is not implemented for P2 advection equation"<<endl; exit(4);
      }
     break;
}	
	
	if(a0x>=0) flux[node][0]=a0x*ul;
	else if(a0x<0) flux[node][0]=a0x*ur;
	
}

	
	
}
