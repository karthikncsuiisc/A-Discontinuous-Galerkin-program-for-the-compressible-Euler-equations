//Function that calculates the rhspo(right hand side) of 1d diffusive equation using DDG method
// and P1 discretization
// The cells are of equal length.
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

void rhs_contrb_burger_p1()
{
	
double u1,u1x,ugp;
double funvalgp1,funvalgp2;
double sourceconv;
	
for(int i=2;i<=nelem+1;i++)
{	
		u1=unkel[i][0][0];
		u1x=unkel[i][0][1];
				
		//Gauss point 1
		ugp=u1+u1x/sqrt(3.0);
		funvalgp1=coefa*ugp+coefb*coefb*ugp*ugp/2.0;
				
		//Gauss point 2
		ugp=u1-u1x/sqrt(3.0);
		funvalgp2=coefa*ugp+coefb*coefb*ugp*ugp/2.0;
				
		sourceconv=(funvalgp1+funvalgp2);
				
		rhspo[i][0][1]=rhspo[i][0][1]+sourceconv;
}			

}
