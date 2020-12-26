//Function that contributes the source terms for burger equation

// The cells are of equal length. The function gives the results as an output file: results_1d_wave.dat
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

void rhs_contrb_burger_p2()
{
	
//	cout<<"entering diff p2"<<endl;

double u1,u1x,u1xx,B2,B3,ugp,sourceadvux,sourceadvuxx;
	
for(int i=2;i<=nelem+1;i++)
{	
	
//-----------------------------------------------------
//	Contribution for the Ux rhspo
//------------------------------------------------------
 
		u1=unkel[i][0][0];
		u1x=unkel[i][0][1];
		u1xx=unkel[i][0][2];				
		
		//Coontribution by advection flux
		//Gauss point 1
		B2=0;
		B3=-1.0/6.0;
		ugp=u1+B2*u1x+B3*u1xx;
		
		sourceadvux=8.0/9.0*(coefa*ugp+coefb*coefb*ugp*ugp/2.0);
		sourceadvuxx=8.0/9.0*B2*(coefa*ugp+coefb*coefb*ugp*ugp/2.0);

		//Gauss point 2
		B2=sqrt(15)/5.0;
		B3=2.0/15.0;
		ugp=u1+B2*u1x+B3*u1xx;				

		sourceadvux=sourceadvux+5.0/9.0*(coefa*ugp+coefb*coefb*ugp*ugp/2.0);
		sourceadvuxx=sourceadvuxx+5.0/9.0*B2*(coefa*ugp+coefb*coefb*ugp*ugp/2.0);

		//Gauss point 3
		B2=-sqrt(15)/5.0;
		B3=2.0/15.0;
		ugp=u1+B2*u1x+B3*u1xx;	

		sourceadvux=sourceadvux+5.0/9.0*(coefa*ugp+coefb*coefb*ugp*ugp/2.0);
		sourceadvuxx=sourceadvuxx+5.0/9.0*B2*(coefa*ugp+coefb*coefb*ugp*ugp/2.0);

		rhspo[i][0][1]=rhspo[i][0][1]+sourceadvux;
		rhspo[i][0][2]=rhspo[i][0][2]+sourceadvuxx;
}			

}
