//Function to right hand side matrix for euler euqations using dgp0 method
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

void dgp0_fun_1d()
{
	
	for(int i=2;i<=nelem+1;i++)
	{
	rhspo[i][0][0]=-(flux[i+1][0]-flux[i][0]);	
	rhspo[i][1][0]=-(flux[i+1][1]-flux[i][1]);	
	rhspo[i][2][0]=-(flux[i+1][2]-flux[i][2]);
	}

}
