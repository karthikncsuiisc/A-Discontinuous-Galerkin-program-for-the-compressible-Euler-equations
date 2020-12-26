//Function to right hand side matrix for euler euqations using dgp1 nethod
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

void dgp1_fun_1d()
{
	
	for(int i=2;i<=nelem+1;i++)
	{
	double u1gp,u2gp,u3gp;
	double rhogp,velgp,presgp;
	double funval1,funval2,funval3;
	double u1,u2,u3;
	double u1x,u2x,u3x;
	
	rhspo[i][0][0]=-(flux[i+1][0]-flux[i][0]);	
	rhspo[i][1][0]=-(flux[i+1][1]-flux[i][1]);	
	rhspo[i][2][0]=-(flux[i+1][2]-flux[i][2]);
	
	u1=unkel[i][0][0];
	u2=unkel[i][1][0];
	u3=unkel[i][2][0];

	u1x=unkel[i][0][1];
	u2x=unkel[i][1][1];
	u3x=unkel[i][2][1];
		
//-----------------Equation 1 ux term-----------------------------	
	rhspo[i][0][1]=2*u2-(flux[i+1][0]+flux[i][0]);	

//-----------------Equation 2 ux term-----------------------------	
	
	
	//Gauss point 1
	u1gp=u1+u1x/sqrt(3);
	u2gp=u2+u2x/sqrt(3);
	u3gp=u3+u3x/sqrt(3);

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval1=(rhogp*velgp*velgp+presgp);

	//Gauss point 2
	u1gp=u1-u1x/sqrt(3);
	u2gp=u2-u2x/sqrt(3);
	u3gp=u3-u3x/sqrt(3);

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval2=(rhogp*velgp*velgp+presgp);

	rhspo[i][1][1]=funval1+funval2-(flux[i+1][1]+flux[i][1]);	
	
//-----------------Equation 3 ux term-----------------------------	
	
	//Gauss point 1
	u1gp=u1;
	u2gp=u2;
	u3gp=u3;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval1=velgp*(u3gp+presgp);

	//Gauss point 2
	u1gp=u1+u1x*sqrt(15.0)/5.0;
	u2gp=u2+u2x*sqrt(15.0)/5.0;
	u3gp=u3+u3x*sqrt(15.0)/5.0;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval2=velgp*(u3gp+presgp);

	//Gauss point 3
	u1gp=u1-u1x*sqrt(15.0)/5.0;
	u2gp=u2-u2x*sqrt(15.0)/5.0;
	u3gp=u3-u3x*sqrt(15.0)/5.0;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval3=velgp*(u3gp+presgp);

	rhspo[i][2][1]=8/9.0*funval1+5.0/9.0*funval2+5.0/9.0*funval3-(flux[i+1][2]+flux[i][2]);
		
	}

}
