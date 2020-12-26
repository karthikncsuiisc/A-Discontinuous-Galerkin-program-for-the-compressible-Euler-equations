//Function to decide which function to choose.

#include "functions.h"
#include "oned_header.h"

void dgp2_fun_1d()
{
	
	for(int i=2;i<=nelem+1;i++)
	{
	double u1gp,u2gp,u3gp;
	double rhogp,velgp,presgp;
	double funval1,funval2,funval3;
	double u1,u2,u3;
	double u1x,u2x,u3x;
	double u1xx,u2xx,u3xx;
	
	
	rhspo[i][0][0]=-(flux[i+1][0]-flux[i][0]);	
	rhspo[i][1][0]=-(flux[i+1][1]-flux[i][1]);	
	rhspo[i][2][0]=-(flux[i+1][2]-flux[i][2]);
	
	u1=unkel[i][0][0];
	u2=unkel[i][1][0];
	u3=unkel[i][2][0];

	u1x=unkel[i][0][1];
	u2x=unkel[i][1][1];
	u3x=unkel[i][2][1];

	u1xx=unkel[i][0][2];
	u2xx=unkel[i][1][2];
	u3xx=unkel[i][2][2];

		
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
	u1gp=u1+u1x*sqrt(15)/5.0;
	u2gp=u2+u2x*sqrt(15)/5.0;
	u3gp=u3+u3x*sqrt(15)/5.0;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval2=velgp*(u3gp+presgp);

	//Gauss point 3
	u1gp=u1-u1x*sqrt(15)/5.0;
	u2gp=u2-u2x*sqrt(15)/5.0;
	u3gp=u3-u3x*sqrt(15)/5.0;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval3=velgp*(u3gp+presgp);

	rhspo[i][2][1]=8.0/9.0*funval1+5.0/9.0*funval2+5.0/9.0*funval3-(flux[i+1][2]+flux[i][2]);
	
	//Calculating uxx variable
	
//-----------------Equation 1 uxx term-----------------------------

	//Gauss point 1
	double funval1_u1,funval2_u1,funval3_u1;
	double funval1_u2,funval2_u2,funval3_u2;
	double funval1_u3,funval2_u3,funval3_u3;
	double B2,B3;
	double w1,w2,w3;

	B2=0;
	B3=-1.0/6.0;
	u1gp=u1+B2*u1x+B3*u1xx;
	u2gp=u2+B2*u2x+B3*u2xx;
	u3gp=u3+B2*u3x+B3*u3xx;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval1_u1=u2gp*B2;
	funval1_u2=(rhogp*velgp*velgp+presgp)*B2;
	funval1_u3=velgp*(u3gp+presgp)*B2;

	//Gauss point 2
	B2=sqrt(15)/5.0;
	B3=2.0/15.0;
	
	u1gp=u1+B2*u1x+B3*u1xx;
	u2gp=u2+B2*u2x+B3*u2xx;
	u3gp=u3+B2*u3x+B3*u3xx;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval2_u1=u2gp*B2;
	funval2_u2=(rhogp*velgp*velgp+presgp)*B2;
	funval2_u3=velgp*(u3gp+presgp)*B2;

	//Gauss point 3
	B2=-sqrt(15)/5.0;
	B3=2.0/15.0;
	
	u1gp=u1+B2*u1x+B3*u1xx;
	u2gp=u2+B2*u2x+B3*u2xx;
	u3gp=u3+B2*u3x+B3*u3xx;

	rhogp=u1gp;
	velgp=u2gp/u1gp;
	presgp=(gama-1)*(u3gp-0.5*pow(u2gp,2)/u1gp);

	funval3_u1=u2gp*B2;
	funval3_u2=(rhogp*velgp*velgp+presgp)*B2;
	funval3_u3=velgp*(u3gp+presgp)*B2;
    
    w1=8.0/9.0;
    w2=5.0/9.0;
    w3=5.0/9.0;
    
	rhspo[i][0][2]=w1*funval1_u1+w2*funval2_u1+w3*funval3_u1-(flux[i+1][0]-flux[i][0])/3.0;
	rhspo[i][1][2]=w1*funval1_u2+w2*funval2_u2+w3*funval3_u2-(flux[i+1][1]-flux[i][1])/3.0;
	rhspo[i][2][2]=w1*funval1_u3+w2*funval2_u3+w3*funval3_u3-(flux[i+1][2]-flux[i][2])/3.0;
	}	

}
