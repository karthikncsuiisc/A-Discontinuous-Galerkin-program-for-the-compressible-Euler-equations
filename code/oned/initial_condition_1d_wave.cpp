//Function to prescribe initial conditions

#include "functions.h"
#include "oned_header.h"

void initial_condition_1d_wave()
{


for(int i=2;i<=nelem+1;i++)
{
double xval=0.5*(coord_1d[i+1]+coord_1d[i]);

	if(xval>=0&&xval<0.6) unkel[i][0][0]=exp(-200*pow(xval-0.3,2));
	else if(xval>=0.6&&xval<=0.8) unkel[i][0][0]=1;
	else unkel[i][0][0]=0;

if(ndegr>=2) unkel[i][0][1]=0;
if(ndegr==3) unkel[i][0][2]=0;

}



}
