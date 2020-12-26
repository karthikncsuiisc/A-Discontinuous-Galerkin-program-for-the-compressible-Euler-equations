//Subroutine for calculating the gradients for reconstruction
#include "functions.h"
#include "twod_header.h"

void recons_coef_calculation()
{
	int ileft, iright;
	double dx,dy,weigh;
	
	rhsel=new double* [nelem+nbface+1];
	for(int ielem=0;ielem<=nelem+nbface;ielem++) rhsel[ielem]=new double[3];
	
	for(int ielem=1;ielem<=nelem+nbface;ielem++)
	{
		rhsel[ielem][0]=0.0;
		rhsel[ielem][1]=0.0;
		rhsel[ielem][2]=0.0;
	}
	
	
	for(int iface=1; iface<=naface;iface++)
	{
		ileft=intface[iface][1];
		iright=intface[iface][2];
		
		dx=centerx[iright]-centerx[ileft];
		dy=centery[iright]-centery[ileft];
//		weigh=pow(sqrt(dx*dx+dy*dy),-1);
		weigh=1;
		
		rhsel[ileft][0]=rhsel[ileft][0]+weigh*weigh*dx*dx;
		rhsel[ileft][1]=rhsel[ileft][1]+weigh*weigh*dy*dy;
		rhsel[ileft][2]=rhsel[ileft][2]+weigh*weigh*dx*dy;
		
		rhsel[iright][0]=rhsel[iright][0]+weigh*weigh*dx*dx;
		rhsel[iright][1]=rhsel[iright][1]+weigh*weigh*dy*dy;
		rhsel[iright][2]=rhsel[iright][2]+weigh*weigh*dx*dy;
	} 	



}


