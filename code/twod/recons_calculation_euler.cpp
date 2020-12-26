//Subroutine for calculating the gradients for reconstruction
#include "functions.h"
#include "twod_header.h"

void recons_calculation_euler()
{

	int ileft, iright;
	double dx,dy,weigh;
	double Axx, Ayy, Axy,det;
	double du1dx,du2dx,du3dx,du4dx;
	double du1dy,du2dy,du3dy,du4dy;
	
	for(int ielem=0;ielem<=nelem+nbface;ielem++)
	for(int ivar=0;ivar<nvar;ivar++)
	for(int i=0; i<2;i++)
	gradu[ielem][ivar][i]=0.0;
	
	for(int iface=1; iface<=naface;iface++)
	{
		ileft=intface[iface][1];
		iright=intface[iface][2];
		
		dx=centerx[iright]-centerx[ileft];
		dy=centery[iright]-centery[ileft];
//		weigh=pow(sqrt(dx*dx+dy*dy),-1);
		weigh=1;
		
		gradu[ileft][0][0]=gradu[ileft][0][0]+weigh*weigh*dx*(unkel[iright][0][0]-unkel[ileft][0][0]);
		gradu[ileft][0][1]=gradu[ileft][0][1]+weigh*weigh*dy*(unkel[iright][0][0]-unkel[ileft][0][0]);

		gradu[ileft][1][0]=gradu[ileft][1][0]+weigh*weigh*dx*(unkel[iright][1][0]-unkel[ileft][1][0]);
		gradu[ileft][1][1]=gradu[ileft][1][1]+weigh*weigh*dy*(unkel[iright][1][0]-unkel[ileft][1][0]);

		gradu[ileft][2][0]=gradu[ileft][2][0]+weigh*weigh*dx*(unkel[iright][2][0]-unkel[ileft][2][0]);
		gradu[ileft][2][1]=gradu[ileft][2][1]+weigh*weigh*dy*(unkel[iright][2][0]-unkel[ileft][2][0]);

		gradu[ileft][3][0]=gradu[ileft][3][0]+weigh*weigh*dx*(unkel[iright][3][0]-unkel[ileft][3][0]);
		gradu[ileft][3][1]=gradu[ileft][3][1]+weigh*weigh*dy*(unkel[iright][3][0]-unkel[ileft][3][0]);

		gradu[iright][0][0]=gradu[iright][0][0]+weigh*weigh*dx*(unkel[iright][0][0]-unkel[ileft][0][0]);
		gradu[iright][0][1]=gradu[iright][0][1]+weigh*weigh*dy*(unkel[iright][0][0]-unkel[ileft][0][0]);

		gradu[iright][1][0]=gradu[iright][1][0]+weigh*weigh*dx*(unkel[iright][1][0]-unkel[ileft][1][0]);
		gradu[iright][1][1]=gradu[iright][1][1]+weigh*weigh*dy*(unkel[iright][1][0]-unkel[ileft][1][0]);

		gradu[iright][2][0]=gradu[iright][2][0]+weigh*weigh*dx*(unkel[iright][2][0]-unkel[ileft][2][0]);
		gradu[iright][2][1]=gradu[iright][2][1]+weigh*weigh*dy*(unkel[iright][2][0]-unkel[ileft][2][0]);

		gradu[iright][3][0]=gradu[iright][3][0]+weigh*weigh*dx*(unkel[iright][3][0]-unkel[ileft][3][0]);
		gradu[iright][3][1]=gradu[iright][3][1]+weigh*weigh*dy*(unkel[iright][3][0]-unkel[ileft][3][0]);
	} 
	

	
	for(int ielem=1;ielem<=nelem;ielem++)
	{
		Axx=rhsel[ielem][0];
		Ayy=rhsel[ielem][1];
		Axy=rhsel[ielem][2];
		
		det=Axx*Ayy-Axy*Axy+1e-16;
		
//		cout<<Axx<<" "<<Ayy<<" "<<Axy<<endl;
//		cout<<"Determinent"<<det<<endl;
		
		du1dx=gradu[ielem][0][0];
		du1dy=gradu[ielem][0][1];

		du2dx=gradu[ielem][1][0];
		du2dy=gradu[ielem][1][1];

		du3dx=gradu[ielem][2][0];
		du3dy=gradu[ielem][2][1];

		du4dx=gradu[ielem][3][0];
		du4dy=gradu[ielem][3][1];
		
		gradu[ielem][0][0]=(du1dx*Ayy-du1dy*Axy)/det;
		gradu[ielem][0][1]=(du1dy*Axx-du1dx*Axy)/det;
		
		gradu[ielem][1][0]=(du2dx*Ayy-du2dy*Axy)/det;
		gradu[ielem][1][1]=(du2dy*Axx-du2dx*Axy)/det;

		gradu[ielem][2][0]=(du3dx*Ayy-du3dy*Axy)/det;
		gradu[ielem][2][1]=(du3dy*Axx-du3dx*Axy)/det;

		gradu[ielem][3][0]=(du4dx*Ayy-du4dy*Axy)/det;
		gradu[ielem][3][1]=(du4dy*Axx-du4dx*Axy)/det;
	}
	


}


