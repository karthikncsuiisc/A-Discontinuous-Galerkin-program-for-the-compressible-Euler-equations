//Subroutine for calculating the fluexes uing vanleer flux splitting schem
#include "functions.h"
#include "twod_header.h"

void flux2D_eul_calc()
{
	
int ileft,iright,ip1,ip2;
double u1,u2x,u2y,u3;
double rho0,velx0,vely0,pres0,c0,vdotn;
double fiplus0,fiplus1x,fiplus1y,fiplus2;
double fiminus0,fiminus1x,fiminus1y,fiminus2;
double Mach;


for(int iface=1;iface<=naface;iface++)
{
	
	ileft=intface[iface][1];
	iright=intface[iface][2];
	ip1=intface[iface][3];
	ip2=intface[iface][4];
	
//--------------------------------------------------------------------	
//   		Calculating rhspo for P0
//--------------------------------------------------------------------	
	
//-------------Fi+(left cell) calcualtion---------------

		u1=unkel[ileft][0][0];
	   u2x=unkel[ileft][1][0];
	   u2y=unkel[ileft][2][0];
		u3=unkel[ileft][3][0];

		
	rho0=u1;
	velx0=u2x/u1;
	vely0=u2y/u1;
//	pres0=max((gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0),1e-16);
//	c0=sqrt(gama*pres0/rho0);

	pres0=(gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0);
	c0=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	
	vdotn=velx0*nx[iface]+vely0*ny[iface];

	Mach=vdotn/c0;
	
	if(Mach<-1.0)
	{
	fiplus0=0;
	fiplus1x=0;
	fiplus1y=0;
	fiplus2=0;
    }
	else if(Mach>=-1.0&&Mach<=1.0)
	{
	double factor;
	
	factor=rho0*c0*(Mach+1)*(Mach+1)*0.25;
	
	fiplus0=factor;
	fiplus1x=factor*(velx0+nx[iface]*(-vdotn+2*c0)/gama);
	fiplus1y=factor*(vely0+ny[iface]*(-vdotn+2*c0)/gama);
	
	fiplus2=0.5*(velx0*velx0+vely0*vely0)-0.5*vdotn*vdotn;
	fiplus2=fiplus2+((gama-1)*vdotn+2.0*c0)*((gama-1)*vdotn+2.0*c0)/(2.0*gama*gama-2.0);
	fiplus2=factor*fiplus2;
    }
	else if(Mach>1.0)
	{
	fiplus0=rho0*vdotn;
	fiplus1x=rho0*vdotn*velx0+pres0*nx[iface];
	fiplus1y=rho0*vdotn*vely0+pres0*ny[iface];
	fiplus2=vdotn*(u3+pres0);
    }
    
//-------------Fj-(right cell) calcualtion---------------


	u1=unkel[iright][0][0];
   u2x=unkel[iright][1][0];
   u2y=unkel[iright][2][0];
	u3=unkel[iright][3][0];
	
	rho0=u1;
	velx0=u2x/u1;
	vely0=u2y/u1;
//	pres0=max((gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0),1e-16);
//	c0=sqrt(gama*pres0/rho0);

	pres0=(gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0);
	c0=sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);

	vdotn=velx0*nx[iface]+vely0*ny[iface];

	Mach=vdotn/c0;
		
	if(Mach>1.0)
	{
	fiminus0=0;
	fiminus1x=0;
	fiminus1y=0;
	fiminus2=0;
    }
	else if(Mach>=-1.0&&Mach<=1.0)
	{
	double factor;
	
	factor=-rho0*c0*(Mach-1)*(Mach-1)*0.25;
	
	fiminus0=factor;
	fiminus1x=factor*(velx0+nx[iface]*(-vdotn-2*c0)/gama);
	fiminus1y=factor*(vely0+ny[iface]*(-vdotn-2*c0)/gama);
	
	fiminus2=0.5*(velx0*velx0+vely0*vely0)-0.5*vdotn*vdotn;
	fiminus2=fiminus2+((gama-1)*vdotn-2.0*c0)*((gama-1)*vdotn-2.0*c0)/(2*gama*gama-2.0);
	fiminus2=factor*fiminus2;
    }
	else if(Mach<-1.0)
	{
	fiminus0=rho0*vdotn;
	fiminus1x=rho0*vdotn*velx0+pres0*nx[iface];
	fiminus1y=rho0*vdotn*vely0+pres0*ny[iface];
	fiminus2=vdotn*(u3+pres0);
    }
    
 //---------------------Calculating rhspo------------------------------
    
   	if(ileft<=nelem)
   	{
    rhspo[ileft][0][0]=rhspo[ileft][0][0]-(fiplus0+fiminus0)*length[iface];
	rhspo[ileft][1][0]=rhspo[ileft][1][0]-(fiplus1x+fiminus1x)*length[iface];
    rhspo[ileft][2][0]=rhspo[ileft][2][0]-(fiplus1y+fiminus1y)*length[iface];
    rhspo[ileft][3][0]=rhspo[ileft][3][0]-(fiplus2+fiminus2)*length[iface];
	}
    
    if(iright<=nelem)
    {
    rhspo[iright][0][0]=rhspo[iright][0][0]+(fiplus0+fiminus0)*length[iface];
    rhspo[iright][1][0]=rhspo[iright][1][0]+(fiplus1x+fiminus1x)*length[iface];
    rhspo[iright][2][0]=rhspo[iright][2][0]+(fiplus1y+fiminus1y)*length[iface];
    rhspo[iright][3][0]=rhspo[iright][3][0]+(fiplus2+fiminus2)*length[iface];
	}
	
}



}


