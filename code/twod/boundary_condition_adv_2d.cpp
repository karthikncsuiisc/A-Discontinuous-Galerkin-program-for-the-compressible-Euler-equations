//Function to implement boundary condition.

#include "functions.h"
#include "twod_header.h"

void boundary_condition_adv_2d()
{

int ileft,iright,jleft;
double adotn;

for(int iface=1;iface<=nbface;iface++)
{
	ileft=intface[iface][1];
	iright=intface[iface][2];
	
	
	if(bound_cond[iface][1]==2.0) //Wall boundary condition
	{
		unkel[iright][0][0]=-unkel[ileft][0][0];
	}
	else if(bound_cond[iface][1]==4.0) //Far field boundary condition
	{
		adotn=a0x*nx[iface]+a0y*ny[iface];
		if(adotn<0) unkel[iright][0][0]=bound_cond[iright][2];
		else if(adotn>0) unkel[iright][0][0]=unkel[ileft][0][0];
	}
	else if(bound_cond[iface][1]==6.0) //Periodic boundary condition
	{
		for(int jface=1;jface<=nbface;jface++)
		{
			if(jface!=iface)
			{
			if(bound_cond[jface][3]==bound_cond[iface][3])
			{
				jleft=intface[jface][1];
				unkel[iright][0][0]=unkel[jleft][0][0];
			}
			}
		}
	}
}

}
