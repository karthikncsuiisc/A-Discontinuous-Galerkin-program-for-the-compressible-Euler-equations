//Subroutine for finding elements surrounding point.

#include "functions.h"
#include "twod_header.h"

void flux2D_adv_calc()
{
int ileft,iright,ip1,ip2;
double adotn,fluxval;

for(int iface=1;iface<=naface;iface++)
{
	
	ileft=intface[iface][1];
	iright=intface[iface][2];
	ip1=intface[iface][3];
	ip2=intface[iface][4];
	
	adotn=a0x*nx[ifcace]+a0y*ny[iface];
	fluxval=0.0;
	if(iright<=nelem)
	{
		for(int igp=1;igp<=2;igp++)
		{		
		if(adotn>=0)
		{
			u1gp=unkel[ileft][0][0]+unkel[ileft][0][1]*Bxline_quad[1][igp]
												+unkel[ileft][0][2]*Byline_quad[1][igp];
			
			fluxval=fluxval+u1gp*adotn*Weline_lin;
		}
		else
		{
			u1gp=(unkel[iright][0][0]+unkel[iright][0][1]*Bxline_quad[2][igp]
													+unkel[iright][0][2]*Byline_quad[2][igp];
		
			 fluxval=fluxval+u1gp*adotn*Weline_lin
		}
	}
	else
	{
		for(int igp=1;igp<=2;igp++)
		{	
		//Applying boundary condiiton
			if(bound_cond[iface][1]==2.0) //Wall boundary condition
			{
			u1gp=unkel[ileft][0][0]+unkel[ileft][0][1]*Bxline_quad[1][igp]
												+unkel[ileft][0][2]*Byline_quad[1][igp];
												
			fluxval=fluxval+u1gp*abs(adotn)*Weline_lin;
			}
			else if(bound_cond[iface][1]==4.0) //Far field boundary condition
			{
				if(adotn<0)
				{
					 fluxval=fluxval+bound_cond[iright][2]*abs(adotn)*Weline_lin;
				 }
				else
				{ 
					u1gp=unkel[ileft][0][0]+unkel[ileft][0][1]*Bxline_quad[1][igp]
														+unkel[ileft][0][2]*Byline_quad[1][igp];
					
					fluxval=fluxval+u1gp*abs(adotn)*Weline_lin
			    }
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
								
							    u1gp=unkel[jleft][0][0]+unkel[jleft][0][1]*Bxline_quad[jface][1][igp]
														            +unkel[ileft][0][2]*Byline_quad[jface][1][igp];
								fluxval=fluxval+u1gp*abs(adotn)*Weline_lin
							}
						}
		           }
		  }
		
	    }
	 }
	 
	 if(iright<=nelem) rhspo[iright][0][0]=rhspo[iright][0][0]+fluxval*length[iface]/2;
	if(ileft<=nelem) rhspo[ileft][0][0]=rhspo[ileft][0][0]-fluxval*length[iface]/2;
}

///Flux calculation

}
