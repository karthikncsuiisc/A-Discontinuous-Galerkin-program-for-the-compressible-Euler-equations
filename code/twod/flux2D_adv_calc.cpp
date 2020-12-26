//Subroutine for finding elements surrounding point.

#include "functions.h"
#include "twod_header.h"

void flux2D_adv_calc()
{
int ileft,iright,ip1,ip2;
double nx,ny,adotn,fluxval;

for(int iface=1;iface<=naface;iface++)
{
	ileft=intface[iface][1];
	iright=intface[iface][2];
	ip1=intface[iface][3];
	ip2=intface[iface][4];
	
	nx=coord[ip2][2]-coord[ip1][2];
	ny=-coord[ip2][1]+coord[ip1][1];
	
	adotn=a0x*nx+a0y*ny;
	
	if(adotn>=0) fluxval=unkel[ileft][0][0]*adotn;
	else fluxval=unkel[iright][0][0]*adotn;
	
	if(iright<=nelem) rhspo[iright][0][0]=rhspo[iright][0][0]+fluxval;
	if(ileft<=nelem) rhspo[ileft][0][0]=rhspo[ileft][0][0]-fluxval;
}

}
