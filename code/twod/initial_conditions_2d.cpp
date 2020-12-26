//Program for applying Initial Conditions

#include "functions.h"
#include "twod_header.h"

void initial_conditions_2d()
{
//----------------------------------------------------------
// Assigning initial conditions
//----------------------------------------------------------

int ip1,ip2,ip3;
double velxref,velyref;
velxref=velref*cos(alpha*3.14159265/180.0);
velyref=velref*sin(alpha*3.14159265/180.0);


if(phytype==0)
{
for(int ielem=1;ielem<=nelem;ielem++)
{
	ip1=inpoel[ielem][1];
	ip2=inpoel[ielem][2];
	ip3=inpoel[ielem][3];
	
	unkel[ielem][0][0]=(ini_cond[ip1][0]+ini_cond[ip2][0]+ini_cond[ip3][0])/3.0;
}
}
else if(phytype==1)
{
for(int ielem=1;ielem<=nelem;ielem++)
{
	unkel[ielem][0][0]=rhoref;
	unkel[ielem][1][0]=rhoref*velxref;
	unkel[ielem][2][0]=rhoref*velyref;
	unkel[ielem][3][0]=presref/(gama-1.0)+0.5*rhoref*(velxref*velxref+velyref*velyref);
}

	
}
}

