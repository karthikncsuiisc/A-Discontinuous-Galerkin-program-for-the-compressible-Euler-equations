//Program contributes diffusive terms to rhspo using BR2 method

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

void rhs_contrb_diff_br2_p1()
{

/*if(disc_order!=1)
{
	cout<<"BR2 method is implemented only for P1 method"<<endl;
	exit(1);
}*/

double ui,uj;
double fluxleft,fluxright;
double M11,M22,ujump,B2,dx;

for(int ielem=2;ielem<=nelem+1;ielem++)
{
	M11=coord_1d[ielem+1]-coord_1d[ielem];
	M22=(coord_1d[ielem+1]-coord_1d[ielem])/3.0;
	
	ui=unkel[ielem][0][0]-unkel[ielem][0][1];
	uj=unkel[ielem-1][0][0]+unkel[ielem-1][0][1];
	ujump=uj-ui;
	B2=-1;

	liftopr[ielem][0][0]=-etabr2*ujump/(2.0*M11);
	liftopr[ielem][0][1]=-etabr2*ujump*B2/(2.0*M22);

	ui=unkel[ielem][0][0]+unkel[ielem][0][1];
	uj=unkel[ielem+1][0][0]-unkel[ielem+1][0][1];
	ujump=uj-ui;
	B2=1;

	liftopr[ielem][1][0]=-etabr2*ujump/(2.0*M11);
	liftopr[ielem][1][1]=-etabr2*ujump*B2/(2.0*M22);
}
//cout<<"-------------------------"<<endl;

for(int ielem=2;ielem<=nelem+1;ielem++)
{
	dx=coord_1d[ielem+1]-coord_1d[ielem];
	
	B2=-1;
	fluxleft=unkel[ielem][0][1]/dx-(liftopr[ielem][0][0]+B2*liftopr[ielem][0][1]);
	
	B2=1;
	fluxright=unkel[ielem][0][1]/dx-(liftopr[ielem][1][0]+B2*liftopr[ielem][1][1]);
	
	rhspo[ielem][0][0]=rhspo[ielem][0][0]+coefc*(fluxright-fluxleft);
	
	B2=-1;
	fluxleft=unkel[ielem][0][1]/dx-(liftopr[ielem][0][0]+B2*liftopr[ielem][0][1]);
	fluxleft=B2*(fluxleft);
	
	B2=1;
	fluxright=unkel[ielem][0][1]/dx-(liftopr[ielem][1][0]+B2*liftopr[ielem][1][1]);
	fluxright=B2*fluxright;
	
	rhspo[ielem][0][1]=rhspo[ielem][0][1]+coefc*(fluxright-fluxleft);
	
	rhspo[ielem][0][1]=rhspo[ielem][0][1]
					  -2.0*coefc*(unkel[ielem][0][1]/dx-liftopr[ielem][0][0]-liftopr[ielem][1][0]);

//	cout<<ielem<<" "<<rhspo[ielem][0][0]<<" "<<rhspo[ielem][0][1]<<endl;
}


}
