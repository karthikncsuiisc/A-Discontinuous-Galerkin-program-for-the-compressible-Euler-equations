//Function to write results for oned
// The cells are of equal length. The function gives the results as an output file: results_1d_wave.dat
//
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

void onedresults()
{
	
cout<<"writing results to the file : results.dat"<<endl;
//Writing results at the nodes or element centers(yesnode=0 on element, 1 on nodes)
int yesnode=1;

if(yesnode==1)
{
ofstream fout;
fout.open("results.dat", ios::out | ios::trunc);

for(int node=3;node<=nelem+1;node++)
{
double xval;
double rho_dum,vel_dum,pres_dum;
double u1,u2,u3;

xval=coord_1d[node];

	if(disc_order==0)
	{
	u1=unkel[node][0][0];
	u2=unkel[node][1][0];
	u3=unkel[node][2][0];
    }
    else if(disc_order==1)
	{
	u1=unkel[node][0][0]-unkel[node][0][1];
	u2=unkel[node][1][0]-unkel[node][1][1];
	u3=unkel[node][2][0]-unkel[node][2][1];
	
    }
    else if(disc_order==2)
	{
	u1=unkel[node][0][0]-unkel[node][0][1]+unkel[node][0][2]/3.0;
	u2=unkel[node][1][0]-unkel[node][1][1]+unkel[node][1][2]/3.0;
	u3=unkel[node][2][0]-unkel[node][2][1]+unkel[node][2][2]/3.0;
    }
    
    rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=max((gama-1)*(u3-0.5*pow(u2,2)/u1),1e-8);

//fout<<xval<<" "<<rho_dum<<" "<<vel_dum<<" "<<pres_dum<<endl;

    if(disc_order==0)
    {
	u1=unkel[node-1][0][0];
	u2=unkel[node-1][1][0];
	u3=unkel[node-1][2][0];
    }
    else if(disc_order==1)
    {
	u1=unkel[node-1][0][0]+unkel[node-1][0][1];
	u2=unkel[node-1][1][0]+unkel[node-1][1][1];
	u3=unkel[node-1][2][0]+unkel[node-1][2][1];
    }
    else if(disc_order==2)
    {
	u1=unkel[node-1][0][0]+unkel[node-1][0][1]+unkel[node-1][0][2]/3.0;
	u2=unkel[node-1][1][0]+unkel[node-1][1][1]+unkel[node-1][1][2]/3.0;
	u3=unkel[node-1][2][0]+unkel[node-1][2][1]+unkel[node-1][2][2]/3.0;
    }	

	rho_dum=rho_dum+u1;
	vel_dum=vel_dum+u2/u1;
	pres_dum=pres_dum+max((gama-1)*(u3-pow(u2,2)/2/u1),1e-8);

fout<<xval<<" "<<rho_dum/2.0<<" "<<vel_dum/2.0<<" "<<pres_dum/2.0<<endl;
}
}
else if(yesnode==0)
{
ofstream fout;
fout.open("results.dat", ios::out | ios::trunc);

for(int i=2;i<=nelem+1;i++)
{
double xval=0.5*(coord_1d[i]+coord_1d[i+1]);
double rho_dum,vel_dum,pres_dum;
double u1,u2,u3;

	u1=unkel[i][0][0];
	u2=unkel[i][1][0];
	u3=unkel[i][2][0];

	rho_dum=u1;
	vel_dum=u2/u1;
	pres_dum=(gama-1)*(u3-pow(u2,2)/2/u1);

fout<<xval<<" "<<rho_dum<<" "<<vel_dum<<" "<<pres_dum<<endl;
}

}

}
