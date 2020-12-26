//This function writes output file containing phi, pressure and velocity vector
//in a VTK format and wrties the x,y,velocity values into a text format for
//comparing the results with exact solution.  

#include "functions.h"
#include "twod_header.h"

void results2D()
{
string filename;

ofstream fout;

if(phytype==1) error_function();

double u1,u2x,u2y,u3,pres0,velx0,vely0,rho0,Mach;

filename=grid_name+"_results.vtk";
fout.open(filename.c_str(), ios::out | ios::trunc);
cout<<"Creating 2D results data file :"<<filename<<endl;
cout<<"writing the geometry information"<<endl;

fout<<"# vtk DataFile Version 2.0"<<endl;
fout<<"2D Unstructured Grid of Linear Triangles"<<endl;
fout<<"ASCII\n"<<endl;
fout<<"DATASET UNSTRUCTURED_GRID"<<endl;
fout<<"POINTS "<<npoint<<" float"<<endl;
for(int j=1;j<=npoint;j++)
{
fout<<coord[j][1]<<" "<<coord[j][2]<<" 0.00"<<endl;	
}

fout<<"\nCells "<<nelem<<" "<<4*(nelem)<<endl;

for(int j=1;j<=nelem;j++)
{
fout<<"3 "<<inpoel[j][1]-1<<" "<<inpoel[j][2]-1<<" "<<inpoel[j][3]-1<<endl;	
}
fout<<"\nCELL_TYPES "<<nelem<<endl;

for(int j=1;j<=nelem;j++)
{
fout<<"5"<<endl;	
}

//-----------------------------------------------------------------------------
// 					Results for Euler equations
//-----------------------------------------------------------------------------

if(phytype==1)
{

double **velnode;
double *Machnode;
double *presnode;
double *rhonode;
double **velelem;
double *Machelem;
double *preselem;
double *rhoelem;

velnode=new double* [npoint+1];
Machnode=new double[npoint+1];
presnode=new double[npoint+1];
rhonode=new double[npoint+1];

velelem=new double* [nelem+1];
Machelem=new double[nelem+1];
preselem=new double[nelem+1];
rhoelem=new double[nelem+1];

for(int ipoint=1;ipoint<=npoint;ipoint++) velnode[ipoint] =new double[2];
for(int ielem=1;ielem<=nelem;ielem++) velelem[ielem] =new double[2];

double we;
double *wetot;
wetot = new double[npoint+1];

int ind,nnode=3;
double xp,yp,dist;

for (int ipoint=1;ipoint<=npoint;ipoint++)
{
	velnode[ipoint][0]=0;
	velnode[ipoint][1]=0;
	Machnode[ipoint]=0;
	presnode[ipoint]=0;
	rhonode[ipoint]=0;
	wetot[ipoint]=1e-16;
}

cout<<"Calculating the elemental and nodal values"<<endl;
//Inner cells
for(int ielem=1;ielem<=nelem;ielem++)
{
    u1=unkel[ielem][0][0];
	u2x=unkel[ielem][1][0];
	u2y=unkel[ielem][2][0];
	u3=unkel[ielem][3][0];

	rho0=u1;
	velx0=u2x/u1;
	vely0=u2y/u1;
	pres0=(gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0);
	Mach=sqrt(velx0*velx0+vely0*vely0)/sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	
	velelem[ielem][0]=velx0;
	velelem[ielem][1]=vely0;
	Machelem[ielem]=Mach;
	preselem[ielem]=pres0;
	rhoelem[ielem]=rho0;
	
	for(int inode=1;inode<=nnode;inode++)
	{
		ind=inpoel[ielem][inode];
		xp=coord[ind][1];
		yp=coord[ind][2];
		if(disc_order==1)
		{
	    u1=unkel[ielem][0][0]+unkel[ielem][0][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][0][2]*(yp-centery[ielem])/deltay[ielem];
							 
   	   u2x=unkel[ielem][1][0]+unkel[ielem][1][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][1][2]*(yp-centery[ielem])/deltay[ielem];
							 
	   u2y=unkel[ielem][2][0]+unkel[ielem][2][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][2][2]*(yp-centery[ielem])/deltay[ielem];
							 
		u3=unkel[ielem][3][0]+unkel[ielem][3][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][3][2]*(yp-centery[ielem])/deltay[ielem];
		}
		else if(disc_order==0)
		{
	    u1=unkel[ielem][0][0];	    					 
   	   u2x=unkel[ielem][1][0]; 
	   u2y=unkel[ielem][2][0];
	   	u3=unkel[ielem][3][0];		
		}

		rho0=u1;
		velx0=u2x/u1;
		vely0=u2y/u1;
		pres0=(gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0);
		Mach=sqrt(velx0*velx0+vely0*vely0)/sqrt(gama*max(pres0,1e-16)/rho0);
	
		dist=sqrt(pow(centerx[ielem]-xp,2)+pow(centery[ielem]-yp,2));
		we=1/(dist);
		
		wetot[ind]=wetot[ind]+we;
		velnode[ind][0]=velnode[ind][0]+velx0*we;
		velnode[ind][1]=velnode[ind][1]+vely0*we;
		
		presnode[ind]=presnode[ind]+pres0*we;
		rhonode[ind]=rhonode[ind]+rho0*we;
		Machnode[ind]=Machnode[ind]+Mach*we;
	}
}
cout<<"Calculating the nodal values from ghost cells"<<endl;

//Ghost cells
if(disc_order==0)
{
for(int iface=1;iface<=nbface;iface++)
{
	if(bound_cond[iface][1]!=4.0)
	{
		
	int ielem=intface[iface][2];
    u1=unkel[ielem][0][0];
	u2x=unkel[ielem][1][0];
	u2y=unkel[ielem][2][0];
	u3=unkel[ielem][3][0];

	rho0=u1;
	velx0=u2x/u1;
	vely0=u2y/u1;
	pres0=max((gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0),1e-8);
	Mach=sqrt(velx0*velx0+vely0*vely0)/sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
	
	for(int inode=1;inode<=2;inode++)
	{
		ind=intface[iface][inode+2];
		xp=coord[ind][1];
		yp=coord[ind][2];
	
		dist=sqrt(pow(centerx[ielem]-xp,2)+pow(centery[ielem]-yp,2));
		we=1/(dist);
		
		wetot[ind]=wetot[ind]+we;
		velnode[ind][0]=velnode[ind][0]+velx0*we;
		velnode[ind][1]=velnode[ind][1]+vely0*we;
		
		presnode[ind]=presnode[ind]+pres0*we;
		rhonode[ind]=rhonode[ind]+rho0*we;
		Machnode[ind]=Machnode[ind]+Mach*we;
	}
	
	}
}
}
else if(disc_order==1)
{
for(int iface=1;iface<=nbface;iface++)
{
	if(bound_cond[iface][1]==2.0)
	{
	int ielem=intface[iface][1];
	
	for(int inode=1;inode<=2;inode++)
	{
		ind=intface[iface][inode+2];
		xp=coord[ind][1];
		yp=coord[ind][2];
	
		dist=sqrt(pow(centerx[ielem]-xp,2)+pow(centery[ielem]-yp,2));
		we=1/(dist);

		u1=unkel[ielem][0][0]+unkel[ielem][0][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][0][2]*(yp-centery[ielem])/deltay[ielem];
		u2x=unkel[ielem][1][0]+unkel[ielem][1][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][1][2]*(yp-centery[ielem])/deltay[ielem];
		u2y=unkel[ielem][2][0]+unkel[ielem][2][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][2][2]*(yp-centery[ielem])/deltay[ielem];
		u3=unkel[ielem][3][0]+unkel[ielem][3][1]*(xp-centerx[ielem])/deltax[ielem]
							 +unkel[ielem][3][2]*(yp-centery[ielem])/deltay[ielem];
		rho0=u1;
		velx0=u2x/u1;
		vely0=u2y/u1;
		pres0=(gama-1)*(u3-0.5*(velx0*velx0+vely0*vely0)*rho0);
	
		velx0=velx0-2*(velx0*nx[iface]+vely0*ny[iface])*nx[iface];
		vely0=vely0-2*(velx0*nx[iface]+vely0*ny[iface])*ny[iface];

		Mach=sqrt(velx0*velx0+vely0*vely0)/sqrt(gama*max(pres0,hrd_lmtr_cin)/rho0);
		
		wetot[ind]=wetot[ind]+we;
		velnode[ind][0]=velnode[ind][0]+velx0*we;
		velnode[ind][1]=velnode[ind][1]+vely0*we;
		
		presnode[ind]=presnode[ind]+pres0*we;
		rhonode[ind]=rhonode[ind]+rho0*we;
		Machnode[ind]=Machnode[ind]+Mach*we;
	}
	
	}
}

	
}

for(int i=1;i<=npoint;i++)
{
rhonode[i]=rhonode[i]/wetot[i];
presnode[i]=presnode[i]/wetot[i];
Machnode[i]=Machnode[i]/wetot[i];

velnode[i][0]=velnode[i][0]/wetot[i];
velnode[i][1]=velnode[i][1]/wetot[i];
}
cout<<"Writing the results at cell centers"<<endl;

					//--------------------------------------------------
					//		  Writing results on cell centers
					//--------------------------------------------------

fout<<"\nCELL_DATA "<<nelem<<endl;
fout<<"SCALARS density float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ielem=1;ielem<=nelem;ielem++)
fout<<rhoelem[ielem]<<endl;
fout<<endl;

fout<<"SCALARS Mach float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ielem=1;ielem<=nelem;ielem++)
fout<<Machelem[ielem]<<endl;
fout<<endl;

fout<<"SCALARS pressure float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ielem=1;ielem<=nelem;ielem++)
fout<<preselem[ielem]<<endl;
fout<<endl;

cout<<"Writing the results at cell centers"<<endl;

fout<<"SCALARS entropy float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ielem=1;ielem<=nelem;ielem++)
fout<<entropy[ielem]<<endl;
fout<<endl;

fout<<"SCALARS entropyprod float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ielem=1;ielem<=nelem;ielem++)
fout<<(entropy[ielem]-presref/pow(rhoref,gama))<<endl;
fout<<endl;

fout<<"VECTORS velocity float"<<endl;
//fout<<"LOOKUP_TABLE default"<<endl;
for(int ielem=1;ielem<=nelem;ielem++)
fout<<velelem[ielem][0]<<" "<<velelem[ielem][1]<<" 0"<<endl;
fout<<endl;
					//--------------------------------------------------
					//	End of writing results on cell centers
					//--------------------------------------------------
cout<<"Writing the results at node points"<<endl;

					//--------------------------------------------------
					//	Writing results on node points
					//--------------------------------------------------
fout<<"\nPOINT_DATA "<<npoint<<endl;
fout<<"SCALARS densitynode float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ipoint=1;ipoint<=npoint;ipoint++)
fout<<rhonode[ipoint]<<endl;
fout<<endl;

fout<<"SCALARS Machnode float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ipoint=1;ipoint<=npoint;ipoint++)
fout<<Machnode[ipoint]<<endl;
fout<<endl;

fout<<"SCALARS pressurenode float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ipoint=1;ipoint<=npoint;ipoint++)
fout<<presnode[ipoint]<<endl;
fout<<endl;

fout<<"VECTORS velocitynode float"<<endl;
//fout<<"LOOKUP_TABLE default"<<endl;
for(int ipoint=1;ipoint<=npoint;ipoint++)
fout<<velnode[ipoint][0]<<" "<<velnode[ipoint][1]<<" 0"<<endl;
fout<<endl;
					//--------------------------------------------------
					//	End of writing results on node points
					//--------------------------------------------------
}
//-----------------------------------------------------------------------------
// 					End of Results for Euler equations
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//				    Results for Advection equation
//-----------------------------------------------------------------------------

else if(phytype==0)
{
//Inner cells
double *unode;
double u0;
unode =new double[npoint+1];

double we,xp,yp,dist;
int ind;
double *wetot;
wetot = new double[npoint+1];

for (int ipoint=1;ipoint<=npoint;ipoint++)
{
	unode[ipoint]=0.0;
	wetot[ipoint]=0.0;
}

for(int ielem=1;ielem<=nelem;ielem++)
{
    u0=unkel[ielem][0][0];
	
	for(int inode=1;inode<=3;inode++)
	{
		ind=inpoel[ielem][inode];
		xp=coord[ind][1];
		yp=coord[ind][2];
	
		dist=sqrt(pow(centerx[ielem]-xp,2)+pow(centery[ielem]-yp,2));
		we=1/(dist);
		
		wetot[ind]=wetot[ind]+we;
		unode[ind]=unode[ind]+u0*we;
	}
}

//Ghost cells
for(int iface=1;iface<=nbface;iface++)
{
	int ielem=intface[iface][2];
    u0=unkel[ielem][0][0];
	
	for(int inode=1;inode<=2;inode++)
	{
		ind=intface[iface][inode+2];
		xp=coord[ind][1];
		yp=coord[ind][2];
	
		dist=sqrt(pow(centerx[ielem]-xp,2)+pow(centery[ielem]-yp,2));
		we=1/(dist);
		wetot[ind]=wetot[ind]+we;
		
		if(bound_cond[iface][1]==2.0) unode[ind]=unode[ind]-u0*we;
		else unode[ind]=unode[ind]+u0*we;
	}
	
}

for(int i=1;i<=npoint;i++) unode[i]=unode[i]/wetot[i];

					//--------------------------------------------------
					//		  Writing results on cell centers
					//--------------------------------------------------

fout<<"\nCELL_DATA "<<nelem<<endl;
fout<<"SCALARS u float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ielem=1;ielem<=nelem;ielem++)
fout<<unkel[ielem][0][0]<<endl;
fout<<endl;

					//--------------------------------------------------
					//	End of writing results on cell centers
					//--------------------------------------------------

					//--------------------------------------------------
					//	Writing results on node points
					//--------------------------------------------------
fout<<"\nPOINT_DATA "<<npoint<<endl;
fout<<"SCALARS unode float"<<endl;
fout<<"LOOKUP_TABLE default"<<endl;
for(int ipoint=1;ipoint<=npoint;ipoint++)
fout<<unode[ipoint]<<endl;
fout<<endl;
					//--------------------------------------------------
					//	End of writing results on node points
					//--------------------------------------------------

}

fout.close();

//Deleting the arrays
	for(int ielem=1;ielem<=nelem+nface; ielem++) for(int ivar=0;ivar<nvar;ivar++) delete [] unkel[ielem][ivar];
	for(int ielem=1;ielem<=nelem+nface; ielem++) delete [] unkel[ielem];
	delete [] unkel;
	
	delete [] centerx;
	delete [] centery;
	delete [] nx;
	delete [] ny;
	delete [] area;
	
	for(int ielem=0;ielem<=nelem; ielem++) delete [] elsuel[ielem];
	delete [] elsuel;
	
	for (int j=1; j <=nelem+nface; j++) delete [] inpoel[j];
	delete [] inpoel;
	
	for (int j=1; j <=npoint+nface; j++) delete [] coord[j];
	delete [] coord;

	for (int j=1; j <=npoint; j++) delete [] ini_cond[j];
	delete [] ini_cond; 
cout<<"End of Creating 2D results data file :"<<filename<<endl;

//------------------------------------------------------------------
//End of writing data file
//-----------------------------------------------------------------
}
