/*Program to create a square grid*/
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <cstdlib>
using namespace std;

int main()
{

int Nx,Ny,nelem,npoint,nface;
double lengthx,lengthy,hx,hy,X0,Y0;

lengthx=1.0;
lengthy=0.5;

cout<<"Enter the value of Nx: ";
cin>>Nx;
cout<<"Enter the value of Ny: ";
cin>>Ny;

hx=lengthx/(double)Nx;
hy=lengthy/(double)Ny;

X0=0;
Y0=0;

string filename ("shocktube_");
ostringstream convert;
convert << Nx;
filename=filename+convert.str();
filename=filename+'_';
ostringstream convert2;
convert2 << Ny;
filename=filename+convert2.str();
filename=filename+".dat";

fstream fout;
fout.open(filename.c_str(), ios::out | ios::trunc);

fout<<"3"<<endl;

fout<<"-----------------------------------------------------"<<endl;
fout<<"mesh in a channel with a circular bump on the bottom"<<endl;
fout<<"----------------------------------------------------"<<endl;                           
fout<<" ndimn ntype"<<endl;
fout<<"     2     3"<<endl;
fout<<" nelem npoin nface           time"<<endl;

nelem=2*Nx*Ny;
npoint=(Nx+1)*(Ny+1);
nface=2*(Nx+Ny);
fout<<"  "<<nelem<<"   "<<npoint<<"   "<<nface<<"     0.0000E+00"<<endl;
fout<<" intmat: nodal points corresponding to each element"<<endl;

double **inpoel;
inpoel=new double* [nelem+1];
for(int ipoint=1;ipoint<=nelem;ipoint++) inpoel[ipoint]=new double[4];

int ielem=1;

for(int j=1;j<=Ny;j++)
for(int i=1;i<=Nx;i++)
{	
	inpoel[ielem][1]=(j-1)*(Nx+1);
	inpoel[ielem][2]=j*(Nx+1)+i+1;
	inpoel[ielem][3]=j*(Nx+1)+i;
		
	fout<<ielem<<" "<<(j-1)*(Nx+1)+i<<" "<<j*(Nx+1)+i+1<<" "<<j*(Nx+1)+i<<"   0  0  0  0  0  0  0  1"<<endl;
	ielem++;

	inpoel[ielem][1]=(j-1)*(Nx+1)+i;
	inpoel[ielem][2]=(j-1)*(Nx+1)+i+1;
	inpoel[ielem][3]=j*(Nx+1)+i+1;

	fout<<ielem<<" "<<(j-1)*(Nx+1)+i<<" "<<(j-1)*(Nx+1)+i+1<<" "<<j*(Nx+1)+i+1<<"   0  0  0  0  0  0  0  1"<<endl;
	ielem++;
}

fout<<"coordinates of the points"<<endl;
double **coord;
coord=new double* [npoint+1];
for(int ipoint=1;ipoint<=npoint;ipoint++) coord[ipoint]=new double[4];
int inode=1;

for(int j=1;j<=Ny+1;j++)
for(int i=1;i<=Nx+1;i++)
{
fout<<inode<<" "<<X0+(i-1)*hx<<" "<<Y0+(j-1)*hy<<endl;
coord[inode][1]=X0+(i-1)*hx;
coord[inode][2]=Y0+(j-1)*hy;
inode++;
}

inode=1;
fout<<"initial values for the unknowns and velocity field"<<endl;
double unkel,x1;
for(int j=1;j<=Ny+1;j++)
for(int i=1;i<=Nx+1;i++)
{
	x1=coord[inode][1];
//	y1=coord[(j-1)*Nx+i][2];
	
	if(x1>=0&&x1<0.6) unkel=exp(-200*pow(x1-0.3,2));
	else if(x1>=0.6&&x1<=0.8) unkel=1;
	else unkel=0;
	
fout<<inode<<" "<<unkel<<"  0.0  0.0  0.0  0.0  0.0"<<endl;	
inode++;
}
fout<<" boundary faces and conditions"<<endl;

int boun_num=1;
for(int i=1;i<=Nx;i++)
{
	fout<<boun_num<<" "<<i<<" "<<i+1<<" 2  0  1"<<endl;
	boun_num++; 
}

int per_bound[Ny+1];
for(int i=0;i<=Ny;i++) per_bound[i]=0;

for(int j=1;j<=Ny;j++)
{
	per_bound[j]=per_bound[j-1]+1;
	fout<<boun_num<<" "<<(Nx+1)*j<<" "<<(Nx+1)*(j+1)<<" 6  0  "<<per_bound[j]<<endl;
	boun_num++;
}

for(int i=1;i<=Nx;i++)
{
	fout<<boun_num<<" "<<(Nx+1)*(Ny+1)-i+1<<" "<<(Nx+1)*(Ny+1)-i<<" 2  0  1"<<endl;
	boun_num++; 
}

for(int j=Ny;j>=1;j--)
{
	fout<<boun_num<<" "<<j*(Nx+1)+1<<" "<<(Nx+1)*(j-1)+1<<" 6  0  "<<per_bound[j]<<endl;
	boun_num++;
}

fout.close();
return(0);
}
