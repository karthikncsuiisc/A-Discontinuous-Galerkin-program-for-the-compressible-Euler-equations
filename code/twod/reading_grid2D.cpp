//Subroutine for reading the mesh file for this particular format
//New subrouinte need to be written for different format

#include "functions.h"
#include "twod_header.h"

void reading_grid2D()
{

int temp=1;	
std::string line;
ifstream fin;

fin.open(grid_name.c_str());
    // Check for successful opening of file
    if(fin.fail())
    {
       cout<<"Error Opening File :"<<grid_name<<endl;
       exit(1);
    }


while(1>0)
{
fin>>line;
if(line=="ndimn") break;
}
getline(fin,line);
fin>>ndimn>>ntype;
getline(fin,line);
getline(fin,line);
fin>>nelem>>npoint>>nface;
getline(fin,line);
getline(fin,line);

cout<<"-----------------------------------------"<<endl;
cout<<"Number of Elements:"<<nelem<<"\nNumber of Points: "<<npoint<<"\nNumber of Boundary faces "<<nface<<endl;

//-----------------------------------------------------------------
//Initializing element connectivity array and reading from the file
inpoel = new int* [nelem+1+nface];
for (int j=1; j <=nelem+nface; j++) inpoel[j] = new int[4];
 
for(int j=1;j<=nelem;j++)
{
fin>>temp>>inpoel[j][1]>>inpoel[j][2]>>inpoel[j][3];
getline(fin,line);
}
getline(fin,line);

//Initializing coordinates of the points array and reading from the file
coord=new double* [npoint+1+nface];
for (int j=1; j <=npoint+nface; j++) coord[j] = new double[3];
for(int j=1;j<=npoint;j++)
{
fin>>temp>>coord[j][1]>>coord[j][2];
}
getline(fin,line);
getline(fin,line);

//End of the coordinate section
//--------------------------------------------------------------------
//Initializing the initial condition array and reading the values from file

ini_cond=new double* [npoint+1]; 
for (int j=1; j <=npoint; j++) ini_cond[j] = new double[4];

for(int j=1;j<=npoint;j++)
{
fin>>temp;
for(int i=0;i<4;i++) fin>>ini_cond[j][i]; 
getline(fin,line);
}
getline(fin,line);
//End of initial condition

//---------------------------------------------------------------------
// Initalizing the boundary faces vector and reading from the file
bface=new double* [nface+1]; 
for (int j=1; j <=nface; j++) bface[j] = new double[6];
for(int j=1;j<=nface;j++)
{
fin>>temp;
for(int i=1;i<=5;i++) 
{
	fin>>bface[j][i];
//	getline(fin,line);
}
getline(fin,line);
}

fin.close();

cout<<"End of reading grid: "<<grid_name<<endl;

}
