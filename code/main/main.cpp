// A C++ program for solving Different physics. 
// The Program is arranged into different subroutines. All the functions
// are declared in the header file "functions.h". A make file is aviable
// to compile all the subroutines.
//

#include "functions.h"
#include "oned_header.h"
#include "twod_header.h"

#include <time.h>
#include <sys/time.h>

//------------Geometry variables----------------------------------
int spacedim; 				//Space Dimenstion type 1D/2D/3D
double leftcoord; 		    //Left Coordinate
double rightcoord; 			//Right Coordinate 
int nelem, nelemtot;			//Number of elements in the mesh
int nfiles;
string grid_names[20];            //Mesh file name for 2d
string grid_name;
int ndimn; 					//Dimension of the problem
int ntype; 					//Type problem

//-----------Physics variables-------------------------------------
int phytype; 				//Physics Type
double a0x,a0y;                     //Advection equation speed
int nvar;
double ***unkel;  			
double **flux;
double ***rhspo;
double **ini_cond;
double rhoref,velref,Machref,gama,alpha,presref;
double hrd_lmtr_cin;
double coefa,coefb,coefc;
double ***liftopr;					
double etabr2;
double *entropy;


//----------Discretization Variables-------------------------------
int disc_type; 				//Discretization type
int disc_order; 				//Discretization order
int recons;
int reconsvariable;
int limiter;
int limitertype;
int ndegr;					//ndegr=disc_order+1
int basis_func; 				//Basis Function type
int flux_type;
int diff_type;				//Indicator of diffusion discretization

//----------Boundary condition variables---------------------
int bcl; 					//Left boundary condition
int bcr;						//Right boundary condition

//---------Study type and time stepping----------------------
int studytype;				//Study type
int steadystate_type;		//Steady state type
int timedep_type;			//Time Dependent type
double t0;					//Start time
double tf;					//Final time
double dt;					//Time Step
int CFL_use;                 //Option for using CFL
double CFL;					//CFL Value
int timedep_method;			//Time Dependent Method
int Niter;					//Maximum Number of iterations
double tol;
int timstp_option;			//Time step option

//--------ONED header file variables
double *coord_1d;
int npoint_1d;

//--------TWOD header file variables
int npoint; 				// Number of nodes in the mesh
int nface; 				// Number of boundary faces
int **inpoel; 			// Node numbers of each element
double **coord; 			// Coordinates of node points
double **bface; 			// Boundary elements, corresponding points and 
double **lhspo; 			// System matrix

//-----ELements surrounding points variables----------------
int *esup1;
int *esup2;

//-----Points surrounding points variables----------------
int *psup1;
int *psup2;

//-----Elements surrounding elements variables----------------
int **elsuel;

//-----Internal facee variables------------------------------
int **intface;
int naface;
int nbface;
double **bound_cond;
double *nx,*ny,*length,*area;
double *deltax,*deltay;
double *centerx,*centery;

//-----Gauss point variables------------------------------
double **Bxquad,**Byquad,*Wequad;
double **Bxcubic,**Bycubic,*Wecubic;
double ***Bxline_quad,***Byline_quad,*Weline_quad;
double ***Bxline_cubic,***Byline_cubic,*Weline_cubic;
double **Bxline_lin,**Byline_lin,Weline_lin;

//---------Reconstruction variables-------------------------
double ***gradu;
double **rhsel;

 
int main()
{
	
double startcputime, endcputime,cpu_time;
startcputime = (double)clock();

read_input();
phymanager();

endcputime=(double)clock();
cpu_time=(endcputime-startcputime)/CLOCKS_PER_SEC;

cout<<"CPU Time(in sec) :"<<cpu_time<<endl;
cout<<"CPU Time(in min) :"<<cpu_time/60.0<<endl;
return 0;
}





