#ifndef HEADER_H
#define HEADER_H
#include <math.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cmath>

using namespace std;

//------------Geometry variables----------------------------------
extern int spacedim; 				//Space Dimenstion type 1D/2D/3D
extern double leftcoord; 		    //Left Coordinate
extern double rightcoord; 			//Right Coordinate 
extern int nelem, nelemtot;			//Number of elements in the mesh
extern int nfiles;
extern string grid_names[];            //Mesh file name for 2d
extern string grid_name;
extern int ndimn; 					//Dimension of the problem
extern int ntype; 					//Type problem

//-----------Physics variables-------------------------------------
extern int phytype; 											//Physics Type
extern double a0x,a0y;                   					    //Advection equation speed
extern int nvar;										    	//Discretization order variables
extern double ***unkel;  										//Unknowns
extern double **flux;											//Flux varaibles
extern double ***rhspo;											//right hand side
extern double rhoref,velref,Machref,gama,alpha,presref;			//Reference values used for euler euqation
extern double hrd_lmtr_cin;										//Hard limiter, used to limit the pressure to positive values when calcualtiong the speed of sound. This may happen for p1, p2 methods
extern double coefa,coefb,coefc;								//Advection equation coefficient
extern double ***liftopr;										//Life operator variables used for diffusion problem
extern double etabr2;
extern double *entropy;
//----------Discretization Variables-------------------------------
extern int disc_type; 				//Discretization type
extern int disc_order; 				//Discretization order
extern int recons;
extern int reconsvariable;
extern int limiter;
extern int limitertype;
extern int ndegr;
extern int basis_func; 				//Basis Function type
extern int flux_type;
extern int diff_type;				//Indicator of diffusion discretization

//----------Boundary condition variables---------------------
extern int bcl; 					//Left boundary condition
extern int bcr;						//Right boundary condition

//---------Study type and time stepping----------------------
extern int studytype;				//Study type
extern int steadystate_type;		//Steady state type
extern int timedep_type;			//Time Dependent type
extern double t0;					//Start time
extern double tf;					//Final time
extern double dt;					//Time Step
extern int CFL_use;                 //Option for using CFL
extern double CFL;					//CFL Value
extern int timedep_method;			//Time Dependent Method
extern int Niter;					//Maximum Number of iterations
extern double tol;
extern int timstp_option;			//Time step option
				 
void read_input();					//Reading the input file input_file.yaml
void phymanager(); 					//Calls corresponding physics
double minmod(double,double);
double MClimiter(double,double,double);
#endif







