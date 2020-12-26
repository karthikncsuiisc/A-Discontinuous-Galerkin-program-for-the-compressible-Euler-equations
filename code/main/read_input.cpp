//Subroutine reads the input file.
//This subroutine reads the parameter from the input files which will 
//be used by the program. To include more parameters in the input file
//change the program apropriately

#include "functions.h"

void read_input()
{
ifstream fin;
string line;
string line1;

fin.open("input_file.yaml");

    // Check for successful opening of file
    if(fin.fail())
    {
       cout<<"Error Opening File :"<<"input_file.yaml"<<endl;
       exit(1);
    }

while(1>0)
{
getline(fin,line);
if(line=="start spacedim") break;	
}
fin>>line1>>spacedim;
cout<<line1<<spacedim<<endl;

while(1>0)
{
getline(fin,line);
if(line=="start geom1d") break;	
}
fin>>line1>>leftcoord;
fin>>line1>>rightcoord;
fin>>line1>>nelem;
cout<<nelem<<endl;

while(1>0)
{
getline(fin,line);
if(line=="start geom2d") break;	
}
fin>>line1>>nfiles;
fin>>line1;

for(int i=0;i<nfiles;i++) fin>>grid_names[i];

while(1>0)
{
getline(fin,line);
if(line=="start phy_type") break;	
}
fin>>line1>>phytype;
fin>>line1>>a0x>>a0y;
fin>>line1>>gama;
fin>>line1>>rhoref;
fin>>line1>>velref;
fin>>line1>>Machref;
fin>>line1>>alpha;
presref=rhoref*velref*velref/(gama*Machref*Machref);
fin>>line1>>coefa>>coefb>>coefc;

if(phytype==0) nvar=1;
if(phytype==1&&spacedim==1) nvar=3;
if(phytype==1&&spacedim==2) nvar=4;

while(1>0)
{
getline(fin,line);
if(line=="start disc") break;	
}
fin>>line1>>disc_type;
fin>>line1>>disc_order;
ndegr=disc_order+1;
fin>>line1>>recons;
fin>>line1>>reconsvariable;
fin>>line1>>limiter;
fin>>line1>>limitertype;
fin>>line1>>basis_func;
fin>>line1>>flux_type;
fin>>line1>>hrd_lmtr_cin;
fin>>line1>>diff_type;
fin>>line1>>etabr2;

while(1>0)
{
getline(fin,line);
if(line=="start bc") break;	
}
fin>>line1>>bcl;
fin>>line1>>bcr;

while(1>0)
{
getline(fin,line);
if(line=="start study") break;	
}
fin>>line1>>studytype;
fin>>line1>>steadystate_type;
fin>>line1>>timedep_type;
fin>>line1>>t0;
fin>>line1>>tf;
fin>>line1>>dt;
fin>>line1>>CFL_use;
fin>>line1>>CFL;
fin>>line1>>timedep_method;
fin>>line1>>Niter;
fin>>line1>>tol;
fin>>line1>>timstp_option;
fin.close();

//Writing the input file  information into another output file
ofstream fout;
fout.open("inputfile_output.yaml", ios::out | ios::trunc);

fout<<"#--------Space Dimension-------- "<<endl;
switch(spacedim)
{
	case 0: fout<<"ODE Equation"<<endl; break;
	case 1:
		fout<<"One dimensional"<<endl;
		fout<<"#--------Geometry and mesh information(only for 1D)--------"<<endl;
		fout<<"Left end coordinate :"<<leftcoord<<endl;
		fout<<"Right end coordinate :"<<rightcoord<<endl;
		fout<<"Number of elements: "<<nelem<<endl;
		fout<<"#--------Boundary Conditions (Implemented only for 1D for now)--------"<<endl;
		fout<<"Left boundary :";
		switch(bcl)
		{
			case 0: fout<<"Dirchilet"<<endl; break;
			case 1: fout<<"Neumann"<<endl; break;
			case 2: fout<<"Periodic"<<endl; break;
		}
		fout<<"Right boundary :";
		switch(bcr)
		{
			case 0: fout<<"Dirchilet"<<endl; break;
			case 1: fout<<"Neumann"<<endl; break;
			case 2: fout<<"Periodic"<<endl; break;
		}
	break;
	case 2: fout<<"Two dimensinal"<<endl; break;
	case 3: fout<<"Three dimensinal"<<endl; break;
}

fout<<"#--------Physics type and Discretization--------"<<endl;
fout<<"Physics type :";
switch(phytype)
{
	case 0: fout<<"Linear advection Equation"<<endl; break;
	case 1: fout<<"Euler Equations"<<endl; break;
}

fout<<"Discretization type :";
switch(disc_type)
{
	case 0: fout<<"Finite volume (DGP0)"<<endl; break;
	case 1: fout<<"Continous Galerkin Finite Element"<<endl; break;
	case 2: fout<<"Discontinous Galerkin Finite Element(DG)"<<endl; break;

}

fout<<"Discretization order :";
switch(disc_order)
{
	case 0: fout<<"Constant(Finite volume, P0)"<<endl; break;
	case 1: fout<<"Linear, P1"<<endl; break;
	case 2: fout<<"Quadratic, P2 "<<endl; break;
}

fout<<"Reconstruction :";
switch(recons)
{
	case 0: fout<<"Not used"<<endl; break;
	case 1: fout<<"Used"<<endl; break;
}

fout<<"Basis function type :";
switch(basis_func)
{
	case 0: fout<<"Lagrangian"<<endl; break;
	case 1: fout<<"Taylor Basis(For DG Methods)"<<endl; break;
    case 2: fout<<"Legendre Basis(For DG Methods)"<<endl; break;
}
fout<<"#--------Study type and Solvers--------"<<endl;
fout<<"Study type :";
switch(studytype)
{
	case 0:
	    fout<<"Steady State solution"<<endl;
		fout<<"Analysis type :";
		switch(steadystate_type)
		{
			case 0: fout<<"Static steady state"<<endl; break;
			case 1: fout<<"Frequency Domain"<<endl; break;
			case 2: fout<<"Eigen Frequency"<<endl; break;
		}
	break;
	case 1:
	    fout<<"Time Dependent Solution"<<endl;
		fout<<"Type :";
		switch(timedep_type)
		{
			case 0: fout<<"Explicit"<<endl; break;
			case 1: fout<<"Implicit"<<endl; break;
		}
		fout<<"Method :";
		switch(timedep_method)
		{
			case 0: fout<<"Forward Euler"<<endl; break;
			case 1: fout<<"Backward Euler"<<endl; break;
			case 2: fout<<"Two stage TVD Range Kutta"<<endl; break;
			case 3: fout<<"Three stage TVD Range Kutta"<<endl; break;
		}
		fout<<"Start Time :"<<t0<<endl;
		fout<<"Final Time :"<<tf<<endl;
		fout<<"Time Step :"<<dt<<endl;
		switch(CFL_use)
		{
			case 0: fout<<"Not Using CFL"<<endl; break;
			case 1: fout<<"Using CFL, "<<"CFL Number :"<<CFL<<endl; break;
		}
	break;
}
	
fout.close();
}

