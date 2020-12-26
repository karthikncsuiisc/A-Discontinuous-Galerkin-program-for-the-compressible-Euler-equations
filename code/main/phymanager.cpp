//Function to decide which function to choose.

#include "functions.h"
#include "oned_header.h"
#include "twod_header.h"

void phymanager()
{
switch(spacedim)
{
	case 0:
		cout<<"ODE problems are not implemented"<<endl; exit(1);
	case 1:
		onedproblem_manager(); break;
	case 2:
		twod_manager(); break;
	case 3:
		cout<<"3D problems are not yet implemented in this package"<<endl; exit(3);
}

}
