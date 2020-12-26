//Function to decide which function to choose.

#include "functions.h"
#include "oned_header.h"

void onedproblem_manager()
{
switch(phytype) 
{
	case 0:
		onedproblem_wave_dg(); break;
	case 1:
		onedproblem_euler_dg(); break;
	case 2:
		onedproblem_burger_dg(); break;
	default: cout<<"wrong choice of the physics seletion"<<endl;
}

}
