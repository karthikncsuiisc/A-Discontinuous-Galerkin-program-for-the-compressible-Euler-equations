//Twodmanger which calls he aproporatiate function based on the input

#include "functions.h"
#include "twod_header.h"

void twod_manager()
{

for(int i=0;i<nfiles;i++)
{
cout<<"--------------------------------------------------------------------------------------"<<endl;
cout<<"Solving for the grid: "<<grid_names[i]<<endl;	
grid_name=grid_names[i];

reading_grid2D();
elements_su_points();
elements_su_elements();		
findfaces();

switch(phytype)
{
	case 0: tvd_rk3_advection(); break;
	case 1: switch(disc_order)
			{
				case 0: tvd_rk3_method(); break;
				case 1: tvd_RK_euler_p1(); break;
				case 2: cout<<"P2 is not yet coded"<<endl; exit(1);
				default: cout<<"Wrong choice of discretization order"<<endl; exit(1);
			}
	break;
	default: cout<<"invalide physics option for 2D problem, check twodmanager folder"<<endl; exit(1);
}

results2D();
}

}
