#ifndef oned_header_h
#define onedheader_h

extern double *coord_1d;
extern int npoint_1d;

void onedproblem_manager();
void onedproblem_euler_dg();               //One Dimensional problem function
void onedproblem_wave_dg();

//Euler equation functions
void initial_condition_1d();
void periodicboundary_1d(); 
void noflux_1d();
void stegar_warming_flux_1d();
void godunov_method_1d(double);
void van_leer_FVS_1d();
void dgp0_fun_1d();
void dgp1_fun_1d();
void dgp2_fun_1d();

//Advection-equations
void noflux_1d_wave();
void periodicboundary_1d_wave();
void initial_condition_1d_wave();
void rieman_solver_oned_wave();
void onedresults();

//Advection burger eqation
void onedproblem_burger_dg();
void initial_condition_1d_burger();
void rhs_flux_contrb_burger();
void rhs_flux_contrb_diff_ddg();
void rhs_contrb_burger_p1();
void rhs_contrb_ddg_p1();
void rhs_contrb_burger_p2();
void rhs_contrb_ddg_p2();
void rhs_contrb_diff_br2_p1();


#endif







