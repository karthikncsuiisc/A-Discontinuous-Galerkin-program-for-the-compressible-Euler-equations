CXX :=g++
CXXFLAGS := -g -ansi -Wall -O3

FILESMAIN := main.o phymanager.o read_input.o minmod.o MClimiter.o

ADD1D := ./../oned
FILES1D := $(ADD1D)/onedproblem_euler_dg.o $(ADD1D)/stegar_warming_flux_1d.o $(ADD1D)/initial_condition_1d.o $(ADD1D)/periodicboundary_1d.o $(ADD1D)/noflux_1d.o $(ADD1D)/godunov_method_1d.o $(ADD1D)/van_leer_FVS_1d.o $(ADD1D)/dgp0_fun_1d.o $(ADD1D)/dgp1_fun_1d.o $(ADD1D)/dgp2_fun_1d.o $(ADD1D)/onedproblem_manager.o $(ADD1D)/onedproblem_wave_dg.o $(ADD1D)/rieman_solver_oned_wave.o $(ADD1D)/noflux_1d_wave.o $(ADD1D)/periodicboundary_1d_wave.o $(ADD1D)/initial_condition_1d_wave.o $(ADD1D)/onedresults.o $(ADD1D)/onedproblem_burger_dg.o $(ADD1D)/rhs_flux_contrb_burger.o $(ADD1D)/rhs_flux_contrb_diff_ddg.o $(ADD1D)/rhs_contrb_burger_p1.o $(ADD1D)/rhs_contrb_ddg_p1.o $(ADD1D)/rhs_contrb_burger_p2.o $(ADD1D)/rhs_contrb_ddg_p2.o $(ADD1D)/initial_condition_1d_burger.o $(ADD1D)/rhs_contrb_diff_br2_p1.o

ADD2D :=./../twod
FILES2D := $(ADD2D)/twod_manager.o $(ADD2D)/reading_grid2D.o $(ADD2D)/elements_su_points.o $(ADD2D)/points_su_points.o $(ADD2D)/elements_su_elements.o $(ADD2D)/findfaces.o $(ADD2D)/tvd_rk3_method.o $(ADD2D)/triangle_area.o $(ADD2D)/flux2D_adv_calc.o $(ADD2D)/flux2D_eul_calc.o $(ADD2D)/boundary_condition_euler_2d.o $(ADD2D)/boundary_condition_adv_2d.o $(ADD2D)/initial_conditions_2d.o $(ADD2D)/gauss_point_values.o $(ADD2D)/results2D.o $(ADD2D)/recons_calculation_euler.o $(ADD2D)/recons_coef_calculation.o $(ADD2D)/flux2D_eul_calc_recons.o $(ADD2D)/flux2D_eul_calc_recons_limiter.o $(ADD2D)/tvd_rk3_advection.o $(ADD2D)/tvd_RK_euler_p1.o $(ADD2D)/flux2D_eul_calc_p1.o $(ADD2D)/error_function.o
 
TOTADD := -I $(ADD1D) -I $(ADD2D)

OBJECTS := $(FILESMAIN) $(FILES1D) $(FILES2D)

all: project

project: $(OBJECTS) 
	 $(CXX) $(OBJECTS) -o project -O3

$(FILESMAIN): %.o: %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@ $(TOTADD) 

$(FILES1D): %.o: %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@ -I.

$(FILES2D): %.o: %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@ -I.

clean:
	rm -rf $(OBJECTS) project

cleano:
	rm -rf $(OBJECTS)
