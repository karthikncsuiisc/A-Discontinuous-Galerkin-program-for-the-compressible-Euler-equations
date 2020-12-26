#ifndef twod_header_h
#define twodheader_h

extern int npoint; 				// Number of nodes in the mesh
extern int nface; 				// Number of boundary faces
extern int **inpoel; 			// Node numbers of each element
extern double **coord; 			// Coordinates of node points
extern double **bface; 			// Boundary elements, corresponding points and 
extern double **lhspo; 			// System matrix
extern double **ini_cond;		//Initial Conditions

void twod_manager();
void reading_grid2D();

//-----ELements surrounding points variables----------------
extern int *esup1;
extern int *esup2;
void elements_su_points();

//-----Points surrounding points variables----------------
extern int *psup1;
extern int *psup2;

//-----Elements surrounding elements variables----------------
extern int **elsuel;
void elements_su_elements();
double triangle_area(double,double,double,double,double,double);

//-----Internal face variables------------------------------
extern int **intface;
extern int naface;              // Number of all  faces
extern int nbface;				// Number of boundary faces
extern double **bound_cond;
extern double *nx,*ny,*length,*area;
extern double *deltax,*deltay;
extern double *centerx,*centery;
void findfaces();

//-----Gauss point variables------------------------------
extern double **Bxquad,**Byquad,*Wequad;
extern double **Bxcubic,**Bycubic,*Wecubic;
extern double ***Bxline_quad,***Byline_quad,*Weline_quad;
extern double ***Bxline_cubic,***Byline_cubic,*Weline_cubic;
extern double **Bxline_lin,**Byline_lin,Weline_lin;

//---------Reconstruction variables-------------------------
extern double ***gradu;
extern double **rhsel;


void gauss_point_values();
void flux2D_adv_calc();
void tvd_rk3_advection();
void boundary_condition_adv_2d();

void tvd_rk3_method();
void flux2D_eul_calc();
void flux2D_eul_calc_recons();
void flux2D_eul_calc_recons_limiter();
void recons_coef_calculation();
void recons_calculation_euler();

void tvd_RK_euler_p1();
void flux2D_eul_calc_p1();

void initial_conditions_2d();
void boundary_condition_euler_2d();
void error_function();
void results2D();

#endif







