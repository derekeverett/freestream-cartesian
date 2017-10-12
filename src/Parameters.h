//IC_ENERGY sets the initial condition on the energy density profile with options
// 1 : Smooth Oblate Gaussian
// 2 : Fluctuating Oblate Gaussian
// 3 : read in from data file in whitespace delimited format : x y z epsilon(x,y,z)
//IC_BARYON sets the initial baryon density profile with options
// 1 : Smooth Oblate Gaussian
// 2 : Fluctuating Oblate Gaussian
// 3 : read in from data file in whitespace delimited format : x y z n_B(x,y,z)
#define BARYON 1 //if true, code will freestream baryon current as well as stress tensor
#define IC_ENERGY 1
#define IC_BARYON 1
#define PI 3.141592654f //question bro?
#define DIM_X 101 //number of grid points in x direction
#define DIM_Y 101 //number of grid points in y direction
#define DIM_Z 51 //number of grid points in z direction
#define DIM (DIM_X * DIM_Y * DIM_Z) //total number of spatial grid points
#define DIM_THETAP 31 //number of grid points in theta_p momentum polar angle
#define DIM_PHIP 31 //number of grid points in phi_p momentum azimuthal angle
#define DX 0.1f //spacing of grid in x direction
#define DY 0.1f //spacing of grid in y direction
#define DZ 0.1f //spacing of grid in eta direction
#define DT 0.5f //free streaming time step size
#define EOS_TYPE 1 // 1 for conformal EOS, 2 for Wuppertal-Budhapest parameterization
