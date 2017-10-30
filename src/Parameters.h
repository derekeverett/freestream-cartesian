//IC_ENERGY sets the initial condition on the energy density profile with options
// 1 : Smooth Oblate Gaussian
// 2 : Fluctuating Oblate Gaussian
// 3 : read in from data file in whitespace delimited format : x y z epsilon(x,y,z)
//IC_BARYON sets the initial baryon density profile with options
// 1 : Smooth Oblate Gaussian
// 2 : Fluctuating Oblate Gaussian
// 3 : read in from data file in whitespace delimited format : x y z n_B(x,y,z)
#define BARYON 0 //if true, code will freestream baryon current as well as stress tensor
#define IC_ENERGY 1
#define IC_BARYON 1
#define PI 3.141592654f //question bro?
#define DIM_X 51 //number of grid points in x direction
#define DIM_Y 51 //number of grid points in y direction
#define DIM_Z 51 //number of grid points in z direction
#define DIM_ETA 51 //number of grid points in eta direction (spacetime rapidity)
#define DIM_T 5 //total number of free-streaming cartesian time steps
#define DIM (DIM_X * DIM_Y * DIM_Z) //total number of spatial grid points in cartesian coords
#define DIM_MILNE (DIM_X * DIM_Y * DIM_ETA) //total number of spatial grid points in milne coords
#define DIM_THETAP 31 //number of grid points in theta_p momentum polar angle
#define DIM_PHIP 31 //number of grid points in phi_p momentum azimuthal angle
#define DX 0.1f //spacing of grid in x direction
#define DY 0.1f //spacing of grid in y direction
#define DZ 0.1f //spacing of grid in z direction
#define DT 0.05f //free streaming time step size
#define T0 0.0 //initial cartesian time at which densities are initialized 
#define TAU 0.5f //final proper time at which Landau Matching is performed in Milne Coordinates
#define EOS_TYPE 1 // 1 for conformal EOS, 2 for Wuppertal-Budhapest parameterization
