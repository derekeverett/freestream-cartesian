#define PI 3.141592654f //question bro?
#define DIM_X 41 //number of grid point in x direction
#define DIM_Y 41 //number of grid point in y direction
#define DIM_Z 41 //number of grid point in z direction
#define DIM (DIM_X * DIM_Y * DIM_Z) //total number of spatial gridpoints
#define DIM_THETAP 31 //number of grid points in theta_p momentum polar angle
#define DIM_PHIP 31 //number of grid points in phi_p momentum azimuthal angle
#define DX 0.1f //spacing of grid in x direction
#define DY 0.1f //spacing of grid in y direction
#define DZ 0.1f //spacing of grid in eta direction
#define DT 0.5f //free streaming time step size
#define EOS_TYPE 2 // 1 for conformal EOS, 2 for Wuppertal-Budhapest parameterization
