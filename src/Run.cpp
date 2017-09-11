// CPU-freestream
// simulates 3 dimensional free streaming of massless particles given an
// initial energy density profile, assuming initial momentum space isotropy in px,py,pz.
// Then performs Landau Matching to find the components of Energy Momentum Tensor.

#include "Parameters.h"
#include "FreeStream.cpp"
#include "InitialConditions.cpp"
#include "LandauMatch.cpp"
#include "Memory.cpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//using namespace std;

int main(void)
{
  printf("Welcome to freestream-cartesian\n");
  printf("Parameters are ...\n");
  printf("(DIM_X, DIM_Y, DIM_Z) = (%d, %d, %d)\n", DIM_X, DIM_Y, DIM_Z);
  printf("(DX, DY, DZ, DT) = (%.2f, %.2f, %.2f, %.2f)\n", DX, DY, DZ, DT);
  //set up the coordinates with center at (x = 0, y = 0, z = 0)
  //this works for an even number of lattice points; include conditional for odd number
  float xmin, xmax, ymin, ymax, zmin, zmax;
  xmax = (DIM_X % 2 == 0) ? (float(DIM_X) / 2.0 * DX) : (float(DIM_X - 1) / 2.0 * DX);
  xmin = (-1.0) * xmax;
  ymax = (DIM_Y % 2 == 0) ? (float(DIM_Y) / 2.0 * DY) : (float(DIM_Y - 1) / 2.0 * DY);
  ymin = (-1.0) * ymax;
  zmax = (DIM_Z % 2 == 0) ? (float(DIM_Z) / 2.0 * DZ) : (float(DIM_Z - 1) / 2.0 * DZ);
  zmin = (-1.0) * zmax;

  //allocate and initialize memory
  //write functions in Memory.cpp to do this
  printf("Allocating memory\n");
  //the initial energy density spatial profile
  float ***initialEnergyDensity;
  initialEnergyDensity = calloc3dArray(initialEnergyDensity, DIM_X, DIM_Y, DIM_Z);

  /*
  initialEnergyDensity = (float ***)calloc(DIM_X, sizeof(float **));
  for (int ix = 0; ix < DIM_X; ix++)
  {
    initialEnergyDensity[ix] = (float **)calloc(DIM_Y, sizeof(float *));
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      initialEnergyDensity[ix][iy] = (float *)calloc(DIM_Z, sizeof(float));
    }
  }
  */
  //the initial density G^(0,0) at time t0
  float ***density;
  density = calloc3dArray(density, DIM_X, DIM_Y, DIM_Z);
  /*
  density = (float ***)calloc(DIM_X, sizeof(float **));
  for (int ix = 0; ix < DIM_X; ix++)
  {
    density[ix] = (float **)calloc(DIM_Y, sizeof(float *));
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      density[ix][iy] = (float *)calloc(DIM_Z, sizeof(float));
    }
  }
  */
  //the shited density profile G^(0,0) at time t
  float *****shiftedDensity;
  shiftedDensity = calloc5dArray(shiftedDensity, DIM_X, DIM_Y, DIM_Z, DIM_THETAP, DIM_PHIP);
  /*
  shiftedDensity = (float *****)calloc(DIM_X, sizeof(float ****));
  for (int ix = 0; ix < DIM_X; ix++)
  {
    shiftedDensity[ix] = (float ****)calloc(DIM_Y, sizeof(float ***));
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      shiftedDensity[ix][iy] = (float ***)calloc(DIM_Z, sizeof(float **));
      for (int iz = 0; iz < DIM_Z; iz++)
      {
        shiftedDensity[ix][iy][iz] = (float **)calloc(DIM_THETAP, sizeof(float *));
        for(int ithetap = 0; ithetap < DIM_THETAP; ithetap++)
        {
          shiftedDensity[ix][iy][iz][ithetap] = (float *)calloc(DIM_PHIP, sizeof(float));
        }
      }
    }
  }
  */
  //the ten independent components of the stress tensor
  float ****stressTensor;
  stressTensor = calloc4dArray(stressTensor, 10, DIM_X, DIM_Y, DIM_Z);
  /*
  stressTensor = (float ****)calloc(10, sizeof(float ***)); //10 independent components in 3+1d
  for (int ivar = 0; ivar < 10; ivar++)
  {
    stressTensor[ivar] = (float ***)calloc(DIM_X, sizeof(float **));
    for (int ix = 0; ix < DIM_X; ix++)
    {
      stressTensor[ivar][ix] = (float **)calloc(DIM_Y, sizeof(float *));
      for (int iy = 0; iy < DIM_Y; iy++)
      {
        stressTensor[ivar][ix][iy] = (float *)calloc(DIM_Z, sizeof(float));
      }
    }
  }
  */
  //a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu) normalized by energy
  float ***trigTable;
  trigTable = calloc3dArray(trigTable, 10, DIM_THETAP, DIM_PHIP);
  /*
  trigTable = (float ***)calloc(10, sizeof(float **));
  for (int ivar = 0; ivar < 10; ivar++)
  {
    trigTable[ivar] = (float **)calloc(DIM_THETAP, sizeof(float *));
    for (int ithetap = 0; ithetap < DIM_THETAP; ithetap++)
    {
      trigTable[ivar][ithetap] = (float *)calloc(DIM_PHIP, sizeof(float));
    }
  }
  */
  //initialize energy density; here we use gaussian for testing
  printf("setting initial conditions on energy density\n");
  initializeGauss(initialEnergyDensity, xmax, ymax, zmax, 1.0);

  //read in the initial energy density profile (from file)
  //readInitialEnergyDensity(initialEnergyDensity);

  //convert the energy density profile into the initial density profile to be streamed - just a normalization
  convertInitialDensity(initialEnergyDensity, density);

  //perform the free streaming time-update step
  printf("performing the free streaming time step\n");
  clock_t t;
  t = clock();
  freeStream(density, shiftedDensity, xmin, ymin, zmin);
  t = clock() - t;
  printf("Free streaming took %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

  //Landau matching to find the components of energy-momentum tensor
  printf("Landau matching to find hydrodynamic variables\n");

  printf("calculating trig table\n");
  t = clock();
  calculateTrigTable(trigTable);
  t = clock() - t;
  printf("calculating trig table took %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

  //calculate the ten independent components of the stress tensor by integrating over momentum angles
  printf("calculating independent components of stress tensor\n");
  t = clock();
  calculateStressTensor(stressTensor, shiftedDensity, trigTable);
  t = clock() - t;
  printf("calculating stress tensor took %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

  //solve the eigenvalue problem for the energy density and flow velocity
  printf("solving eigenvalue problem for energy density and flow velocity\n");
  t = clock();
  solveEigenSystem(stressTensor);
  t = clock() - t;
  printf("solving eigenvalue problem took %f seconds\n", ((float)t)/CLOCKS_PER_SEC);
  //printf("writing hydro variables to file\n");

  //free the memory
  free3dArray(initialEnergyDensity, DIM_X, DIM_Y);
  free3dArray(density, DIM_X, DIM_Y);
  free5dArray(shiftedDensity, DIM_X, DIM_Y, DIM_Z, DIM_THETAP);
  free4dArray(stressTensor, 10, DIM_X, DIM_Y);
  free3dArray(trigTable, 10, DIM_THETAP);
  printf("Done... Goodbye!\n");
}
