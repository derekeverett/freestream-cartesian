// CPU-freestream
// simulates 3 dimensional free streaming of massless particles given an
// initial energy density profile, assuming initial momentum space isotropy in px,py,pz.
// Then performs Landau Matching to find the components of Energy Momentum Tensor and baryon current
// and performs standard viscous hydro decomposition.

#include "Parameters.h"
#include "FreeStream.cpp"
#include "InitialConditions.cpp"
#include "LandauMatch.cpp"
#include "EquationOfState.cpp"
#include "Memory.cpp"
#include "FileIO.cpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
//using namespace std;

int main(void)
{
  printf("Welcome to freestream-cartesian\n");
  printf("Parameters are ...\n");
  printf("(DIM_X, DIM_Y, DIM_Z) = (%d, %d, %d)\n", DIM_X, DIM_Y, DIM_Z);
  printf("(DX, DY, DZ, DT) = (%.2f, %.2f, %.2f, %.2f)\n", DX, DY, DZ, DT);

  //allocate and initialize memory
  printf("Allocating memory\n");
  //the initial energy density spatial profile
  float *initialEnergyDensity;
  initialEnergyDensity = (float *)calloc(DIM, sizeof(float));
  //the initial baryon density spatial profile
  float *initialChargeDensity;
  initialChargeDensity = (float *)calloc(DIM, sizeof(float));
  //the initial density G^(0,0) at time t0
  float *density;
  density = (float *)calloc(DIM, sizeof(float));
  //the initial density J^0 at time t0
  float *chargeDensity;
  chargeDensity = (float *)calloc(DIM, sizeof(float));
  //the shifted density profile G^(0,0) at time t
  float ***shiftedDensity;
  shiftedDensity = calloc3dArray(shiftedDensity, DIM, DIM_THETAP, DIM_PHIP);
  //the shifted baryon density profile J^0 at time t
  float ***shiftedChargeDensity;
  shiftedChargeDensity = calloc3dArray(shiftedChargeDensity, DIM, DIM_THETAP, DIM_PHIP);
  //the ten independent components of the stress tensor
  float **stressTensor;
  stressTensor = calloc2dArray(stressTensor, 10, DIM);
  //four components of baryon number current four-vector
  float **baryonCurrent;
  baryonCurrent = calloc2dArray(baryonCurrent, 4, DIM);
  //a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu) normalized by energy
  float ***trigTable;
  trigTable = calloc3dArray(trigTable, 11, DIM_THETAP, DIM_PHIP);

  //variables to store the hydrodynamic variables after the Landau matching is performed
  //the energy density
  float *energyDensity;
  energyDensity = (float *)calloc(DIM, sizeof(float));
  //the baryon density
  float *baryonDensity;
  baryonDensity = (float *)calloc(DIM, sizeof(float));
  //the flow velocity
  float **flowVelocity;
  flowVelocity = calloc2dArray(flowVelocity, 4, DIM);
  //the pressure
  float *pressure;
  pressure = (float *)calloc(DIM, sizeof(float));
  float *bulkPressure;
  bulkPressure = (float *)calloc(DIM, sizeof(float));
  float **shearTensor;
  shearTensor = calloc2dArray(shearTensor, 10, DIM); //calculate 10 components; can check tracelessness
  // and orthogonality later as a consistency check
  //the baryon diffusion current vector
  float **baryonDiffusion;
  baryonDiffusion = calloc2dArray(baryonDiffusion, 4, DIM);
  //initialize energy density
  printf("setting initial conditions on energy density : ");
  if (IC_ENERGY == 1)
  {
    initializeEllipticalGauss(initialEnergyDensity, 0.5, 1.0, 0.5);
    printf("Smooth Oblate Gaussian \n");
  }
  else if (IC_ENERGY == 2)
  {
    initializeEllipticalMCGauss(initialEnergyDensity, 0.5, 1.0, 0.5);
    printf("Fluctuating Oblate Gaussian \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return 0;
  }

  if (BARYON)
  {
    //initialize baryon density
    printf("setting initial conditions on baryon density : ");
    if (IC_BARYON == 1)
    {
      initializeEllipticalGauss(initialChargeDensity, 0.5, 1.0, 0.5);
      printf("Smooth Oblate Gaussian \n");
    }
    else if (IC_BARYON == 2)
    {
      initializeEllipticalMCGauss(initialChargeDensity, 0.5, 1.0, 0.5);
      printf("Fluctuating Oblate Gaussian \n");
    }
    else
    {
      printf("Not a valid initial Condition... Goodbye\n");
      return 0;
    }
  }
  //write initial profile to file
  writeScalarToFile(initialEnergyDensity, "initial_e");
  if (BARYON) writeScalarToFile(initialChargeDensity, "initial_nB");
  writeScalarToFileProjection(initialEnergyDensity, "initial_e_projection");
  if (BARYON) writeScalarToFileProjection(initialChargeDensity, "initial_nB_projection");

  //convert the energy density profile into the initial density profile to be streamed - just a normalization
  convertInitialDensity(initialEnergyDensity, density);
  //convert baryon density profile into initial profile J^0 to be streamed - just a normalization
  if (BARYON) convertInitialDensity(initialChargeDensity, chargeDensity);
  //perform the free streaming time-update step
  //pretabulate trig functions before this step to save time?
  printf("performing the free streaming time step\n");
  double sec;
  sec = omp_get_wtime();
  freeStream(density, shiftedDensity);
  if (BARYON) freeStream(chargeDensity, shiftedChargeDensity);
  sec = omp_get_wtime() - sec;
  printf("Free streaming took %f seconds\n", sec);

  //Landau matching to find the components of energy-momentum tensor and baryon current
  printf("Landau matching to find hydrodynamic variables\n");

  calculateTrigTable(trigTable);

  //calculate the ten independent components of the stress tensor by integrating over momentum angles
  printf("calculating independent components of stress tensor\n");
  sec = omp_get_wtime();
  calculateStressTensor(stressTensor, shiftedDensity, trigTable);
  sec = omp_get_wtime() - sec;
  printf("calculating stress tensor took %f seconds\n", sec);


  //solve the eigenvalue problem for the energy density and flow velocity
  printf("solving eigenvalue problem for energy density and flow velocity\n");
  sec = omp_get_wtime();
  solveEigenSystem(stressTensor, energyDensity, flowVelocity);
  sec = omp_get_wtime() - sec;
  printf("solving eigenvalue problem took %f seconds\n", sec);

  calculatePressure(energyDensity, pressure);
  calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure);
  calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor);

  if (BARYON)
  {
    printf("calculating independent components of baryon current\n");
    sec = omp_get_wtime();
    calculateBaryonCurrent(baryonCurrent, shiftedChargeDensity, trigTable);
    sec = omp_get_wtime() - sec;
    printf("calculating baryon current took %f seconds\n", sec);

    calculateBaryonDensity(baryonDensity, baryonCurrent, flowVelocity);
    calculateBaryonDiffusion(baryonDiffusion, baryonCurrent, baryonDensity, flowVelocity);
  }

  printf("writing hydro variables to file\n");
  writeScalarToFile(energyDensity, "e");
  writeScalarToFile(pressure, "p");
  writeScalarToFile(bulkPressure, "bulk_PI");
  writeScalarToFileProjection(energyDensity, "e_projection");
  writeVectorToFile(flowVelocity, "u_x", 1);
  writeVectorToFile(flowVelocity, "u_y", 2);
  writeVectorToFile(flowVelocity, "u_z", 3);

  if (BARYON)
  {
    writeScalarToFile(baryonDensity, "nB");
    writeScalarToFileProjection(baryonDensity, "nB_projection");
  }
  //free the memory
  free(initialEnergyDensity);
  free(density);
  free3dArray(shiftedDensity, DIM, DIM_THETAP);
  free2dArray(stressTensor, 10);
  free3dArray(trigTable, 11, DIM_THETAP);
  free(energyDensity);
  free2dArray(flowVelocity, 4);
  free(pressure);
  free(bulkPressure);
  free2dArray(shearTensor, 10);

  if (BARYON)
  {
    free(initialChargeDensity);
    free(chargeDensity);
    free3dArray(shiftedChargeDensity, DIM, DIM_THETAP);
    free2dArray(baryonCurrent, 4);
    free(baryonDensity);
    free2dArray(baryonDiffusion, 4);
  }

  printf("Done... Goodbye!\n");
}
