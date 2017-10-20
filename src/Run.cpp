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

  //the initial energy density spatial profile in cartesian coords
  float *initialEnergyDensity;
  initialEnergyDensity = (float *)calloc(DIM, sizeof(float));

  //the initial baryon density spatial profile in cartesian coords
  float *initialChargeDensity;
  initialChargeDensity = (float *)calloc(DIM, sizeof(float));

  //the initial density G^(0,0) at time t0 in cartesian coords
  float *density;
  density = (float *)calloc(DIM, sizeof(float));

  //the initial density J^0 at time t0 in cartesian coords
  float *chargeDensity;
  chargeDensity = (float *)calloc(DIM, sizeof(float));

  //the shifted density profile G^(0,0) at time t in cartesian coords
  float ***shiftedDensity;
  shiftedDensity = calloc3dArray(shiftedDensity, DIM, DIM_THETAP, DIM_PHIP);

  //the shifted baryon density profile J^0 at time t in cartesian coords
  float ***shiftedChargeDensity;
  shiftedChargeDensity = calloc3dArray(shiftedChargeDensity, DIM, DIM_THETAP, DIM_PHIP);

  //the stress tensor as a function of cartesian time and z
  float ***timeDependentStressTensor;
  timeDependentStressTensor = calloc3dArray(timeDependentStressTensor, 10, DIM_T, DIM);
  //the baryon current as a function of cartesian time and z
  float ***timeDependentBaryonCurrent;
  timeDependentBaryonCurrent = calloc3dArray(timeDependentBaryonCurrent, 4, DIM_T, DIM);

  //the ten components of the final stress tensor in milne coords
  //must be constructed by interpolating T(t,z) on a regular grid in t,z to a
  //regular grid in eta, evaluated at a fixed proper time.
  float **stressTensor;
  stressTensor = calloc2dArray(stressTensor, 10, DIM_MILNE);

  //four components of the final baryon number current four-vector in milne coords
  float **baryonCurrent;
  baryonCurrent = calloc2dArray(baryonCurrent, 4, DIM_MILNE);
  //a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu) normalized by energy
  float ***trigTable;
  trigTable = calloc3dArray(trigTable, 11, DIM_THETAP, DIM_PHIP);

  //variables to store the hydrodynamic variables after the Landau matching is performed

  //the energy density in milne coords
  float *energyDensity;
  energyDensity = (float *)calloc(DIM_MILNE, sizeof(float));

  //the baryon density in milne coords
  float *baryonDensity;
  baryonDensity = (float *)calloc(DIM_MILNE, sizeof(float));

  //the flow velocity in milne coords
  float **flowVelocity;
  flowVelocity = calloc2dArray(flowVelocity, 4, DIM_MILNE);

  //the pressure in milne coords
  float *pressure;
  pressure = (float *)calloc(DIM_MILNE, sizeof(float));

  //the bulk pressure in milne coords
  float *bulkPressure;
  bulkPressure = (float *)calloc(DIM_MILNE, sizeof(float));

  //the shear viscous tensor in milne coords
  float **shearTensor;
  shearTensor = calloc2dArray(shearTensor, 10, DIM_MILNE); //calculate 10 components; can check tracelessness
  // and orthogonality later as a consistency check

  //the baryon diffusion current vector in milne coords
  float **baryonDiffusion;
  baryonDiffusion = calloc2dArray(baryonDiffusion, 4, DIM_MILNE);

  //initialize energy density in cartesian coordinates
  printf("setting initial conditions on energy density : ");
  if (IC_ENERGY == 1)
  {
    initializeEllipticalGauss(initialEnergyDensity, 1.0, 1.0, 1.0);
    printf("Smooth Oblate Gaussian \n");
  }
  else if (IC_ENERGY == 2)
  {
    initializeEllipticalMCGauss(initialEnergyDensity, 1.0, 1.0, 1.0);
    printf("Fluctuating Oblate Gaussian \n");
  }
  else if (IC_ENERGY == 3)
  {
    readDensityFile(initialEnergyDensity, "initial_profiles/e.dat");
    printf("Reading from energy density file in initial_profiles/ \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return 0;
  }

  if (BARYON)
  {
    //initialize baryon density in cartesian coordinates
    printf("setting initial conditions on baryon density : ");
    if (IC_BARYON == 1)
    {
      initializeEllipticalGauss(initialChargeDensity, 1.0, 1.0, 1.0);
      printf("Smooth Oblate Gaussian \n");
    }
    else if (IC_BARYON == 2)
    {
      initializeEllipticalMCGauss(initialChargeDensity, 1.0, 1.0, 1.0);
      printf("Fluctuating Oblate Gaussian \n");
    }
    else if (IC_BARYON == 3)
    {
      readDensityFile(initialChargeDensity, "initial_profiles/nB.dat");
      printf("Reading from baryon density file in initial_profiles/ \n");
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

  //calculate trig table to speed up calculation of stress tensor and baryon current
  calculateTrigTable(trigTable);

  double sec;
  printf("performing the free streaming and calculating stress tensor (and baryon current) \n");
  sec = omp_get_wtime();

  //main time step loop for free-streaming
  for (int itime = 0; itime < DIM_T; itime++)
  {
    //perform the free streaming time-update steps
    float dt = (itime + 1) * DT;  //the streaming time
    freeStream(density, shiftedDensity, dt);
    if (BARYON) freeStream(chargeDensity, shiftedChargeDensity, dt);
    //calculate the time dependent stress tensor and baryon current
    calculateStressTensor(timeDependentStressTensor, shiftedDensity, trigTable, itime);
    if (BARYON) calculateBaryonCurrent(timeDependentBaryonCurrent, shiftedChargeDensity, trigTable, itime);
    printf("time step %d finished\n", itime);
  }
  sec = omp_get_wtime() - sec;
  printf("Free streaming took and calculating stress tensor (and baryon current) took %f seconds\n", sec);

  //interpolate T^(mu,nu) and j^(mu) along a regular grid in spacetime rapidity for a fixed longitudinal proper time
  //and perform appropriate transformations to obtain both quantities in milne coords (e.g. j^(mu) = (J^(tau), j^(x),j^(y),j^(eta)))

  //solve the eigenvalue problem for the energy density and flow velocity
  printf("solving eigenvalue problem for energy density and flow velocity\n");
  sec = omp_get_wtime();
  solveEigenSystem(stressTensor, energyDensity, flowVelocity);
  sec = omp_get_wtime() - sec;
  printf("solving eigenvalue problem took %f seconds\n", sec);

  if (BARYON)
  {
    printf("calculating independent components of baryon current\n");
    sec = omp_get_wtime();

    sec = omp_get_wtime() - sec;
    printf("calculating baryon current took %f seconds\n", sec);

    calculateBaryonDensity(baryonDensity, baryonCurrent, flowVelocity);
    calculateBaryonDiffusion(baryonDiffusion, baryonCurrent, baryonDensity, flowVelocity);
  }

  //need EoS which includes baryon number!!! fix this
  calculatePressure(energyDensity, pressure);
  calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure);
  calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor);

  printf("writing hydro variables to file\n");
  writeScalarToFile(energyDensity, "e");
  writeScalarToFile(pressure, "p");
  writeScalarToFile(bulkPressure, "bulk_PI");
  writeScalarToFileProjection(energyDensity, "e_projection");
  writeVectorToFile(flowVelocity, "u_x", 1);
  writeVectorToFile(flowVelocity, "u_y", 2);
  writeVectorToFile(flowVelocity, "u_z", 3);
  writeVectorToFileProjection(flowVelocity, "u_x_projection", 1);
  writeVectorToFileProjection(flowVelocity, "u_y_projection", 2);
  writeVectorToFileProjection(flowVelocity, "u_z_projection", 3);
  writeVectorToFile1DProjection(flowVelocity, "u_x_1Dprojection", 1);

  if (BARYON)
  {
    writeScalarToFile(baryonDensity, "nB");
    writeScalarToFileProjection(baryonDensity, "nB_projection");
    writeVectorToFile(baryonDiffusion, "V_x", 1);
    writeVectorToFile(baryonDiffusion, "V_y", 2);
    writeVectorToFile(baryonDiffusion, "V_z", 3);
  }
  //free the memory
  printf("freeing memory \n");
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
