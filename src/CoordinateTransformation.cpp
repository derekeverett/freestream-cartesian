#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "Memory.cpp"

interpolateToMilneGrid(float ***timeDependentQuantity, float **quantityAtFixedTau, int dimVar)
{

  //fill the coordinate arrays with coordinate values
  double zmin = (-1.0) * ((double)(DIM_Z-1) / 2.0) * DZ;
  double tArray[DIM_T];
  double zArray[DIM_Z];
  for (int itime = 0; itime < DIM_T; itime++)
  {
    tArray[itime] = T0 + (double)itime * DT; //check that this matches free streaming
  }
  for (int iz = 0; iz < DIM_Z; iz++)
  {
    zArray[iz] = (double)iz * DZ  + zmin;
  }
  //for all indices of conserved current ...
  for (int ivar = 0; ivar < dimVar; ivar++)
  {
    for (int ix = 0; ix < DIM_X; ix++)
    {
      for (int iy = 0; iy < DIM_Y; iy++) //...must interpolate from t,z to tau,eta at every point in transverse plane
      {
        float timeDependentArray[DIM_T * DIM_Z];
        for (int itime = 0; itime < DIM_T; itime++)
        {
          for (int iz = 0; iz < DIM_Z; iz++)
          {
            int is = (DIM_Y * DIM_Z * ix) + (DIM_Z * iy) + iz;
            //set the values of the conserved current at the fixed coordinate points
            //see https://www.gnu.org/software/gsl/doc/html/interp.html
            timeDependentArray[itime + DIM_T * iz] = (double)timeDependentQuantity[ivar][itime][is];
          }
        }
        //now set up the interpolation
        const gsl_interp2d_type *T = gsl_interp2d_bicubic;
        gsl_spline2d *spline = gsl_spline2d_alloc(T, sizeof(tArray) / sizeof(double), sizeof(zArray) / sizeof(double));
        gsl_spline2d_init(spline, timeArray, zArray, timeDependentArray, sizeof(tArray) / sizeof(double), sizeof(zArray) / sizeof(double));
        gsl_interp_accel *tacc = gsl_interp_accel_alloc();
        gsl_interp_accel *zacc = gsl_interp_accel_alloc();

        //try interpolating a value near a known value
        double tnew = tArray[0] + 0.0001
        double znew = zArray[0] + 0.0001
        double value = gsl_spline2d_eval(spline, tnew, znew, tacc, zacc);
        printf("interpolated value is %f \n", value);
      }
    }
  }
}

transformTensorToMilne(stressTensor, stressTensor)
{
  
}

transformVectorToMilne(baryonCurrent, baryonCurrent)
{

}
