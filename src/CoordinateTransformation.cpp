#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <stdio.h>

void interpolateToMilneGrid(float ***timeDependentQuantity, float **quantityAtFixedTau, int dimVar)
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
        double timeDependentArray[DIM_T * DIM_Z];
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
        gsl_spline2d_init(spline, tArray, zArray, timeDependentArray, sizeof(tArray) / sizeof(double), sizeof(zArray) / sizeof(double));
        gsl_interp_accel *tacc = gsl_interp_accel_alloc();
        gsl_interp_accel *zacc = gsl_interp_accel_alloc();

        //evaluate interpolation at regular spaced points in eta and a fixed proper time
        double etamin = (-1.0) * ((double)(DIM_ETA-1) / 2.0) * DETA;
        for (int ieta = 0; ieta < DIM_ETA; ieta++)
        {
          //find the values of t,z corresponging to regularly spaced eta and fixed tau
          double eta = (double)ieta * DETA  + etamin;
          double tnew = TAU * cosh(eta);
          double znew = TAU * sinh(eta);

          //fill the new variable defined on an eta grid with interpolated values
          double value = gsl_spline2d_eval(spline, tnew, znew, tacc, zacc);
          int isnew = (DIM_Y * DIM_ETA * ix) + (DIM_ETA * iy) + ieta;
          quantityAtFixedTau[ivar][isnew] = (float)value;
        }
      }
    }
  }
}

void transformTensorToMilne(float **Tensor, float **newTensor)
{

}

void transformVectorToMilne(float **vector, float **newVector)
{

}
