#include <math.h>
void initializeZero(float *density)
{
  //#pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    density[is] = 0.0;
  }
}
//this doesnt seems to center the gaussian profile! whats wrong ?
void initializeGauss(float *density, float b) // b is the variance
{
  //#pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_Z);
    int iy = (is - (DIM_Y * DIM_Z * ix))/ DIM_Z;
    int iz = is - (DIM_Y * DIM_Z * ix) - (DIM_Z * iy);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float z = (float)iz * DZ  - ((float)(DIM_Z-1)) / 2.0 * DZ;

    density[is] = exp(-(1.0 / b) * ((x * x) + (y * y) + (z * z)));
  }
}
