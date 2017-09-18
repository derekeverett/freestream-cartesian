#include <math.h>
void initializeZero(float *density)
{
  //#pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    density[is] = 0.0;
  }
}
void initializeGauss(float *density, float xmax, float ymax, float zmax, float b)
{
  //#pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_Z);
    int iy = (is - (DIM_Y * DIM_Z * ix))/ DIM_Z;
    int iz = is - (DIM_Y * DIM_Z * ix) - (DIM_Z * iy);

    float x = ((float)ix * DX) - (xmax / 2.0);
    float y = ((float)iy * DY) - (ymax / 2.0);
    float z = ((float)iz * DZ) - (zmax / 2.0);
    density[is] = exp(-(1.0 / b) * ((x * x) + (y * y) + (z * z)));
  }
}
