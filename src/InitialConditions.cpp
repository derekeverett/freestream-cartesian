#include <math.h>
void initializeZero(float ***density)
{
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int iz = 0; iz < DIM_Z; iz++)
      {
        density[ix][iy][iz] = 0.0;
      }
    }
  }
}
void initializeGauss(float ***density, float xmax, float ymax, float zmax, float b)
{
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int iz = 0; iz < DIM_Z; iz++)
      {
        float x = ((float)ix * DX) - (xmax / 2.0);
        float y = ((float)iy * DY) - (ymax / 2.0);
        float z = ((float)iz * DZ) - (zmax / 2.0);
        density[ix][iy][iz] = exp(-(1.0 / b) * ((x * x) + (y * y) + (z * z)));
      }
    }
  }
}
