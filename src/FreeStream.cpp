#include <math.h>

void freeStream(float *density, float ***shiftedDensity, float dt)
{
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    for (int ithetap = 0; ithetap < DIM_THETAP; ithetap++)
    {
      for (int iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        int ix = is / (DIM_Y * DIM_Z);
        int iy = (is - (DIM_Y * DIM_Z * ix))/ DIM_Z;
        int iz = is - (DIM_Y * DIM_Z * ix) - (DIM_Z * iy);

        float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
        float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
        float zmin = (-1.0) * ((float)(DIM_Z-1) / 2.0) * DZ;

        float x = (float)ix * DX  + xmin;
        float y = (float)iy * DY  + ymin;
        float z = (float)iz * DZ  + zmin;

        float thetap = float(ithetap) * (PI) / float(DIM_THETAP);
        float phip = float(iphip) * (2.0 * PI) / float(DIM_PHIP);

        //these trig functions could be tabulated ahead of time!
        float x_new = x - (dt * sin(thetap) * cos(phip));
        float y_new = y - (dt * sin(thetap) * sin(phip));
        float z_new = z - (dt * cos(thetap));

        int ix_new = (int)round((x_new - xmin) / DX);
        int iy_new = (int)round((y_new - ymin) / DY);
        int iz_new = (int)round((z_new - zmin) / DZ);

        int is_new = (DIM_Y * DIM_Z * ix_new) + (DIM_Z * iy_new) + iz_new;

        //prevent from going out of array bounds
        //note this may be causing problems! what happens when it goes out of array bounds?
        if ((ix_new >= 0) && (ix_new < DIM_X) && (iy_new >= 0) && (iy_new < DIM_Y) && (iz_new >= 0) && (iz_new < DIM_Z))
        {
          shiftedDensity[is][ithetap][iphip] = density[is_new];
        }
      }
    }
  }
}
void convertInitialDensity(float *initialDensity, float *density)
{
  float norm_factor = 1.0 / (4.0 * PI); //the normalization constant
  //relating the intial energy /baryon density to the intial density profile G^(0,0) or J^0
  for (int is = 0; is < DIM; is++)
  {
    density[is] = initialDensity[is] * norm_factor;
  }
}
