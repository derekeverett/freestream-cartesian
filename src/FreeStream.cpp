#include <math.h>

void freeStream(float ***density, float *****shiftedDensity, float xmin, float ymin, float zmin)
{
for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int iz = 0; iz < DIM_Z; iz++)
      {
        for (int ithetap = 0; ithetap < DIM_THETAP; ithetap++)
        {
          for (int iphip = 0; iphip < DIM_PHIP; iphip++)
          {
            float x, y, z, x_new, y_new, z_new, thetap, phip;
            int ix_new, iy_new, iz_new;

            x = xmin + (float(ix) * DX);
            y = ymin + (float(iy) * DY);
            z = zmin + (float(iz) * DZ);
            thetap = float(ithetap) * (PI) / float(DIM_THETAP);
            phip = float(iphip) * (2.0 * PI) / float(DIM_PHIP);

            //could instead tabulate sin(thetap) and cos(thetap) ahead of time
            // to reduce runtime but program is already fairly fast - may be uncessary
            x_new = x - (DT * sin(thetap) * cos(phip));
            y_new = y - (DT * sin(thetap) * sin(phip));
            z_new = z - (DT * cos(thetap));

            ix_new = (int)round((x_new - xmin) / DX);
            iy_new = (int)round((y_new - ymin) / DY);
            iz_new = (int)round((z_new - zmin) / DZ);

            //prevent from going out of array bounds
            if ((ix_new >= 0 && ix_new < DIM_X) &&  (iy_new >= 0 && iy_new < DIM_Y) && (iz_new >= 0 && iz_new < DIM_Z))
            {
              shiftedDensity[ix][iy][iz][ithetap][iphip] = density[ix_new][iy_new][iz_new];
            }
          }
        }
      }
    }
  }
}
void convertInitialDensity(float ***initialEnergyDensity, float ***density)
{
  float norm_factor = 1.0; //the normalization constant relating the intial energy density to the intial density profile H^(0,0)
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int iz = 0; iz < DIM_Z; iz++)
      {
        density[ix][iy][iz] = initialEnergyDensity[ix][iy][iz] / norm_factor;
      }
    }
  }
}
