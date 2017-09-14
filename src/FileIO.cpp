#include <unistd.h>
#include <stdio.h>
#include <fstream>
void writeVarToFile(float ***var, char name[255])
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int iz = 0; iz < DIM_Z; iz++)
      {
        float x = (float)ix * DX  - (((float)(DIM_X-1)) / 2.0 * DX);
        float y = (float)iy * DY  - (((float)(DIM_Y-1)) / 2.0 * DY);
        float z = (float)iz * DZ  - (((float)(DIM_Z-1)) / 2.0 * DZ);
        myfile << x << "\t" << y << "\t" << z << "\t" << var[ix][iy][iz] << "\n";
      }
    }
  }
  myfile.close();
}
