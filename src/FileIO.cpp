#include <unistd.h>
#include <stdio.h>
#include <fstream>

void writeScalarToFile(float *var, char name[255])
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int iz = 0; iz < DIM_Z; iz++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ix = 0; ix < DIM_X; ix++)
      {
        float x = (float)ix * DX  - (((float)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX);
        float y = (float)iy * DY  - (((float)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        float z = (float)iz * DZ  - (((float)(DIM_Z-1)) / 2.0 * DZ);
        z = DZ * roundf(z / DZ);

        int is = (DIM_Y * DIM_Z) * ix + (DIM_Z) * iy + iz; //the column packed index spanning x, y, z

        myfile << x << "\t" << y << "\t" << z << "\t" << var[is] << "\n";
      }
    }
  }
  myfile.close();
}

void writeVectorToFile(float **var, char name[255], int idx)
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int iz = 0; iz < DIM_Z; iz++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ix = 0; ix < DIM_X; ix++)
      {
        float x = (float)ix * DX  - (((float)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX); //rounding for regularly spaced values
        float y = (float)iy * DY  - (((float)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        float z = (float)iz * DZ  - (((float)(DIM_Z-1)) / 2.0 * DZ);
        z = DZ * roundf(z / DZ);

        int is = (DIM_Y * DIM_Z) * ix + (DIM_Z) * iy + iz; //the column packed index spanning x, y, z

        myfile << x << "\t" << y << "\t" << z << "\t" << var[idx][is] << "\n";
      }
    }
  }
  myfile.close();
}

//this function writes the transverse density of a variable at z = 0
// as regularly spaced values
void writeScalarToFileProjection(float *var, char name[255])
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
      int iz = (DIM_Z - 1) / 2; // at z = 0
      int is = (DIM_Y * DIM_Z) * ix + (DIM_Z) * iy + iz; //the column packed index spanning x, y, z
      myfile << var[is] << "\t"; //different columns for y values
    }
    myfile << "\n"; // different rows correspond to different x values
  }
  myfile.close();
}

void writeVectorToFileProjection(float **var, char name[255], int idx)
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
      int iz = (DIM_Z - 1) / 2; //at z = 0
      int is = (DIM_Y * DIM_Z) * ix + (DIM_Z) * iy + iz; //the column packed index spanning x, y, z
      myfile << var[idx][is] << "\t"; //different columns for y values
    }
    myfile << "\n"; // different rows correspond to different x values
  }
  myfile.close();
}
