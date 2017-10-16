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

        myfile << x << " " << y << " " << z << " " << var[is] << "\n";
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

        myfile << x << " " << y << " " << z << " " << var[idx][is] << "\n";
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
  for (int iy = 0; iy < DIM_Y; iy++)
  {
    for (int ix = 0; ix < DIM_X; ix++)
    {
      int iz = (DIM_Z - 1) / 2; // at z = 0
      int is = (DIM_Y * DIM_Z) * ix + (DIM_Z) * iy + iz; //the column packed index spanning x, y, z
      myfile << var[is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
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
  for (int iy = 0; iy < DIM_Y; iy++)
  {
    for (int ix = 0; ix < DIM_X; ix++)
    {
      int iz = (DIM_Z - 1) / 2; //at z = 0
      int is = (DIM_Y * DIM_Z) * ix + (DIM_Z) * iy + iz; //the column packed index spanning x, y, z
      myfile << var[idx][is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
  }
  myfile.close();
}
//this function writes a variable as a function of x at y = z = 0.
void writeVectorToFile1DProjection(float **var, char name[255], int idx)
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = 0; ix < DIM_X; ix++)
  {
    int iz = (DIM_Z - 1) / 2; //at z = 0
    int iy = (DIM_Y - 1) / 2; //at y = 0
    int is = (DIM_Y * DIM_Z) * ix + (DIM_Z) * iy + iz; //the column packed index spanning x, y, z
    myfile << var[idx][is] << "\n";
  }
  myfile.close();
}


//use this function to read in initial energy density or baryon density profile
//format should be whitespace delimited : x   y   z   value(x,y,z)
void readDensityFile(float *density, char name[255])
{
  float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
  float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
  float zmin = (-1.0) * ((float)(DIM_Z-1) / 2.0) * DZ;

  char filename[255] = "";
  sprintf(filename, "%s.dat", name);
  std::ifstream infile;
  infile.open(filename);
  for (int irow = 0; irow < DIM; irow++)
  {
    float x;
    float y;
    float z;
    float value;

    infile >> x;
    infile >> y;
    infile >> z;
    infile >> value; //can be either the energy density or baryon density

    int ix = (int)round((x - xmin) / DX);
    int iy = (int)round((y - ymin) / DY);
    int iz = (int)round((z - zmin) / DZ);
    int is = (DIM_Y * DIM_Z * ix) + (DIM_Z * iy) + iz;

    density[is] = value;
  }
  infile.close();
}
