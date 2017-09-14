//trigTable is a table with 10 rows for ten combinations or p_(mu)p_(nu) normalized by the energy
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
void calculateTrigTable(float ***trigTable)
{
  for (int ithetap = 0; ithetap < DIM_THETAP; ithetap++)
  {
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float thetap = float(ithetap) * (PI) / float(DIM_THETAP);
      float phip = float(iphip) * (2.0 * PI) / float(DIM_PHIP);
      trigTable[0][ithetap][iphip] = 1.0; //p_t, p_t component
      trigTable[1][ithetap][iphip] = sin(thetap) * cos(phip); //p_t, p_x
      trigTable[2][ithetap][iphip] = sin(thetap) * sin(phip); //p_t, p_y
      trigTable[3][ithetap][iphip] = cos(thetap); //p_t, p_z
      trigTable[4][ithetap][iphip] = sin(thetap) * cos(phip) * sin(thetap) * cos(phip); //p_x, p_x
      trigTable[5][ithetap][iphip] = sin(thetap) * cos(phip) * sin(thetap) * sin(phip); //p_x, p_y
      trigTable[6][ithetap][iphip] = sin(thetap) * cos(phip) * cos(thetap); //p_x, p_z
      trigTable[7][ithetap][iphip] = sin(thetap) * sin(phip) * sin(thetap) * sin(phip); //p_y, p_y
      trigTable[8][ithetap][iphip] = sin(thetap) * sin(phip) * cos(thetap); //p_y, p_z
      trigTable[9][ithetap][iphip] = cos(thetap) * cos(thetap); //p_z, p_z
    }
  }
}
void calculateStressTensor(float ****stressTensor, float *****shiftedDensity, float ***trigTable)
{
  float d_thetap = PI / float(DIM_THETAP);
  float d_phip = (2.0 * PI) / float(DIM_PHIP);
  for (int ivar = 0; ivar < 10; ivar++)
  {
    for (int ix = 0; ix < DIM_X; ix++)
    {
      for (int iy = 0; iy < DIM_Y; iy++)
      {
        for (int iz = 0; iz < DIM_Z; iz++)
        {
          for (int ithetap = 0; ithetap < DIM_THETAP; ithetap++)
          {
            float thetap = float(ithetap) * (PI) / float(DIM_THETAP);
            for (int iphip = 0; iphip < DIM_PHIP; iphip++)
            {
              //rather than gauss quadrature, just doing a elementary Riemann sum here; check convergence!
              stressTensor[ivar][ix][iy][iz] += shiftedDensity[ix][iy][iz][ithetap][iphip] * trigTable[ivar][ithetap][iphip] * sin(thetap) * d_thetap * d_phip;
            }
          }
        }
      }
    }
  }
}
void solveEigenSystem(float ****stressTensor, float ***energyDensity, float ****flowVelocity)
{
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int iz = 0; iz < DIM_Z; iz++)
      {
        gsl_matrix *Tmunu; //T^(mu,nu) with two contravariant indices; we need to lower an index
        //using the metric to find the eigenvectors of T^(mu)_(nu) with one contravariant and one contravariant index
        Tmunu = gsl_matrix_alloc(4,4);
        gsl_matrix *gmunu;
        gmunu = gsl_matrix_alloc(4,4);
        gsl_matrix_complex *eigen_vectors;
        eigen_vectors = gsl_matrix_complex_alloc(4,4);
        gsl_vector_complex *eigen_values;
        eigen_values = gsl_vector_complex_alloc(4);
        //set the values of the energy momentum tensor
        gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][ix][iy][iz]); //tt
        gsl_matrix_set(Tmunu, 0, 1, stressTensor[1][ix][iy][iz]); //tx
        gsl_matrix_set(Tmunu, 0, 2, stressTensor[2][ix][iy][iz]); //ty
        gsl_matrix_set(Tmunu, 0, 3, stressTensor[3][ix][iy][iz]); //tz
        gsl_matrix_set(Tmunu, 1, 1, stressTensor[4][ix][iy][iz]); //xx
        gsl_matrix_set(Tmunu, 1, 2, stressTensor[5][ix][iy][iz]); //xy
        gsl_matrix_set(Tmunu, 1, 3, stressTensor[6][ix][iy][iz]); //xz
        gsl_matrix_set(Tmunu, 2, 2, stressTensor[7][ix][iy][iz]); //yy
        gsl_matrix_set(Tmunu, 2, 3, stressTensor[8][ix][iy][iz]); //yz
        gsl_matrix_set(Tmunu, 3, 3, stressTensor[9][ix][iy][iz]); //zz
        gsl_matrix_set(Tmunu, 1, 0, stressTensor[1][ix][iy][iz]); //xt
        gsl_matrix_set(Tmunu, 2, 0, stressTensor[2][ix][iy][iz]); //yt
        gsl_matrix_set(Tmunu, 3, 0, stressTensor[3][ix][iy][iz]); //zt
        gsl_matrix_set(Tmunu, 2, 1, stressTensor[5][ix][iy][iz]); //yx
        gsl_matrix_set(Tmunu, 3, 1, stressTensor[6][ix][iy][iz]); //zx
        gsl_matrix_set(Tmunu, 3, 2, stressTensor[8][ix][iy][iz]); //zy

        //set the values of the "metric"; not really the metric, but the numerical constants
        //which are multiplied by the elements of T^(mu,nu) to get the values of T^(mu)_(nu)
        gsl_matrix_set(gmunu, 0, 0, 1.0); //tt
        gsl_matrix_set(gmunu, 0, 1, -1.0); //tx
        gsl_matrix_set(gmunu, 0, 2, -1.0); //ty
        gsl_matrix_set(gmunu, 0, 3, -1.0); //tz
        gsl_matrix_set(gmunu, 1, 0, 1.0); //xt
        gsl_matrix_set(gmunu, 1, 1, -1.0); //xx
        gsl_matrix_set(gmunu, 1, 2, -1.0); //xy
        gsl_matrix_set(gmunu, 1, 3, -1.0); //xz
        gsl_matrix_set(gmunu, 2, 0, 1.0); //yt
        gsl_matrix_set(gmunu, 2, 1, -1.0); //yx
        gsl_matrix_set(gmunu, 2, 2, -1.0); //yy
        gsl_matrix_set(gmunu, 2, 3, -1.0); //yz
        gsl_matrix_set(gmunu, 3, 0, 1.0); //zt
        gsl_matrix_set(gmunu, 3, 1, -1.0); //zx
        gsl_matrix_set(gmunu, 3, 2, -1.0); //zy
        gsl_matrix_set(gmunu, 3, 3, -1.0); //zz
        //lower one index of the stress tensor; save it to the same matrix to save memory
        gsl_matrix_mul_elements(Tmunu, gmunu); //result stored in Tmunu !this multiplies element-wise, not ordinary matrix multiplication!
        gsl_eigen_nonsymmv_workspace *eigen_workspace;
        eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
        gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
        gsl_eigen_nonsymmv_free(eigen_workspace);
        for (int i = 0; i < 4; i++)
        {
          gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);

          if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0.0)
          {
            double v0 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 0));
            double v1 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 1));
            double v2 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 2));
            double v3 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 3));
            //double euclideanLength = v0*v0 + v1*v1 + v2*v2 + v3*v3; //gsl normalizes eigenvectors to euclideanLength = 1; this is just a check
            double minkowskiLength = v0*v0 - (v1*v1 + v2*v2 + v3*v3); //we want to flow velocity normalized s.t. minkowskiLength = 1
            double scaleFactor = 1.0 / sqrt(minkowskiLength); //so we need to scale all the elements of the eigenvector by scaleFactor
            v0 = scaleFactor * v0;
            v1 = scaleFactor * v1;
            v2 = scaleFactor * v2;
            v3 = scaleFactor * v3;
            //printf("scaled eigenvector %d is (%f ,%f , %f, %f) and eigenvalue %d is %f\n", i, v0, v1, v2, v3, i, GSL_REAL(eigenvalue));
            //set values of energy density and flow velocity
            energyDensity[ix][iy][iz] = GSL_REAL(eigenvalue);
            flowVelocity[0][ix][iy][iz] = v0;
            flowVelocity[1][ix][iy][iz] = v1;
            flowVelocity[2][ix][iy][iz] = v2;
            flowVelocity[3][ix][iy][iz] = v3;
          }
        }
      }
    }
  }
}