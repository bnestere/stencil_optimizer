/*
 * Reference for ISO one-dimensional 2nd order accurate (spatial) explicit method 
 * with variable coefficients
 *
 * Original equation: U_t = U_xx
 * Solved by: u(t+1,x) = 
 *    c[0][x] * u(t,x)
 *  + c[1][x] * u(t,x-1) 
 *  + c[2][x] * u(t,x+1)
 *
 *
 * @author Brandon Nesterenko (bnestere@uccs.edu)
 * @date 5-19-2020
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stencil_config.h>

#if defined(_OPENMP)
#	include <omp.h>
#endif

#define N_COEFFS 4


/**
 * Get current time in seconds.
 */
double seconds ()
{
    struct timeval tv;
    gettimeofday (&tv, NULL);
    return ((double) tv.tv_sec) + 1e-6 * tv.tv_usec;
}

/**
 * Do the calculation.
 */
int main(int argc, char** argv)
{
    int x_max, y_max, z_max;
    int i, j, k, t;
    int x, y, z;
    int w;
    int T_MAX;
    double tim1, tim2, nFlops;
    double alpha, beta;

	if (argc != 5)
	{
		printf ("Wrong number of parameters.\n", argv[0]);
		exit (-1);
	}
	
	x_max = atoi (argv[1]);
	y_max = atoi (argv[2]);
	z_max = atoi (argv[3]);
  T_MAX = atoi (argv[4]);
 
    /* allocate memory */
    double (*u_0_0)[x_max][y_max][z_max] = (double*) malloc (2 * x_max * y_max * z_max * sizeof (double));

    double (*c)[x_max][y_max][z_max] = (double*) malloc (N_COEFFS * x_max * y_max * z_max * sizeof (double));

    alpha = 1.f / (double) x_max;
    beta = 2.f / (double) y_max;

    /* initialize the first timesteps */
    for (i = 0; i < x_max; i++)
    {
      for (j = 0; j < y_max; j++)
      {
        for (k = 0; k < z_max; k++)
        {
          u_0_0[0][i][j][k] = 1. + i*0.1 + j*0.01 + k*0.001;

          // Initialize constants
          for(w = 0; w < N_COEFFS; w++) {
            c[w][i][j][k] = w + i*0.1 + j*0.01 + k*0.001;
          }
        }
      }
    }
	
    /* do the calculation */ 
	//tim1 = seconds();
  benchInit();
  benchBeginStencil();
#pragma stencil
	for (t = 0; t < T_MAX; t++)
	{
    int idp0 = t % 2;
    int idp1 = 1-idp0;
    for (x = 1; x < x_max - 1; x++)
    {
      for (y = 1; y < y_max - 1; y++)
      {
        for (z = 1; z < z_max - 1; z++)
        {
          u_0_0[idp1][x][y][z] = 
              c[0][x][y][z] * u_0_0[idp0][x][y][z]
            + c[1][x][y][z] * (u_0_0[idp0][x-1][y][z] + u_0_0[idp0][x+1][y][z])
            + c[2][x][y][z] * (u_0_0[idp0][x][y-1][z] + u_0_0[idp0][x][y+1][z])
            + c[3][x][y][z] * (u_0_0[idp0][x][y][z-1] + u_0_0[idp0][x][y][z+1]);
        }
      }
    }
	}

  benchEndStencil();
  benchSetEnv();
  benchSetDomain(3, 2, 1);
  benchSetProblemSize(x_max, y_max, z_max, T_MAX);
  benchSetArithProps(6,0,4,0);
  benchSetMemProps(10, 1);
  benchSetFpSize(sizeof(double));
  benchSetMatProps(N_COEFFS+1, 1);
  benchFinalize();

    /* clean up */
	free (u_0_0);
  free (c);
	
	return EXIT_SUCCESS;
}


