/*
 * Reference for laplacian two-dimensional 4th order accurate (spatial) explicit method
 *
 * Original equation: U_xx + U_yy = 0
 * Solved by: u(t+1,x,y) = 
 *    (1/12u(t,x-2,y) - 4/3u(t,x-1,y) + 5/2(t,x,y) - 4/3u(t,x+1,y) + 1/12u(t,x+2,y))
 *    + (1/12u(t,x,y-2) - 4/3u(t,x,y-1) + 5/2(t,x,y) - 4/3u(t,x,y+1) + 1/12u(t,x,y+2))
 *
 * @author Brandon Nesterenko (bnestere@uccs.edu)
 * @date 9-16-2018
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
  int w, x, y, z;
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

  double e = 2.9e-5;
  double p = 5e-7;
  double ie = 1/e;


  /* allocate memory */
  double (*u_0_0)[x_max][y_max] = (double*) malloc (2 * x_max * y_max * sizeof (double));
  double (*c)[x_max][y_max] = (double*) malloc (10 * x_max * y_max * sizeof (double));

  alpha = 1.f / (double) x_max;
  beta = 2.f / (double) x_max;
  //beta = 2.f / (double) y_max;


  /* initialize the first timesteps */
  for (i = 0; i < x_max; i++)
  {
    for (i = 0; i < x_max; i++)
    {
      u_0_0[0][i][j] = 1. + i*0.1 + j*0.01;

      for(w = 0; w < 10; w++) {
        c[w][i][j] = w + i*0.1 + j*0.01;
      }
    }
  }

  /* do the calculation */ 
  benchInit();
  benchBeginStencil();
#pragma scop
  for (t = 0; t < T_MAX; t++)
  {
    int t0 = t%2;
    int tp1 = 1-t0;

    for (x = 2; x < x_max - 2; x++)
    {
      for (y = 2; y < y_max - 2; y++)
      {
        u_0_0[tp1][x][y] = (-p * ie)
          * ((c[0][x][y]*u_0_0[t0][x-2][y] + c[1][x][y]*u_0_0[t0][x-1][y] + c[2][x][y] * u_0_0[t0][x][y] + c[3][x][y] * u_0_0[t0][x+1][y] + c[4][x][y] * u_0_0[t0][x+2][y])
              + (c[5][x][y]*u_0_0[t0][x][y-2] + c[6][x][y]*u_0_0[t0][x][y-1] + c[7][x][y] * u_0_0[t0][x][y] + c[8][x][y] * u_0_0[t0][x][y+1] + c[9][x][y] *u_0_0[t0][x][y+2]));
      }
    }

  }
#pragma endscop

  benchEndStencil();
  benchSetEnv();
  benchSetDomain(2, 4, 0);
  benchSetProblemSize(x_max, y_max, z_max, T_MAX);
  benchSetArithProps(9,0,12,0);
  benchSetMatProps(1, 1);
  benchSetMemProps(11, 1);
  benchSetFpSize(sizeof(double));
  benchFinalize();

  /* clean up */
  free (u_0_0);

  return EXIT_SUCCESS;
}


