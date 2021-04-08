/*
 * Reference for laplacian one-dimensional 4th order accurate (spatial) explicit method
 *
 * Original equation: U_xx = 0
 * Solved by: u(t+1,x) = 
 *    (1/12u(t,x-2) - 4/3u(t,x-1) + 5/2(t,x) - 4/3u(t,x+1) + 1/12u(t,x+2))
 *
 * @author Brandon Nesterenko (bnestere@uccs.edu)
 * @date 8-26-2018
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


  /* allocate memory */
  double (*u_0_0)[x_max] = (double*) malloc (2 * x_max * sizeof (double));
  double (*c)[x_max] = (double*) malloc (5 * x_max * sizeof (double));

  alpha = 1.f / (double) x_max;
  beta = 2.f / (double) x_max;
  //beta = 2.f / (double) y_max;

  /* initialize the first timesteps */
  for (i = 0; i < x_max; i++)
  {
    u_0_0[0][i] = 1. + i*0.1;

    for(w = 0; w < 5; w++) {
      c[w][i] = w + i*0.1;
    }
  }

  /* do the calculation */ 
  benchInit();
  benchBeginStencil();

#pragma stencil
  for (t = 0; t < T_MAX; t++)
  {
    int t0 = t%2;
    int tp1 = 1-t0;

    for (x = 2; x < x_max - 2; x++)
    {
      u_0_0[tp1][x] = (c[0][x]*u_0_0[t0][x-2] -  c[1][x]*u_0_0[t0][x-1] + c[2][x]* u_0_0[t0][x] - c[3][x]*u_0_0[t0][x+1] + c[4][x]*u_0_0[t0][x+2]);
    }
  }

  benchEndStencil();
  benchSetEnv();
  benchSetDomain(1, 4, 1);
  benchSetProblemSize(x_max, y_max, z_max, T_MAX);
  benchSetArithProps(2,2,5,0);
  benchSetMatProps(6, 1);
  benchSetMemProps(2, 1);
  benchSetFpSize(sizeof(double));
  benchFinalize();

  /* clean up */
  free (u_0_0);

  return EXIT_SUCCESS;
}


