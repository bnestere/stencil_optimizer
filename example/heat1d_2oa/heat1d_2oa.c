/*
 * Reference for heat one-dimensional 2nd order accurate (spatial) explicit method
 *
 * Original equation: U_t = U_xx
 * Solved by: u(t+1,x) = u(t,x)
 *  + 0.125 * (u(t,x-1) - 2(t,x) + u(t,x+1))
 *
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

int validate_results(int x_max, int y_max, int z_max, int T_MAX, void **_opt_res) {

  double tim1, tim2, nFlops;
  int i, j, k, t;
  int x, y, z;

  double alpha, beta;

  /* allocate memory */
  double (*u_0_0)[x_max] = (double*) malloc (2 * x_max * sizeof (double));

  alpha = 1.f / (double) x_max;
  beta = 2.f / (double) y_max;

  /* initialize the first timesteps */
  for (i = 0; i < x_max; i++)
  {
    u_0_0[0][i] = 1. + i*0.1;
  }

  /* do the calculation */ 
  tim1 = seconds();
  for (t = 0; t < T_MAX; t++)
  {
    for (x = 1; x < x_max - 1; x++)
    {
      u_0_0[1-t%2][x] = u_0_0[t%2][x] +
        0.125 * (u_0_0[t%2][x+1] - 2*u_0_0[t%2][x] + u_0_0[t%2][x-1]);
    }
  }
  tim2 = seconds();

  double (*opt_res)[x_max] = _opt_res;
  for (x = 1; x < x_max - 1; x++)
  {
    if(u_0_0[1][x] != opt_res[1][x]) {
      fprintf( stderr, "Index %d differs by %.16lf", x, u_0_0[1][x]-opt_res[1][x]); 
      exit(EXIT_FAILURE);
    }
  }

  printf("%s\n%s\n%s\n",
      "---------------------------------------------------------",
      "Opt stencil successfully validated against normal version",
      "---------------------------------------------------------");

  nFlops = (double) (x_max-2) * T_MAX * 5.0;
  printf ("ORIGINAL FLOPs in stencil code:      %e\n", nFlops);    
  printf ("ORIGINAL Time spent in stencil code: %f\n", tim2 - tim1);
  printf ("ORIGINAL Performance in GFlop/s:     %f\n", nFlops / (1e9 * (tim2 -tim1)));

  free(u_0_0);

  return 0;
}

/**
 * Do the calculation.
 */
int main(int argc, char** argv)
{
  int x_max, y_max, z_max;
  int i, j, k, t;
  int x, y, z;
  int T_MAX;
  double tim1, tim2, nFlops;

  double alpha, beta;

  if (argc != 5)
  {
    printf ("Wrong number of parameters, <xmax> <ymax> <zmax> <timesteps>.\n", argv[0]);
    exit (-1);
  }

  x_max = atoi (argv[1]);
  y_max = atoi (argv[2]);
  z_max = atoi (argv[3]);
  T_MAX = atoi (argv[4]);

  /* allocate memory */
  double (*u_0_0)[x_max] = (double*) malloc (2 * x_max * sizeof (double));

  alpha = 1.f / (double) x_max;
  beta = 2.f / (double) y_max;

  /* initialize the first timesteps */
  for (i = 0; i < x_max; i++)
  {
    u_0_0[0][i] = 1. + i*0.1;
  }

  /* do the calculation */ 
  tim1 = seconds();
#pragma stencil
  for (t = 0; t < T_MAX; t++)
  {
    for (x = 1; x < x_max - 1; x++)
    {
      u_0_0[1-t%2][x] = u_0_0[t%2][x] +
        0.125 * (u_0_0[t%2][x+1] - 2*u_0_0[t%2][x] + u_0_0[t%2][x-1]);
    }
  }
  tim2 = seconds ();

  validate_results(x_max, y_max, z_max, T_MAX, u_0_0);

  /* print statistics */    
  nFlops = (double) (x_max-2) * T_MAX * 5.0;
  printf ("\nOPTIMIZED FLOPs in stencil code:      %e\n", nFlops);    
  printf ("OPTIMIZED Time spent in stencil code: %f\n", tim2 - tim1);
  printf ("OPTIMIZED Performance in GFlop/s:     %f\n", nFlops / (1e9 * (tim2 -tim1)));

  /* clean up */
  free (u_0_0);

  return EXIT_SUCCESS;
}


