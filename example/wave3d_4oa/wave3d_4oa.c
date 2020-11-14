/*
 * Reference for wave three-dimensional 4th order accurate (spatial) explicit method
 *
 * Original equation: U_tt = U_xx + U_yy + U_zz
 * Solved by: u(t+1,x,y) = 2u(t,x,y,z) - u(t-1,x,y,z)
 *  + (-1/12u(t,x-2,y,z) + 4/3u(t,x-1,y,z) - 5/2(t,x,y,z) + 4/3u(t,x+1,y,z) - 1/12u(t,x+2,y,z))  
 *  + (-1/12u(t,x,y-2,z) + 4/3u(t,x,y-1,z) - 5/2(t,x,y,z) + 4/3u(t,x,y+1,z) - 1/12u(t,x,y+2,z))  
 *  + (-1/12u(t,x,y,z-2) + 4/3u(t,x,y,z-1) - 5/2(t,x,y,z) + 4/3u(t,x,y,z+1) - 1/12u(t,x,y,z+2))  
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

#define NONUMA

#ifndef M_PI
#	define M_PI 3.14159265358979323846
#endif

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


int malloc_error (const char* err)
{
  fprintf (stderr, "Failed to allocate the field %s.\n", err);
  return EXIT_FAILURE;
}

int validate_results(int x_max, int y_max, int z_max, int T_MAX, void ****_opt_res) {

  double tim1, tim2, nFlops;
  int i, j, k, t;
  int x, y, z;
  
  const double MIN = -1.f; const double MAX = 1.f;
  const double DX = (MAX - MIN) / (x_max - 3);
  const double DT = DX / 2.0f;

  const double DT_DX_SQUARE = DT * DT / (DX * DX);

  /* allocate memory */
  double (*u)[x_max][y_max][z_max] = (double*) malloc (3 * x_max * y_max * z_max * sizeof (double));
  if (u == NULL)
  {
    return malloc_error ("u");
  }

  /* initialize the first timesteps */
  memset (u, 0,3* x_max * sizeof (double));

  for (k = 2; k < z_max - 2; k++)
  {
    for (j = 2; j < y_max - 2; j++)
    {
      for (i = 2; i < x_max - 2; i++)
      {
        double x = (i - 1) * DX + MIN;
        double y = (j - 1) * DX + MIN;
        double z = (k - 1) * DX + MIN;

        if (k == 2)
        {
          u[0][i][j][0] = 0;
          u[0][i][j][1] = 0;
          u[1][i][j][0] = 0;
          u[1][i][j][1] = 0;
        }
        if (k == z_max - 3)
        {
          u[0][i][j][z_max - 2] = 0;
          u[0][i][j][z_max - 1] = 0;
          u[1][i][j][z_max - 2] = 0;
          u[1][i][j][z_max - 1] = 0;
        }
        if (j == 2)
        {
          u[0][i][0][k] = 0;
          u[0][i][1][k] = 0;
          u[1][i][0][k] = 0;
          u[1][i][1][k] = 0;
        }
        if (j == y_max - 3)
        {
          u[0][i][y_max - 2][k] = 0;
          u[0][i][y_max - 1][k] = 0;
          u[1][i][y_max - 2][k] = 0;
          u[1][i][y_max - 1][k] = 0;
        }
        if (i == 2)
        {
          u[0][0][j][k] = 0;
          u[0][1][j][k] = 0;
          u[1][0][j][k] = 0;
          u[1][1][j][k] = 0;
        }
        if (i == x_max - 3)
        {
          u[0][x_max - 2][j][k] = 0;
          u[0][x_max - 1][j][k] = 0;
          u[1][x_max - 2][j][k] = 0;
          u[1][x_max - 1][j][k] = 0;
        }

        u[1][i][j][k] = (double) (sin (2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z));
        u[0][i][j][k] = u[1][i][j][k];
      }
    }
  }

  double sc1 = 1.0/12;
  double sc2 = 4.0/3.0;
  double sc3 = 5.0/2.0;

  /* do the calculation */ 
  tim1 = seconds();
  for (t = 0; t < T_MAX; t++)
  {
    for (i = 2; i < x_max - 2; i++) {
      for (j = 2; j < y_max - 2; j++) {
        for (k = 2; k < z_max - 2; k++) {
          u[(t+2)%3][i][j][k] = 
            (-sc1*u[(t+1)%3][i-2][j][k] +  sc2*u[(t+1)%3][i-1][j][k] - sc3* u[(t+1)%3][i][j][k] + sc2*u[(t+1)%3][i+1][j][k] - sc1*u[(t+1)%3][i+2][j][k])
            + (-sc1*u[(t+1)%3][i][j-2][k] +  sc2*u[(t+1)%3][i][j-1][k] - sc3* u[(t+1)%3][i][j][k] + sc2*u[(t+1)%3][i][j+1][k] - sc1*u[(t+1)%3][i][j+2][k])
            + (-sc1*u[(t+1)%3][i][j][k-2] +  sc2*u[(t+1)%3][i][j][k-1] - sc3* u[(t+1)%3][i][j][k] + sc2*u[(t+1)%3][i][j][k+1] - sc1*u[(t+1)%3][i][j][k+2])
            + (2*u[(t+1)%3][i][j][k]) - u[(t%3)][i][j][k];
        }
      }
    }
  }
  tim2 = seconds();

  double (*opt_res)[x_max][y_max][z_max] = _opt_res;
  for (x = 1; x < x_max - 1; x++)
  {
    for (y = 1; y < y_max - 1; y++)
    {
      for (z = 1; z < z_max - 1; z++)
      {
        if(u[2][x][y][z] != opt_res[2][x][y][z]) {
          fprintf( stderr, "Index (%d, %d, %d) differs by %.16lf", x, y, z, u[2][x][y][z]-opt_res[2][x][y][z]); 
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  printf("%s\n%s\n%s\n",
      "---------------------------------------------------------",
      "Opt stencil successfully validated against normal version",
      "---------------------------------------------------------");

  /* print statistics */    
  nFlops = (double) (x_max-4) * (double) (y_max-4)* (double)(z_max-4) * T_MAX * 35.0;
  printf ("ORIGINAL FLOPs in stencil code:      %e\n", nFlops);    
  printf ("ORIGINAL Time spent in stencil code: %f\n", tim2 - tim1);
  printf ("ORIGINAL Performance in GFlop/s:     %f\n", nFlops / (1e9 * (tim2 -tim1)));

  /* clean up */
  free (u);
  return 0;
}

/**
 * Do the calculation.
 */
int main(int argc, char** argv)
{
  int i=0, j=0, k=0, t;
  double tim1, tim2, nFlops;

  if (argc != 5)
  {
    printf ("Wrong number of parameters, <xmax> <ymax> <zmax> <timesteps>.\n", argv[0]);
    exit (-1);
  }

  int x_max = atoi (argv[1]) + 4;
  int y_max = atoi (argv[2]) + 4;
  int z_max = atoi (argv[3]) + 4;
  int T_MAX = atoi (argv[4]);

  const double MIN = -1.f; const double MAX = 1.f;
  const double DX = (MAX - MIN) / (x_max - 3);
  const double DT = DX / 2.0f;

  const double DT_DX_SQUARE = DT * DT / (DX * DX);


  /* allocate memory */
  double (*u)[x_max][y_max][z_max] = (double*) malloc (3 * x_max * y_max * z_max * sizeof (double));
  if (u == NULL)
  {
    return malloc_error ("u");
  }

  /* initialize the first timesteps */
  memset (u, 0,3* x_max * sizeof (double));

  for (k = 2; k < z_max - 2; k++)
  {
    for (j = 2; j < y_max - 2; j++)
    {
      for (i = 2; i < x_max - 2; i++)
      {
        double x = (i - 1) * DX + MIN;
        double y = (j - 1) * DX + MIN;
        double z = (k - 1) * DX + MIN;

        if (k == 2)
        {
          u[0][i][j][0] = 0;
          u[0][i][j][1] = 0;
          u[1][i][j][0] = 0;
          u[1][i][j][1] = 0;
        }
        if (k == z_max - 3)
        {
          u[0][i][j][z_max - 2] = 0;
          u[0][i][j][z_max - 1] = 0;
          u[1][i][j][z_max - 2] = 0;
          u[1][i][j][z_max - 1] = 0;
        }
        if (j == 2)
        {
          u[0][i][0][k] = 0;
          u[0][i][1][k] = 0;
          u[1][i][0][k] = 0;
          u[1][i][1][k] = 0;
        }
        if (j == y_max - 3)
        {
          u[0][i][y_max - 2][k] = 0;
          u[0][i][y_max - 1][k] = 0;
          u[1][i][y_max - 2][k] = 0;
          u[1][i][y_max - 1][k] = 0;
        }
        if (i == 2)
        {
          u[0][0][j][k] = 0;
          u[0][1][j][k] = 0;
          u[1][0][j][k] = 0;
          u[1][1][j][k] = 0;
        }
        if (i == x_max - 3)
        {
          u[0][x_max - 2][j][k] = 0;
          u[0][x_max - 1][j][k] = 0;
          u[1][x_max - 2][j][k] = 0;
          u[1][x_max - 1][j][k] = 0;
        }

        u[1][i][j][k] = (double) (sin (2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z));
        u[0][i][j][k] = u[1][i][j][k];
      }
    }
  }

  double sc1 = 1.0/12;
  double sc2 = 4.0/3.0;
  double sc3 = 5.0/2.0;

  /* do the calculation */ 
  tim1 = seconds();

#pragma stencil
  for (t = 0; t < T_MAX; t++) {
    for (i = 2; i < x_max - 2; i++) {
      for (j = 2; j < y_max - 2; j++) {
        for (k = 2; k < z_max - 2; k++) {
          u[(t+2)%3][i][j][k] = 
            (-sc1*u[(t+1)%3][i-2][j][k] +  sc2*u[(t+1)%3][i-1][j][k] - sc3* u[(t+1)%3][i][j][k] + sc2*u[(t+1)%3][i+1][j][k] - sc1*u[(t+1)%3][i+2][j][k])
            + (-sc1*u[(t+1)%3][i][j-2][k] +  sc2*u[(t+1)%3][i][j-1][k] - sc3* u[(t+1)%3][i][j][k] + sc2*u[(t+1)%3][i][j+1][k] - sc1*u[(t+1)%3][i][j+2][k])
            + (-sc1*u[(t+1)%3][i][j][k-2] +  sc2*u[(t+1)%3][i][j][k-1] - sc3* u[(t+1)%3][i][j][k] + sc2*u[(t+1)%3][i][j][k+1] - sc1*u[(t+1)%3][i][j][k+2])
            + (2*u[(t+1)%3][i][j][k]) - u[(t%3)][i][j][k];
        }
      }
    }
  }

  tim2 = seconds();

  validate_results(x_max, y_max, z_max, T_MAX, u);

  /* print statistics */    
  nFlops = (double) (x_max-4) * (double) (y_max-4)* (double)(z_max-4) * T_MAX * 35.0;
  printf ("\nOPTIMIZED FLOPs in stencil code:      %e\n", nFlops);    
  printf ("OPTIMIZED Time spent in stencil code: %f\n", tim2 - tim1);
  printf ("OPTIMIZED Performance in GFlop/s:     %f\n", nFlops / (1e9 * (tim2 -tim1)));

  /* clean up */
  free (u);

  return EXIT_SUCCESS;
}


