/*
 * Reference for wave two-dimensional 4th order accurate (spatial) explicit method
 *
 * Original equation: U_tt = U_xx + U_yy
 * Solved by: u(t+1,x,y) = 2u(t,x,y) - u(t-1,x,y)
 *  + (-1/12u(t,x-2,y) + 4/3u(t,x-1,y) - 5/2(t,x,y) + 4/3u(t,x+1,y) - 1/12u(t,x+2,y))  
 *  + (-1/12u(t,x,y-2) + 4/3u(t,x,y-1) - 5/2(t,x,y) + 4/3u(t,x,y+1) - 1/12u(t,x,y+2))  
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

/**
 * Do the calculation.
 */
int main(int argc, char** argv)
{
    int i=0, j=0, k=0, t;
    double tim1, tim2, nFlops;


	if (argc != 5)
	{
		printf ("Wrong number of parameters. Syntax:\n%s <x_max> <y_max> <z_max>\n", argv[0]);
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
  double (*u)[x_max][y_max] = (double*) malloc(3 * x_max * y_max * sizeof(double));
    if (u== NULL)
    {
        free (u);
        return malloc_error ("u");
    }

    /* initialize the first timesteps */
#ifdef NONUMA
	memset (u, 0,3 * x_max * sizeof (double));
#endif

  for (j = 2; j < y_max - 2; j++)
  {
    for (i = 2; i < x_max - 2; i++)
    {
      double x = (i - 1) * DX + MIN;
      double y = (j - 1) * DX + MIN;

#ifndef NONUMA
      if (j == 2)
      {
        u[0][i][0] = 0;
        u[0][i][1] = 0;
        u[1][i][0] = 0;
        u[1][i][1] = 0;
      }
      if (j == y_max - 3)
      {
        u[0][i][y_max - 2] = 0;
        u[0][i][y_max - 1] = 0;
        u[1][i][y_max - 2] = 0;
        u[1][i][y_max - 1] = 0;
      }
      if (i == 2)
      {
        u[0][0][j] = 0;
        u[0][1][j] = 0;
        u[1][0][j] = 0;
        u[1][1][j] = 0;
      }
      if (i == x_max - 3)
      {
        u[0][x_max - 2][j] = 0;
        u[0][x_max - 1][j] = 0;
        u[1][x_max - 2][j] = 0;
        u[1][x_max - 1][j] = 0;
      }
#endif

      u[1][i][j] = (double) (sin (2 * M_PI * x) * sin(2 * M_PI * y));
      u[0][i][j] = u[1][i][j];
    }
  }
	


    double sc1 = 1.0/12;
    double sc2 = 4.0/3.0;
    double sc3 = 5.0/2.0;

    /* do the calculation */ 
    //tim1 = seconds();
    benchInit();
    benchBeginStencil();
#pragma stencil
    for (t = 0; t < T_MAX; t++)
    {
      int tnew = (t+2)%3;
      int tm1 = (t+1)%3;
      int tm2 = (t)%3;
      for (i = 2; i < x_max - 2; i++)
      {
        for (j = 2; j < y_max - 2; j++)
        {
          u[tnew][i][j] = 
            (-sc1*u[tm1][i-2][j] +  sc2*u[tm1][i-1][j] - sc3* u[tm1][i][j] + sc2*u[tm1][i+1][j] - sc1*u[tm1][i+2][j])
            + (-sc1*u[tm1][i][j-2] +  sc2*u[tm1][i][j-1] - sc3* u[tm1][i][j] + sc2*u[tm1][i][j+1] - sc1*u[tm1][i][j+2])
            + (2*u[tm1][i][j]) - u[tm2][i][j];
        }
    }


	}

    benchEndStencil();
    benchSetEnv();
    benchSetDomain(2, 4, 2);
    benchSetProblemSize(x_max, y_max, z_max, T_MAX);
    benchSetArithProps(6,7,11,0);
    benchSetMemProps(7, 1);
    benchSetMatProps(2, 1);
    benchSetFpSize(sizeof(double));
    benchFinalize();

	free (u);
	
	return EXIT_SUCCESS;
}


