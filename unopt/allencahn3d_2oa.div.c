/*
 * Reference for Allen Cahn two-dimensional 2nd order accurate (spatial) explicit method
 *
 * TODO: Fix this for Allen-Cahn
 * Original equation: U_t = U_xx
 * Solved by: u(t+1,x) = u(t,x)
 *  + 0.125 * (u(t,x-1) - 2(t,x) + u(t,x+1))
 *
 *
 * @author Brandon Nesterenko (bnestere@uccs.edu)
 * @date 1-22-2020
 */                                                                                                                  
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include <stencil_config.h>

#if defined(_OPENMP)
# include <omp.h>
#endif

#include <math.h>

/**
 * Get current time in seconds.
 */
double seconds ()
{
    struct timeval tv;
    gettimeofday (&tv, NULL);
    return ((double) tv.tv_sec) + 1e-6 * tv.tv_usec;
}

int main(int argc, char **argv) {

  int i,j,k,t;
  double nFlops;
  double tim1, tim2;
  double r_nuclei, r;
  int x_max, y_max, z_max, t_max;

  double dx = 0.5e-6;
  double dy = dx;
  double dz = dx;
  double idx = 1/dx;
  double idy = 1/dy;
  double idz = 1/dz;
  double eee = 5.0e+5;
  double sigma = 1.0;
  double delta = 4.0*dx;
  double amobi = 4.0e-14;
  double ram = 0.1;
  double bbb = 2.0*log((1.-(1.*ram))/(1.-(1.-2.*ram)))/2.;

// Phase-field parameters
  double aaa = sqrt(3.0*delta*sigma/bbb);
  double www = 6.0*sigma*bbb/delta;
  double pmobi = amobi*sqrt(2.0*www)/(6.0*aaa);
  double iwe = 1./(2.*www);

// time increment and time steps
  double dt = (dx*dx)/(5.0*pmobi*aaa*aaa)/2.0;

  x_max = atoi (argv[1]);
  y_max = atoi (argv[2]);
  z_max = atoi (argv[3]);
  t_max = atoi (argv[4]);
  /* allocate memory */
  double (*p)[x_max][y_max][z_max] = (double*) malloc (2 * x_max * y_max * z_max * sizeof (double));

  r_nuclei = 5.*dx;
  for(i = 0; i < x_max; i++) {
    for(j = 0; j < y_max; j++) {
      for(k = 0; k < z_max; k++) {
        r = sqrt(pow((i*dx),2) + pow((j*dy),2) + pow((k*dz),2)) - r_nuclei;
        p[0][i][j][k] = 0.5*(1.-tanh(sqrt(2.*www)/(2.*aaa)*r));
      }
    }
  }

//  tim1 = seconds();

  benchInit();
  benchBeginStencil();
#pragma scop
  for(t = 0; t < t_max; t++) {

    int t0 = t%2;
    int t1 = 1-t0;

    for(i = 1; i < x_max -1; i++) {
      for(j = 1; j < y_max -1; j++) {
        for(k = 1; k < z_max -1; k++) {
            //* (4.*www*p[t0][i][j][k] * (1.-p[t0][i][j][k])*(p[t0][i][j][k]-0.5+3.*iwe*eee)
          p[t1][i][j][k] = p[t0][i][j][k] + pmobi 
            * (4.*www*p[t0][i][j][k] * (1.-p[t0][i][j][k])*(p[t0][i][j][k]-0.5+3./(2.*www)*eee)
                + aaa*aaa
                *(  (p[t0][i+1][j][k] - 2*p[t0][i][j][k] + p[t0][i-1][j][k])/dx/dx
                  + (p[t0][i][j+1][k] - 2*p[t0][i][j][k] + p[t0][i][j-1][k])/dy/dy
                  + (p[t0][i][j][k+1] - 2*p[t0][i][j][k] + p[t0][i][j][k-1])/dz/dz) ) * dt;
        }
      }
    }
  }
#pragma endscop

 // tim2 = seconds();
  benchEndStencil();
  benchSetEnv();
  benchSetDomain(3, 2, 1);
  benchSetProblemSize(x_max, y_max, z_max, t_max);
  benchSetArithProps(8, 5, 20, 0);
  benchSetMemProps(6, 1);
  benchSetMatProps(2, 1);
  benchSetFpSize(sizeof(double));
  benchFinalize();

  /* print statistics */
//  nFlops = (double) (x_max-2) * (double) (y_max-2) * (double) (z_max-2) * t_max * 33.0;
//  printf ("FLOPs in stencil code:      %e\n", nFlops);
//  printf ("Time spent in stencil code: %f\n", tim2 - tim1);
//  printf ("Performance in GFlop/s:     %f\n", nFlops / (1e9 * (tim2 -tim1)));

  free(p);
  return EXIT_SUCCESS;
}
