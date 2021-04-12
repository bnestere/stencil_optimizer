#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <stencil_config.h>

double GetWallTime(void)
   {
      struct timeval tp;
      static long start=0, startu;
      if (!start)
      {
         gettimeofday(&tp, NULL);
         start = tp.tv_sec;
         startu = tp.tv_usec;
         return(0.0);
      }
      gettimeofday(&tp, NULL);
      return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
   }

#ifndef RANDSEED
#define RANDSEED 1
#endif

int main(int argc, char ** argv) 
{
  if (argc != 5) {
    printf ("Wrong number of parameters.\n", argv[0]);
    exit (-1);
  }

  long N0 = atoi (argv[1]); 
  long N1 = atoi (argv[2]);
  long N2 = atoi (argv[3]);
  long TS = atoi (argv[4]);

  double (*u)[N0]=(double*)malloc(2*(N0*sizeof (double)));
  double (*c)[N0]=(double*)malloc(5*(N0*sizeof (double)));;
#ifdef VALIDATE
  long tlast = 0;
  double (*u_test)[N0]=(double*)malloc(2*(N0*sizeof (double)));
  double (*c)[N0]=(double*)malloc(5*(N0*sizeof (double)));;
#endif
  
  long i0,i,t,l;
  long tnew;long tm1; 

  

  srand(RANDSEED);

  for (i0=0; i0<N0; i0+=1)
    {
       
       u[0][i0] = rand();
       #ifdef VALIDATE
       u_test[0][i0] = u[0][i0];
       #endif
       u[1][i0] = rand();
       #ifdef VALIDATE
       u_test[1][i0] = u[1][i0];
       #endif
       for (l=0; l<5; l+=1)
         {
            c[l][i0] = (double)fmod((double)rand(),2);
         }
    }

  benchInit();
  benchBeginStencil();

#pragma stencil
  for (t = 0; t < TS; t++) {
    
    tnew = (t+(2-1)) % 2;
    tm1 = (t+(2-1-1)) % 2; 

    for (i0=2; i0<-2+N0; i0+=1)
      {
         u[tnew][i0] = c[2][i0]*u[tm1][-2+i0]+(c[1][i0]*u[tm1][-1+i0])+(c[0][i0]*u[tm1][i0])+(c[3][i0]*u[tm1][1+i0])+(c[4][i0]*u[tm1][2+i0]);
      }
  }

  benchEndStencil();
  benchSetEnv();
  benchSetDomain(1,4,0);
  benchSetProblemSize(N0, 0, 0, TS);
  benchSetArithProps(4, 0, 5, 0);
  benchSetMemProps(6, 1);
  benchSetMatProps(1, 1); //Fix for (long distance mats, immediate mats)
  benchSetFpSize(sizeof(double));
  benchFinalize();


#ifdef VALIDATE
  for (t = 0; t < TS; t++) {
    
    tnew = (t+(2-1)) % 2;
    tm1 = (t+(2-1-1)) % 2; 

    for (i0=2; i0<-2+N0; i0+=1)
      {
         u_test[tnew][i0] = c[2][i0]*u_test[tm1][-2+i0]+(c[1][i0]*u_test[tm1][-1+i0])+(c[0][i0]*u_test[tm1][i0])+(c[3][i0]*u_test[tm1][1+i0])+(c[4][i0]*u_test[tm1][2+i0]);
      }

    tlast = tnew;
  }

  for (i0=0; i0<N0; i0+=1)
    {
       
       if (u[tlast][i0]!=u_test[tlast][i0]) 
         {
            fprintf(stderr,"Comparison Failed at index %d (%lf!=%lf)\n",i0,u[tlast][i0],u_test[tlast][i0]);
            exit(EXIT_FAILURE);
         }
    }
  free(u_test);
#endif

  free(u);
  free(c);

  return EXIT_SUCCESS;
}

