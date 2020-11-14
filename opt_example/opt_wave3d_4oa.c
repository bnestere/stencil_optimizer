#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stencil_config.h>
#define NONUMA
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#if defined(_OPENMP)
#endif
double seconds() {
   struct timeval tv;
   gettimeofday(&tv,NULL);
   return (double)(tv.tv_sec)+1e-6*tv.tv_usec;
}
int malloc_error(const char* err) {
   fprintf(stderr,"Failed to allocate the field %s.\n",err);
   return EXIT_FAILURE;
}
int validate_results(int x_max,int y_max,int z_max,int T_MAX,void**** _opt_res) {
   double tim1;double tim2;double nFlops;
   int i;int j;int k;int t;
   int x;int y;int z;
   const double MIN=-1.f;
   const double MAX=1.f;
   const double DX=(MAX-MIN)/(x_max-3);
   const double DT=DX/2.0f;
   const double DT_DX_SQUARE=DT*DT/(DX*DX);
   double(*u)[x_max][y_max][z_max] = (double*)malloc(3*x_max*y_max*z_max*sizeof (double));
   if (u==NULL) 
     {
        return malloc_error("u");
     }
   memset(u,0,3*x_max*sizeof (double));
   for (k=2; k<z_max-2; k+=1)
     {
        for (j=2; j<y_max-2; j+=1)
          {
             for (i=2; i<x_max-2; i+=1)
               {
                  double x=(i-1)*DX+MIN;
                  double y=(j-1)*DX+MIN;
                  double z=(k-1)*DX+MIN;
                  if (k==2) 
                    {
                       u[0][i][j][0] = 0;
                       u[0][i][j][1] = 0;
                       u[1][i][j][0] = 0;
                       u[1][i][j][1] = 0;
                    }
                  if (k==z_max-3) 
                    {
                       u[0][i][j][z_max-2] = 0;
                       u[0][i][j][z_max-1] = 0;
                       u[1][i][j][z_max-2] = 0;
                       u[1][i][j][z_max-1] = 0;
                    }
                  if (j==2) 
                    {
                       u[0][i][0][k] = 0;
                       u[0][i][1][k] = 0;
                       u[1][i][0][k] = 0;
                       u[1][i][1][k] = 0;
                    }
                  if (j==y_max-3) 
                    {
                       u[0][i][y_max-2][k] = 0;
                       u[0][i][y_max-1][k] = 0;
                       u[1][i][y_max-2][k] = 0;
                       u[1][i][y_max-1][k] = 0;
                    }
                  if (i==2) 
                    {
                       u[0][0][j][k] = 0;
                       u[0][1][j][k] = 0;
                       u[1][0][j][k] = 0;
                       u[1][1][j][k] = 0;
                    }
                  if (i==x_max-3) 
                    {
                       u[0][x_max-2][j][k] = 0;
                       u[0][x_max-1][j][k] = 0;
                       u[1][x_max-2][j][k] = 0;
                       u[1][x_max-1][j][k] = 0;
                    }
                  u[1][i][j][k] = (double)(sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));
                  u[0][i][j][k] = u[1][i][j][k];
               }
          }
     }
   double sc1=1.0/12;
   double sc2=4.0/3.0;
   double sc3=5.0/2.0;
   tim1 = seconds();
   for (t=0; t<T_MAX; t+=1)
     {
        for (i=2; i<x_max-2; i+=1)
          {
             for (j=2; j<y_max-2; j+=1)
               {
                  for (k=2; k<z_max-2; k+=1)
                    {
                       u[(t+2)%3][i][j][k] = -sc1*u[(t+1)%3][i-2][j][k]+sc2*u[(t+1)%3][i-1][j][k]-sc3*u[(t+1)%3][i][j][k]+sc2*u[(t+1)%3][i+1][j][k]-sc1*u[(t+1)%3][i+2][j][k]+(-sc1*u[(t+1)%3][i][j-2][k]+sc2*u[(t+1)%3][i][j-1][k]-sc3*u[(t+1)%3][i][j][k]+sc2*u[(t+1)%3][i][j+1][k]-sc1*u[(t+1)%3][i][j+2][k])+(-sc1*u[(t+1)%3][i][j][k-2]+sc2*u[(t+1)%3][i][j][k-1]-sc3*u[(t+1)%3][i][j][k]+sc2*u[(t+1)%3][i][j][k+1]-sc1*u[(t+1)%3][i][j][k+2])+2*u[(t+1)%3][i][j][k]-u[t%3][i][j][k];
                    }
               }
          }
     }
   tim2 = seconds();
   double(*opt_res)[x_max][y_max][z_max] = _opt_res;
   for (x=1; x<x_max-1; x+=1)
     {
        for (y=1; y<y_max-1; y+=1)
          {
             for (z=1; z<z_max-1; z+=1)
               {
                  if (u[2][x][y][z]!=opt_res[2][x][y][z]) 
                    {
                       fprintf(stderr,"Index (%d, %d, %d) differs by %.16lf",x,y,z,u[2][x][y][z]-opt_res[2][x][y][z]);
                       exit(EXIT_FAILURE);
                    }
               }
          }
     }
   printf("%s\n%s\n%s\n","---------------------------------------------------------","Opt stencil successfully validated against normal version","---------------------------------------------------------");
   nFlops = (double)((x_max-4)*(double)((y_max-4)*(double)((z_max-4)*T_MAX*35.0)));
   printf("ORIGINAL FLOPs in stencil code:      %e\n",nFlops);
   printf("ORIGINAL Time spent in stencil code: %f\n",tim2-tim1);
   printf("ORIGINAL Performance in GFlop/s:     %f\n",nFlops/(1e9*(tim2-tim1)));
   free(u);
   return 0;
}
int main(int argc,char** argv) {
   int i=0;int j=0;int k=0;int t;
   double tim1;double tim2;double nFlops;
   if (argc!=5) 
     {
        printf("Wrong number of parameters, <xmax> <ymax> <zmax> <timesteps>.\n",argv[0]);
        exit(-1);
     }
   int x_max=atoi(argv[1])+4;
   int y_max=atoi(argv[2])+4;
   int z_max=atoi(argv[3])+4;
   int T_MAX=atoi(argv[4]);
   const double MIN=-1.f;
   const double MAX=1.f;
   const double DX=(MAX-MIN)/(x_max-3);
   const double DT=DX/2.0f;
   const double DT_DX_SQUARE=DT*DT/(DX*DX);
   double(*u)[x_max][y_max][z_max] = (double*)malloc(3*x_max*y_max*z_max*sizeof (double));
   if (u==NULL) 
     {
        return malloc_error("u");
     }
   memset(u,0,3*x_max*sizeof (double));
   for (k=2; k<z_max-2; k+=1)
     {
        for (j=2; j<y_max-2; j+=1)
          {
             for (i=2; i<x_max-2; i+=1)
               {
                  double x=(i-1)*DX+MIN;
                  double y=(j-1)*DX+MIN;
                  double z=(k-1)*DX+MIN;
                  if (k==2) 
                    {
                       u[0][i][j][0] = 0;
                       u[0][i][j][1] = 0;
                       u[1][i][j][0] = 0;
                       u[1][i][j][1] = 0;
                    }
                  if (k==z_max-3) 
                    {
                       u[0][i][j][z_max-2] = 0;
                       u[0][i][j][z_max-1] = 0;
                       u[1][i][j][z_max-2] = 0;
                       u[1][i][j][z_max-1] = 0;
                    }
                  if (j==2) 
                    {
                       u[0][i][0][k] = 0;
                       u[0][i][1][k] = 0;
                       u[1][i][0][k] = 0;
                       u[1][i][1][k] = 0;
                    }
                  if (j==y_max-3) 
                    {
                       u[0][i][y_max-2][k] = 0;
                       u[0][i][y_max-1][k] = 0;
                       u[1][i][y_max-2][k] = 0;
                       u[1][i][y_max-1][k] = 0;
                    }
                  if (i==2) 
                    {
                       u[0][0][j][k] = 0;
                       u[0][1][j][k] = 0;
                       u[1][0][j][k] = 0;
                       u[1][1][j][k] = 0;
                    }
                  if (i==x_max-3) 
                    {
                       u[0][x_max-2][j][k] = 0;
                       u[0][x_max-1][j][k] = 0;
                       u[1][x_max-2][j][k] = 0;
                       u[1][x_max-1][j][k] = 0;
                    }
                  u[1][i][j][k] = (double)(sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));
                  u[0][i][j][k] = u[1][i][j][k];
               }
          }
     }
   double sc1=1.0/12;
   double sc2=4.0/3.0;
   double sc3=5.0/2.0;
   tim1 = seconds();
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((T_MAX >= 1) && (T_MAX <= 2147483646) && (x_max >= 5) && (y_max >= 5) && (z_max >= 5)) {
  for (t1=-1;t1<=floord(T_MAX-1,8);t1++) {
    lbp=max(ceild(t1,2),ceild(16*t1-T_MAX+2,16));
    ubp=min(floord(16*t1+x_max+12,32),floord(x_max+2*T_MAX-5,32));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(0,ceild(t1-1,2)),ceild(32*t2-x_max-25,32));t3<=min(min(min(floord(16*t1+y_max+27,32),floord(32*t2+y_max+25,32)),floord(y_max+2*T_MAX-5,32)),floord(32*t1-32*t2+y_max+x_max+25,32));t3++) {
        for (t4=max(max(max(0,ceild(t1-1,2)),ceild(32*t2-x_max-25,32)),ceild(32*t3-y_max-25,32));t4<=min(min(min(min(floord(16*t1+z_max+27,32),floord(32*t2+z_max+25,32)),floord(32*t3+z_max+25,32)),floord(z_max+2*T_MAX-5,32)),floord(32*t1-32*t2+z_max+x_max+25,32));t4++) {
          for (t5=max(max(max(max(max(0,ceild(32*t2-x_max+3,2)),ceild(32*t3-y_max+3,2)),ceild(32*t4-z_max+3,2)),8*t1),16*t1-16*t2+1);t5<=min(min(min(min(min(floord(32*t1-32*t2+x_max+28,2),T_MAX-1),8*t1+15),16*t2+14),16*t3+14),16*t4+14);t5++) {
            for (t6=max(max(32*t2,2*t5+2),-32*t1+32*t2+4*t5-31);t6<=min(min(32*t2+31,-32*t1+32*t2+4*t5),2*t5+x_max-3);t6++) {
              for (t7=max(32*t3,2*t5+2);t7<=min(32*t3+31,2*t5+y_max-3);t7++) {
                lbv=max(32*t4,2*t5+2);
                ubv=min(32*t4+31,2*t5+z_max-3);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  u[(t5 + 2) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = ((((((((((-sc1) * u[(t5 + 1) % 3][(-2*t5+t6) - 2][(-2*t5+t7)][(-2*t5+t8)]) + (sc2 * u[(t5 + 1) % 3][(-2*t5+t6) - 1][(-2*t5+t7)][(-2*t5+t8)])) - (sc3 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)])) + (sc2 * u[(t5 + 1) % 3][(-2*t5+t6) + 1][(-2*t5+t7)][(-2*t5+t8)])) - (sc1 * u[(t5 + 1) % 3][(-2*t5+t6) + 2][(-2*t5+t7)][(-2*t5+t8)])) + ((((((-sc1) * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7) - 2][(-2*t5+t8)]) + (sc2 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7) - 1][(-2*t5+t8)])) - (sc3 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)])) + (sc2 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7) + 1][(-2*t5+t8)])) - (sc1 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7) + 2][(-2*t5+t8)]))) + ((((((-sc1) * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8) - 2]) + (sc2 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8) - 1])) - (sc3 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)])) + (sc2 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8) + 1])) - (sc1 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8) + 2]))) + (2 * u[(t5 + 1) % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)])) - u[t5 % 3][(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]);;
                }
              }
            }
          }
        }
      }
    }
  }
}
   tim2 = seconds();
   validate_results(x_max,y_max,z_max,T_MAX,u);
   nFlops = (double)((x_max-4)*(double)((y_max-4)*(double)((z_max-4)*T_MAX*35.0)));
   printf("\nOPTIMIZED FLOPs in stencil code:      %e\n",nFlops);
   printf("OPTIMIZED Time spent in stencil code: %f\n",tim2-tim1);
   printf("OPTIMIZED Performance in GFlop/s:     %f\n",nFlops/(1e9*(tim2-tim1)));
   free(u);
   return EXIT_SUCCESS;
}
